// PNC.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "Modulation.h"
#include "Channel.h"
#include "Demapping.h"
int main()
{
	ofstream log;
	int type =CoMP;
	int modtype = QPSK;
	switch (type)
	{
	case CoMP:
		log.open("comp_result.txt", ios::trunc | ios::out);
		break;
	case PNC:
		log.open("pnc_result.txt", ios::trunc | ios::out);
		break;
	default:
		break;
	}
	
	int block_num = 100;
	int msg_len = 10000;
	double Es = 1;
	int total_bit = 0;
	int error = 0;
	Modulation mod;
	Channel channel;
	Demapping dp;
	VectorXi EsN0dB=VectorXi::LinSpaced(21,0,20);
	VectorXd EsN0(21);
	for (int i = 0; i < 21; i++)
	{
		EsN0(i) = pow(10.0,0.1*EsN0dB(i));
	}
	VectorXd N0 = Es*EsN0.array().pow(-1);
	VectorXd sqrt_N0 = N0.array().sqrt();
	VectorXd sigma2 = N0 / 2;
	VectorXd sigma = sigma2.array().sqrt();
	VectorXi msg_u1(msg_len), msg_u2(msg_len);
	MatrixXi A(2, 2);

	A << 1, 1,
		0, 1;
	cout << "#EsN0dB       #err         SER" << endl;
	//cout << sigma << endl;
	//log << "#EsN0dB     SER" << endl;
	for (int k = 0; k < EsN0dB.size(); k++)
	{
		error = 0;
		total_bit = 0;
		for (int i = 0; i < block_num; i++)
		{
			//======================================
			// generate message bits
			//======================================
			VectorXd msg_u1_tmp = VectorXd::Random(msg_len) + VectorXd::Constant(msg_len, 1.0);
			VectorXd msg_u2_tmp = VectorXd::Random(msg_len) + VectorXd::Constant(msg_len, 1.0);
			for (int j = 0; j< msg_len; j++)
			{
				msg_u1(j) = (int)msg_u1_tmp(j);
				msg_u2(j) = (int)msg_u2_tmp(j);
			}
			//======================================
			// modulation
			//======================================
			Matrix<complex<double>, Dynamic, Dynamic> cv_txsig_u1, cv_txsig_u2;
			switch (modtype)
			{
			case BPSK:
				// cv_txsug_u1 = mod.mod_BPSK(msg_u1);
				// cv_txsug_u2 = mod.mod_BPSK(msg_u2);
				break;
			case QPSK:
				cv_txsig_u1 = mod.mod_QPSK(msg_u1);
				cv_txsig_u2 = mod.mod_QPSK(msg_u2);
				break;
			default:
				break;
			}
			//cout << msg_u1 << endl;
			//cout << "--------------" << endl;
			//cout << cv_txsig_u1 << endl;
			//======================================
			// Channel
			//======================================
			Matrix<complex<double>, Dynamic, Dynamic> Rx_sig = channel.AWGN(sigma(k), cv_txsig_u1, cv_txsig_u2, A);
			//======================================
			// PNC Demapping
			//======================================
			MatrixXd res;
			if (type == CoMP)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_comp_bpsk(Rx_sig, A);
					break;
				case QPSK:
					res = dp.dp_comp_qpsk(Rx_sig, A);
					break;
				default:
					break;
				}
				
			}
			else if (type == PNC)
			{
				switch (modtype)
				{
				case BPSK:
					//res = dp.dp_pnc_bpsk(Rx_sig);
					break;
				case QPSK:
					  res = dp.dp_pnc_qpsk(Rx_sig);
				default:
					break;
				}
			}
			//cout << msg_u1 << endl;
			//cout << "--------------" << endl;
			//cout << msg_u2 << endl;
			//cout << "--------------" << endl;
			//cout << res << endl;
			
			VectorXd res_u1 = res.col(0);
			VectorXd res_u2 = res.col(1);
			
			//======================================
			// BER calculation
			//======================================
			for (int j = 0; j < res_u1.size(); j++)
			{
				if (res_u1(j) != msg_u1(j))
					error++;
				if (res_u2(j) != msg_u2(j))
					error++;
			}
			total_bit += 2 * msg_len;
			//cout << "--------------" << endl;
			//cout << error << endl;
		}

		std::printf(" %3i        %6i      %1.3e\n", EsN0dB(k), error, ((double)error) / total_bit);
		log << EsN0dB(k)<<"        " << ((double)error) / total_bit << endl;
	}
	log.close();
	
    return 0;
}
