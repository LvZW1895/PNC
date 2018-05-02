#include "stdafx.h"
#include "Demapping.h"
#include <map>

Demapping::Demapping()
{
}


Demapping::~Demapping()
{
}


MatrixXi Demapping::dp_pnc_bpsk(MatrixXd Rx_sig)
{
	MatrixXi res(Rx_sig.rows(), Rx_sig.cols());
	VectorXi s1xors2(Rx_sig.rows());
	VectorXi s1(Rx_sig.rows());
	VectorXi s2(Rx_sig.rows());
	//======================================
	// R1 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 0) > -1.0&&Rx_sig(i,0)<=1.0)
		{
			s1xors2(i) = 1;
		}
		else
		{
			s1xors2(i) = 0;
		}
	}
	//======================================
	// R2 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 1) >= 0.0)
		{
			s2(i) = 1;
		}
		else
		{
			s2(i) = 0;
		}
	}
	for (int i = 0; i < s1.size(); i++)
	{
		s1(i) = s1xors2(i) ^ s2(i);
	}
	
	res << s1, s2;
	//cout << s1<< endl;
	//cout << "--------" << endl;
	return res;
}


MatrixXi Demapping::dp_comp_bpsk(MatrixXd Rx_sig, MatrixXi A)
{
	MatrixXi res(Rx_sig.rows(), Rx_sig.cols());
	MatrixXi bit(2, 4);
	MatrixXi bit_bpsk(2, 4);
	bit << 0, 0, 1, 1,
		   0, 1, 0, 1;	  
	bit_bpsk << -1, -1, 1, 1,
	         	-1, 1, -1, 1;
	MatrixXi res_tmp=A*bit_bpsk;
	
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE*1.0;
		double a1 = Rx_sig(i, 0);
		double a2 = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(a1 - res_tmp(0, j), 2) + pow(a2 - res_tmp(1, j), 2);
			if (dis < min)
			{
				min = dis;
				res.row(i) = bit.col(j);
			}
		}
	}
	//cout << Rx_sig<<endl;
	//cout << "-------------" << endl;
	//cout << res<<endl;
	return res;
}

MatrixXd Demapping::dp_comp_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	MatrixXd res(Rx_sig.rows()*2, 2);
	MatrixXd bit(4, 16);
	Matrix<complex<double>, Dynamic, Dynamic> bit_bpsk(2,16);
	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
		   0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
	
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit_bpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		        c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;

	/*Matrix<complex<double>, Dynamic, Dynamic> A_temp(2, 2);
	A_temp << A(0, 0), A(0, 1),
		      A(1, 0), A(1, 1);*/
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_bpsk;
	//cout << res_tmp << endl;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE*1.0;
		complex<double> a1 = Rx_sig(i, 0);
		complex<double> a2 = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(a1.real() - res_tmp(0, j).real(), 2) 
				+ pow(a1.imag() - res_tmp(0, j).imag(), 2)
				+ pow(a2.real() - res_tmp(1, j).real(), 2)
				+ pow(a2.imag() - res_tmp(1, j).imag(), 2);
			/*double dis = abs(a1.real() - res_tmp(0, j).real())
				+ abs(a1.imag() - res_tmp(0, j).imag())
				+ abs(a2.real() - res_tmp(1, j).real())
				+ abs(a2.imag() - res_tmp(1, j).imag());*/
			if (dis < min)
			{
				min = dis;
				res(2*i,0) = bit(0,j);
				res(2*i+1,0) = bit(1,j);
				res(2 * i, 1) = bit(2, j);
				res(2 * i + 1, 1) = bit(3, j);
			}
		}
	}
	//cout << Rx_sig<<endl;
	//cout << "-------------" << endl;
	//cout << res<<endl;
	return res;
}


MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	VectorXd s1xors2(Rx_sig.rows()*2);
	VectorXd s1(Rx_sig.rows()*2);
	VectorXd s2(Rx_sig.rows()*2);
	//======================================
	// R1 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 0).real() <= 1 && Rx_sig(i, 0).real() > -1)
		{
			s1xors2(2 * i) = 1;
		}
		else
		{
			s1xors2(2 * i) = 0;
		}
		if (Rx_sig(i, 0).imag() <= 1 && Rx_sig(i, 0).imag() > -1)
		{
			s1xors2(2 * i+1) = 1;
		}
		else
		{
			s1xors2(2 * i+1) = 0;
		}
	}
	//======================================
	// R2 Demapping
	//======================================
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		if (Rx_sig(i, 1).real() <= 0)
		{
			s2(2 * i) = 0;
		}
		else
		{
			s2(2 * i) = 1;
		}
		if (Rx_sig(i, 1).imag() <= 0)
		{
			s2(2 * i+1) = 0;
		}
		else
		{
			s2(2 * i+1) = 1;
		}

	}
	for (int i = 0; i < s1.size(); i++)
	{
		s1(i) = (int)s1xors2(i) ^(int)s2(i);
	}

	res << s1, s2;
	//cout << res<< endl;
	//cout << "--------" << endl;
	return res;

}

//MatrixXd Demapping::dp_llr_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A,double sigma,int bits)
//{
//	MatrixXd res(Rx_sig.rows() * 2, 2);
//	MatrixXd llr1(Rx_sig.rows() * 2, 2),llr2(Rx_sig.rows() * 2, 2),llr(Rx_sig.rows() * 2, 2);
//	Matrix<complex<double>, Dynamic, Dynamic> bit11_qpsk(2, 16), bit12_qpsk(2, 16), bit21_qpsk(2, 16), bit22_qpsk(2, 16);
//	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
//	bit11_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
//		          c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
//	bit12_qpsk << c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4,
//		          c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
//	bit21_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
//		          c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4;
//	bit22_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
//		          c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp11 = A*bit11_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp12 = A*bit12_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp21 = A*bit21_qpsk;
//	Matrix<complex<double>, Dynamic, Dynamic> res_tmp22 = A*bit22_qpsk;
//	for (int i = 0; i < Rx_sig.rows(); i++)
//	{
//		complex<double> a1 = Rx_sig(i, 0);
//		complex<double> a2 = Rx_sig(i, 1);
//		//==========
//		//LLR1
//		//==========
//		double sum1 = 0, sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2) / 2 * abs(pow(a1.real() - res_tmp11(0, j).real(), 2)+pow(a1.imag() - res_tmp11(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2) / 2 * abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp12(0, j).real(), 2)+pow(a1.imag() - res_tmp12(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i + 1, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp21(0, j).real(), 2)+ pow(a1.imag() - res_tmp21(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i, 1) = log(sum1 / sum0);
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp22(0, j).real(), 2)+pow(a1.imag() - res_tmp22(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr1(2 * i + 1, 1) = log(sum1 / sum0);
//		//======
//		//LLR2
//		//======
//		sum1 = 0, sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i + 1, 0) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i, 1) = log(sum1 / sum0);
//		sum1 = 0;
//		sum0 = 0;
//		for (int j = 0; j < 8; j++)
//		{
//			sum0 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		for (int j = 8; j < 16; j++)
//		{
//			sum1 += exp(-1 / pow(sigma, 2)/2*abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2))) / 2 / PI / pow(sigma, 2);
//		}
//		llr2(2 * i + 1, 1) = log(sum1 / sum0);
//
//	}
//	//
//	//quantize llr1
//	//
//		MatrixXd llr1_q = llr1;
//		int qnum = pow(2, bits);
//		double maxvalue = 400.0, minvalue = -400.0;
//		double deta = (maxvalue - minvalue) / qnum;
//		for (int i = 0; i < llr1.rows(); i++)
//		{
//			llr1_q(i, 0) = llr1_q(i, 0) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 0);
//			llr1_q(i, 0) = llr1_q(i, 0) < minvalue ? minvalue + deta / 2 : llr1_q(i, 0);
//			llr1_q(i, 1) = llr1_q(i, 1) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 1);
//			llr1_q(i, 1) = llr1_q(i, 1) < minvalue ? minvalue + deta / 2 : llr1_q(i, 1);
//			for (int j = 1; j <= qnum; j++)
//			{
//				if (llr1_q(i, 0) >(minvalue + (j - 1)*deta) && llr1_q(i, 0) <= (minvalue + j*deta))
//				{
//					llr1_q(i, 0) = minvalue + j*deta - deta / 2;
//					break;
//				}
//			}
//			for (int j = 1; j <= qnum; j++)
//			{
//				if (llr1_q(i, 1) > (minvalue + (j - 1)*deta) && llr1_q(i, 1) <= (minvalue + j*deta))
//				{
//					llr1_q(i, 1) = minvalue + j*deta - deta / 2;
//					break;
//				}
//			}
//		}
//	//
//	//decode 
//	//
//	llr = llr1_q+llr2;
//	for (int i = 0; i < llr.rows(); i++)
//	{
//		res(i, 0) = llr(i, 0) >= 0 ? 1 : 0;
//		res(i, 1) = llr(i, 1) >= 0 ? 1 : 0;
//	}
//	
//	//cout << llr1_q << endl;
//	//cout << "---------"<< endl;
//	return res;
//}
MatrixXd Demapping::dp_llr_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A, double sigma, int bits)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	MatrixXd llr1(Rx_sig.rows() * 2, 2), llr2(Rx_sig.rows() * 2, 2), llr(Rx_sig.rows() * 2, 2);
	Matrix<complex<double>, Dynamic, Dynamic> bit11_qpsk(2, 16), bit12_qpsk(2, 16), bit21_qpsk(2, 16), bit22_qpsk(2, 16);
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit11_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	bit12_qpsk << c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4,
		c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	bit21_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
		c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4;
	bit22_qpsk << c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4,
		c1, c1, c1, c1, c3, c3, c3, c3, c2, c2, c2, c2, c4, c4, c4, c4;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp11 = A*bit11_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp12 = A*bit12_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp21 = A*bit21_qpsk;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp22 = A*bit22_qpsk;
	//double maxvalue = MIN_VALUE, minvalue = MAX_VALUE;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		complex<double> a1 = Rx_sig(i, 0);
		complex<double> a2 = Rx_sig(i, 1);
		//==========
		//LLR1
		//==========
		double sum1 = 0, sum0 = 0,min=MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s= abs(pow(a1.real() - res_tmp11(0, j).real(), 2) + pow(a1.imag() - res_tmp11(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i, 0) = (sum0-sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp12(0, j).real(), 2) + pow(a1.imag() - res_tmp12(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i + 1, 0) = (sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp21(0, j).real(), 2) + pow(a1.imag() - res_tmp21(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i , 1) = (sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a1.real() - res_tmp22(0, j).real(), 2) + pow(a1.imag() - res_tmp22(0, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr1(2 * i + 1, 1) = (sum0 - sum1);

		//======
		//LLR2
		//======
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s =abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp11(1, j).real(), 2) + pow(a2.imag() - res_tmp11(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i, 0) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp12(1, j).real(), 2) + pow(a2.imag() - res_tmp12(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i + 1, 0) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp21(1, j).real(), 2) + pow(a2.imag() - res_tmp21(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i, 1) = 1 / sigma*(sum0 - sum1);
		sum1 = 0;
		sum0 = 0;
		min = MAX_VALUE;
		for (int j = 0; j < 8; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum0 = tmp_s;
				min = tmp_s;
			}
		}
		min = MAX_VALUE;
		for (int j = 8; j < 16; j++)
		{
			double tmp_s = abs(pow(a2.real() - res_tmp22(1, j).real(), 2) + pow(a2.imag() - res_tmp22(1, j).imag(), 2));
			if (tmp_s < min)
			{
				sum1 = tmp_s;
				min = tmp_s;
			}
		}
		llr2(2 * i + 1, 1) = 1 / sigma*(sum0 - sum1);

	}
	//
	//quantizer llr1
	//
	MatrixXd llr1_q = llr1;
	int qnum = pow(2, bits);
	double maxvalue = 5.0, minvalue = -5.0;
	double deta = (maxvalue - minvalue) / qnum;
	for (int i = 0; i < llr1.rows(); i++)
	{
		llr1_q(i, 0) = llr1_q(i, 0) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 0);
		llr1_q(i, 0) = llr1_q(i, 0) < minvalue ? minvalue + deta / 2 : llr1_q(i, 0);
		llr1_q(i, 1) = llr1_q(i, 1) > maxvalue ? maxvalue - deta / 2 : llr1_q(i, 1);
		llr1_q(i, 1) = llr1_q(i, 1) < minvalue ? minvalue + deta / 2 : llr1_q(i, 1);
		for (int j = 1; j <= qnum; j++)
		{
			if (llr1_q(i, 0) >(minvalue + (j - 1)*deta) && llr1_q(i, 0) <= (minvalue + j*deta))
			{
				llr1_q(i, 0) = minvalue + j*deta - deta / 2;
				break;
			}
		}
		for (int j = 1; j <= qnum; j++)
		{
			if (llr1_q(i, 1) > (minvalue + (j - 1)*deta) && llr1_q(i, 1) <= (minvalue + j*deta))
			{
				llr1_q(i, 1) = minvalue + j*deta - deta / 2;
				break;
			}
		}
	}






	//
	//decode
	//
	llr = llr1_q+llr2;
	for (int i = 0; i < llr.rows(); i++)
	{
		res(i, 0) = llr(i, 0) >= 0 ? 1 : 0;
		res(i, 1) = llr(i, 1) >= 0 ? 1 : 0;
	}

	return res;
}

MatrixXd Demapping::dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	MatrixXd res(Rx_sig.rows() * 2, 2);
	VectorXd r1_s1(Rx_sig.rows() * 2);
	VectorXd r1_s2(Rx_sig.rows() * 2);
	VectorXd r2_s1(Rx_sig.rows() * 2);
	VectorXd r2_s2(Rx_sig.rows() * 2);
	MatrixXd bit(4, 16);
	Matrix<complex<double>, Dynamic, Dynamic> bit_qpsk(2, 16);
	MatrixXd ncm(9, 4);
	bit << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
		   0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
		   0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
	complex<double> c1(-1, -1), c2(-1, 1), c3(1, -1), c4(1, 1);
	bit_qpsk << c1, c1, c1, c1, c2, c2, c2, c2, c3, c3, c3, c3, c4, c4, c4, c4,
		        c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4, c1, c2, c3, c4;
	Matrix<complex<double>, Dynamic, Dynamic> res_tmp = A*bit_qpsk;


	//
	//Computational network coding matrix
	//
	ncm << 1, 0, 0, 0,
		   0, 1, 0, 0,
		   0, 0, 1, 0,
		   0, 0, 0, 1,
		   1, 0, 1, 0,
		   0, 1, 0, 1,
		   1, 0, 0, 1,
		   0, 1, 1, 0,
		   1, 1, 0, 1;
	Matrix<double, 2, 4> nc[72];
	int n = 0;
	for (int i = 0; i < ncm.rows(); i++)
	{
		for (int j = 0; j < ncm.rows(); j++)
		{
			if (i != j)
			{
				nc[n].row(0) = ncm.row(i);
				nc[n].row(1) = ncm.row(j);
				n++;
			}
		}
	}
	//r1 matrix
	double maxdis = 0;
	int index1 = 0, index2 = 0;
	for (int k = 0; k < 72; k++)
	{
		double min = MAX_VALUE;
		MatrixXd pnc = nc[k] * bit;
		for (int i = 0; i < pnc.rows(); i++)
		{
			for (int j = 0; j < pnc.cols(); j++)
			{
				
					pnc(i, j) = (int)pnc(i,j)%2;
				
			}
		}
		for (int i = 0; i < 16; i++)
		{
			for (int j = 0; j < 16; j++)
			{
				if (i != j)
				{
					double dic = pow(res_tmp(0, i).real() - res_tmp(0, j).real(), 2) + pow(res_tmp(0, i).imag() - res_tmp(0, j).imag(), 2);
					if (pnc.col(i) != pnc.col(j))
					{
						if (dic < min)
						{
							min = dic;
						}
					}
				}
			}
		}
		if (maxdis < min)
		{
			maxdis = min;
			index1 = k;
		}
		
	}
	//r2 matrix
	maxdis = 0;
	for (int k = 0; k < 72; k++)
	{
		double min = MAX_VALUE;
		MatrixXd pnc = nc[k] * bit;
		for (int i = 0; i < pnc.rows(); i++)
		{
			for (int j = 0; j < pnc.cols(); j++)
			{

				pnc(i, j) = (int)pnc(i, j) % 2;

			}
		}
		for (int i = 0; i < 16; i++)
		{
			for (int j = 0; j < 16; j++)
			{
				if (i != j)
				{
					double dic = pow(res_tmp(1, i).real() - res_tmp(1, j).real(), 2) + pow(res_tmp(1, i).imag() - res_tmp(1, j).imag(), 2);
					if (pnc.col(i) != pnc.col(j))
					{
						if (dic < min)
						{
							min = dic;
						}
					}
				}
			}
		}
		if (maxdis < min)
		{
			maxdis = min;
			index2 = k;
		}

	}
	//
	//decode s1,s2
	//
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 0);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(0, j).real() - a.real(), 2) + pow(res_tmp(0, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r1_s1(i * 2) = bit(0, j);
				r1_s1(i * 2 + 1) = bit(1, j);
				r1_s2(i * 2) = bit(2, j);
				r1_s2(i * 2 + 1) = bit(3, j);
			}
		}
	}
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double min = MAX_VALUE;
		complex<double> a = Rx_sig(i, 1);
		for (int j = 0; j < res_tmp.cols(); j++)
		{
			double dis = pow(res_tmp(1, j).real() - a.real(), 2) + pow(res_tmp(1, j).imag() - a.imag(), 2);
			if (dis < min)
			{
				min = dis;
				r2_s1(i * 2) = bit(0, j);
				r2_s1(i * 2 + 1) = bit(1, j);
				r2_s2(i * 2) = bit(2, j);
				r2_s2(i * 2 + 1) = bit(3, j);
			}
		}
	}

	//
	//network coding
	//
	MatrixXd pnc_a(4, 4);
	pnc_a << nc[index1],
		     nc[index2];
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		VectorXd r1(4), r2(4);
		r1(0) = r1_s1(2 * i);
		r1(1) = r1_s1(2 * i + 1);
		r1(2) = r1_s2(2 * i);
		r1(3) = r1_s2(2 * i + 1);
		r2(0) = r2_s1(2 * i);
		r2(1) = r2_s1(2 * i + 1);
		r2(2) = r2_s2(2 * i);
		r2(3) = r2_s2(2 * i + 1);
		VectorXd tmp1(2), tmp2(2), tmp(4);
		tmp1 = nc[index1] * r1;
		tmp2 = nc[index2] * r2;
		tmp1(0) = (int)tmp1(0) % 2;
		tmp1(1) = (int)tmp1(1) % 2;
		tmp2(0) = (int)tmp2(0) % 2;
		tmp2(1) = (int)tmp2(1) % 2;
		tmp << tmp1,
			   tmp2;
		VectorXd res_v = pnc_a*tmp;
		res(2 * i, 0) = (int)res_v(0)%2;
		res(2 * i + 1, 0) = (int)res_v(1)%2;
		res(2 * i, 1) = (int)res_v(2)%2;
		res(2 * i + 1, 1) = (int)res_v(3)%2;
	}
	
	
	return res;
}


