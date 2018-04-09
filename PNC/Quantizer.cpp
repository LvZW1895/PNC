#include "stdafx.h"
#include "Quantizer.h"


Quantizer::Quantizer()
{
}


Quantizer::~Quantizer()
{
}


Matrix<complex<double>, Dynamic, Dynamic> Quantizer::quantize(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig,int bits,double maxvalue,double minvalue)
{
	Matrix<complex<double>, Dynamic, Dynamic> res(Rx_sig.rows(),Rx_sig.cols());
	int qnum = pow(2, bits / 2);
	double deta = (maxvalue - minvalue) / qnum;
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
			double tmpreal = Rx_sig(i, 0).real();
			double tmpimag = Rx_sig(i, 0).imag();
			tmpreal = tmpreal > maxvalue ? maxvalue-deta/2 : tmpreal;
			tmpreal = tmpreal < minvalue ? minvalue+deta/2 : tmpreal;
			tmpimag = tmpimag > maxvalue ? maxvalue-deta/2 : tmpimag;
			tmpimag = tmpimag < minvalue ? minvalue+deta/2 : tmpimag;
			for (int j = 1; j <= qnum; j++)
			{
				if (tmpreal > (minvalue + (j - 1)*deta) && tmpreal <= (minvalue + j*deta))
				{
					tmpreal = minvalue + j*deta - deta / 2;
					break;
				}
			}
			for (int j = 1; j <= qnum; j++)
			{
				if (tmpimag > (minvalue + (j - 1)*deta) && tmpimag <= (minvalue + j*deta))
				{
					tmpimag = minvalue + j*deta - deta / 2;
					break;
				}
			}
			res(i, 0) = complex<double>(tmpreal, tmpimag);		
			//res(i, 1) = Rx_sig(i, 1);
	}
	for (int i = 0; i < Rx_sig.rows(); i++)
	{
		double tmpreal = Rx_sig(i, 1).real();
		double tmpimag = Rx_sig(i, 1).imag();
		tmpreal = tmpreal > maxvalue ? maxvalue - deta / 2 : tmpreal;
		tmpreal = tmpreal < minvalue ? minvalue + deta / 2 : tmpreal;
		tmpimag = tmpimag > maxvalue ? maxvalue - deta / 2 : tmpimag;
		tmpimag = tmpimag < minvalue ? minvalue + deta / 2 : tmpimag;
		for (int j = 1; j <= qnum; j++)
		{
			if (tmpreal >(minvalue + (j - 1)*deta) && tmpreal <= (minvalue + j*deta))
			{
				tmpreal = minvalue + j*deta - deta / 2;
				break;
			}
		}
		for (int j = 1; j <= qnum; j++)
		{
			if (tmpimag > (minvalue + (j - 1)*deta) && tmpimag <= (minvalue + j*deta))
			{
				tmpimag = minvalue + j*deta - deta / 2;
				break;
			}
		}
		res(i, 1) = complex<double>(tmpreal, tmpimag);
		//res(i, 1) = Rx_sig(i, 1);
	}
	/*cout << Rx_sig << endl;
	cout << "------------" << endl;*/
	//cout << res << endl;
	return res;
}
