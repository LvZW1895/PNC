#include "stdafx.h"
#include "Channel.h"


Channel::Channel()
{
}


Channel::~Channel()
{
}

double RandFloat()
{
	// Generate a uniform distributed varable in the interval (0,1).   
	return (double)(rand() + 1) / (double)(RAND_MAX + 2);
}

/*MatrixXd Channel::AWGN(double sigma, VectorXi Tx_sig1, VectorXi Tx_sig2, MatrixXi A)
{
	MatrixXd rx(Tx_sig1.size(), 2);
	MatrixXd noise(Tx_sig1.size(), 2);
	VectorXd temp1(Tx_sig1.size()), temp2(Tx_sig2.size());
	for (int i = 0; i < Tx_sig1.size(); i++)
	{
		noise(i,0)= sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		noise(i,1) =sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	}
	temp1 = A(0, 0)*Tx_sig1.cast<double>() + A(0, 1)*Tx_sig2.cast<double>();
	temp2 = A(1, 0)*Tx_sig1.cast<double>() + A(1, 1)*Tx_sig2.cast<double>();
	rx << temp1, temp2;
	rx = rx + noise;
	//cout << rx<< endl;
	//cout << "----------" << endl;
	return rx;
}*/

Matrix<complex<double>, Dynamic, Dynamic> Channel::AWGN(double sigma, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig1, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig2, MatrixXi A)
{
	Matrix<complex<double>, Dynamic, Dynamic> rx(Tx_sig1.rows(), 2);
	Matrix<complex<double>, Dynamic, Dynamic> noise(Tx_sig1.rows(),2);
	Matrix<complex<double>, Dynamic, Dynamic> temp1(Tx_sig1.rows(), Tx_sig1.cols()), temp2(Tx_sig2.rows(), Tx_sig2.cols());
	for (int i = 0; i < Tx_sig1.rows(); i++)
	{
		double n0= sqrt(2) / 2 * sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		double n1= sqrt(2) / 2 * sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		complex<double> n0_tmp(n0, n0), n1_tmp(n1, n1);
		noise(i, 0) = n0_tmp;
		noise(i, 1) = n1_tmp;
	}

	temp1 = A(0, 0)*Tx_sig1 + A(0, 1)*Tx_sig2;
	temp2 = A(1, 0)*Tx_sig1 + A(1, 1)*Tx_sig2;
	rx << temp1, temp2;
	rx = rx + noise;
	//cout << rx<< endl;
	//cout << "----------" << endl;
	return rx;
}
