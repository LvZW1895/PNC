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
	//srand((unsigned)time(NULL));
	return (double)(rand() + 1) / (double)(RAND_MAX + 2);
}


double randomRayleigh(double sigma)
{
	double pv = sigma*sqrt(-2 * log(RandFloat()));
	return pv;
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

Matrix<complex<double>, Dynamic, Dynamic> Channel::AWGN(double sigma, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig1, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig2, Matrix<complex<double>, Dynamic, Dynamic> A)
{
	Matrix<complex<double>, Dynamic, Dynamic> rx(Tx_sig1.rows(), 2);
	Matrix<complex<double>, Dynamic, Dynamic> noise(Tx_sig1.rows(),2);
	Matrix<complex<double>, Dynamic, Dynamic> temp1(Tx_sig1.rows(), Tx_sig1.cols()), temp2(Tx_sig2.rows(), Tx_sig2.cols());
	for (int i = 0; i < Tx_sig1.rows(); i++)
	{
		double n11 =  sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		double n12 =  sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		double n21 =  sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		double n22 =  sigma * sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
		complex<double> n0_tmp(n11, n12), n1_tmp(n21, n22);
		noise(i, 0) = n0_tmp;
		noise(i, 1) = n1_tmp;
	}

	temp1 = A(0, 0)*Tx_sig1 + A(0, 1)*Tx_sig2;
	temp2 = A(1, 0)*Tx_sig1 + A(1, 1)*Tx_sig2;
	rx << temp1, temp2;
	rx = rx + noise;
	//cout << noise<< endl;
	//cout << "----------" << endl;
	return rx;
}


Matrix<complex<double>, Dynamic, Dynamic> Channel::CreatH()
{
	Matrix<complex<double>, Dynamic, Dynamic> H(2, 2);
	double h11r, h12r, h21r, h22r,h11i,h12i,h21i,h22i;
	double tmp = 0;
	h11r = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h21r = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h12r = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h22r = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h11i = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h21i = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h12i = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	h22i = sqrt(-2.0 * log(RandFloat())) * cos(6.283185307 * RandFloat());
	tmp = sqrt(pow(h11r, 2) + pow(h11i, 2));
	if (tmp > 1)
	{
		h11r = h11r / tmp;
		h11i = h11i / tmp;
	}
	tmp = sqrt(pow(h12r, 2) + pow(h12i, 2));
	if (tmp > 1)
	{
		h12r = h12r / tmp;
		h12i = h12i / tmp;
	}
	tmp = sqrt(pow(h21r, 2) + pow(h21i, 2));
	if (tmp > 1)
	{
		h21r = h21r / tmp;
		h21i = h21i / tmp;
	}

	tmp = sqrt(pow(h22r, 2) + pow(h22i, 2));
	if (tmp > 1)
	{
		h22r = h22r / tmp;
		h22i = h22i / tmp;
	}



	H << complex<double>(h11r, h11i), complex<double>(h12r, h12i),
		complex<double>(h21r, h21i), complex<double>(h22r, h22i);

	return H;
}
