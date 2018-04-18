#pragma once
class Channel
{
public:
	Channel();
	~Channel();
	//MatrixXd AWGN(double sigma, VectorXi Tx_sig1, VectorXi Tx_sig2,MatrixXi A);
	Matrix<complex<double>, Dynamic, Dynamic> AWGN(double sigma, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig1, Matrix<complex<double>, Dynamic, Dynamic> Tx_sig2, Matrix<complex<double>, Dynamic, Dynamic> A);
	Matrix<complex<double>, Dynamic, Dynamic> CreatH();
};

