#pragma once
class Quantizer
{
public:
	Quantizer();
	~Quantizer();
	Matrix<complex<double>, Dynamic, Dynamic> quantize(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig,int bits,double maxvalue, double minvalue);
};

