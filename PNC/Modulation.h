#pragma once
class Modulation
{
public:
	int modulationtype;
public:
	Modulation();
	~Modulation();
	VectorXi mod_BPSK(VectorXi msg_bits);
	Matrix<complex<double>, Dynamic, Dynamic> mod_QPSK(VectorXi msg_bits);
};

