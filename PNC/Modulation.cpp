#include "stdafx.h"
#include "Modulation.h"


Modulation::Modulation()
{
}




Modulation::~Modulation()
{
}


VectorXi Modulation::mod_BPSK(VectorXi msg_bits)
{
	int len = msg_bits.size();
	VectorXi res(len);
	for (int i = 0; i < len; i++)
	{
		if (msg_bits(i) == 0)
			res(i) = -1;
		else if (msg_bits(i) == 1)
			res(i) = 1;
	}
	return res;
}


Matrix<complex<double>, Dynamic, Dynamic> Modulation::mod_QPSK(VectorXi msg_bits)
{
	int len = msg_bits.size()/2;
	Matrix<complex<double>, Dynamic, Dynamic> res(len,1);
	for (int i = 0; i < len; i++)
	{
		int tmp1 = msg_bits(2*i) == 0 ? -1.0: 1.0;
		int tmp2 = msg_bits(2 * i + 1) == 0 ? -1.0 : 1.0;
		complex<int> c(tmp1, tmp2);
		res(i, 0) = c;

	}
	return res;
}
