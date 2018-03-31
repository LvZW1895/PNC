#pragma once
class Demapping
{
public:
	Demapping();
	~Demapping();
	MatrixXi dp_pnc_bpsk(MatrixXd Rx_sig);
	MatrixXi dp_comp_bpsk(MatrixXd Rx_sig, MatrixXi A);
	MatrixXd dp_comp_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig, MatrixXi A);
	MatrixXd dp_pnc_qpsk(Matrix<complex<double>, Dynamic, Dynamic> Rx_sig);
};

