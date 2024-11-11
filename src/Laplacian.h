#ifndef _LAPLACIAN_H

#include "Function.h"

class Laplacian
{
private:
	DataFile* _df;
	Function* _fct;
public:
	// Constructeur
	Laplacian(Function* fct, DataFile* df);
	std::vector<double> MatVecProd(const std::vector<double> &U);
	std::vector<double> RHS(const double t, const std::vector<double> &U);
	void InitialCondition(std::vector<double> &U);
	std::vector<double> ExactSol(const double t);
};

#define _LAPLACIAN_H
#endif
