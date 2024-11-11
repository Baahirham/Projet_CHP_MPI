#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
private:
   // Pointeur de la classe DataFile pour récupérer toutes les
   // valeurs de paramètres
   const DataFile* _df;
   // Diffusion coefficient

   public: // Méthodes et opérateurs de la classe
   Function(DataFile* df);
   double Exact_solution(const double x, const double y, const double t) const;
   double Initial_condition(const double x, const double y) const;
   double Source(const double x, const double y, const double t) const;
   double Dirichlet_Gamma_0(const double x, const double y, const double t) const;
   double Dirichlet_Gamma_1(const double x, const double y, const double t) const;
};

double DotProduct(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> MultiplyBy(const std::vector<double> &a, const double &b);
std::vector<double> AddVector(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> SubVector(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> AbsVector(const std::vector<double> &a);
std::vector<double> TruncVector(const std::vector<double> &a, const int Nx);
void charge(int me, int n, int np, int *iBeg, int *iEnd);

#define _FUNCTION_H
#endif
