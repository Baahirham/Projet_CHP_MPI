#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* df) :
_df(df)
{
   
}

double Function::Initial_condition(const double x, const double y) const
{
   if (this->_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return 0.0;
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Exact_solution(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      double pi(std::acos(-1.0));
      double xmax(_df->Get_xmax()), ymax(_df->Get_ymax()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
      return sin(t)*sin(2.0*pi*x/(_df->Get_xmax()-_df->Get_xmin()))*sin(2.0*pi*y/(_df->Get_ymax()-_df->Get_ymin()));
   }
   else if (this->_df->Get_cas() == 1)
   {
      return x*(1-x)*y*(1-y);
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      std::cout << "Pas de solution exacte" << std::endl;
      std::exit(0);
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Source(const double x, const double y, const double t) const
{
   double pi(std::acos(-1.0));
   if (_df->Get_cas() == 0)
   {
      double xmax(_df->Get_xmax()), ymax(_df->Get_ymax()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
      return sin(2.0*pi*x/(_df->Get_xmax()-_df->Get_xmin()))*sin(2.0*pi*y/(_df->Get_ymax()-_df->Get_ymin()))*(cos(t) + sin(t)*_df->Get_D()*4.0*pi*pi*(1.0/((xmax - xmin)*(xmax - xmin)) + 1.0/((ymax - ymin)*(ymax - ymin))));
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 2.0*(x - x*x + y - y*y);
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      double Lx(_df->Get_xmax()-_df->Get_xmin()), Ly(_df->Get_ymax()-_df->Get_ymin());
      return exp(-(x-Lx/2.0)*(x-Lx/2.0))*exp(-(y-Ly/2.0)*(y-Ly/2.0))*cos((pi/2.0)*t);
   }
   else
   {
      return 0.0;
   }
}

double Function::Dirichlet_Gamma_0(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Dirichlet_Gamma_1(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}

double DotProduct(const std::vector<double> &a, const std::vector<double> &b){
    double result(0.0);
    for(int i = 0; i < a.size(); i++){
        result += a[i]*b[i];
    }
    return result;
}

std::vector<double> MultiplyBy(const std::vector<double> &a, const double &b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = b*a[i];
    }
    return result;
}

std::vector<double> AddVector(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> SubVector(const std::vector<double> &a, const std::vector<double> &b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> AbsVector(const std::vector<double> &a){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = std::fabs(a[i]);
    }
    return result;
}

std::vector<double> TruncVector(const std::vector<double> &a, const int Nx){
   int Np, Me;
   MPI_Comm_rank(MPI_COMM_WORLD, &Me);
   MPI_Comm_size(MPI_COMM_WORLD, &Np);

   std::vector<double> a_trunc(a);
   
   if(Me == 0){
      a_trunc.erase(a_trunc.begin(), a_trunc.end()-(Nx-1));
   }
   else if (Me == Np-1){
      a_trunc.resize(Nx-1);
   }
   else {
      a_trunc.erase(a_trunc.begin()+(Nx-1), a_trunc.end()-(Nx-1));
   }

   return a_trunc;

}

void charge(int me, int n, int np, int *iBeg, int *iEnd){
    // Fonction charge
    if (me < n%np){
        *iBeg = (n/np+1)*me;
        *iEnd = (n/np+1)*(me+1) - 1;
    }else{
        *iBeg = n%np*(n/np+1) + (me - n%np)*(n/np);
        *iEnd = *iBeg + n/np - 1;
    }
}


#define _FUNCTION_CPP
#endif
