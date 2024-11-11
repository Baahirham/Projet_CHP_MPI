#ifndef _TIME_SCHEME_H

#include "Laplacian.h"

class TimeScheme
{
    protected:
    DataFile* _df;
    Laplacian* _lap;

    public:
    TimeScheme(DataFile* df, Laplacian* lap);
    virtual void Integrate(double &t, std::vector<double> &U) = 0;
    void SaveSol(const std::vector<double> &U, std::string n_sol, int n);
};

class ImplicitScheme : public TimeScheme
{
    public:
    ImplicitScheme(DataFile* dd, Laplacian* lap);
    void Integrate(double &t, std::vector<double> &U);
    std::vector<double> Jacobi(const std::vector<double> &U, const std::vector<double> &F);
    std::vector<double> CG(const std::vector<double> &U, const std::vector<double> &F);
    std::vector<double> BiCGstab(std::vector<double> &U, const std::vector<double> &F);
};

#define _TIME_SCHEME_H
#endif