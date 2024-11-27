#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"

TimeScheme::TimeScheme(DataFile* df, Laplacian* lap) :
_df(df), _lap(lap)
{

}

void TimeScheme::SaveSol(const std::vector<double> &U, std::string n_sol, int n){
    std::string n_file = "../res/" + n_sol + "." + std::to_string(n) + ".vtk";
    std::ofstream monflux;
    int Np, Me, iBeg, iEnd, jBeg, jEnd, Nx(_df->Get_Nx());
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    int r(_df->Get_r());

    charge(Me, (_df->Get_Ny()-1), Np, &iBeg, &iEnd);

    if (Me == 0){
        jBeg = iBeg;
        jEnd = iEnd + r/2 + r%2;
    }
    else if (Me == Np-1){
        jBeg = iBeg - r/2;
        jEnd = iEnd;
    }
    else{
        jBeg = iBeg - r/2;
        jEnd = iEnd + r/2 + r%2;
    }

    monflux.open(n_file, std::ios::out);
    monflux << "# vtk DataFile Version 2.0" << std::endl;
    monflux << "U & mpi_rank" << std::endl;
    monflux << "ASCII" << std::endl;
    monflux << "DATASET STRUCTURED_GRID" << std::endl;
    monflux << "DIMENSIONS" << " " << Nx-1 << " " << jEnd-jBeg+1 << " " << 1 << std::endl;
    monflux << "POINTS" << " " << (jEnd-jBeg+1)*(Nx-1) << " " << "double" << std::endl;
    for(int j = jBeg+1; j <= jEnd+1; j++){
        for(int i = 1; i < _df->Get_Nx(); i++){
            monflux << _df->Get_xmin() + i*_df->Get_dx() << " " << _df->Get_ymin() + j*_df->Get_dy() << " " << 0.0 << std::endl; 
        }
    }
    monflux << "POINT_DATA" << " " << (jEnd-jBeg+1)*(Nx-1) << std::endl;
    monflux << "SCALARS U double" << std::endl;
    monflux << "LOOKUP_TABLE default" << std::endl;
    int k(0);
    for(int j = jBeg+1; j <= jEnd+1; j++){
        for(int i = 1; i < _df->Get_Nx(); i++){
            monflux << U[k] << std::endl; 
            ++k;
        }
    }
    monflux << "SCALARS mpi_rank int" << std::endl;
    monflux << "LOOKUP_TABLE default" << std::endl;
    for(int j = jBeg+1; j <= jEnd+1; j++){
        for(int i = 1; i < _df->Get_Nx(); i++){
            monflux << Me << std::endl;
        }
    }
    monflux.close();
}

ImplicitScheme::ImplicitScheme(DataFile* df, Laplacian* lap) : 
TimeScheme(df, lap)
{

}

std::vector<double> ImplicitScheme::Jacobi(const std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> x(U.size()), r(U.size());
    int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    double dx(_df->Get_dx()), dy(_df->Get_dy()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
    double dt(_df->Get_dt()), D(_df->Get_D());
    double a(1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy)));
    int Nmax(10000), it(0);
    x = U;
    r = SubVector(F,_lap->MatVecProd(x));
    while ((it < Nmax) && (std::sqrt(DotProduct(r,r))/std::sqrt(DotProduct(AddVector(U,MultiplyBy(F,dt)),AddVector(U,MultiplyBy(F,dt)))) > 1e-12)){
        x = AddVector(x,MultiplyBy(r,(1.0/a)));
        r = SubVector(AddVector(U,MultiplyBy(F,dt)),_lap->MatVecProd(x));
        ++it;
    }

    if (it >= Nmax){
        printf("Pas de convergence (Jacobi)\n");
        std::exit(0);
    }

    return x;
}

std::vector<double> ImplicitScheme::CG(const std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> x(U.size()), r(U.size()), z(U.size()), p(U.size()), q(U.size());
    double rho, rho0, alpha, delta, gamma;
    int Nmax(10000), k(0);
    double dt(_df->Get_dt());
    x = U;
    r = SubVector(F,_lap->MatVecProd(x));
    while ((k < Nmax) && (std::sqrt(DotProduct(r,r))/std::sqrt(DotProduct(AddVector(U,MultiplyBy(F,dt)),AddVector(U,MultiplyBy(F,dt)))) > 1e-12)){
        z = r;
        rho = DotProduct(r,z);
        if (rho == 0.0){
            printf("Pas de convergence (CG)\n");
            std::exit(0);
        }
        if (k == 0){
            p = z;
        }
        else{
            gamma = rho/rho0;
            p = AddVector(MultiplyBy(p,gamma),z);
        }
        q = _lap->MatVecProd(p);
        delta = DotProduct(p,q);
        if (delta == 0.0){
            printf("Pas de convergence (CG)\n");
            std::exit(0);
        }
        alpha = rho/delta;
        x = AddVector(x,MultiplyBy(p,alpha));
        r = SubVector(r,MultiplyBy(q,alpha));
        rho0 = rho;
        ++k;
    }

    if (k >= Nmax){
        printf("Pas de convergence (CG)\n");
        std::exit(0);
    }

    return x;
}

std::vector<double> ImplicitScheme::BiCGstab(std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> r(U.size()), r_old(U.size()), r_tilde(U.size()), p(U.size()), nu(U.size()), h(U.size()), s(U.size());
    std::vector<double> t1(U.size()), x(U.size()), p_old(U.size());
    double rho(0.0), rho_old(0.0), alpha(0.0), omega(0.0), beta(0.0);
    int Nmax(10000), k(0);
    double dt(_df->Get_dt());

    int Np, Me;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    x = U;
    r = SubVector(F,_lap->MatVecProd(x));
    r_tilde = r;
    rho = DotProduct(r_tilde,r);
    p = r;
    while (k <= Nmax){
        nu = _lap->MatVecProd(p);
        alpha = rho/DotProduct(r_tilde,nu);
        h = AddVector(x,MultiplyBy(p,alpha));
        s = SubVector(r,MultiplyBy(nu,alpha));
        if (sqrt(DotProduct(s,s)) <= 1.e-12){
            x = h;
            return x;
        }
        t1 = _lap->MatVecProd(s);
        omega = DotProduct(t1,s)/DotProduct(t1,t1);
        x = AddVector(h,MultiplyBy(s,omega));
        r_old = r;
        r = SubVector(s,MultiplyBy(t1,omega));
        if (sqrt(DotProduct(r,r)) <= 1.e-12) {
            return x;
        }
        rho_old = rho;
        rho = DotProduct(r_tilde,r);
        beta = (rho/rho_old)*(alpha/omega);
        p_old = p;
        p = AddVector(r,MultiplyBy(SubVector(p_old,MultiplyBy(nu,omega)),beta));
        k++;
    }

    if (k >= Nmax){
        printf("Pas de convergence (BiCGstab)\n");
        std::exit(0);
    }

    return x;
}

void ImplicitScheme::Integrate(double &t, std::vector<double> &U){
    t += _df->Get_dt();
    std::vector<double> U_old, U_diff, U_Me, U_Me_old;
    double max(1e12), max_loc(1e12);
    int Nx(_df->Get_Nx());
    int Np, Me, k(0), kmax(1000);
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    while ((k < kmax)&&(max > 1e-12)){
        U_old = U;
        U_Me_old = TruncVector(U_old,_df->Get_Nx());
        if (_df->Get_Solver() == "Jacobi"){
            U = Jacobi(U, _lap->RHS(t,U));
        } 
        if (_df->Get_Solver() == "BiCGstab"){
            U = BiCGstab(U, _lap->RHS(t,U));
        }
        else if (_df->Get_Solver() == "CG"){
            U = CG(U, _lap->RHS(t,U));
        }
        else{
            printf("Pas de solveur\n");
            std::exit(0);
        }
        U_Me = TruncVector(U,_df->Get_Nx());
        U_diff = AbsVector(SubVector(U_Me,U_Me_old));
        max_loc = *std::max_element(U_diff.begin(), U_diff.end());
        MPI_Allreduce(&max_loc, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        ++k;
    }

    if (k >= kmax){
        printf("Pas de convergence (Schwarz)\n");
        std::exit(0);
    }
}

#define _TIME_SCHEME_CPP
#endif