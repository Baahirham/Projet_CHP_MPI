#ifndef _LAPLACIAN_CPP

#include "Laplacian.h"

using namespace std;

// Constructeur
Laplacian::Laplacian(Function* fct, DataFile* df) :
_fct(fct), _df(df)
{
    
}

void Laplacian::InitialCondition(std::vector<double> &U){
    int k(0);
    int Np, Me, iBeg, iEnd, jBeg, jEnd;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    int r(_df->Get_r());

    charge(Me, (_df->Get_Ny()-1), Np, &iBeg, &iEnd);

    if (Me == 0){
        U.resize((iEnd-iBeg+1 + r/2 + r%2)*(_df->Get_Nx()-1),0.0);
        jBeg = iBeg;
        jEnd = iEnd + r/2 + r%2;
    }
    else if (Me == Np-1){
        U.resize((iEnd-iBeg+1 + r/2)*(_df->Get_Nx()-1),0.0);
        jBeg = iBeg - r/2;
        jEnd = iEnd;
    }
    else{
        U.resize((iEnd-iBeg+1 + r)*(_df->Get_Nx()-1),0.0);
        jBeg = iBeg - r/2;
        jEnd = iEnd + r/2 + r%2;
    }


    for (int j = jBeg+1; j <= jEnd+1; ++j){
        for (int i = 1; i<_df->Get_Nx(); ++i){
            U[k] = _fct->Initial_condition(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy());
            ++k;
        }
    }
}

std::vector<double> Laplacian::MatVecProd(const std::vector<double> &U){

    int k(0), N(U.size());
    std::vector<double> X(N);
    int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    double dx(_df->Get_dx()), dy(_df->Get_dy()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
    double dt(_df->Get_dt()), D(_df->Get_D()), alpha(_df->Get_alpha()), beta(_df->Get_beta());
    double a, b, c, d((2.0*D*dt*beta)/(alpha*dy));
    int i, j,r(_df->Get_r());
    int Np, Me, iBeg, iEnd;

    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    charge(Me, (_df->Get_Ny()-1), Np, &iBeg, &iEnd);
    
    if ((Me != 0) && (Me != Np-1)){
        for(int k = 0; k < N; ++k){
            if ((Nx-2 < k) && (k < N-Nx+1)) {
        
                a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy));
                b = -D*dt/(dx*dx);
                c = -D*dt/(dy*dy);
                
            }
            else{
                if(alpha == 0.0)
                {
                    a = 1.0/(dx*dx) + 1.0/(dy*dy);
                    b = 0;
                    c = 0;
                }
                else
                {
                    a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy)) + d;
                    b = -D*dt/(dx*dx);
                    c = -2.0*D*dt/(dy*dy);
                }
                
            }
            i = k%(Nx-1);
            j = k/(Nx-1);
            X[k] = a * U[k];
            if (i > 0) X[k] += b * U[k - 1];
            if (i < Nx-2) X[k] += b * U[k + 1];
            if (j > 0) X[k] += c * U[k - Nx+1];
            if (j < (iEnd-iBeg + 1 + r - 1)) X[k] += c * U[k + Nx-1];
        }
    }
    if ((Me == 0)){
        for(int k = 0; k < N; ++k){
            if (k < N-Nx+1) {
                a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy));
                b = -D*dt/(dx*dx);
                c = -D*dt/(dy*dy);
            }
            else{
                if(alpha == 0.0)
                {
                    a = 1.0/(dx*dx) + 1.0/(dy*dy);
                    b = 0;
                    c = 0;
                }
                else
                {
                    a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy)) + d;
                    b = -D*dt/(dx*dx);
                    c = -2.0*D*dt/(dy*dy);
                }
            }
            i = k%(Nx-1);
            j = k/(Nx-1);
            X[k] = a * U[k];
            if (i > 0) X[k] += b * U[k - 1];
            if (i < Nx-2) X[k] += b * U[k + 1];
            if (j > 0) X[k] += c * U[k - Nx+1];
            if (j < (iEnd-iBeg + r/2 + r%2)) X[k] += c * U[k + Nx-1];
        }
    }
    if (Me == Np-1){
        for(int k = 0; k < N; ++k){
            if (Nx-2 < k) {
                a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy));
                b = -D*dt/(dx*dx);
                c = -D*dt/(dy*dy);
            }
            else{
                if(alpha == 0.0)
                {
                    a = 1.0/(dx*dx) + 1.0/(dy*dy);
                    b = 0;
                    c = 0;
                }
                else
                {
                    a = 1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy)) + d;
                    b = -D*dt/(dx*dx);
                    c = -2.0*D*dt/(dy*dy);
                }
            }
            i = k%(Nx-1);
            j = k/(Nx-1);
            X[k] = a * U[k];
            if (i > 0) X[k] += b * U[k - 1];
            if (i < Nx-2) X[k] += b * U[k + 1];
            if (j > 0) X[k] += c * U[k - Nx+1];
            if (j < (iEnd-iBeg + r/2)) X[k] += c * U[k + Nx-1];
        }
    }

    return X;
}

std::vector<double> Laplacian::RHS(const double t, const std::vector<double> &U){
    const int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    std::vector<double> rhs(U.size());
    int k(0);
    double dx(_df->Get_dx()), dy(_df->Get_dy()), alpha(_df->Get_alpha()), beta(_df->Get_beta()), dt(_df->Get_dt()), D(_df->Get_D());
    double x_i, xmin(_df->Get_xmin());
    double y_j, ymin(_df->Get_ymin());
    int Np, Me, r(_df->Get_r()), N(U.size()), iBeg, jBeg, iEnd, jEnd, r1, r2;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    int Tag11(11), Tag12(12), Mep1, Mem1,tmessage;

    if (alpha == 0)
    {
        tmessage = 1;
        r1 = r-1;
        r2 = r+1;
    }
    else
    {
        tmessage = 3;
        r1 = r;
        r2 = r;
    }
    std::vector<double> U_Mem1(tmessage*(Nx-1),0.0), U_Mep1(tmessage*(Nx-1),0.0);

    charge(Me, (_df->Get_Ny()-1), Np, &iBeg, &iEnd);

    if (Me == 0){
        jBeg = iBeg;
        jEnd = iEnd + r/2 + r%2;
        Mep1 = Me+1;
        Mem1 = MPI_PROC_NULL;
    }
    if (Me == Np-1)
    {
        jBeg = iBeg - r/2;
        jEnd = iEnd;
        Mep1 = MPI_PROC_NULL;
        Mem1 = Me-1;
    }
    if ((Me != 0) && (Me != Np-1)){
        jBeg = iBeg - r/2;
        jEnd = iEnd + r/2 + r%2;
        Mep1 = Me+1;
        Mem1 = Me-1;
    }

    MPI_Send(&U[N-(r1+1)*(Nx-1)], tmessage*(Nx-1), MPI_DOUBLE, Mep1, Tag12, MPI_COMM_WORLD);
    MPI_Recv(&U_Mem1[0], tmessage*(Nx-1), MPI_DOUBLE, Mem1, Tag12, MPI_COMM_WORLD, &status);
    MPI_Send(&U[(r2-2)*(Nx-1)], tmessage*(Nx-1), MPI_DOUBLE, Mem1, Tag11, MPI_COMM_WORLD);
    MPI_Recv(&U_Mep1[0], tmessage*(Nx-1), MPI_DOUBLE, Mep1, Tag11, MPI_COMM_WORLD, &status);

    for(int j = jBeg+1; j <= jEnd+1; ++j){
        for(int i = 1; i < Nx; ++i){
            x_i = xmin + i*dx;
            y_j = ymin + j*dy;
            rhs[k] = U[k] + dt*_fct->Source(x_i,y_j,t); 
            if(i == 1){
                rhs[k] += dt*D*_fct->Dirichlet_Gamma_1(x_i - dx,y_j,t)/(dx*dx);
            }        
            if(i == Nx-1){
                rhs[k] += dt*D*_fct->Dirichlet_Gamma_1(x_i + dx,y_j,t)/(dx*dx);
            } 
            if(j == jBeg+1){
                if(alpha == 0)
                {
                    if (Me == 0) rhs[k] += dt*D*_fct->Dirichlet_Gamma_0(x_i,y_j - dy,t)/(dy*dy);
                    else rhs[k] = U_Mem1[i-1]*(1.0/(dx*dx) + 1.0/(dy*dy));
                }
                else
                {
                    if (Me == 0) rhs[k] += dt*D*_fct->Dirichlet_Gamma_0(x_i,y_j - dy,t)/(dy*dy);
                    else rhs[k] += dt*(D/(dy*dy))*(U_Mem1[i-1] - U_Mem1[i-1 + 2*(Nx-1)] + (2*beta*dy/alpha)*U_Mem1[i-1 + Nx-1]);
                }
            } 
            if(j == jEnd+1){
                if(alpha == 0)
                {
                    if (Me == Np-1) rhs[k] += dt*D*_fct->Dirichlet_Gamma_0(x_i,y_j + dy,t)/(dy*dy);
                    else rhs[k] = U_Mep1[i-1]*(1.0/(dx*dx) + 1.0/(dy*dy));
                }
                else
                {
                    if (Me == Np-1) rhs[k] += dt*D*_fct->Dirichlet_Gamma_0(x_i,y_j + dy,t)/(dy*dy);
                    else rhs[k] += dt*(D/(dy*dy))*(U_Mep1[i-1 + 2*(Nx-1)] - U_Mep1[i-1] + (2*beta*dy/alpha)*U_Mep1[i-1 + Nx-1]);
                }
            } 
            ++k;
        }
    }
    return rhs;
}

std::vector<double> Laplacian::ExactSol(const double t){
    std::vector<double> ExactSol(_df->Get_Nx()*_df->Get_Ny());
    int k(0), Np, Me, iBeg, iEnd, jBeg, jEnd;
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

    for (int j = jBeg+1; j <= jEnd+1; ++j){
        for (int i = 1; i<_df->Get_Nx(); ++i){
            ExactSol[k] = _fct->Exact_solution(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy(),t);
            ++k;
        }
    }
    return ExactSol;
}

#define _LAPLACIAN_CPP
#endif
