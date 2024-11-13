#include "TimeScheme.h"

int main(int argc, char** argv)
{

   if (argc < 2)
   {
      printf("Please, enter the name of your data file.\n");
      exit(0);
   }

      MPI_Init(&argc, &argv);

      const std::string data_file_name = argv[1];
      DataFile* df = new DataFile(data_file_name);
      Function* fct = new Function(df);
      Laplacian* lap = new Laplacian(fct, df);
      TimeScheme* ts = new ImplicitScheme(df, lap);
      std::vector<double> U;
      double t1, t2;
      int Np, Me;
      
      MPI_Comm_rank(MPI_COMM_WORLD, &Me);
      MPI_Comm_size(MPI_COMM_WORLD, &Np);

      t1 = MPI_Wtime();

      lap->InitialCondition(U);
      int it(0);
      ts->SaveSol(U,"sol_" + std::to_string(Me),it);
      double t(0.0);

      while (t < df->Get_Tf()){
         ts->Integrate(t, U);
         ++it;
         ts->SaveSol(U,"sol_" + std::to_string(Me),it);
      }

      t2 = MPI_Wtime();

      printf("Me = %d, temps = %lf\n", Me, t2-t1);

      // ------------------- Pour valider l'ordre du schéma ----------------

      if (df->Get_cas() != 3){

         double Erreur(0.0), Normalise(0.0), ErreurNorm(0.0);
         std::vector<double> ExacteSol(U.size());

         ExacteSol = lap->ExactSol(t);

         for(int k = 0; k < (U.size()); ++k){
         Erreur += (U[k] - ExacteSol[k])*(U[k] - ExacteSol[k]);
         Normalise += ExacteSol[k]*ExacteSol[k];
         }

         ErreurNorm = sqrt(Erreur/Normalise);

         printf("Me = %d, temps = %lf, ln(dx) = %lf, ln(Erreur) = %lf\n", Me, t2-t1, log(df->Get_dx()), log(ErreurNorm));

      }
      else{
         printf("Pas de solution exacte\n");
      }

      // ----------------------------------------------------------------------

      delete df, delete fct, delete lap, delete ts;

      MPI_Finalize();

   return 0;
}
