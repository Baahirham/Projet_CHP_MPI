# Résolution d'une équation de diffusion à par décompostion de Schwarz additive

## Prérequis
- [MPI](https://www.open-mpi.org/) installé sur votre système
- [Paraview](https://www.paraview.org/) pour la visualisation des résultats
- Compilateur C/C++ compatible (gcc, clang, etc.)
- Make pour la compilation

## Lancer une simulation

1. **Préparation du fichier de configuration**  
   - Créez un fichier `data.toml` dans le dossier `data/` avec les paramètres suivants :  
     ```toml
     cas = 1   # Choisir 1, 2 ou 3
     Nx = 100  # Nombre de points en x
     Ny = 100  # Nombre de points en y
     r = 2     # Recouvrement (minimum 2)
     ```

2. **Compilation**  
   - Ouvrez un terminal dans le dossier `src/` et exécutez :  
     ```sh
     make clean && make
     ```

3. **Exécution**  
   - Lancez la simulation avec MPI (ici avec 4 processus) :  
     ```sh
     mpirun -n 4 ./run ../data/data.toml
     ```

4. **Récupération des résultats**  
   - Les fichiers de sortie (`*.vtk`) seront générés dans le dossier `res/`.  
   - Utilisez Paraview pour les visualiser :  
     ```sh
     paraview res/*.vtk
     ```

## Structure des dossiers

.
├── data/       # Contient data.toml
├── res/        # Fichiers résultats VTK
└── src/        # Code source


## Notes
- Ajustez le nombre de processus (`-n 4`) en fonction de votre système.
- Assurez-vous que le recouvrement `r` ≥ 2 dans `data.toml`.

