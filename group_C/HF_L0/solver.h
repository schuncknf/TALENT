// solves the Hartree-Fock according to the matrix elements in Vacbd
// return value will contain the resulting energies and eigenvectors
eig_t solve_HF(Vab_t **Vacbd, int N_dim, int N_occ);
