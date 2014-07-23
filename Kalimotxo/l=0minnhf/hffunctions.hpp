void spinit(imat& spbasis,int nmax);
void tpbasisinit(imat& spbasis,imat& tpbasis,int nmax);
void onebodymatinit(mat& onebodymat,imat& spbasis,double data[]);
void twobodymatinit(mat& twobodymat,imat& spbasis, imat& tpbasis,double data[]);
double enotexpect(mat& densities,mat& onebodymat,mat& twobodymat,imat& spbasis,imat& tpbasis);
void fockinit(mat& fock,mat& onebodymat,mat& twobodymat,mat& densities,imat& spbasis,imat& tpbasis);
void densityupdate(mat& fock,mat& densities,imat& spbasis,int nmaxoccupied);
