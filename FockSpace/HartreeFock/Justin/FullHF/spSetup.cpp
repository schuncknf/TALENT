#include <iostream>
#include <fstream>

using namespace std;

ofstream spFile;

int main()
{
  spFile.open("spFile.dat");
  
  // Nshell = 2*nshell + l
  int maxN = 4; // l = 0 for now.
  int particleNum = 2;
  int l = 0;
  double h_barOmega = 10; // MeV
  double spEnergy;

  spFile << "# number of particles" << endl;
  spFile << particleNum << endl;
  spFile << "# number of single-particle states" << endl;
  spFile << 2*(maxN+1) << endl;
  spFile << "# legend for single-particles states j=spin" << endl;
  spFile << "# nr n  l  2j 2j_z  energy" << endl;
  for(int n = 0; n <= maxN; n++)
    {
      spEnergy = (2*n + 1.5)*h_barOmega;
      spFile << "  " << 2*n << "  " << n << "  " << l << "  " << 1 << "  " << -1 << "    " << spEnergy << endl;
      spFile << "  " << 2*n+1 << "  " << n << "  " << l << "  " << 1 << "  " << 1 << "     " << spEnergy << endl;
    }


  cout << "meow mix please deliver!" << endl;

  spFile.close();

  return 0;
}
