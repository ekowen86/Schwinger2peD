#include "utils.h"
#include "io.h"

void writeGauge(field3D<Complex> *gauge, string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  int Nz = gauge->p.Nz;  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int z=0; z<Nz; z++)
	for(int mu=0; mu<2; mu++)
	  outPutFile << setprecision(12) <<  setw(20) << arg(gauge->read(x,y,z,mu)) << endl;
  
  outPutFile.close();
  return;  
}

void readGauge(field3D<Complex> *gauge, string name)
{
  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  int Nz = gauge->p.Nz;
  
  for(int x=0; x<Nx; x++) {
    for(int y=0; y<Ny; y++) {
      for(int z=0; z<Nz; z++) {
	for(int mu=0; mu<2; mu++) {
	  getline(inPutFile, val);
	  gauge->write(x, y, z, mu, polar(1.0, stod(val)));
	}
      }
    }
  }
  
  return;
}
