#include "utils.h"
#include "io.h"

using namespace std;

void writeGaugeText(field3D<Complex> *gauge, string name){

  fstream outPutFile;
  outPutFile.open(name, ios::out | ios::trunc);
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

void readGaugeText(field3D<Complex> *gauge, string name)
{
  fstream inPutFile;
  inPutFile.open(name, ios::in);
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

void writeGaugeBinary(field3D<Complex>& gauge, string name){

    ofstream file(name, ios::out | ios::trunc | ios::binary);

    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;
    int Nz = gauge.p.Nz;
    Complex value;
    double valueR;
    double valueI;
    char* pR = reinterpret_cast<char*>(&valueR);
    char* pI = reinterpret_cast<char*>(&valueI);
    for (int i = 0; i < gauge.data.size(); i++) {
        value = gauge.data[i];
        valueR = real(value);
        valueI = imag(value);
        file.write(pR, sizeof(double));
        file.write(pI, sizeof(double));
    }

    file.close();
}

void readGaugeBinary(field3D<Complex>& gauge, string name){

    ifstream file(name, ios::in | ios::binary);

    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;
    int Nz = gauge.p.Nz;
    double valueR;
    double valueI;
    char* pR = reinterpret_cast<char*>(&valueR);
    char* pI = reinterpret_cast<char*>(&valueI);
    for (int i = 0; i < gauge.data.size(); i++) {
        file.read(pR, sizeof(double));
        file.read(pI, sizeof(double));
        gauge.data[i] = Complex(valueR,valueI);
    }

    file.close();
}

void writeRngState(field3D<Complex>& gauge, string name){

    ofstream file(name, ios::out | ios::trunc);
    file.setf(ios::dec | ios::left);
    file.fill(' ');
    file << *gauge.p.gen;
    file.close();
}

void readRngState(field3D<Complex>& gauge, string name){

    ifstream file(name, ios::in);
    file.setf(ios::dec);
    file >> *gauge.p.gen;
    file.close();
}
