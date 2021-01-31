#pragma once

#include "schwinger2peD_internal.h"
#include "blas.h"
// #include "utils.h"


template<typename T> class field {

public:

  ~field(){
    data.resize(0);
  }

  std::vector<T> data;
  param_t p;

  field(std::vector<T> &data, param_t p) : data(data), p(p)
  {
  }

  field(param_t p) : p(p)
  {
    data.resize(p.Nx * p.Ny * 2);
    for(unsigned int i=0; i<data.size(); i++) data[i] = 0.0;
  }

  T read(int x, int y, int mu) const {
    return data[2*(x + p.Nx * y) + mu];
  }

  void write(int x, int y, int mu, const T elem) {
    data[2*(x + p.Nx * y) + mu] = elem;
  }

  void copy(const field<T> *in){
    blas::copy(data, in->data);
  }

  unsigned int size() { return data.size(); }

  void print() {
    for(int x=0; x<p.Nx; x++) {
      for(int y=0; y<p.Ny; y++) {
	for(int mu=0; mu<2; mu++) {
	  cout << "elem("<<x<<","<<y<<":" << mu << ") = " << data[2*(x + p.Nx * y) + mu] << endl;
	}
      }
    }
  }

  void print(int n) {
    cout << "elem " << n << " = " << data[n] << endl;
  }

    double rand_double(double min = 0.0, double max = 1.0) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(p.gen);
    }

    double rand_int(int min, int max) {
        std::uniform_int_distribution<double> dist(min, max);
        return dist(p.gen);
    }

    double rand_normal(double mean = 0.0, double stdev = 1.0) {
        std::normal_distribution<double> dist(mean, stdev);
        return dist(p.gen);
    }

    Complex rand_u1() {
        double arg = rand_double(-M_PI, M_PI);
        return polar(1.0, arg);
    }
};


template<typename T> class field3D {

public:

  ~field3D(){
    data.resize(0);
  }

  std::vector<T> data;
  param_t p;

  field3D(std::vector<T> &data, param_t p) : data(data), p(p)
  {
  }

  field3D(param_t p) : p(p)
  {
    data.resize(p.Nx * p.Ny * p.Nz * 2);
    for(unsigned int i=0; i<data.size(); i++) data[i] = 0.0;
  }

  T read(int x, int y, int z, int mu) const {
    return data[2 * (x + p.Nx * y + p.Nx * p.Ny * z) + mu];
  }

  void write(int x, int y, int z, int mu, const T elem) {
    data[2 * (x + p.Nx * y + p.Nx * p.Ny * z) + mu] = elem;
  }

  void copy(const field3D<T> *in){
    blas::copy(data, in->data);
  }

  unsigned int size() { return data.size(); }

  void print() {
    for(int x=0; x<p.Nx; x++) {
      for(int y=0; y<p.Ny; y++) {
	for(int z=0; z<p.Nz; z++) {
	  for(int mu=0; mu<2; mu++) {
	    cout << "elem("<<x<<","<<y<<","<<z<<":"<< mu << ") = " << data[2 * (x + p.Nx * y + p.Nx * p.Ny * z) + mu] << endl;
	  }
	}
      }
    }
  }

  void print(int n) {
    cout << "elem " << n << " = " << data[n] << endl;
  }

    double rand_double(double min = 0.0, double max = 1.0) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(p.gen);
    }

    double rand_int(int min, int max) {
        std::uniform_int_distribution<double> dist(min, max);
        return dist(p.gen);
    }

    double rand_normal(double mean = 0.0, double stdev = 1.0) {
        std::normal_distribution<double> dist(mean, stdev);
        return dist(p.gen);
    }

    Complex rand_u1() {
        double arg = rand_double(-M_PI, M_PI);
        return polar(1.0, arg);
    }
};
