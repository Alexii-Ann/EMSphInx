#include "sht/square_sht.hpp"
#include <iostream>
#include <vector>


int main(int argc, char *argv[]){
  const size_t bw = atoi(argv[1]);

  std::ofstream csv("lambert.csv");
  csv << "x,y,z\n";
  std::vector<double> lam = emsphinx::square::normals<double>(bw, emsphinx::square::Layout::Lambert);
  for(size_t j = 0; j < bw*bw; j++) {
      csv << lam[3*j+0] << ',' << lam[3*j+1] << ',' << lam[3*j+2] << '\n';
  }


  std::ofstream csv2("legendre.csv");
  csv2 << "x,y,z\n";
  std::vector<double> leg = emsphinx::square::normals<double>(bw, emsphinx::square::Layout::Legendre);
  for(size_t j = 0; j < bw*bw; j++) {
      csv2 << leg[3*j+0] << ',' << leg[3*j+1] << ',' << leg[3*j+2] << '\n';
  }
  
