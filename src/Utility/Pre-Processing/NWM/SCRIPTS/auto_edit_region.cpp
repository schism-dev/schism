//
//  Note that division by zero is avoided because the division is protected
//  by the "if" clause which surrounds it.
//
//  g++ -o auto_edit_region auto_edit_region.cpp
//
#include <iostream>
#include <fstream>
#include <new>
#include <string>
#include <stdlib.h>   
#include <math.h>  

int polyCorners;
double * constant;
double * multiple;
double * polyX;
double * polyY;
double * gridX;
double * gridY;
double * gridZ;



void precalc_values() {

int   i, j=polyCorners-1 ;

for(i=0; i<polyCorners; i++) {
  if(polyY[j]==polyY[i]) {
    constant[i]=polyX[i];
    multiple[i]=0; }
  else {
    constant[i]=polyX[i]-(polyY[i]*polyX[j])/(polyY[j]-polyY[i])+(polyY[i]*polyX[i])/(polyY[j]-polyY[i]);
    multiple[i]=(polyX[j]-polyX[i])/(polyY[j]-polyY[i]); }
    j=i; }}



bool pointInPolygon(double x, double y) {

int i, j=polyCorners-1;
bool oddNodes=false;

for (i=0; i<polyCorners; i++) {
  if ((polyY[i]< y && polyY[j]>=y
  ||   polyY[j]< y && polyY[i]>=y)) {
    oddNodes^=(y*multiple[i]+constant[i]<x); }
    j=i; }

return oddNodes; }

bool is_digits(const std::string &str) {
  return str.find_first_not_of("0123456789") == std::string::npos;
}

int main(int argc, char *argv[]) {
  using namespace std; 
  bool in_poly=false;
  string line;
  double x=0,y=0,dummy;
  double in_value=0.0, out_value=0.0;
  int np,ne,iMaxMin;
  bool i_default = true;

  if ( argc < 5 || argc > 6 || !is_digits(argv[1]) ) {
    cout << "wrong number of arguments" << endl;
    cout << "usage: ptNpoly 0 input.reg hgrid.gr3 in_value [out_value]" <<endl;
    cout << "the default outside-region value is the original z in the input grid" << endl;
    cout << "options:" << endl;
    cout << "first argument (0): inside region: z=value" << endl;
    cout << "first argument (1): inside region: z=max(z,value)" << endl;
    cout << "first argument (2): inside region: z=min(z,value)" << endl;
    return -1;
  } else {
    iMaxMin = atoi(argv[1]); 
    if (iMaxMin < 0 || iMaxMin > 2) {
      cout << "wrong option, must choose from (0, 1, or 2)" << endl;
      return -1;
    } else {
      in_value = atof(argv[4]);
      if (argc==6) {
        i_default = false;
        out_value = atof(argv[5]);
      }
    }
  }

  cout << "option = " << iMaxMin << endl;
  cout << "in_value = " << in_value << endl;
  cout << "out_value = " << out_value << endl;


  //read polygon from *.reg
  ifstream fin (argv[2]);
  if (!fin.is_open()) {
    cout << "cannot open " << argv[2] << " for reading" << endl;
    return -1;
  }
  getline (fin,line);
  getline (fin,line);
  fin >> polyCorners >> dummy;
  cout<< "polyCorners " << polyCorners << endl;

  constant = new double[polyCorners];
  multiple = new double[polyCorners];
  polyX = new double[polyCorners];
  polyY = new double[polyCorners];

  cout.precision(16);
  for (int i=0; i<polyCorners; i++){
    fin >> polyX[i] >> polyY[i];
    cout << polyX[i] << " " << polyY[i] << endl;
  }
  
 
  //read polygon from *.reg
  ifstream f_grid (argv[3]);
  if (!f_grid.is_open()) {
    cout << "cannot open " << argv[2] << " for reading" << endl;
    return -1;
  }
  getline (f_grid,line);
  f_grid >> ne >> np;
  cout<< "hgrid, ne/np:  " << ne <<"/" << np << endl;
  gridX = new double[np];
  gridY = new double[np];
  gridZ = new double[np];
  for (int i=0; i<np; i++){
    f_grid >> dummy >> gridX[i] >> gridY[i] >> gridZ[i];
    if (i==0 || i==np-1) {
      cout << "node " << i+1 << "; dummy = " << dummy << endl;
      cout << gridX[i] << " " << gridY[i] << " " << gridZ[i] << endl;
    }
  }

  precalc_values();

  //options:
  //0
  for (int i=0; i<np; i++) {
    if (pointInPolygon(gridX[i],gridY[i])) {
      switch (iMaxMin) {
        case 0:
          gridZ[i] = in_value; break;
        case 1:
          gridZ[i] = fmax(gridZ[i],in_value); break;
        case 2:
          gridZ[i] = fmin(gridZ[i],in_value); break;
        default:
          cout << "wrong option" << endl;
      }
    }else if (!i_default){
      gridZ[i] = out_value;
    }
  }
  
  ofstream fout ("out.gr3");
  if (!fout.is_open()) {
    cout << "cannot open " << argv[2] << " for writing" << endl;
    return -1;
  }
  fout.precision(16);
  fout << argv[2] << ", in: " << in_value << ", out: " << out_value << endl;
  fout << ne << " " << np << endl;
  for (int i=0; i<np; i++) {
    fout << i+1 << " " <<  gridX[i] << " " << gridY[i] << " " << gridZ[i] << endl;
  }

  //elements
  getline (f_grid,line);
  for (int i=0; i<ne; i++) {
    getline (f_grid,line);
    fout << line << endl;
  }

  delete [] constant; delete [] multiple;
  delete [] polyX; delete [] polyY;
  delete [] gridX; delete [] gridY;

  return 0;
}

