//
//  Note that division by zero is avoided because the division is protected
//  by the "if" clause which surrounds it.
//
//  g++ -o auto_edit_prop auto_edit_prop.cpp
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
  int i34, nds[4];
  double * ele_prop;
  double * default_prop;

  if ( argc !=5  ) {
    cout << "wrong arguments" << endl;
    cout << "usage: ptNpoly input.reg hgrid.gr3 in_value out_value" <<endl;
    cout << "options:" << endl;
    return -1;
  } else {
    in_value = atof(argv[3]);
    out_value = atof(argv[4]);
    if (int(out_value) == -9999) {
      i_default = true;
    }else{
      i_default = false;
    }
  }

  cout << "i_default = " << i_default << endl;
  cout << "in_value = " << in_value << endl;
  cout << "out_value = " << out_value << endl;


  //read polygon from *.reg
  ifstream fin (argv[1]);
  if (!fin.is_open()) {
    cout << "cannot open " << argv[1] << " for reading" << endl;
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
  
 
  //read hgrid
  ifstream f_grid (argv[2]);
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
  ele_prop = new double[ne];
  default_prop = new double[ne];
  for (int i=0; i<np; i++){
    f_grid >> dummy >> gridX[i] >> gridY[i] >> gridZ[i];
    if (i==0 || i==np-1) {
      cout << "node " << i+1 << "; dummy = " << dummy << endl;
      cout << gridX[i] << " " << gridY[i] << " " << gridZ[i] << endl;
    }
  }

  precalc_values();


  //read polygon from default.prop
  if (i_default) {
    ifstream f_default_prop ("default.prop");
    if (!f_default_prop.is_open()) {
      cout << "cannot open default prop for reading" << endl;
      return -1;
    }
    for (int i=0; i<ne; i++) {
      f_default_prop >> dummy >> default_prop[i];
    }
  }


  //continue reading hgrid for elements
  getline (f_grid,line);
  for (int i=0; i<ne; i++) {
    x = 0; y = 0;
    f_grid >> dummy >> i34;
    for (int j=0; j<i34; j++) {
      f_grid >> nds[j];
      x += gridX[nds[j]-1]; 
      y += gridY[nds[j]-1]; 
    }
    x/=i34; y/=i34;

    if (i==0 || i==ne-1) {
      cout << "ele " << i+1 << "; dummy = " << dummy << endl;
      cout << dummy <<" " << i34 << " " << nds[0] << " " << nds[1] << " " << nds[2] << " " << nds[3] << " " << endl;
      cout << x << " " << y << " " << endl;
    }

    if (pointInPolygon(x,y)) {
      ele_prop[i] = in_value;
    }else{
      ele_prop[i] = i_default ? int(default_prop[i]) : int(out_value); 
    }
  }
  

  ofstream fout ("out.prop");
  if (!fout.is_open()) {
    cout << "cannot open out.prop for writing" << endl;
    return -1;
  }
  fout.precision(16);
  for (int i=0; i<ne; i++) {
    fout << i+1 << " " << ele_prop[i] << endl;
  }


  delete [] constant; delete [] multiple;
  delete [] polyX; delete [] polyY;
  delete [] gridX; delete [] gridY;
  delete [] ele_prop; delete [] default_prop;

  return 0;
}

