/* Copyright 2007, 2010 K.L. Lee, J.W. Shang, W.K. Chua, S.Y. Looi and B.-G. Englert
// This file is part of SOMIM.
SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <ctime>

using namespace std;

typedef complex<double> Complex;

class SOMIM{
 private:
  int povmDimension, directGradChance, maxDimen, maxMatrix;
  double MI, new_MI, alpha, phi;
  double tol_MI, tol_alpha, tol_S;
  vector< vector<double> > prob_table;
  int N, J, K;
  Complex zero;
  vector<Complex> rho2; /* 3rd dimension of rho */
  vector< vector<Complex> > rho1; /* 2nd dimension of rho */
  vector< vector< vector<Complex> > > rho; /* 1st dimension of rho. Dimension of rho is J x N x N */
  vector<Complex> povm2; /* 3rd dimension of povm */
  vector< vector<Complex> > povm1; /* 2nd dimension of povm */
  vector< vector< vector<Complex> > > povm; /* 1st dimension of povm. Dimension of povm is K x N x N */
						
  // various functions that deal with actual calculation
  double calc_beta(vector< vector<Complex> >&, vector< vector<Complex> >&);
  void calc_delta_A_dagger_CG(double, int, vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  void calc_grad_MI(vector< vector<Complex> >&, vector< vector< vector<Complex> > >&, vector< vector<Complex> >&);
  void calc_L_inverse(vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  inline double calc_MI(double some_alpha, vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  double calc_MI_from_prob_table(vector< vector<double> >&);
  void calc_new_A_dagger(double some_alpha, vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  void calc_prob_table(vector< vector<Complex> >&);
  void calc_R(vector< vector< vector<Complex> > >&);
  void find_alpha_GSS(vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  void gen_A_dagger(vector< vector<Complex> >&, vector< vector<Complex> >&);
		
  // functions initializing 2D and 3D matrices
  vector< vector<Complex> > Init2DMatrix(int X, int Y);
  vector< vector<double> > Init2DRealMatrix(int X, int Y);
  vector< vector< vector<Complex> > > Init3DMatrix(int X, int Y, int Z);
		
  // functions multiplying scalars, vectors and matrices
  Complex inner_product(vector<Complex>&, vector<Complex>&);
  double LOG2(double);
  void matrix_add_matrix(vector< vector<Complex> >&, vector< vector<Complex> >&);
  void matrix_times_matrix(vector< vector<Complex> >&, vector< vector<Complex> >&, vector< vector<Complex> >&);
  void matrix_times_vector(vector<Complex>&, vector< vector<Complex> >&, vector<Complex>&);
  void outer_product_add(vector< vector<Complex> >&, vector<Complex>&, vector<Complex>&);
  void scalar_times_vector(vector<Complex>&, double, vector<Complex>&);
  void swap2D(vector< vector<Complex> >&, vector< vector<Complex> >&);
  void transpose_conjugate(vector< vector<Complex> >&);		

 public:
  SOMIM(int n, int j, int k);		
  void setRho(int j, int n, int m, Complex z);
  void setPOVM(int j, int n, int m, Complex z);
  Complex getRho(int j, int n, int m);
  Complex getPOVM(int j, int n, int m);
  int getN();
  int getJ();
  int getK();
  void setDirectGrad(int g); 
  void setTolerance(double t); 
  double getMI();
};
