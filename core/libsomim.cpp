/* Copyright 2007, 2010 K.L. Lee, J.W. Shang, W.K. Chua, S.Y. Looi and B.-G. Englert
// This file is part of SOMIM.
SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

#include "libsomim.h"

typedef complex<double> Complex;

SOMIM::SOMIM(int n, int j, int k)
:tol_MI(1.0e-15),tol_alpha(1E-13),tol_S(1E-15),N(n),J(j),K(k),zero(0.0,0.0),povmDimension(n),directGradChance(0)
{
  for(int s=0;s<N;s++) rho2.push_back(zero);
  for(int s=0;s<N;s++) rho1.push_back(rho2);
  for(int s=0;s<J;s++) rho.push_back(rho1);
	
  for(int s=0;s<N;s++) povm2.push_back(zero);
  for(int s=0;s<N;s++) povm1.push_back(povm2);
  for(int s=0;s<K;s++) povm.push_back(povm1);

  phi = 0.5+(0.5*sqrt(5.0));
}


void SOMIM::setRho(int j, int n, int m, Complex z) {
  rho[j][n][m] = z;
}

void SOMIM::setPOVM(int j, int n, int m, Complex z) {
  povm[j][n][m] = z;
}

Complex SOMIM::getRho(int j, int n, int m) {
  return rho[j][n][m];
}

Complex SOMIM::getPOVM(int j, int n, int m) {
  return povm[j][n][m];
}

int SOMIM::getN() {
  return N;
}

int SOMIM::getJ() {
  return J;
}

int SOMIM::getK() {
  return K;
}

void SOMIM::setDirectGrad(int g) {
  directGradChance = g;
}

void SOMIM::setTolerance(double t) {
  tol_MI = t;
}

double SOMIM::getMI() {
  int i, j, k, l, cmpPOVM=0;
  Complex ratio1, ratio2;
  Complex temp(0,0);
  int restarted=0;
  int counter=0;
  double beta=0.0;
  int beta_count=0;
  int seed=(unsigned int)time(NULL);
  new_MI=0.0; MI=0.0; alpha=1.0;
  srand(seed);
  int test_import=0;
    
  // to determine the machine_epsilon
  double float_radix = 2.0;
  double inverse_radix = 1.0/float_radix;
  double machine_epsilon = 1.0;
  double temp_value = 1.0 + machine_epsilon;

  while(temp_value != 1.0) {
    machine_epsilon *= inverse_radix;
    temp_value = 1.0 + machine_epsilon;
  }
		
  // clear the probabilities stored previously
  prob_table.clear();
		
  vector< vector<Complex> > A_dagger = Init2DMatrix(K,N);
  vector< vector<Complex> > new_A_dagger = Init2DMatrix(K,N);
  vector< vector<Complex> > S = Init2DMatrix(N,N);
  vector< vector<Complex> > grad_MI = Init2DMatrix(K,N);
  vector< vector<Complex> >	delta_A_dagger = Init2DMatrix(K,N);
  vector< vector<Complex> > old_delta_A_dagger = Init2DMatrix(K,N);
  vector< vector<Complex> > old_grad_MI = Init2DMatrix(K,N);
  old_grad_MI[0][0] = 1.0;
  vector< vector< vector<Complex> > > R = Init3DMatrix(K,N,N);
	
  // // check whether to randomise intial POVMs or not
  // for(k=0; k<K; k++){
  //   for(i=0; i<N; i++){
  //     for(j=0; j<N; j++){
  // 	if(povm[k][i][j]!=temp){test_import = 1; break;}
  //     }
  //   }
  // }

  // if(test_import==1){ // import the manually input initial POVMs to S
  //   for(k=0; k<K; k++){
  //     matrix_add_matrix(S, povm[k]);} // S is not zero matrix anymore
  // } // otherwise do nothing, to randomly generate the initial POVMs
    
  gen_A_dagger(A_dagger, S); // initialise A_dagger

 restart:
  if (restarted==1)
    {
      counter = 0;
      seed = seed+1;
      srand(seed);
      gen_A_dagger(A_dagger, S);
    }
  restarted = 1;
  calc_prob_table(A_dagger); // for the specified rho, calculate the probability table using POVM A_dagger
  MI = calc_MI_from_prob_table(prob_table); // calculate the mutual information from the probability table
  if ((isnan(MI) || MI<0.0)) goto restart;

  while(1)
    {
      counter++;
      if (counter>10000) goto restart;
      calc_R(R);
      calc_grad_MI(grad_MI, R, A_dagger);
			
      // evaluate the chance to use direct gradient in the search for steepest ascent
      if((rand()%2)*100<directGradChance){
	for (i=0; i<K; i++){
	  for (j=0; j<N; j++) {
	    delta_A_dagger[i][j] = grad_MI[i][j];
	    old_delta_A_dagger[i][j] = grad_MI[i][j];
	  }
	}
      }else{
	calc_delta_A_dagger_CG(beta, beta_count, grad_MI, old_grad_MI, delta_A_dagger, old_delta_A_dagger);
      }
      // to find the improved rho and povm
      find_alpha_GSS(new_A_dagger, A_dagger, delta_A_dagger);
      new_MI = calc_MI(alpha, new_A_dagger, A_dagger, delta_A_dagger);

      if ((isnan(new_MI) || new_MI<0.0)) goto restart;
      if (2.0*(new_MI - MI) <= tol_MI*(new_MI+MI)+machine_epsilon && counter > 10) break; // termination criteria
      swap2D(new_A_dagger, A_dagger);
      MI = new_MI;
    }

  // to form the POVM before elimination
  for(k=0; k<K; k++){
    outer_product_add(povm[k], A_dagger[k], A_dagger[k]);
  }

  //  eliminate the extra povms
  calc_prob_table(A_dagger);
  double test, temp1, temp2, temp3;
  for(k=K-1; k>-1; k--){
    for(i=0; i<k; i++){
      temp1=0;
      temp2=0;
      temp3=0;
      for(j=0; j<J; j++){
	temp1 += prob_table[j][k]*prob_table[j][i];
	temp2 += prob_table[j][k]*prob_table[j][k];
	temp3 += prob_table[j][i]*prob_table[j][i];
      }

      test = 1.0 - temp1*temp1/(temp2*temp3);
      if(test<0.001){
	matrix_add_matrix(povm[i], povm[k]);
	swap2D(povm[k], povm[K-1]);
	K--;
	break;
      }
    }
  }

  return MI;
}

//------------------------------------------------------------------------------
//	functions dealing with actual calculation
//------------------------------------------------------------------------------
// The calculated value of beta equals to gamma_i defined in Eq. (10.6.7) in p. 422 of NRC, or beta_k in Eq.(2.7.32) in p. 84 of NRC
double SOMIM::calc_beta(vector< vector<Complex> >& grad_MI, vector< vector<Complex> >& old_grad_MI)
{
  double beta;
  int k;
  Complex tempnum = 0.0;
  Complex tempden = 0.0;
  for(k=0; k<K; k++){
    tempnum += inner_product(grad_MI[k], grad_MI[k])-inner_product(grad_MI[k], old_grad_MI[k]);
    tempden += inner_product(old_grad_MI[k], old_grad_MI[k]);
  }
  beta = real(tempnum)/real(tempden);
  return beta;
}

// This subroutine computes the conjugate gradient in a single iteration following the method in p. 83-86, or p. 420-423 of NRC
void SOMIM::calc_delta_A_dagger_CG(double beta, int beta_count, vector< vector<Complex> >& grad_MI, vector< vector<Complex> >& old_grad_MI, vector< vector<Complex> >& delta_A_dagger, vector< vector<Complex> >& old_delta_A_dagger)
{
  /* Conjugate Gradient (CG) algorithm */
  int k, n;
  beta = calc_beta(grad_MI, old_grad_MI);
  Complex tempsum = 0.0;
  if((beta_count < K*N*2) && (beta > 0)){
    for(k=0; k<K; k++){
      for(n=0; n<N; n++){
	// this line correspond to Eq. (10.6.2) in NRC
	// delta_A_dagger corresponds to h_{i+1} in Eq. (10.6.2), or p_{k+1} in Eq. (2.7.32) in NRC
	// grad_MI corresponds to g_{i+1} in Eq. (10.6.2), or r_{k+1} in Eq. (2.7.32) in NRC
	delta_A_dagger[k][n] = grad_MI[k][n] + beta * old_delta_A_dagger[k][n];
      }
    }
		
    for(k=0; k<K; k++){
      tempsum += inner_product(delta_A_dagger[k], grad_MI[k]);
    }
    // check that inner_product(delta_A_dagger, grad_MI) is positive
    if(real(tempsum) > 0 ){
      beta_count++;
    }else{
      swap2D(delta_A_dagger, grad_MI);
      beta_count=0;
    }
  }
  else{
    swap2D(delta_A_dagger, grad_MI);
    beta_count=0;
  }

  // backup grad_MI to old_grad_MI; after that it doesn't matter what happens to grad_MI so we can just swap
  swap2D(grad_MI, old_grad_MI);
	
  // backup delta_A_dagger to old_delta_A_dagger; cannot swap because delta_A_dagger is required for calculating new_A_dagger
  for(k=0;k<K;k++){
    for(n=0;n<N;n++){
      old_delta_A_dagger[k][n] = delta_A_dagger[k][n];
    }
  }
}

void SOMIM::calc_grad_MI(vector< vector<Complex> >& m, vector< vector< vector<Complex> > >& a, vector< vector<Complex> >& b)
{
  int i, j, k;
  vector< vector<Complex> > GAMMA = Init2DMatrix(N, N);
  vector< vector<Complex> > R_minus_GAMMA = Init2DMatrix(N, N);
  vector<Complex> RA_dagger(N, 0);
	
  for(k=0; k<K; k++){
    matrix_times_vector(RA_dagger, a[k], b[k]);
    outer_product_add(GAMMA, RA_dagger, b[k]);
  }
  transpose_conjugate(GAMMA);

  for(k=0; k<K; k++){
    // subtract GAMMA from R
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
	R_minus_GAMMA[i][j] = a[k][i][j] - GAMMA[i][j];
      }
    }
    matrix_times_vector(m[k], R_minus_GAMMA, b[k]);
  }
  // when this routine is executed till this line, grad_MI[k] = (R_k - \sum_l \pi_l R_l)A_dagger_k
}

void SOMIM::calc_L_inverse(vector< vector<Complex> >& L, vector< vector<Complex> >& L_inv, vector< vector<Complex> >& S)
{
  // calculate L where S = L x L_dagger, where L is a lower triangular matrix, then find L_inverse by back substitution
  // first calculate L
  // L_ij = 1/L_jj (S_ij - \sum_{k=0}^{j-1}L_ik L*_jk), i>j
  // L_ii = \sqrt{S_ii - \sum_{k=0}^{i-1} L_{ik}L*_{jk}}
  int row, col, k;
  Complex sum;
  // This routine computes the LU decomposition of S. Since S is hermitian, U=dagger of L.
  // The computation is well explained in section 2.3 of NRC, p. 43-49.
  for(row=0; row<N; row++)
    {
      for(col=0; col<row+1; col++)
	{
	  // first obtain summation \sum_{k=0}^{j-1}L_ik L*_jk
	  sum = 0.0;
	  for(k=0; k<col; k++)
	    {
	      sum = sum + L[row][k] * conj(L[col][k]);
	    }

	  if(col < row)
	    {
	      L[row][col] = (S[row][col] - sum) / conj(L[col][col]);
	    }
	  else
	    {
	      L[row][col] = sqrt(S[row][col] - sum);
	    }
	}
    }
	
  // now calculate L_inverse, assuming that L_inverse is a lower triangular matrix as well
  // for i!=j, \sum_{k=j}^{k=i-1} L_ik L^{-1}_kj = - L_ii L^{-1}_ij
  for(col=0; col<N; col++)
    {
      L_inv[col][col] = 1.0 / L[col][col];
      for(row=col+1; row<N; row++)
	{
	  sum = 0;
	  for(k=col; k<row; k++)
	    {
	      sum = sum + L[row][k] * L_inv[k][col];
	    }
	  L_inv[row][col] = -sum / L[row][row];
	}
    }
}

inline double SOMIM::calc_MI(double some_alpha, vector< vector<Complex> >& new_A_dagger, vector< vector<Complex> >& A_dagger, vector< vector<Complex> >& delta_A_dagger)
{
  calc_new_A_dagger(some_alpha, new_A_dagger, A_dagger, delta_A_dagger);
  calc_prob_table(new_A_dagger);
  return calc_MI_from_prob_table(prob_table);
}

double SOMIM::calc_MI_from_prob_table(vector< vector<double> >& a)
{
  int j, k;
  double temp = 0.0;
  for(j=0; j<J; j++)
    {
      for(k=0; k<K; k++)
	{
	  temp = temp + a[j][k] * ( LOG2(a[j][k]) - LOG2(a[j][K]) - LOG2(a[J][k]) );
	}
    }
  return temp;
}

void SOMIM::calc_new_A_dagger(double some_alpha, vector< vector<Complex> >& new_A_dagger, vector< vector<Complex> >& A_dagger, vector< vector<Complex> >& delta_A_dagger)
{
  int k, n;
  vector< vector<Complex> > S = Init2DMatrix(N,N);
  vector< vector<Complex> > tmp_S_inverse_sqrt = Init2DMatrix(N,N);
  vector< vector<Complex> > S_inverse_sqrt = Init2DMatrix(N,N);
  vector< vector<Complex> > tmp_A_dagger = Init2DMatrix(K,N);
	
  for(k=0; k<K; k++){
    for(n=0; n<N; n++){
      new_A_dagger[k][n] = A_dagger[k][n] + some_alpha * delta_A_dagger[k][n];
    }
  }

  for(k=0; k<K; k++){
    outer_product_add(S, new_A_dagger[k], new_A_dagger[k]);
  }
	
  // calc_S_inverse_sqrt();
  calc_L_inverse(tmp_S_inverse_sqrt, S_inverse_sqrt, S);

  for(k=0; k<K; k++){
    matrix_times_vector(tmp_A_dagger[k], S_inverse_sqrt, new_A_dagger[k]);
    scalar_times_vector(new_A_dagger[k], 1.0, tmp_A_dagger[k]);
  }
}

void SOMIM::calc_prob_table(vector< vector<Complex> >& any_A_dagger)
{
  int j, k;
  double temp;
  vector< vector<Complex> > tmp_A_dagger = Init2DMatrix(K,N);
  prob_table = Init2DRealMatrix(J+1, K+1);
	
  for(j=0; j<J; j++)
    {
      temp = 0.0;
      for(k=0; k<K; k++)
	{
	  matrix_times_vector(tmp_A_dagger[k], rho[j], any_A_dagger[k]);
	  prob_table[j][k]=real(inner_product(any_A_dagger[k], tmp_A_dagger[k]));
	  if (prob_table[j][k] < 0.0) prob_table[j][k] = 0.0;
	  temp += prob_table[j][k];
	}
      prob_table[j][K] = temp;
    }
  for(k=0; k<K; k++)
    {
      temp = 0.0;
      for(j=0; j<J; j++)
	{
	  temp += prob_table[j][k];
	}
      prob_table[J][k] = temp;
    }
}

void SOMIM::calc_R(vector< vector< vector<Complex> > >& a) // R depends on rho and prob_table
{
  int j, k, row, col;
  for(k=0; k<K; k++){
    for(row=0; row<N; row++){
      for(col=0; col<N; col++){
	a[k][row][col]=0.0;
	for(j=0; j<J; j++){
	  a[k][row][col] += rho[j][row][col] * LOG2( prob_table[j][k] / (prob_table[j][K] * prob_table[J][k]) );
	}
      }
    }
  }
}

// Given the direction of change, i.e. delta_A_dagger, this routine attempts to find the size of step to take in the direction of change, i.e. alpha.
// The computation of alpha is carried out following the routine golden() in p. 401-402 of NRC. The idea of the algorithm is explained in section 10.2 from p. 397-400 NRC.
void SOMIM::find_alpha_GSS(vector< vector<Complex> >& new_A_dagger, vector< vector<Complex> >& A_dagger, vector< vector<Complex> >& delta_A_dagger)
{
  double alphaR, alphaL, MIR, MIL, interval, current_max;
  alpha = 1.0;
  new_MI = calc_MI(alpha, new_A_dagger, A_dagger, delta_A_dagger);
	
  /* find an upper bound for the optimal alpha */
  if(new_MI > MI)
    {
      current_max = new_MI;
      while(1)
	{
	  /* increase alpha until MI(alpha) < MI */
	  alpha = alpha + 1.0;
	  new_MI = calc_MI(alpha, new_A_dagger, A_dagger, delta_A_dagger);
	  if(isnan(new_MI) || new_MI < current_max)
	    {
	      alpha = alpha - 1.0; break; // p. 399 in NRC
	    }
	  else current_max = new_MI;
	  if(alpha > 100.0) break; // limit alpha
	}
    }
  else
    {
      while(1)
	{
	  /* half alpha until MI(alpha) > MI */
	  alpha = alpha * 0.5;
	  if (alpha < tol_alpha) break;
	  new_MI = calc_MI(alpha, new_A_dagger, A_dagger, delta_A_dagger);
	  if(new_MI > MI) {alpha = alpha * 2.0; break;}
	}
    }

  /* apply Golden Section Search */
  alphaR = phi * alpha;
  alphaL = (phi-1.0) * alpha;
  MIR = calc_MI(alphaR, new_A_dagger, A_dagger, delta_A_dagger);
  MIL = calc_MI(alphaL, new_A_dagger, A_dagger, delta_A_dagger);
  interval = alpha / (2.0*phi + 1.0);
  while(interval > tol_alpha * (alphaR + alphaL))
    {
      interval = interval * (phi-1.0);
      if (MIR > MIL)
	{
	  alphaL = alphaR;
	  MIL = MIR;
	  alphaR = alphaR + interval;
	  MIR = calc_MI(alphaR, new_A_dagger, A_dagger, delta_A_dagger);
	}
      else
	{
	  alphaR = alphaL;
	  MIR = MIL;
	  alphaL = alphaL - interval;
	  MIL = calc_MI(alphaL, new_A_dagger, A_dagger, delta_A_dagger);
	}
    }
  alpha = (alphaR + alphaL) * 0.5;
}

void SOMIM::gen_A_dagger(vector< vector<Complex> >& A_dagger, vector< vector<Complex> >& S)
{
  int i, j, k, n;
  vector< vector<Complex> > tmp_A_dagger = Init2DMatrix(K, N);
  vector< vector<Complex> > tmp_S_inverse_sqrt = Init2DMatrix(N, N);
  vector< vector<Complex> > S_inverse_sqrt = Init2DMatrix(N, N);
	
  for(n=0; n<N; n++){
    A_dagger[n][n] = 1.0;
  }
	
  for(k=0; k<K; k++){
    for(n=0; n<N; n++){
      Complex random = complex<double>(0.002*((rand()/(double)RAND_MAX) - 0.5), 0.002*((rand()/(double)RAND_MAX)-0.5));
      A_dagger[k][n] = A_dagger[k][n] + random;
    }
  }
	
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      S[i][j] = 0.0;
    }
  }

  int test = 0;
  for(k=0; k<K; k++){
    outer_product_add(S, A_dagger[k], A_dagger[k]);
  }
	
  calc_L_inverse(tmp_S_inverse_sqrt, S_inverse_sqrt, S);

  for(k=0; k<K; k++){
    matrix_times_vector(tmp_A_dagger[k], S_inverse_sqrt, A_dagger[k]);
    scalar_times_vector(A_dagger[k], 1.0, tmp_A_dagger[k]);
  }
}


//-------------------------------------------------------------------------------------
//	functions to initialize 2D and 3D zero matrices, with all entries set to 0
//-------------------------------------------------------------------------------------
vector< vector<Complex> > SOMIM::Init2DMatrix(int X, int Y)
{
  vector<Complex> b(Y, 0); /* 2nd dimension */
  vector< vector<Complex> > a(X, b);
	
  return a;
}

vector< vector<double> > SOMIM::Init2DRealMatrix(int X, int Y)
{
  vector<double> b(Y, 0); /* 2nd dimension */
  vector< vector<double> > a(X, b);
	
  return a;
}

vector< vector< vector<Complex> > > SOMIM::Init3DMatrix(int X, int Y, int Z)
{
  vector<Complex> c(Z, 0); /* 3rd dimension */
  vector< vector<Complex> > b(Y, c); /* 2nd dimension */
  vector< vector< vector<Complex> > > a(X, b);
	
  return a;
}
//------------------------------------------------------------------------------
//	functions to multiply scalar, vectors and matrices
//------------------------------------------------------------------------------

Complex SOMIM::inner_product(vector<Complex>& v1, vector<Complex>& v2)
{
  int i;
  Complex temp = 0.0;
  for(i=0; i<N; i++)
    {
      temp += conj(v1[i]) * v2[i];
    }
  return temp;
}

double SOMIM::LOG2(double x)
{
  if(x==0.0) return 0.0;
  return log(x)/log(2.0);
}

// add two N by N matrix m2 and m3 and assign to m2;
void SOMIM::matrix_add_matrix(vector< vector<Complex> >& m2, vector< vector<Complex> >& m3)
{
  int i, j;
  for(i=0; i<N; i++)
    {
      for(j=0; j<N; j++)
	m2[i][j] += m3[i][j];
    }
}

// multiply two N by N matrix m2 and m3 and assign the product to m1
void SOMIM::matrix_times_matrix(vector< vector<Complex> >& m1, vector< vector<Complex> >& m2, vector< vector<Complex> >& m3)
{
  int i, j, k;
  for(i=0; i<N; i++)
    {
      for(k=0; k<N; k++)
	{
	  m1[i][k] = 0;
	  for(j=0; j<N; j++)
	    {
	      m1[i][k] += m2[i][j] * m3[j][k];
	    }
	}
    }
}

// multiply an N by N matrix m with an N-element vector v2 and assign the product to v1
void SOMIM::matrix_times_vector(vector<Complex>& v1, vector< vector<Complex> >& m, vector<Complex>& v2)
{
  int i, j;
  for(i=0; i<N; i++)
    {
      v1[i] = 0;
      for(j=0; j<N; j++)
	{
	  v1[i] += m[i][j] * v2[j];
	}
    }
}

void SOMIM::outer_product_add(vector< vector<Complex> >& m, vector<Complex>& v1, vector<Complex>& v2)
{
  int row, col;
  Complex temp;
  for(col=0; col<N; col++)
    {
      temp = conj(v2[col]);
      for(row=0; row<N; row++)
	m[row][col] += v1[row] * temp;
    }
}

// multiply an N-element vector v2 with a scalar s and assign the product to an N-element vector v1
void SOMIM::scalar_times_vector(vector<Complex>& v1, double s, vector<Complex>& v2)
{
  int i;
  for(i=0; i<N; i++)
    {
      v1[i] = s * v2[i];
    }
}

void SOMIM::swap2D(vector< vector<Complex> >& m1, vector< vector<Complex> >& m2)
{
  m1.swap(m2);
}

void SOMIM::transpose_conjugate(vector< vector<Complex> >& m)
{
  int row, col;
  Complex temp;
  for(row=0; row<N; row++){
    for(col=0; col <= row; col++){
      temp = m[row][col];
      m[row][col] = conj(m[col][row]);
      m[col][row] = conj(temp);
    }
  }
}


