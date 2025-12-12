/* Copyright 2013 Michele Dall'Arno
// This file is part of SOMIM.
SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

#include "accinfo.h"

void setRho(ComplexNDArray* E, SOMIM* somim) {
  int N = E->dims()(0);
  int J = E->dims()(2);
  for (int j = 0; j < J; j++)
    for (int n = 0; n < N; n++)
      for (int m = 0; m < N; m++)
	somim->setRho(j,n,m, E->elem(n,m,j));
}

void setPOVM(ComplexNDArray* P, SOMIM* somim) {
  for (int k = 0; k < P->dims()(2); k++)
    for (int n = 0; n < P->dims()(0); n++)
      for (int m = 0; m < P->dims()(1); m++)
	somim->setPOVM(k,n,m, P->elem(n,m,k));
}

void getPOVM(ComplexNDArray* P, SOMIM* somim) {
  for (int k = 0; k < P->dims()(2); k++)
    for (int n = 0; n < P->dims()(0); n++)
      for (int m = 0; m < P->dims()(1); m++)
	P->elem(n,m,k) = somim->getPOVM(k,n,m);
}

DEFUN_DLD (accinfo, args, , USAGE)
{
  double A;
  int N, M, J, K;
  ComplexNDArray E, P;
  int nargin = args.length ();
  octave_value_list retval;

  if (nargin == 0) {print_usage (); return octave_value_list (retval);}

  if (nargin >= 1) {
    E = args(0).complex_array_value ();
    N = E.dims()(0); M = E.dims()(1); J = E.dims()(2); K = N*N;
  }

  if (nargin >= 2) {
    P = args(1).complex_array_value ();
    if (P.dims()(0) == 1) {
      K = (int) P.elem(0).real();
      if (K < N || K > N*N) {print_usage (); return octave_value_list (retval);}
    }
    else K = P.dims()(2);
  }

  SOMIM* somim = new SOMIM(N, J, K);

  if (nargin >= 3) somim->setDirectGrad(args(2).matrix_value().elem(0));
  if (nargin >= 4) somim->setTolerance(args(3).matrix_value().elem(0));

  setRho(&E, somim);
  if (P.dims()(0) > 1) setPOVM(&P, somim);
  
  A = somim->getMI();
  K = somim->getK();

  dim_vector dv(N, N, K);
  ComplexNDArray Q(dv);
  getPOVM(&Q, somim);

  retval(0) = octave_value (A);
  retval(1) = octave_value (Q);

  return octave_value_list (retval);
}

