/* Copyright 2013 Michele Dall'Arno
// This file is part of SOMIM.
SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

#include "libsomim.h"
#include <octave/oct.h>

const string USAGE = "\
 -- Loadable Function: [A, Q] = accinfo (E)\n\
 -- Loadable Function: [A, Q] = accinfo (E, P)\n\
 -- Loadable Function: [A, Q] = accinfo (E, P, GRAD)\n\
 -- Loadable Function: [A, Q] = accinfo (E, P, GRAD, TOL)\n\
\n\
Evaluates the accessible information A and the optimal POVM\n\
Q for quantum ensemble E.\n\
\n\
Copyright 2013 Michele Dall'Arno (michele.dallarno@gmail.com).\n\
\n\
Based on SOMIM (http://www.quantumlah.org/publications/software/SOMIM).\n\
\n\
\n\
\n";
