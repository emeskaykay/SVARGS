#include "mex.h"
#include "inv_stdnormal_cdf.c"

/* Copyright 2021, Sandeep K Mody, Govindan Rangarajan
This file is part of the SVARGS program.
SVARGS is free software: you can redistribute it and/or modify it under the terms
of the GNU Lesser General Public License as published by the Free Software Foundation,
either version 2 of the License, or (at your option) any later version.

SVARGS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with
SVARGS. If not, see <https://www.gnu.org/licenses/>.
*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nrhs != 1)
	{ mexErrMsgIdAndTxt("mexFunction:inv_std_normal", "Exactly one input required. "); }
	
	if (nlhs > 1)
	{ mexErrMsgIdAndTxt("mexFunction:inv_std_normal", "Only one output required. "); }
	
	if (!mxIsDouble(prhs[0]))
	{ mexErrMsgIdAndTxt("mexFunction:inv_std_normal", "Argument must be a double precision."); }
	
	const double* p_conf = (const double*)mxGetData(prhs[0]);
	double conf = *p_conf;
	if ( !((conf >= 0) && (conf <= 1)) )
	{ mexErrMsgIdAndTxt("mexFunction:inv_std_normal", "Argument must be in the interval [0,1]."); }
	
	double x = ltqnorm(conf); //inv_stdnormal_cdf.c (in folder __external)
	plhs[0] = mxCreateDoubleScalar(x);
}
