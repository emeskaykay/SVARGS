#include "mex.h"
#include <algorithm> //for std::copy, std::max_element, std::min_element
#include<cmath> //for std::llround
#include<vector>

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
		  
// Helper function for the main mexFunction.
// Finds the half-width of the peak at location r.
// The location r must be the index of a peak.
double halfwidth__(int m, const double* p, int r)
{
	int left = r - r%m;
	int right = left + m;
	double halfmax = p[r]/2;
	
	int l = r;
	for ( ; (p[l] > halfmax) && (p[l-1] <= p[l])  && (l > left); --l);

	int h = r;
	for ( ; (p[h] > halfmax) && (p[h+1] <= p[h]) && (h < right); ++h);
	
	return (double)(h - l);
}

// Find the turning points (local maxima and minma) along the first dimension of the input array.
// INPUTS
// prhs[0]: Numeric array of type double.
//
// OUTPUTS
// plhs[0]: Boolean array indicating locations of local maxima in the input array.
// plhs[1]: Boolean array indicating locations of local minima in the input array.
// plhs[2]: struct with two fields:
//					widths_max: maximum half-bandwidth for each line vector.
//					widths_min: minimum half-bandwidth for each line vector.
//
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nrhs != 1)
	{
		mexErrMsgIdAndTxt("mexFunction:turningpts: ", "Exactly one input required. ");
	}
	
	const mxArray* p_array = prhs[0];
	if ( !mxIsDouble(p_array) )
	{
		mexErrMsgIdAndTxt("mexFunction:turningpts: ", "Input must be numeric double array. ");
	}
	
	if ((nlhs < 1) || (nlhs > 3))
	{
		mexErrMsgIdAndTxt("mexFunction:turningpts: ", "Number of outputs must be between 1 and 3. ");
	}
	
	mwSize ndims = mxGetNumberOfDimensions(p_array);
	const mwSize* psize = mxGetDimensions(p_array);
	int m = psize[0];
	int M = mxGetNumberOfElements(p_array);
	
	plhs[0] = mxCreateLogicalArray(ndims, psize);
	bool* pbmax = (bool*)mxGetData(plhs[0]);
	
	bool* pbmin = nullptr;
	if (nlhs >= 2)
	{
		plhs[1] = mxCreateLogicalArray(ndims, psize);
		pbmin = (bool*)mxGetData(plhs[1]);
	}
	
	double* pwa = nullptr;
	double* pwb = nullptr;
	if (nlhs == 3)
	{
		mxArray* pwidths_max = mxCreateNumericArray(ndims-1, psize+1, mxDOUBLE_CLASS, mxREAL);
		mxArray* pwidths_min = mxCreateNumericArray(ndims-1, psize+1, mxDOUBLE_CLASS, mxREAL);
		const char* fnames[] = {"widths_max", "widths_min"};
		plhs[2] = mxCreateStructMatrix(1, 1, 2, fnames);
		mxSetField(plhs[2], 0, "widths_max", pwidths_max);
		mxSetField(plhs[2], 0, "widths_min", pwidths_min);
		pwa = mxGetPr(pwidths_max);
		pwb = mxGetPr(pwidths_min);
	}
	
	const double* p = mxGetPr(p_array);
	double dr = 0;
	double ds = 0;
	double mid = 0;
	bool plat = false;
	std::vector<double> hwidths;
	for (int i = 1; i <= M; ++i)
	{
		if (i%m == 0)
		{
			dr = 0;
			ds = 0;
			plat = false;
			if ( (pwa) && (!hwidths.empty()) )
			{
				pwa[(i-1)/m] = *max_element(hwidths.begin(), hwidths.end());
				pwb[(i-1)/m] = *min_element(hwidths.begin(), hwidths.end());
				hwidths.clear();
			}
			continue;
		}
		
		double dt = p[i] - p[i-1];
		if (plat)
		{
			if (dt == 0)
			{
				mid += 0.5;
			}
			else
			{
				plat = false;
				if ((dr > 0) && (dt < 0))
				{
					int loc = std::llround(mid);
					pbmax[loc] = true;
					if (pwa)
					{
						hwidths.push_back(halfwidth__(m, p, loc));
					}
				}
				else if (pbmin && (dr < 0) && (dt > 0))
				{
					pbmin[std::llround(mid)] = true;
				}
			}
		}
		else
		{
			mid = i-1;
			if (dt == 0)
			{
				plat = true;
				dr = ds;
				mid += 0.5;
			}
			else if ((ds > 0) && (dt < 0))
			{
				int loc = std::llround(mid);
				pbmax[loc] = true;
				if (pwa)
				{
					hwidths.push_back(halfwidth__(m, p, loc));
				}
			}
			else if (pbmin && (ds < 0) && (dt > 0))
			{
				pbmin[std::llround(mid)] = true;
			}
		}
		ds = dt;
	}
}
