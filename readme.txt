ReadMe file for the SVARGS (Sparse Vector Auto-Regressive Greedy Search) implementation in Matlab
_________________________________________________________________________________________________
1. All folders should be added to the Matlab path.

2. The file test.m can be run to give a sample result and model assessment.
	The number of variables (k) can be varied on the first line of the file test.m.
	number of samples (nsamples) can be varied on the second line of the file test.m.

3. Before running test.m, two .cpp files ("turningpts.cpp" and "inv_std_normal.cpp") need to be compiled.
	This can be done using the following Matlab commands:
		mex turningpts.cpp
		mex inv_std_normal.cpp

4. Updated versions of this SVARGS implementation can be obtained from the repository:
https://github.com/emeskaykay/SVARGS.git

For questions about the code or the algorithm contact the authors at
modysk@gmail.com
