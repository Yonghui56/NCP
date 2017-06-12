/* File : example.i */
%module ogsPython

%include "std_vector.i"

%include "typemaps.i"
%template(VecI)  std::vector<int>;
%template(VecD)  std::vector<double>;

%apply std::vector<int> *INOUT { std::vector<int> *rangevec };
%apply std::vector<double> *OUTPUT { std::vector<double> *x };
%apply std::vector<double> *OUTPUT { std::vector<double> *results };

%{
/*	#define SWIG_FILE_WITH_INIT */
	#include "example.h"
%}

/* %include "numpy.i" */
/* %include "carrays.i" */

/*
%init %{
     import_array();
%}

%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* rangevec, int n)}
*/



/* Let's just grab the original header file here */
%include "example.h"

