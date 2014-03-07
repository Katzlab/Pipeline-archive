// $Id: definitions.h 10679 2012-05-29 19:04:27Z cohenofi $

#ifndef ___DEFINITIONS_H
#define ___DEFINITIONS_H

#ifdef _MSC_VER
#define LIMITS_WORKING
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4786)
#pragma warning (disable: 4267)
#pragma warning (disable: 4018)
#pragma warning (disable: 4305) //truncation from 'double' to 'float'
#endif
 

#include <vector>
#include <string>

#ifdef LIMITS_WORKING
	#include <limits>
#endif
using namespace std;

#define MDOUBLE double
//#define MDOUBLE float

// Contants
#define PI  (3.1415926535897932384626433832795028841971693993751058)

typedef vector<MDOUBLE> Vdouble;
typedef vector<int> Vint;
typedef vector<Vint> VVint;
typedef vector<VVint> VVVint;
typedef vector<char> Vchar;
typedef vector<Vdouble> VVdouble;
typedef vector<VVdouble> VVVdouble;
typedef vector<VVVdouble> VVVVdouble;
typedef vector<VVVVdouble> VVVVVdouble;
typedef vector<string> Vstring;

#ifdef LIMITS_WORKING
	const MDOUBLE VERYBIG = numeric_limits<MDOUBLE>::max();
	const MDOUBLE VERYSMALL = -VERYBIG;
	const MDOUBLE EPSILON = numeric_limits<MDOUBLE>::epsilon();
#else
// IF <limits> is not recognized, and MDOUBLE is double.
	const MDOUBLE VERYBIG = 1.79769e+308;
	const MDOUBLE VERYSMALL = -VERYBIG;
	const MDOUBLE EPSILON = 2.22045e-016;
#endif

//The maximum value for type float is:  3.40282e+038
//The maximum value for type double is:  1.79769e+308
//::epsilon() returns the difference between 1 and the smallest value greater than 1 that is representable for the data type.
//epsilon float 1.19209e-007
//epsilon double 2.22045e-016

#ifdef LOGREP
	class logRep;
	typedef vector<logRep> VlogRep;
	typedef vector <vector<logRep> > VVlogRep;
	typedef vector< vector <vector<logRep> > > VVVlogRep;
	typedef logRep doubleRep;
	typedef VlogRep VdoubleRep;
	typedef VVlogRep VVdoubleRep;
	typedef VVVlogRep VVVdoubleRep;
	#include "logRep.h"
#elif defined (DOUBLEREP)
	class doubleRepMantisa;
	typedef vector<doubleRepMantisa> VdoubleRepMantisa;
	typedef vector <vector<doubleRepMantisa> > VVdoubleRepMantisa;
	typedef vector <VVdoubleRepMantisa > VVVdoubleRepMantisa;
	typedef vector <VVVdoubleRepMantisa > VVVVdoubleRepMantisa;
	typedef doubleRepMantisa doubleRep;
	typedef VdoubleRepMantisa VdoubleRep;
	typedef VVdoubleRepMantisa VVdoubleRep;
	typedef VVVdoubleRepMantisa VVVdoubleRep;
	typedef VVVVdoubleRepMantisa VVVVdoubleRep;
	#include "doubleRep.h"
#else
    typedef MDOUBLE  doubleRep;
	typedef Vdouble  VdoubleRep;
	typedef VVdouble VVdoubleRep;
	typedef VVVdouble VVVdoubleRep;
	typedef VVVVdouble VVVVdoubleRep;
	inline MDOUBLE convert (MDOUBLE d) {return (d);}
#endif

#endif
  

