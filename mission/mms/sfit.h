//#ifndef _SFIT_H
#define _SFIT_H 1

// Maximum of terms used for fit is 9;
#define MAXTERMS_FIT 9

// FillVal for NaN
#define NaN -159e7

void spinfit(const int maxIt, const int minPts, const int nTerms, const long t0, const int nSegments,
	const int nData, const double te[], const double az[], const double pha[], const double fitInterv, const double fitEvery,
	double ts[], double sfit[], double sdev[], double iter[], double nout[]);

int onesfit (const int nTerms, const int maxIter, int nIter, int flim, 
	const int nData, const double phaseArray[], const double dataArray[], 
	int &nBadPoints, double x[], double &sigma);
	
int solve(double A[MAXTERMS_FIT][MAXTERMS_FIT+1], const int nTerms,
	double X[MAXTERMS_FIT]);
//#endif
