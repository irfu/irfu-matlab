//
//  Spin fit routine for MMS
//  based on similar code for Cluster (c_efw_spinfit_mx.cpp)
//
//  Modified to allow for overlapping intervals and allow for period to be argument.
//


#include "mex.h"
#include "sfit.h"

#include "cmath"



///////////////////////
// Sub function "solve"
/////////////////////// 

int solve(double A[MAXTERMS_FIT][MAXTERMS_FIT+1], const int nTerms,
	double X[MAXTERMS_FIT])
{
/*
  Equation solver.

  Interface:
  A 	 - equation system
  nTerms - number of equations
  X 	 - result

  Returns:
  error indicator
*/

	bool notFound = true;

	//Find non zero entry in jth column
	for ( int j=0; j<nTerms-1; j++ ) {
		int mm;
		for ( int k=j; k<nTerms; k++ ) {
			mm = k;
			if ( std::abs(A[k][j]) > 1.E-12 ) {
				notFound = false;
				break;
			}
		}
		if (notFound)
			return j;
			
		if (mm != j) {
			//Interchange mmth row with jth row
			for ( int k=j; k<=nTerms; k++) {
				double work = A[mm][k]; 
				A[mm][k] = A[j][k];
				A[j][k] = work;
			}
		}
	
		//Subtract jth row from subsequent rows
		for ( int k=j+1; k<nTerms; k++ ) {
			double y = A[k][j]/A[j][j];
			for ( int l=j; l<=nTerms; l++ )
				A[k][l] -= y*A[j][l];
		}
	}

	//Now solve
	for ( int j=0; j<nTerms; j++ ) {
		int m = nTerms -1 -j;
		X[m] = A[m][nTerms]/A[m][m];
		for ( int k=0; k<m; k++ )
			A[k][nTerms] -= X[m]*A[k][m];
	}

	return (int) 0;
} // END of subfunction "solve"


////////////////////////
// Subfunction "onesfit"
////////////////////////

int onesfit (const int nTerms, const int maxIter, int nIter, int flim, 
	const int nData, const double phaseArray[], const double dataArray[], 
	int &nBadPoints, double x[MAXTERMS_FIT], double &sigma) {
/*
Function name: ONESFIT

Description:
Fit x(1)+x(2)*cos(pha)+x(3)*sin(pha)
         +x(4)*cos(2*pha)+x(5)*sin(2*pha)+... to data

Input:
	nTerms		- number of terms to fit
	maxIter		- maximum number of iterations
	nData		- number of data points
	phaseArray	- phase
	dataArray	- data to fit
Output:
	nIter		- number of iterations performed
	flim 		- ?
	nBadPoints	- number of points disregarded from the fit
	x			- array for resulting coefficients from fit
	sigma		- output value
Returns:
	error indicator (0-success)
*/

	const double cnst0 = 1.4;	// XXX: move to header
	const double dcnst = 0.4;	// XXX: move to header
	bool badPoint[nData];		// Array of bad points
	double s[MAXTERMS_FIT][MAXTERMS_FIT+1];
	double q[MAXTERMS_FIT][MAXTERMS_FIT+1];
	double cnst = cnst0;
	int ier=-1; // Define returning error value.
	double adiff[nData];
	
	nBadPoints = 0;
	
	// Verify inputs, nData > nTerms, nTerms < maxterm, 
    // nTerms odd (ie x(1)+x(2)sin(w*phase)+x(3)cos(w*phase)).
	if ( nData<nTerms+1 || nTerms>MAXTERMS_FIT || nTerms%2==0 )
		return ier;

	// Build normal equations system
	for ( int row=0; row<nTerms; row++ )
	  for ( int col=0; col < nTerms+1; col++ )
	    s[row][col] = 0.0;

	// Add to normal equations
	for ( int i=0; i<nData; i++) {
		double w[MAXTERMS_FIT+1];
		w[0] = 1;
		for ( int row=2; row<nTerms; row+=2 ) {
			double arg = (double)(row/2) * phaseArray[i];
			w[row-1] = cos(arg);
			w[row] = sin(arg);
		}
		w[nTerms] = dataArray[i];

		for ( int row=0; row<nTerms; row++ ) {
			for ( int col=row; col<=nTerms; col++ )
				s[row][col] += w[row]*w[col];
		}
		// Assume it is not a badPoint..
		badPoint[i] = false;
	} // End of for loop, i.
	
	flim = nData;
	
	// Start of iteration loop    
	for ( int iter=1; iter<=maxIter; iter++) {
		// Store the number of Iterations used.
		nIter = iter;
		// Solve normal equations
		if (flim < nTerms+1) {
			ier = -1;
			break;
		}
		
		for ( int row=0; row<nTerms; row++) {
			// diag(q) = diag(s)
			q[row][row] = s[row][row];
			
	    	for ( int col = row+1; col<nTerms; col++ ){
				q[row][col] = q[col][row] = s[row][col];             
			}
	    	q[row][nTerms] = s[row][nTerms];
		}

		// Solve
	  	ier = solve (q,nTerms,x);
	  	if ( ier != 0)
	  		break;		
				
		// Compute sigma
	  	sigma = 0.0;
		for ( int i=0; i<nData; i++ ) {
			if ( badPoint[i] )
				continue;
				
			double w[MAXTERMS_FIT+1];
			w[0] = 1;
			for ( int row=2; row<=nTerms; row+=2 ) {
				double arg = (double)(row/2) * phaseArray[i];
				w[row-1] = cos(arg);
				w[row] = sin(arg);
			}
			double y = 0.0;
			for ( int row=0; row<nTerms; row++ )
				y += x[row] * w[row];
			double diff = dataArray[i] - y;
			adiff[i] = diff;
			sigma += diff*diff;  
		}
		
        sigma = sqrt(sigma/double(flim-1));

	  	if ( nIter<maxIter) {
			double ref = cnst*sigma;

			// Search badPoint points
			bool flagChanged = false;
			for ( int i=0; i<nData; i++) {
				if ( !badPoint[i] && std::abs(adiff[i])>ref ) {
					// Subtract from normal equations
					double w[MAXTERMS_FIT+1];
					w[0] = 1;
					for ( int row=2; row<=nTerms; row+=2 ) {
						double arg = (double)(row/2) * phaseArray[i];
						w[row-1] = cos(arg);
						w[row] = sin(arg);
					}
					w[nTerms] = dataArray[i];
					for ( int row=0; row<nTerms; row++ ) {
						for ( int col=0; col<=nTerms; col++ )
							s[row][col] += w[row]*w[col];
					}
					flim = flim - 1;
					badPoint[i] = true;
					flagChanged = true;
				}
			}
	    	if ( !flagChanged || flim<=1 )
				break;
			cnst = cnst + dcnst;
		}
	}
	
	for ( int i=0; i<nData; i++)
		if ( badPoint[i] )
			nBadPoints++;

	return ier;
} // End of subfunction "onesfit"



////////////////////////
// Subfunction "spinfit"
////////////////////////

void spinfit(const int maxIt, const int minPts, const int nTerms, const double t0, const double tEnd, const int nSegments,
	const int nData, const double te[], const double az[], const double pha[], const double fitInterv, const double fitEvery,
	double ts[], double sfit[], double sdev[], double iter[], double nout[])
{
	// Fill output time, ts, as each fitEvery interval and default other outputs to NaN.
	for ( int i=0; i<nSegments; i++){
		ts[i] = (double)(i)*fitEvery + (double)t0;
		// NaN = default for all other output
		sdev[i] = NaN;
		iter[i] = NaN;
		nout[i] = NaN;
		for (int j=0; j<nTerms; j++){
			sfit[j*nSegments +i] = NaN;
		} // End of for loop, j
	} // End of for loop, i

	// Check if we have enough data for at least one fit.
	if (nData < minPts)
		return;
	// If nSegments are not enough for on fit, exit with fillVal only.
	if (nSegments < 1)
		return;

	int idx = 0;
	for ( int i=0; i<nSegments; i++)
	{
		double startT = 0;
		//mexPrintf("\nSegment: %i, t0: %f, tEnd: %f, ts[i]: %f, fitInterv: %f.", i,t0,tEnd,ts[i],fitInterv );
		if ( te[0] >= (ts[i] - fitInterv/2.0) ) {
			startT = te[0];
		} else if ( ((ts[i] - fitInterv/2.0)>=te[0]) && (ts[i]+fitInterv/2.0<=tEnd) ){
			startT = ts[i] - fitInterv/2.0;
		} else {
			startT = tEnd - fitInterv;
		}

		bool foundStart = false;
		int idxs = 0; // Index of start, for spinfit, i.
		int idxe = 0; // Index of end, for spinfit, i.

		//mexPrintf("\nSegment: %i, startT: %f.", i,startT);

		while (idxe == 0)
		{
			//mexPrintf("\nSegment: %i, startT: %f, idx: %i, idxs: %i, idxe: %i.", i,startT,idx,idxs,idxe);

			// first point alredy more then one spin later then the first one <-- Wait, WHAT??
			if ( (idxs == 0) && (te[idx] >= (startT+fitInterv) ) && !foundStart ) {
				idxe = -1; // Provides exit out of while loop. No points to do spinfit on.
			}
			else
			{
				if (idxs == 0 && te[idx] >= startT && !foundStart){
					// idx is the first data point in spinStart interval (spinStart to spinStart+fitInterv)
					idxs = idx;
					foundStart = true;
				}
				if (idx==nData-1 || te[idx+1] > startT+fitInterv){
					// idx is the last data point in data or just outside interval (spinStart+fitInterv)
					if(foundStart) {
						idxe = idx; // Provides exit out of while loop.
						idx = idxs-1; // Resume next iteration at the start point of previous start point.
					}
					else
					{
						// Gap in time series, no start was found but next value was well after fitInterv
						idxe = -1; // Provides exit out of while loop.
					}
				}
				// Else check next index idx.
				idx++;
			}
		} // End of while loop, idxe==0
		
		//mexPrintf("\nEnded while loop with values:\nSegment i: %i, idx: %i, startT: %f, idxs: %i, idxe: %i.", i,idx,startT,idxs,idxe);
		//mexPrintf("\n te[idxs]: %f, te[idxe]: %f, fitInterv: %f",te[idxs],te[idxe],fitInterv);

		// check number of data points in interval idxs to idxe, and verify they are at least minPts.
		int nn = idxe-idxs+1;

		//mexPrintf("\n Got nn: %i, while min was: %i",nn,minPts);
		
		if ( (nn > minPts) && foundStart )
		{
			double lim, x[MAXTERMS_FIT];
			int nIter, nBad, ierr;
			ierr = onesfit(nTerms,maxIt, nIter, lim,nn,&pha[idxs],&az[idxs], nBad, x, sdev[i]);
			if (ierr == 0)
			{
				for (int j=0; j<nTerms; j++){
					// Store each term.
					sfit[i*nTerms+j] = x[j];
				} // End of for loop, j
				iter[i] = (double)nIter;
				nout[i] = (double)nBad;
			}
		}
	} // End of for loop, i
} // End of subfunction "spinfit"


/////////////////////////
// ENTRY point for MATLAB
/////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
 	/* Check for proper number of arguments. */
	if ( nrhs != 9) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:rhs",
		"This function requires 9 input arguments.");
	}
	if ( nlhs != 5) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:lhs",
		"This function requires 5 output arguments.");
	}

	//	maxIt		argument #1
	// Maximum number of iterations to go through fitting
	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
	    mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
	    	mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:maxitNotScalar",
		"Input MAXIT must be a scalar.");
	}
	int maxIt = int(mxGetScalar(prhs[0]));
	if (maxIt < 1) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:maxitZeroNegative",
		"Input MAXIT must be a positive nonzero number.");
	}

	//	nTerms		argument #3
	// Number of terms to compute fit to, A + Bcos(w) + Csin(w) + Dcos(2w) + Esin(2w) etc..
	if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || 
	    mxGetN(prhs[2])*mxGetM(prhs[2]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:ntermsNotScalar",
		"Input NTERMS must be a scalar.");
	}
	int nTerms = int(mxGetScalar(prhs[2]));
	if ( nTerms <= 1 || nTerms > MAXTERMS_FIT || nTerms%2 != 1) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:ntermsNotOdd",
		"Input NTERMS must be one of 3,5,7.");
	}

	//	minPts		argument #2
	// Minimum of points used for fit of one spin rev.
	if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
	     mxGetN(prhs[1])*mxGetM(prhs[1]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:ntermsNotScalar",
		"Input NTERMS must be a scalar.");
	}
	int minPts = int(mxGetScalar(prhs[1]));
	if ( minPts <= nTerms) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:minptsNotLargerNTerms",
		"Input MINPTS must be larger than NTERMS.");
	}
	
	//	te		argument #4
	// Timestamp for each point of data/phase
	int nData = mxGetN(prhs[3]);
	if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || 
	     mxGetM(prhs[3]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:teNotAVector",
		"Input TE must be n x 1 vector.");
	}
	double *te = mxGetPr(prhs[3]);
	
	//	data		argument #5
	// Measurement data
	if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || 
	     mxGetM(prhs[4]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:dataNotAVector",
		"Input DATA must be n x 1 vector.");
	}
	if ( mxGetN(prhs[4]) != nData ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:dataNotNData",
		"Inputs TE and DATA must be of the same length.");
	}
	double *data = mxGetPr(prhs[4]);
	
	//	phase		argument #6
	// Phase corresponding to each data and timestamp te
	if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || 
	     mxGetM(prhs[5]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:phaseNotAVector",
		"Input PHASE must be n x 1 vector.");
	}
	if ( mxGetN(prhs[5]) != nData ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:phaseNotNData",
		"Inputs TE and PHASE must be of the same length.");
	}
	double *pha = mxGetPr(prhs[5]);
	
	// fitEvery		argument #7
	// Perform a fit every "fitEvery":th second.
	if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
		mxGetN(prhs[6])*mxGetM(prhs[6]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:fitEveryNotScalar",
		"Input FITEVERY must be a scalar.");
	}
	double fitEvery = mxGetScalar(prhs[6]);

	// fitInterv	argument #8
	// Perform a fit over interval "fitInterv" seconds.
	if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
		mxGetN(prhs[7])*mxGetM(prhs[7]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:fitIntervNotScalar",
		"Input FITINTERV must be a scalar.");
	}
	double fitInterv = mxGetScalar(prhs[7]);
	
    // Verify fitEvery <= fitInterv, (don't create gaps in time series).
    // Equal corresponds to no overlap.
    if ( fitEvery > fitInterv) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:fitEveryLargerFitInterv",
		"Input FITINTERV must be larger than equal to FITEVERY.");
	}
    
	// t00		argument #9
	if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) ||
		mxGetN(prhs[8])*mxGetM(prhs[8]) != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:mms_spinfit_mx:t00NotScalar",
		"Input t00 must be a scalar.");
	}
	double t0 = mxGetScalar(prhs[8]);

	//mexPrintf("\n Input arguments... fitInterv: %f, fitEvery: %f, te[0]: %f, nData: %i",fitInterv,fitEvery,te[0],nData);
	
    
	// Get the number of complete segments from start of data to the end and first timestamp.
	const double tEnd = te[nData-1];
	int nSegments = (int)(floor((te[nData-1]-t0)/fitEvery)+1);

	//mexPrintf("\n Calculated... nSegments: %i, t0: %f, nData: %i, tEnd: %f",nSegments,t0,nData,tEnd);

	// Pre allocate double matricies of required size.
	plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)nSegments, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)nTerms, (mwSize)nSegments, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)1, (mwSize)nSegments, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((mwSize)1, (mwSize)nSegments, mxREAL);
	plhs[4] = mxCreateDoubleMatrix((mwSize)1, (mwSize)nSegments, mxREAL);
	// Get pointers, ts = timestamp, sfit=spinfit, sdev=?, iter=iterations used, nout=?
	double *ts = mxGetPr(plhs[0]);
	double *sfit = mxGetPr(plhs[1]);
	double *sdev  = mxGetPr(plhs[2]);
	double *iter = mxGetPr(plhs[3]);
	double *nout = mxGetPr(plhs[4]);


	// Call the actual spinfit, arguments on first line here are inputs, second line output arguments.
	spinfit( maxIt,minPts,nTerms,t0,tEnd,nSegments,nData,te,data,pha,fitInterv,fitEvery,
		ts,sfit,sdev,iter,nout);
} // End of Matlab interface sub function
