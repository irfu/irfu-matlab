/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <yuri@irfu.se> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
 * ----------------------------------------------------------------------------
 *
 * irf_average_mx.c  MEX function to do averages
 *
 * RES = IRF_AVERAGE_MX(X, Y, DT2);
 *
 * Resample X to timeline of Y, using half-window of DT2
 *
 * Compile with:
 *   mex -v irf_average_mx.c CFLAGS='$CFLAGS -O2 -mtune=opteron -funroll-loops'
 *
 * $Id$
 */

#include <limits.h>
#include "mex.h"

/*
 * This typedef is needed for MATLAB < 7.3
 */
#ifndef MWSIZE_MAX
typedef int mwSize;
#endif

void mexFunction(
		 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{
    mwSize ndata, ncomp, ntref, i, comp, start = 0, stop = 0;
    double *data, *tref, *res, dt2;
	double NaN = mxGetNaN();
    
    
    /* Check for proper number of input and output arguments */    
    if ( nrhs != 3 )
		mexErrMsgTxt("Three input arguments required.");

    if ( nlhs > 1 )
		mexErrMsgTxt("Too many output arguments.");

	/* Check data type of input argument  */
    if ( !(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])) || 
			!(mxIsDouble(prhs[2])) )
	{
		mexErrMsgTxt("Input arguments must be of type double.");
    }	

    if( mxIsEmpty(prhs[0]) ) 
		mexErrMsgTxt("First input argument is empty\n");

	if( mxGetNumberOfDimensions(prhs[0]) != 2 )
		mexErrMsgTxt("First input arguments must be a 2D matrix.");

	if(mxIsEmpty(prhs[1]))
		mexErrMsgTxt("Second input argument is empty");

	if( mxGetNumberOfDimensions(prhs[1]) != 2 )
		mexErrMsgTxt("Second input arguments must be a 2D matrix.");
	
	if(mxIsEmpty(prhs[2]))
		mexErrMsgTxt("Third input argument is empty");

    data = mxGetPr(prhs[0]);
    ndata = mxGetM(prhs[0]);
	ncomp = mxGetN(prhs[0]);
	tref = mxGetPr(prhs[1]);
	ntref = mxGetM(prhs[1]);
	dt2 = mxGetScalar(prhs[2]);
	
    /* Create output array */
    plhs[0] = mxCreateDoubleMatrix(ntref,ncomp,0);
    res = mxGetPr(plhs[0]);
    for (i=0; i < ntref; i++)
		res[i] = tref[i];
	
	/* check is there is a total interval mismatch */
	if ( (data[0] > res[ntref-1] + dt2) || (data[ndata-1] <= res[0] - dt2) )
	{
		mexWarnMsgTxt("interval mismatch\n");
		for (i=0; i<ndata*ncomp; i++) 
			res[i] = NaN;
		return;
	}
	
	for (i=0; i < ntref; i++)
	{
		/* check if we have been through all the data 
		 * or that the data starts after the current interval */
		while ( (start<ndata) && (data[start] < res[i] - dt2) ) 
			start++;
		
		/*
		if ( start > stop )
			printf("interval(%d) : skipping %d points\n",i,start-stop);
		 */
		
		if ( (start==ndata) || (data[start] > res[i] + dt2) )
		{
			for (comp=1; comp<ncomp; comp++) 
				res[i+comp*ntref] = NaN;
			/*
			printf("interval(%d) : no data\n",i);
			 */
			continue;
		}
		
		for (comp=1; comp<ncomp; comp++)
		{
			mwSize cur = start, nav = 0;
			while ( (data[cur] <= res[i] + dt2) && (cur < ndata) )
			{
				if ( mxIsNaN(data[cur+comp*ndata]) )
				{
					res[i+comp*ntref] = NaN;
					nav = 0;
					break;
				}
				else
				{
					res[i+comp*ntref] += data[cur+comp*ndata];
					nav++;
				}
				cur++;
			}
			if ( nav )
				res[i+comp*ntref] = res[i+comp*ntref]/(double)nav;
			
			if ( cur>stop )
				stop = cur;
		}
		start = stop;
		
	}
}
