/* $Id$ */
/*=================================================================
 * cefprint_mx.c 
 * MEX function to write data part of CEF files 
 * and GZIP the resulting file
 * 
 * Usage:
 *   STATUS = cefprint_mx(FILENAME,DATA)
 *
 *   STATUS = 0 means everything went OK
 *=============================================================*/
#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* converts ISDAT epoch to ISO time string */
void epoch2iso(double *epoch, char *str)
{
	time_t sec;
	int ms;
	struct tm *t;

	sec = floor(*epoch);
	ms = round((*epoch - (double)sec)*1000000);
	
	t = gmtime(&sec);
	
	sprintf(str, "%.4d-%.2d-%.2d%c%.2d:%.2d:%.2d.%.6d%c",
		t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, 'T',
		t->tm_hour, t->tm_min, t->tm_sec, ms, 'Z');
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	char *f_name;
	char tmp_s[BUFSIZ], tmp_s1[64], buf[BUFSIZ*64], unix_c[BUFSIZ];
	double *data, *res;
	int   buflen,status,t_mrows,t_ncols,d_mrows,d_ncols,i,j;
	FILE *fp;
	
	/* check for proper number of arguments */
	if(nrhs!=2) 
		mexErrMsgTxt("Two inputs required.");
	else if(nlhs > 2) 
		mexErrMsgTxt("Too many output arguments.");
	if(nlhs!=1)
		mexErrMsgTxt("One output required.");
	
	/* check input*/
	if ( mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("First input must be a string.");
	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("Input must be a row vector.");
	if ( mxIsDouble(prhs[1])!= 1)
		mexErrMsgTxt("Second input must be a double.");
	
	/* filename */
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	f_name = mxCalloc(buflen, sizeof(char));
	status = mxGetString(prhs[0], f_name, buflen);
	if(status != 0) 
		mexWarnMsgTxt("Not enough space. String is truncated.");
	/* printf("Filename : %s\n",f_name); */
    
	/* data */
	d_mrows = mxGetM(prhs[1]);
	d_ncols = mxGetN(prhs[1]);
	/* printf("Data     : %dx%d\n",d_mrows,d_ncols); */
	if ( d_ncols < 2 )
		mexErrMsgTxt("Input must have at least two columns.");
	data = mxGetPr(prhs[1]);
	
	/*  set the output pointer to the output matrix */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/*  create a C pointer to a copy of the output matrix */
	res = mxGetPr(plhs[0]);
	*res = 0;
	
	if ((fp = fopen(f_name,"a")) == NULL) {
		mexWarnMsgTxt("Cannot open output file");
		*res = 1;
	} else {
		setbuf(fp,buf);
		status = 0;

		for (i=0; i<d_mrows; i++){
			epoch2iso(data+i, tmp_s);
			sprintf(tmp_s,"%s, %8.3f",tmp_s,*(data +d_mrows +i));
			
			if ( d_ncols > 2)
			for ( j=2; j<d_ncols; j++ ){
				sprintf(tmp_s1,", %8.3f",*(data +d_mrows*j +i));
				strcat(tmp_s,tmp_s1);
			}
			if ( status=fprintf(fp,"%s $\n",tmp_s) < 0 ) break;
		}
		if ( status < 0 ){
			mexWarnMsgTxt("Error writing to output file");
			*res = 1;
		}
		fclose(fp);
	}
	
	strcpy(unix_c,"/usr/bin/gzip ");
	strcat(unix_c,f_name);
	if ( system(unix_c) < 0)
	{
		mexWarnMsgTxt("Error gzipping output file");
		*res = 1;
	}
	return;
}

