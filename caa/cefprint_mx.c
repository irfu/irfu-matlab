/*
 *  $Id$
 *
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <yuri@irfu.se> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
 * ----------------------------------------------------------------------------
 *
 * 
 * MEX function to write data part of CEF files and GZIP the resulting file.
 * CEF header is written fromMatlab by caa_export.m
 * 
 * Usage:
 *   STATUS = cefprint_mx(FILENAME,DATA,FORMAT_STRING)
 *
 *   FORMAT_STRING is optional, and should be an n-column array, where n is the number of data columns in DATA
 *   FORMAT_STRING should be transposed from a n-row string array to account for the different row-column ordering in C
 *   e.g.: FORMAT_STRING=['%8.3f'; '%8.3f'; '%8.3f'; '%2.0f'; '%7.0f']';
 *
 *   STATUS = 0 means everything went OK
 *   STATUS = 1 means error on the data writing stange
 *   STATUS = 2 means error on the compression stange
 */

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

/* converts ISDAT epoch to ISO time string */
void epoch2iso(double *epoch, char *str)
{
	time_t sec;
	int ms;
	struct tm *t;

	sec = floor(*epoch);
	ms = round((*epoch - (double)sec)*1000000);

	/* Check if we round up to a whole second */
	if ( ms >= 1000000 )
	{
		ms = ms - 1000000;
		sec++;
	}
	
	t = gmtime(&sec);
	
	sprintf(str, "%.4d-%.2d-%.2d%c%.2d:%.2d:%.2d.%.6d%c",
		t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, 'T',
		t->tm_hour, t->tm_min, t->tm_sec, ms, 'Z');
}

void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[])
{
	char *f_name, *format, *formats=NULL;
	char tmp_s[BUFSIZ], tmp_s1[128], buf[BUFSIZ*64], unix_c[BUFSIZ];
	double *data, *res;
	int   buflen,status,t_mrows,t_ncols,d_mrows,d_ncols,i,j,formatlen;
	FILE *fp;
    
	/* check for proper number of arguments */
	if(nlhs!=1)
		mexErrMsgTxt("One output required.");
	if(nrhs<2)
		mexErrMsgTxt("At least two inputs required.");
	else if(nrhs > 3)
		mexErrMsgTxt("Too many input arguments.");
    
	/* check input*/
	if ( mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("First input must be a string.");
	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("First input must be a row vector.");
	if ( mxIsDouble(prhs[1])!= 1)
		mexErrMsgTxt("Second input must be a double.");
	if ( (nrhs == 3) && ( mxIsChar(prhs[2])!= 1) )
		mexErrMsgTxt("Third input must be a string.");
    
	/* data */
	d_mrows = mxGetM(prhs[1]);
	d_ncols = mxGetN(prhs[1]);
	/* printf("Data     : %dx%d\n",d_mrows,d_ncols); */
	if ( d_ncols < 2 )
		mexErrMsgTxt("Input must have at least two columns.");
	data = mxGetPr(prhs[1]);
    
	/* set up formatting string(s) */
	if ( nrhs == 3 ) {
		if (mxGetN(prhs[2])!=d_ncols-1)
			mexErrMsgTxt("Third input must have the same number of rows as the number of data columns.");
		formatlen=mxGetM(prhs[2]);
		format=mxCalloc(formatlen+8, sizeof(char));
		formats = mxArrayToString(prhs[2]);
		if(formats == NULL)
			mexErrMsgTxt("Could not convert third input to string.");
	} else {
		format = mxCalloc(16, sizeof(char));
		strcpy(format,"%s, %8.3f");
	}
    
	/* filename */
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	f_name = mxCalloc(buflen, sizeof(mxChar));
	status = mxGetString(prhs[0], f_name, buflen);
	if(status != 0)
		mexWarnMsgTxt("Not enough space. String is truncated.");
	printf("Filename : %s\n",f_name);
    
	/*  set the output pointer to the output matrix */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/*  create a C pointer to a copy of the output matrix */
	res = mxGetPr(plhs[0]);
	*res = 0;
	
	/* first remove the old file if it is there */
	strcpy(unix_c,"/bin/rm -f ");
	strcat(unix_c,f_name);
	strcat(unix_c,".gz");
	system(unix_c);
	
	if ( (fp = fopen(f_name,"a")) == NULL ) {
		mexWarnMsgTxt("Cannot open output file");
		*res = 1;
	} else {
/*		setbuf(fp,buf);*/
		setbuffer(fp,buf,BUFSIZ*64-1);
		status = 0;

		for (i=0; i<d_mrows; i++){
			epoch2iso(data+i, tmp_s);
			if ( nrhs == 3 ) {
				strcpy(format,"%s, ");
				strncat(format,formats,formatlen);                
			}
			sprintf(tmp_s,format,tmp_s,*(data +d_mrows +i));
			
			if ( d_ncols > 2)
				for ( j=2; j<d_ncols; j++ ){
					if ( nrhs == 3) {
						strcpy(format,", ");
						strncat(format,formats+formatlen*(j-1),formatlen);
						sprintf(tmp_s1,format,*(data +d_mrows*j +i));
					} else
						sprintf(tmp_s1,", %8.3f",*(data +d_mrows*j +i));
					strcat(tmp_s,tmp_s1);
				}

			if ( (status = fprintf(fp,"%s $\n",tmp_s)) < 0 )
				break;
		}
		if ( status < 0 ){
			mexWarnMsgTxt("Error writing to output file");
			*res = 1;
		}
		
		status = fprintf(fp,"END_OF_DATA\n");
		if ( status < 0 ){
			mexWarnMsgTxt("Error writing to output file");
			*res = 1;
		}
		
		fclose(fp);
	}

	/* Free allocated dynamic memory*/
	mxFree(format);
	mxFree(formats);
    
	if (*res) {
		mxFree(f_name);
		return;
	}
	
	/* gzip the output */
	if ( access("/usr/bin/gzip",X_OK) == 0 )
		strcpy(unix_c,"/usr/bin/gzip "); /* FreeBSD */
	else if ( access("/bin/gzip",X_OK) == 0 )
		strcpy(unix_c,"/bin/gzip "); /* Linux */
	else {
		mexWarnMsgTxt("gzip not found");
		mxFree(f_name);
		*res = 2;
		return;
	}

	strcat(unix_c,f_name);
	if ( system(unix_c) < 0)
	{
		mexWarnMsgTxt("Error gzipping output file");
		*res = 2;

		/* remove the corrupt .gz file if it is there */
		strcpy(unix_c,"/bin/rm -f ");
		strcat(unix_c,f_name);
		strcat(unix_c,".gz");
		system(unix_c);
	}
   
	mxFree(f_name);
	return;
}

