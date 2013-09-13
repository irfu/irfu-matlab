#include "matrix.h"

int		cef_read (char * filename);

int		cef_close (void);

void		cef_verbosity (int level);

mxArray *	cef_metanames (void);

mxArray *	cef_meta (char * meta);

mxArray *	cef_gattributes (void);

mxArray *	cef_vattributes (char * varname);

mxArray *	cef_gattr (char * attribute);

mxArray *	cef_vattr (char * varname, char * attribute);

mxArray *	cef_var (char * varname);

mxArray *	cef_varnames (void);

mxArray *	cef_depends (char * varname);

mxArray *	milli_to_isotime (mxArray * var, int digits);
