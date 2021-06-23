================================================================================
CESR CEFLIB INSTALLATION AND USER MANUAL
================================================================================

INTRODUCTION
============

The CESR CEFLIB library is a set of tools, written in C langage, which can be used 
to read efficiently Ascii files in CEF format (Cluster Exchange Format).

This delivery contains source code to generate the main CESR CEFLIB library,
and additional software interface to call it from different langages.

The delivery contains 4 sub-directories :

- LIB : general purpose tools used in all CESR softwares, handling times for example

- CEFLIB : the main C library used to parse CEF files

- IDL : source code used to write a DLM (Dynamically Loadable Module),
  which will allow to call this library from IDL

- MATLAB : contains source code to generate a dynamic library which will allow
  to call the CEFLIB from Matlab

WHAT'S NEW
==========


NEW IN VERSION 1.8.0 (2016/03/19)
---------------------------
Updated Matlab API
^^^^^^^^^^^^^^^^^^

 * Minor internal changes on Matlab API to ensure compatibility with Matlab 2016a.

 * This release is mandatory if you are using latest Mac OSX "El Capitan".

 * You will have also to download the latest Matlab 2016a version

 * And probably upgrade to Xcode 7.2.

See Building Matlab API on MacOSX El Capitan

Updated C and Python API
^^^^^^^^^^^^^^^^^^^^^^^^

As needed by Thomasz KLOS, member of CAA teem, the C and Python API have been updated to give the ability of reading only the CEF metadata.

A new function has been added in the Python interface :

ceflib.read_metadata (filename.cef)

This function parse the content of CEF header, to allow further metadata queries, without reading data records, to improve performances.

NEW IN VERSION 1.7.3 (2014/02/20)
---------------------------------

Minor changes, due to internal usage :

 * Adding a new fonction to LIB/src/LISTS.c :

    char *Join_list (t_list *, char * sep);

 * Modification of the generation process of IDL CEFLIB DLM :

    CEFLIB/IDL/bin/idl_cef.dlm

NEW IN VERSION 1.7.2 (2013/09/27)
---------------------------------

A bug was discovered by Nataliya Bourrel, when building the CEFLIB C API in a Fedora 19 computer, that embed the zlib.1.2.7 package.

Changes were necessary in the zlib interface (module ZZLIB).

 * LIB/inc/ZZLIB.h

 * LIB/src/ZZLIB.c

 * C/src/CEF.c

NEW IN VERSION 1.7.1 (2013/09/04)
---------------------------------
Bugs solved on Mac OSX

If the library was extracted on a case insensitive file system, a conflict between local and system header files occurs, and make the building process fail.

Some header files had to be renamed.

Another problem should occur while building the library.

You will need to update manually the mexopts.sh scripts delivered with MATLAB, as indicaed in Known bugs

CEFLIB Python Interface (Beta testing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It’s now possible to call the CEFLIB from Python, using numpy package to handle multi-dimensionnal arrays.

Plots can be done using matplotlib package, that gives plot capabilities similar to Matlab.

This is a preliminar release for this language, and the interface should evolve in the future to become more Pythonic.

Handling BYTE data type
^^^^^^^^^^^^^^^^^^^^^^^

The library was not able to handle variables with VALUE_TYPE = BYTE.

This bug was found by A. Masson while readind a CEF file, generated with Qtran translator from CSDS CDF Prime Parameter file.

The BYTE data type was defined in the "Cluster Exchange Format" document, but not present in the "CAA Metadata dictionnary".

The library has been updated to take in account this datatype.

Handling CEF files generated on Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The library was modified in order to allow reading of CEF files generated on Windows,
with CR + LF (0x0C + 0x0A) end of line terminators.

It’s actually not possible to read CEF files generated on old MacOS systems,
using CR (0x0C) end of line terminators.

Other
^^^^^

An older bug corrected in CEFLIB/Makefile but not yet redelivered.

NEW FEATURES IN VERSION 1.6
---------------------------

Solving bug on 64 bits architecture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The LIB/Makefile and CEFLIB/Makefile have been modified to solve a bug found when building
libraries on a 64 bits architecture.

The flag -fPIC has been added while compiling C source files, to create relocable functions

Converting ISO_TIME variables into ISO_TIME strings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function Milli_to_isotime() has been added in both IDL and MATLAB interface.

This function accepts 2 parameters :

- an ISO_TIME of ISO_TIME_RANGE variable, scalar or array (coded in milli-seconds from 1958)

- an integer in range 0..6, giving the number of fractionnal seconds digits expected 

The return value is an ISO_TIME string (or array of strings) 

NEW FEATURES IN VERSION 1.5
---------------------------

Multi-valued attributes
^^^^^^^^^^^^^^^^^^^^^^^

Variable attributes with multiple values are now allowed in the CEFLIB.

The functions *cef_gattr (global-attribute)* and *cef_vattr (variable, attribute)* have been updated,
in both IDL and MATLAB interface, and now return :

- a string value if the attribute has an unique value
- a list of string values in the case where attribute has multiple ones

Examples : ::

	cef_vattr ("b_vec_xyg_gse", "VALUE_TYPE") => "FLOAT"
	cef_vattr ("b_vec_xyz_gse", "REPRESENTATION_1") => [ "x", "y", "z" ]

The functions *cef_gattributes (global-attribute)* and *cef_vattributes (variable)* have also been update,
for both IDL and MATLAB interface:

- Previously, they were returning a list of (attribute, value), for each attribute of a given variable.

- Now, they just return the list of the attributes names 

ISO_TIME_RANGE
^^^^^^^^^^^^^^

The ISO_TIME_RANGE datatype is now handled by the CEFLIB, and data is returned as a two fields 
structure { START, STOP }

As other time-tags values, ISO_TIME_RANGE are internally coded in milli-seconds from 1958,
and can be converted into others IDL or MATLAB base times values using :

- MILLI_TO_JULIAN () : returns fractional Julian days (IDL interface)
- MILLI_TO_EPOCH () : returns CDF Epoch milli-seconds (IDL interface)
- CEF_DATE () : returns fractionnal days from 1-Jan-0000 (MATLAB interface)

Commas now allowed in string variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

String CEF variables and metadata attributes can now contain commas.

Example ::

	CATDESC = "Quality flag : 1 = Known problems, 2 = Survey data, 3 = Good, 4 = Excellent"

It's also allowed to use " in ' delimited strings, and ' in " delimited strings.

Examples ::

	ENTRY = "9th Workshop in King's College"
	ENTRY = '500 Gb 2.5" disk'

NEW FEATURES IN VERSION 1.4
---------------------------

Include files
^^^^^^^^^^^^^

The possibility to handle include files while reading CEF files has been added, 
to allow parsing files containing lines as below : ::

	INCLUDE = CL_CH_MISSION.ceh

However, it's necessary to specify one or several directories where to find
such include files.

This can be specified by using the CEFPATH environment variable 
(as for other CAA's tools).

Syntax : ::

	# ksh, bash environment
	$ export CEFPATH=/home/CAA/INCLUDE:/home/CIS/INCLUDE:/home/FGM/INCLUDE

	# csh, tcsh environment
	$ setenv CEFPATH /home/CAA/INCLUDE:/home/CIS/INCLUDE:/home/FGM/INCLUDE
	
Include files can be nested: for example a CEF header file describing the content of a dataset, 
included in each dataset data file, can use another INCLUDE directive to include mission 
or instrument metadata description.

Unfortunately, include files cannot actually be zipped files.

Handling missing data type
^^^^^^^^^^^^^^^^^^^^^^^^^^

The DOUBLE data type was not allowed previous releases of the CEFLIB.

This functionality has been added.

KNOWN LIMITATIONS
-----------------

Examples files
^^^^^^^^^^^^^^

Examples script files for IDL and MATLAB have not been updated specifically for CEFLIB V1.5.

Some of them should work whereas other will fail if they are using the updated functions :

- cef_gattributes ()
- cef_gattr ()
- cef_vattributes ()
- cef_vattr ()

Some CEF examples files used in these scripts are also delivered, whereas others are not.

Variable names
^^^^^^^^^^^^^^

Actually, the variable names are truncated before the __<dataset>.

If you desagree with this usage, you can edit the CEFLIB/src/CEF.c file 
and replace the following line ::
	   
#define	TRUNCATE

with ::

#undef	TRUNCATE

The C interface functions used to access to a given variable are case insensitive,
and consequently the IDL and MATLAB interface have the same behaviour.

GENERATION OF THE LIBRARIES
===========================

Each library comes with a generic Makefile, which needs some external variables declarations, 
like compilation options, directories where to find libraries, etc...

These declarations are given in the build script contained in the root directory.

Uncompress and extract sources
------------------------------

Choose a directory to install the software, copy the tar.gz file and extract the sources::

	$ cd <install-directory>
	$ gzip -d CESR_CEFLIB-Vx.y.tar.gz
	$ tar xvf CESR_CEFLIB-Vx.y.tar

Update the build script
-----------------------

In the install directory, you will find a ksh/bash script called build.

This script contains various definitions, depending on your environment (OS, compiler).

It has been tested for both Linux, Solaris and Mac OS (Darwin).

If your are working under other operating systems, try to make it work for your environment,
or contact us in case of problem.

In all cases, you will have to modify the INSTALL_DIR variable defined in the line ::

	INSTALL_DIR=<directory>

Where INSTALL_DIR must be the path of the directory containing the build script.

IDL specific instructions
^^^^^^^^^^^^^^^^^^^^^^^^^

In order to build the IDL interface, you will have to set the IDL_DIR environment variable,
which is done generally under Linux by a setup script file provided with IDL.

For example, if IDL is intalled in the /usr/local directory, enter one of the next commands ::

	$ source /usr/local/rsi/idl/bin/idl_setup		# csh or tcsh shell
	$ . /usr/local/rsi/idl/bin/idl_setup.ksh		# ksh or bash shell

**Bug with IDL 6.0**

	There is a problem between this IDL version and the usual zlib library, used to read
	compressed CEF files (.gz).

	You will have to build the CEFLIB IDL version by linking with the IDL MGPEG library
	which is implementing its own version of zlib.

	Please refer to the instructions given in the <install-dir>/IDL/Makefile, 
	to select the expected zlib library, by commenting or not a makefile declaration.

MATLAB specific instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to build the Matlab interface, you will also have to modify the following line ::

	MATLAB_DIR=<matlab_root_directory>

The Makefile in the MATLAB directory has to know the Matlab root directory, to be able to 
find some C include files located in the <matlab_root_directory>/extern/include/ directory.

Building the libraries
----------------------

You can build the whole libraries, but you will need to have both IDL and
MATLAB installed previously on your system.

In the other case, you can just build the LIB and CEFLIB libraries, 
and then choose to build either the IDL or MATLAB interface.

Recursive build of the 4 directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To generate the whole librairies, enter the following commands::

	$ cd <install-directory>
	$ ./build clean raz all

This command will recursively generate binaries for LIB CEFLIB IDL and MATLAB subdirectories.

Separate build of libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can generate only the libraries for your environment, starting from LIB, then CEFLIB, 
and finishing with the IDL or MATLAB interface ::

	$ cd <install-dir>
	$ cd LIB
	$ ../build clean raz all	// or ../build -B under linux
	$ cd ..
	$ cd CEFLIB
	$ ../build clean raz all	// or ../build -B under Linux
	$ cd ..

Then, you can choose to build the IDL interface ::

	$ cd IDL 
	$ ../build clean raz all	// or ../build -B under Linux

Or the Matlab interface ::

	$ cd MATLAB
	$ ../build clean raz all	// or ../build -B under Linux


IDL USER MANUAL
===============

How to use the CEFLIB from IDL
------------------------------

First, check the building results ::

	$ cd <install-dir>
	$ cd IDL

After binaries generation, you must obtain two files:

- bin/idl_cef.dlm : Descriptive file of DLM module
- bin/idl_cef.so : Shared library

In order to make this library available to IDL, these 2 files must be :

- in the current directory before running IDL

- in the directory where IDL search for it's dynamic libraries ::

	/usr/local/rsi/idl/bin/bin.Linux.x86	; in our computer

- in a alternate directory, specified using the Unix IDL_DLM_PATH environment variable::

	$ export IDL_DLM_PATH=<install-dir>/IDL/bin

The IDL/profile script file gives you the syntax for bash environment.

You can obtain the IDL's dynamic librairies directory by entering::

	$ idl
	IDL> print, !DLM_PATH
	/usr/local/rsi/idl_5.4/bin/bin.linux.x86

Then, you can test the library with 2 small IDL scripts by entering ::

	$ cd pro
	$ idl fgm.pro		; plot FGM data
	$ idl cis.pro		; plot CIS density and velocity
	$ idl cis_whisper.pro	; plot CIS H+, O+ and WHISPER density
	$ idl
	IDL>info, "<cef-file>"	; print all metadata and attributes

You can also use separate commands like ::

	IDL> nrec = CEF_READ (<filename>)
	IDL> print, CEF_VARNAMES()
	IDL> print, CEF_METANAMES()
	IDL> print, CEF_META ("file_name")
	IDL> print, CEF_GATTRIBUTES()

Library reference	
-----------------

At this time, the IDL interface to the CEFLIB library implements the following functions :

*nrec = CEF_READ (<filename>)*

	reads a CEF file, loads all data and metadata in memory.

	returns the number of records in the file

*CEF_CLOSE*

	will release all memory used to handle the file and its data

*CEF_VERBOSITY, <level>*

	will change the verbosity level of the CEFLIB, which can be usefull to debug
	problems while reading CEF files

	Defaut <level> is 1, and can be increased up to 5

*CEF_METANAMES ()*

	returns the list of metadata sections

*CEF_META (<section>)*

	returns the list of string entries for this metadata section

*CEF_GATTRIBUTES ()*

	returns the list of global attributes names

*CEF_GATTR (<attribute>)*

	returns value of a global attribute:
	- result should be a string if the attribute has an unique value
	- or a list of strings for multi-valued attributes

*CEF_VARNAMES ()*

	return the list of CEF variables found in the file

*CEF_VAR (<cef-var-name>)*

	creates an IDL variable, with the content of a given CEF-variable

*CEF_VATTRIBUTES (<cef-varname>)*

	Returns the list of attributes names for a given variable

*CEF_VATTR (<cef-varname>, <attribute>)*

	returns attribute value for a given CEF variable attribute:
	- result should be a string if the attribute has an unique value
	- or a list of strings for multi-valued attributes

*CEF_DEPENDS (<cef-varname>)*

	Returns a list of DEPEND_i variables names for a given CEF variable

	The first one corresponds to the DEPEND_0, then DEPEND_1, DEPEND_2, ...

	Values can be empty strings if DEPEND_i attribute is not set. 

*MILLI_TO_JULIAN (<time-tag-variable>)*

	Time-tags are internally coded in milli-seconds from 1958 in the CEFLIB

	This function converts these time-tags into corresponding Julian days

	Input can be a scalar or an array, of ISO_TIME or ISO_TIME_RANGE data

*MILLI_TO_EPOCH (<time-tag-variable>)*

	Time-tags are internally coded in milli-seconds from 1958 in the CEFLIB

	This function converts these time-tags into corresponding CEF Epoch time

	Input can be a scalar or an array, of ISO_TIME or ISO_TIME_RANGE data

MATLAB CEFLIB USER MANUAL
=========================

How to use the MATLAB CEFLIB interface
--------------------------------------

By looking at the <install-dir>/MATLAB directory you will find :

- lib/libcef.so		dynamic library used to access to the CEFLIB
- lib/libcef.h		include file necessary for MATLAB to access to the binaries
- lib/cef_init.m	small Matlab script file to load the library and facilitate functions call

You will have to update the lib/cef_init.m file, and modify the following line ::

	install_dir = <directory>

Set this variable to the path of the directory which contains libcef.so and libcef.h

You will also have to modify the Matlab path, to make this script callable from Matlab.

You can set this path by using the "File>Set path" menu, or in the Matlab command line with addpath()


Library reference
-----------------

*cef_init ()*

	This function must be called before the other ones, in order to :

	- load the shared library libcef.so,
	- create function handlers to facilitate the call to the functions of this library.

*cef_read ('cef-filename')*

	This function will read the given CEF file, load all data and metadata in memory 

	Return an error code (O : OK, > 100 : error)

*cef_close ()*

	Release memory allocated for the CEF variables and metadata.

*cef_verbosity (level)*

	This function should be used to increase the level of verbosity of the CEFLIB,
	in order to debug why the CEFLIB is not able to read a given CEF file. 

	Default value for level is 1, and could be increased up to 2, 3, 4, 5.

	CAUTION : as Matlab is not using the usual stdout stream when launched with its
	Java interface, you will have to launch it in console mode to saw the CEFLIB trace messages.
	   
	Use the following command lines to perform such debug ::
	   
		$ matlab -nojvm
		$ cef_init ()
		$ cef_read ('cef-filename')

*cef_metanames ()*

	Returns a list of strings (Matlab cell-array), giving the metadata section names.

*cef_meta ('metadata-section')*

	Returns the list of string entries for this metadata section

*cef_gattributes ()*

	Returns the list of global attributes names

*cef_gattr ('attribute')*

	Returns the value of a global attribute:

	- result should be a string if the attribute has an unique value
	- or a list of strings for multi-valued attributes

*cef_vattributes ('varname')*

	Returns the list of attributes names for a given variable

*cef_vattr ('varname', 'attribute')*

	Returns the value of a given CEF variable attribute:

	- result should be a string if the attribute has an unique value
	- or a list of strings for multi-valued attributes

*cef_varnames ()*

	Returns the list of CEF variables found in the file

*cef_var ('varname')*

	Reads the content of a given CEF variable, and returns :

	- a Matlab cell-array for string variables (CEF CHAR)
	- a n-D matrix for numerical variables (CEF ISO_TIME, FLOAT, INTEGER)

*cef_depends ('varname')*

	Returns a list of DEPEND_i variable names for a given CEF variable

	The first one correspond to the DEPEND_0, then DEPEND_1, DEPEND_2, ...

	Values can be empty strings if DEPEND_i attribute is not set. 

*cef_date ()*

	Time-tags are internally coded in milli-seconds from 1958 in the CEFLIB,
	but you can use this function to translate them in fractional days used by Matlab.
	Input can be a scalar or an array

	e.g.:  tt = cef_date (cef_var ('time_tags'))


Examples of use
---------------

You will find in the MATLAB/samples directory some Matlab scripts and CEF files,
to give you an idea of how to use this library.


C LANGUAGE USER MANUAL
======================

It is of course possible to call the CEFLIB directy from C language, but this part of the documentation
is not yet available.

Please contact us if you need some help.

CONTACTS
========

If you need some help while using or generating this library, please contact :
   
alain.barthe@cesr.fr
   
Feedbacks are of course welcome.
