/*
Filename: irfu_cdfwrite_usv.c 

Purpose: To be used with MatLab to write CDF files according to specifications for MMS.
This file should only be used for writing data for quick look dce, as it make huge assumptions
on the input arguments (which data is sent where).

Sources: 
1. CDF C Ref. manual, for v.3.5
http://cdaweb.gsfc.nasa.gov/pub/software/cdf/doc/cdf350/cdf350crm.pdf
2. MMS CDF specificaitons.
3. SDP_Data_products_guide_v01_ttt.docx (dropbox)
4. FIELDSProcessing_20131025.docx

Author: T. Nilsson
Swedish Institute of Space Physics, Uppsala

Date of latest mod: 2014/01/22

Maj.Rev:
0.1 (2013/12/21) Created first outline.

Note:
Use a modern compiler or there might be issues with FillVal for CDF_TIME_TT2000 as it is a rather large number. 
Tested on Xubuntu 12.04 LTS 64bit, MatLab 2013b 64bit, gcc/g++ version 4.7.
*/



// Include needed files
#include "cdf.h"    // For CDFLib & CDF commands, is then including a lot more such as "math.h" etc..
#include "stdio.h"  // For printing messages to screen. Possibly remove this later on.
#include "string.h" // For string operations such as strlen.
#include "mex.h"    // For MatLab mex c
#include "stdint.h" // Datatype & lengths such as int64_t (TT2000).

void UserStatusHandler(int status)
{
// FIXME: This should have a proper handling of error messages, see manual but for now:
// Dont care about actual error handling at this time.. Just print status.. Most likely caused by existing file with identical name.
 mexPrintf ("Error found as: %d\n",status);
//  return;
}

// ENTRENCE for MatLab function call
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

// Input: prhs[0] - filename (Excluding extension and NO directory just the filename).
//        prhs[1] - SC_id (1, 2, 3 or 4).
//        prhs[2] - Epoch (int64_t CDF_TT2000).
//        prhs[3] - ESCP Estimated s/c potential (float, CDF_Float), variable of epoch.
//        prhs[4] - PSP Probe-to-s/c potential (float, CDF_Float), variable of epoch.
//        prhs[5] - Delta, where ESCP = - PSP*(static shortening factor) + Delta, variable of epoch.
//        prhs[6] - PSP_individual (float, M*6, CDF_Float[6]), variable of epoch.
//        prhs[7] - bitmask, variable of epoch.

if(nrhs!=8)
{
  mexErrMsgIdAndTxt("MATLAB:irfu_cdfwrite:inputarg",
     "Eight input arguments required. 'Filename', 'SC_id', 'Epoch', 'mmsX_escp', 'mmsX_psp', 'mmsX_delta', 'mmsX_probes_to_sc', 'bitmask'"); 
}

  char *filename;
  int8_t sc_id;
 
// Seventh parameters mmsX_probes_to_sc (CDF_Float[6]) are read as a single column vector in Mex, re-shape it to matrix to write with CDF lib.
  float buffer7[mxGetM(prhs[6])][mxGetN(prhs[6])];


// First input argument, filename

  if(mxGetClassID(prhs[0])==4)
  {
    filename = mxArrayToString(prhs[0]);
  } 
  else
  {
    mexErrMsgIdAndTxt("MATLAB:irfu_cdfwrite:inputarg", "First input was identified as ClassID: %s\n",mxGetClassName(prhs[0])); 
  }

// Second input argument, SC_id
  if((mxGetClassID(prhs[1])==8))
  {
    sc_id = mxGetScalar(prhs[1]);
    if((sc_id>4)||(sc_id<1))
    {
      mexErrMsgIdAndTxt("MATLAB:irfu_cdfwrite:inputarg","SC id incorrectly found as: %d.\nAllowed values are only 1, 2, 3 and 4.",sc_id);
    }
  }
  else
  {
    mexErrMsgIdAndTxt("MATLAB:irfu_cdfwrite:inputarg","Second input was identified as ClassID: %s\nShould be only be int8().",mxGetClassName(prhs[1]));
  }

// Third input argument, Epoch times

 // Debug message 
  // M = mxGetM(prhs[2]);
  // N = mxGetN(prhs[2]);
  // mexPrintf("M found as: %d and N found as %d\n",M,N);
  
  // Debug message
   // mexPrintf("Buffer3[0]: %ld\n", buffer3[0]);
   // mexPrintf("Buffer3[1]: %ld\n", buffer3[1]);

// Seventh input argument, mmsX_probes_to_sc (N*6, cdf_float[6])
  float *xValues;

  xValues = (float *)mxGetData(prhs[6]);

  // REMEMBER: Matlab begins with index 1, C & Mex begin with index 0. 
  for(int i=0; i<mxGetN(prhs[6]); i++)
  {
    for(int j=0; j<mxGetM(prhs[6]); j++)
    {
       buffer7[j][i] = *xValues++;
    }
  }
 // Debug message
  // mexPrintf("Buffer7[1][0]: %g\n", buffer7[1][0]);
  // mexPrintf("Buffer7[1][1]: %g\n", buffer7[1][1]);
  // mexPrintf("Buffer7[M-1][N-1]: %g\n", buffer7[M-1][N-1]);

// Sixth input, bitmask
  

//////////////////////////////
// Create temp. variables in C
//////////////////////////////

CDFid id;                                 /* CDF identifier. */
CDFstatus status;                         /* Returned status code. */

long EPOCHrecVary = {VARY};           /* EPOCH record variance. */
//long TIMETAGrecVary = {VARY};       /* TIMETAG record variance. */
//long SAMPLERATErecVary = {VARY};      /* SAMPLERATE record variance. */
long LABELrecVary = {NOVARY};           /* LON record variance. */
long SENSORrecVary = {VARY};             /* TMP record variance. */
long SENSORINDrecVary = {VARY};
long BITMASKrecVary = {VARY};
long EPOCHdimVarys[1] = {NOVARY};     /* EPOCH dimension variances. */
//long TIMETAGdimVarys[1] = {NOVARY}; /* TIMETAG dimension variances. */
//long SAMPLERATEdimVarys[1] = {NOVARY};    /* SAMPLERATE dimension variances. */
long LABELdimVarys[1] = {VARY};    /* LON dimension variances. */
long SENSORdimVarys[1] = {VARY};    /* TMP dimension variances. */
long SENSORINDdimVarys[1] = {VARY};
long BITMASKdimVarys[1] = {VARY};

long EPOCHvarNum;                            /* EPOCH zVariable number. */
//long TIMETAGvarNum;                          /* TIMETAG zVariable number. */
//long SAMPLERATEvarNum;                       /* SAMPLERATE zVariable number. */
long LABELvarNum;                              /* LON zVariable number. */
long SENSORvarNum;                              /* TMP zVariable number. */
long SENSORINDvarNum;
long BITMASKvarNum;
long ESCPvarNum;
long PSPvarNum;
           //FIXME: Pissibly dimSizes should be = {1} for epoch?
long EPOCHdimSizes[1] = {3};          /* EPOCH dimension sizes. */
//long TIMETAGdimSizes[1] = {3};      /* TIMETAG dimension sizes. */
//long SAMPLERATEdimSizes[1] = {3};     /* SAMPLERATE dimension sizes.*/
long LABELdimSizes[1] = {6};          /* LON dimension sizes. */
long SENSORdimSizes[1] = {1}; 
long SENSORINDdimSizes[1] = {6};            /* TMP dimension sizes. */
long BITMASKdimSizes[1] = {1};

//////////////////////////
// Create the new CDF file
//////////////////////////
status = CDFcreateCDF (filename, &id);
  if (status != CDF_OK) UserStatusHandler (status);

// Set file encoding and data majority 
long encoding; /* Encoding. */
long majority; /* Majority. */
encoding = NETWORK_ENCODING;
majority = COLUMN_MAJOR;
status = CDFsetEncoding(id, encoding);
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFsetMajority(id, majority);
  if (status != CDF_OK) UserStatusHandler (status);


///////////////////////////////////////////////
// Create needed zVariables in the new CDF file
///////////////////////////////////////////////

// Temporary variables to handle spacecraft numbering later on.
char sc[6];
char scid[2]; sprintf(scid, "%d", sc_id);
sprintf(sc, "mms%d_", sc_id); // Store variable sc as: mmsX_, were X = 1, 2, 3 or 4.
char tmp_string[250]; // Theoretical Max for CDF variable names etc is 255.



status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"), CDF_TIME_TT2000, 1, 0L, EPOCHdimSizes, EPOCHrecVary, EPOCHdimVarys, &EPOCHvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

//status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_timetag_dce"), CDF_TIME_TT2000, 1, 0L, EPOCHdimSizes, EPOCHrecVary, EPOCHdimVarys, &TIMETAGvarNum);
//  if (status != CDF_OK) UserStatusHandler (status);

/* Possibly FIXME: Possible wrong format, according to skeleton CDF_UINT4 but according to Ref2 it should only be CDF_INT2 or CDF_REAL4. */ 
//status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_samplerate_dce"), CDF_UINT4, 1, 0L, SAMPLERATEdimSizes, SAMPLERATErecVary, SAMPLERATEdimVarys, &SAMPLERATEvarNum);
//  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, "DCE_LABL_1", CDF_CHAR, 4, 1L, LABELdimSizes, LABELrecVary, LABELdimVarys, &LABELvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_escp_dcv"), CDF_REAL4, 1, 1L, SENSORdimSizes, SENSORrecVary, SENSORdimVarys, &ESCPvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_psp_dcv"), CDF_REAL4, 1, 1L, SENSORdimSizes, SENSORrecVary, SENSORdimVarys, &PSPvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_delta_dcv"), CDF_REAL4, 1, 1L, SENSORdimSizes, SENSORrecVary, SENSORdimVarys, &SENSORvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_psp_probes_dcv"), CDF_REAL4, 1, 1L, SENSORINDdimSizes, SENSORINDrecVary, SENSORINDdimVarys, &SENSORINDvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFcreatezVar (id, strcat(strcpy(tmp_string,sc),"sdp_bitmask_dcv"), CDF_UINT4, 1, 1L, BITMASKdimSizes, BITMASKrecVary, BITMASKdimVarys, &BITMASKvarNum);
  if (status != CDF_OK) UserStatusHandler (status);

  //printf("SENSOR var num found as: %ld\n",SENSORvarNum);



//////////////////////////////////////
// Create Global Attributes.
//////////////////////////////////////

long TEMPattrNum;                         /*"temporary" attribute number. */

status = CDFattrCreate (id, "Project", GLOBAL_SCOPE, &TEMPattrNum);                     //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Discipline", GLOBAL_SCOPE, &TEMPattrNum);                  //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Validity", GLOBAL_SCOPE, &TEMPattrNum);                    //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Validator", GLOBAL_SCOPE, &TEMPattrNum);                   //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Caveats", GLOBAL_SCOPE, &TEMPattrNum);                     //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Source_name", GLOBAL_SCOPE, &TEMPattrNum);                 //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Data_type", GLOBAL_SCOPE, &TEMPattrNum);                   //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Descriptor", GLOBAL_SCOPE, &TEMPattrNum);                  //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Data_version", GLOBAL_SCOPE, &TEMPattrNum);                //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "TITLE", GLOBAL_SCOPE, &TEMPattrNum);                       //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Logical_file_id", GLOBAL_SCOPE, &TEMPattrNum);             //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Logical_source", GLOBAL_SCOPE, &TEMPattrNum);              //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Logical_source_description", GLOBAL_SCOPE, &TEMPattrNum);  //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Mission_group", GLOBAL_SCOPE, &TEMPattrNum);               //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "PI_name", GLOBAL_SCOPE, &TEMPattrNum);                     //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "PI_affiliation", GLOBAL_SCOPE, &TEMPattrNum);              //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Acknowledgement", GLOBAL_SCOPE, &TEMPattrNum);             //NOT REQ, but REC.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Generated_by", GLOBAL_SCOPE, &TEMPattrNum);                //NOT REQ, but REC.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Generation_date", GLOBAL_SCOPE, &TEMPattrNum);             //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Rules_of_use", GLOBAL_SCOPE, &TEMPattrNum);                //NOT REQ, but OPT.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Skeleton_version", GLOBAL_SCOPE, &TEMPattrNum);            //NOT REQ, but OPT.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Software_version", GLOBAL_SCOPE, &TEMPattrNum);            //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Validate", GLOBAL_SCOPE, &TEMPattrNum);                    //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SC_Eng_id", GLOBAL_SCOPE, &TEMPattrNum);                   //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "File_naming_convention", GLOBAL_SCOPE, &TEMPattrNum);      //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Instrument_type", GLOBAL_SCOPE, &TEMPattrNum);             //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "LINK_TITLE", GLOBAL_SCOPE, &TEMPattrNum);                  //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "HTTP_LINK", GLOBAL_SCOPE, &TEMPattrNum);                   //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "LINK_TEXT", GLOBAL_SCOPE, &TEMPattrNum);                   //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Time_resolution", GLOBAL_SCOPE, &TEMPattrNum);             //NOT REQ, but OPT.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "TEXT", GLOBAL_SCOPE, &TEMPattrNum);                        //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "MODS", GLOBAL_SCOPE, &TEMPattrNum);                        //REQ.
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "ADID_ref", GLOBAL_SCOPE, &TEMPattrNum);                    //NOT REQ, but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Parents", GLOBAL_SCOPE, &TEMPattrNum);                     //NOT REQ, but OPT.
  if (status != CDF_OK) UserStatusHandler (status);



////////////////////////////
// Write Global Attributes. 
////////////////////////////

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Project"), 0, CDF_CHAR, strlen("STP>Solar-Terrestrial Physics"), "STP>Solar-Terrestrial Physics");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Discipline"), 0, CDF_CHAR, strlen("Space Physics>Magnetospheric Science"), "Space Physics>Magnetospheric Science");
  if (status != CDF_OK) UserStatusHandler (status);
 
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Validity"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Validator"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Caveats"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

char tmp[250]; 
sprintf(tmp,"MMS%d>MMS Satellite Number %d",sc_id,sc_id);
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Source_name"), 0, CDF_CHAR, strlen(tmp), tmp);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Data_type"), 0, CDF_CHAR, strlen("DCE>DC Double Probe Electric Field"), "DCE>DC Double Probe Electric Field");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Descriptor"), 0, CDF_CHAR, strlen("ADP-SDP>Axial Double Probe- Spin Plane Double Probe"), "ADP-SDP>Axial Double Probe- Spin Plane Double Probe");
  if (status != CDF_OK) UserStatusHandler (status);
 
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Data_version"), 0, CDF_CHAR, strlen("v.0.0.0"), "v.0.0.0");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"TITLE"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Logical_file_id"), 0, CDF_CHAR, strlen(filename), filename);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Logical_source"), 0, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler (status);

//FIXME: Should be written in words but for now use the value from skeleton.
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Logical_source_description"), 0, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_16c_l1a_dce")), strcat(strcpy(tmp_string,sc),"sdp_16c_l1a_dce"));
  if (status != CDF_OK) UserStatusHandler (status);
 
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Mission_group"), 0, CDF_CHAR, strlen("MMS"), "MMS");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"PI_name"), 0, CDF_CHAR, strlen("Burch, J., Ergun, R., Lindqvist, P."), "Burch, J., Ergun, R., Lindqvist, P.");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"PI_affiliation"), 0, CDF_CHAR, strlen("SwRI, LASP, KTH"), "SwRI, LASP, KTH");
  if (status != CDF_OK) UserStatusHandler (status);

//static char Acknowledgement[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Acknowledgement"), 0, CDF_CHAR, strlen(Acknowledgement), Acknowledgement);
//  if (status != CDF_OK) UserStatusHandler (status);

//static char Generated_by[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Generated_by"), 0, CDF_CHAR, strlen(Generated_by), Generated_by);
//  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Generation_date"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);
 
//static char Rules_of_use[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Rules_of_use"), 0, CDF_CHAR, strlen(Rules_of_use), Rules_of_use);
//  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Skeleton_version"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Software_version"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Validate"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"SC_Eng_id"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"File_naming_convention"), 0, CDF_CHAR, strlen("source_datatype_descriptor"), "source_datatype_descriptor");
  if (status != CDF_OK) UserStatusHandler (status);
 
status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Instrument_type"), 0, CDF_CHAR, strlen("Electric Fields (space)"), "Electric Fields (space)");
  if (status != CDF_OK) UserStatusHandler (status);

//static char LINK_TITLE[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"LINK_TITLE"), 0, CDF_CHAR, strlen(LINK_TITLE), LINK_TITLE);
//  if (status != CDF_OK) UserStatusHandler (status);

//static char HTTP_LINK[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"HTTP_LINK"), 0, CDF_CHAR, strlen(HTTP_LINK), HTTP_LINK);
//  if (status != CDF_OK) UserStatusHandler (status);

//static char LINK_TEXT[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"LINK_TEXT"), 0, CDF_CHAR, strlen(LINK_TEXT), LINK_TEXT);
//  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Time_resolution"), 0, CDF_CHAR, strlen("Configurable"), "Configurable");
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"TEXT"), 0, CDF_CHAR, strlen("L1A DC Electric Field"), "L1A DC Electric Field");
  if (status != CDF_OK) UserStatusHandler (status);
 
//static char MODS[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"MODS"), 0, CDF_CHAR, strlen(MODS), MODS);
//  if (status != CDF_OK) UserStatusHandler (status);

//static char ADID_ref[] = {" "};
//status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"ADID_ref"), 0, CDF_CHAR, strlen(ADID_ref), ADID_ref);
//  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputAttrgEntry (id, CDFgetAttrNum(id,"Parents"), 0, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler (status);



//////////////////////////////
// Create Variable attributes.
//////////////////////////////

status = CDFattrCreate (id, "FIELDNAM", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "VALIDMIN", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for data, support data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "VALIDMAX", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for data, support data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SCALEMIN", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SCALEMAX", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "LABLAXIS", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. (or LABL_PTR_1) [for data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "LABL_PTR_1", VARIABLE_SCOPE, &TEMPattrNum);                //REQ. (or LABLAXIS) [for data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "UNITS", VARIABLE_SCOPE, &TEMPattrNum);                     //REQ. (or UNITS_PTR) [for data, support data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "UNIT_PTR", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. (or UNITS) [for data, support data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "FORMAT", VARIABLE_SCOPE, &TEMPattrNum);                    //REQ. (or FORM_PTR) [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "FORM_PTR", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. (or FORMAT) [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "FILLVAL", VARIABLE_SCOPE, &TEMPattrNum);                   //REQ. [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "VAR_TYPE", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DICT_KEY", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SCALETYP", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "MONOTON", VARIABLE_SCOPE, &TEMPattrNum);                   //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "AVG_TYPE", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "CATDESC", VARIABLE_SCOPE, &TEMPattrNum);                   //REQ. [for data, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DELTA_PLUS_VAR", VARIABLE_SCOPE, &TEMPattrNum);            //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DELTA_MINUS_VAR", VARIABLE_SCOPE, &TEMPattrNum);           //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DEPEND_0", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for datai, support data, meta data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DEPEND_1", VARIABLE_SCOPE, &TEMPattrNum);                  //REQ. [for data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Calib_software", VARIABLE_SCOPE, &TEMPattrNum);            //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Calib_input", VARIABLE_SCOPE, &TEMPattrNum);               //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Frame", VARIABLE_SCOPE, &TEMPattrNum);                     //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SI_conversion", VARIABLE_SCOPE, &TEMPattrNum);             //REQ. [for data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SI_conversion_ptr", VARIABLE_SCOPE, &TEMPattrNum);         //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SC_id", VARIABLE_SCOPE, &TEMPattrNum);                     //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "Sig_digits", VARIABLE_SCOPE, &TEMPattrNum);                //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "DISPLAY_TYPE", VARIABLE_SCOPE, &TEMPattrNum);              //REQ. [for data]
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "VAR_NOTES", VARIABLE_SCOPE, &TEMPattrNum);                 //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);
status = CDFattrCreate (id, "SCAL_PTR", VARIABLE_SCOPE, &TEMPattrNum);                  //NOT REQ. but listed in skeleton
  if (status != CDF_OK) UserStatusHandler (status);



//////////////////////////////
// Write Variable attributes.
//////////////////////////////

/////
// First write variable attributes to mms1_epoch_variable
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), EPOCHvarNum, CDF_CHAR, strlen("Time tags"), "Time tags");
  if (status != CDF_OK) UserStatusHandler(status);

signed long long ValidMin_1[1] = { -431358160000000 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), EPOCHvarNum, CDF_TIME_TT2000, 1, ValidMin_1);
  if (status != CDF_OK) UserStatusHandler(status);

signed long long ValidMax_1[1] = { 946728067183999999 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), EPOCHvarNum, CDF_TIME_TT2000, 1, ValidMax_1);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), EPOCHvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

signed long long int Fillval_1[1] = { -9223372036854775808 }; // NOTE: Use a modern C compiler..
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), EPOCHvarNum, CDF_TIME_TT2000, 1, Fillval_1);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), EPOCHvarNum, CDF_CHAR, strlen("support_data"), "support_data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DICT_KEY"), EPOCHvarNum, CDF_CHAR, strlen("time>TT2000"), "time>TT2000");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), EPOCHvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"MONOTON"), EPOCHvarNum, CDF_CHAR, strlen("INCREASE"), "INCREASE");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), EPOCHvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_software"), EPOCHvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_input"), EPOCHvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Frame"), EPOCHvarNum, CDF_CHAR, strlen("scalar>na"), "scalar>na");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SI_conversion"), EPOCHvarNum, CDF_CHAR, strlen("1.0e-3>s"), "1.0e-3>s");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SC_id"), EPOCHvarNum, CDF_CHAR, strlen(scid), scid);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Sig_digits"), EPOCHvarNum, CDF_CHAR, strlen("14"), "14");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DISPLAY_TYPE"), EPOCHvarNum, CDF_CHAR, strlen("time_series"), "time_series");
  if (status != CDF_OK) UserStatusHandler(status);


/////
// Second write variable attributes to mms1_timestamp_variable
/////
/*
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), TIMETAGvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_timetag_dcv")), strcat(strcpy(tmp_string,sc),"sdp_timetag_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

signed long long ValidMin_2[1] = { -431358160000000 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), TIMETAGvarNum, CDF_TIME_TT2000, 1, ValidMin_2);
  if (status != CDF_OK) UserStatusHandler(status);

signed long long ValidMax_2[1] = { 946728067183999999 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), TIMETAGvarNum, CDF_TIME_TT2000, 1, ValidMax_2);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), TIMETAGvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_timetag_dcv")), strcat(strcpy(tmp_string,sc),"sdp_timetag_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), TIMETAGvarNum, CDF_CHAR, strlen("ns"), "ns");
  if (status != CDF_OK) UserStatusHandler(status);

signed long long int Fillval_2[1] = { -9223372036854775808 }; // NOTE: This value will give some issues if not using a modern C compiler..
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), TIMETAGvarNum, CDF_TIME_TT2000, 1, Fillval_2);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), TIMETAGvarNum, CDF_CHAR, strlen("support_data"), "support_data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), TIMETAGvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"MONOTON"), TIMETAGvarNum, CDF_CHAR, strlen("INCREASE"), "INCREASE");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), TIMETAGvarNum, CDF_CHAR, strlen("packet time tag"), "packet time tag");
  if (status != CDF_OK) UserStatusHandler(status);

/////
// Third write variable attributes to mms1_sdp_samplerate_variable
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), SAMPLERATEvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_samplerate_dce")), strcat(strcpy(tmp_string,sc),"sdp_samplerate_dce"));
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int ValidMin_3[1] = { 1 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), SAMPLERATEvarNum, CDF_UINT4, 1, ValidMin_3);
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int ValidMax_3[1] = { 262144 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), SAMPLERATEvarNum, CDF_UINT4, 1, ValidMax_3);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), SAMPLERATEvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_samplerate_dce")), strcat(strcpy(tmp_string,sc),"sdp_samplerate_dce"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), SAMPLERATEvarNum, CDF_CHAR, strlen("samples per second"), "samples per second");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), SAMPLERATEvarNum, CDF_CHAR, strlen("I7"), "I7");
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int Fillval_3[1] = { 4294967295 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), SAMPLERATEvarNum, CDF_UINT4, 1, Fillval_3);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), SAMPLERATEvarNum, CDF_CHAR, strlen("support_data"), "support_data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), SAMPLERATEvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), SAMPLERATEvarNum, CDF_CHAR, strlen("sample rate"), "sample rate");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), SAMPLERATEvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_timetag_dce")), strcat(strcpy(tmp_string,sc),"sdp_timetag_dce"));
  if (status != CDF_OK) UserStatusHandler(status);
*/
/////
// Fourth (actually second) write variable attributes to DCE_LABL_1_variable
/////

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), LABELvarNum, CDF_CHAR, strlen("DCE_LABL_1"), "DCE_LABL_1");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), LABELvarNum, CDF_CHAR, strlen("A23"), "A23");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), LABELvarNum, CDF_CHAR, strlen("meta_data"), "meta_data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), LABELvarNum, CDF_CHAR, strlen("DCE_LABL_1"), "DCE_LABL_1");
  if (status != CDF_OK) UserStatusHandler(status);

/////
// Fifth (actually third) write variable attributes to mmsX_escp_dcv (named spd_escp_dcv
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), ESCPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_escp_dcv")), strcat(strcpy(tmp_string,sc),"sdp_escp_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMin_5[1] = { -500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), ESCPvarNum, CDF_REAL4, 1, ValidMin_5);
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMax_5[1] = { 500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), ESCPvarNum, CDF_REAL4, 1, ValidMax_5);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), ESCPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_escp_dcv")), strcat(strcpy(tmp_string,sc),"sdp_escp_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), ESCPvarNum, CDF_CHAR, strlen("V"), "V");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), ESCPvarNum, CDF_CHAR, strlen("F8.3"), "F8.3");
  if (status != CDF_OK) UserStatusHandler(status);

float Fillval_5[1] = { -1.0e+31 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), ESCPvarNum, CDF_REAL4, 1, Fillval_5);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), ESCPvarNum, CDF_CHAR, strlen("data"), "data");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: SCALAR? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DICT_KEY"), ESCPvarNum, CDF_CHAR, strlen("dc_electric_field>vector"), "dc_electric_field>vector");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), ESCPvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"AVG_TYPE"), ESCPvarNum, CDF_CHAR, strlen("standard"), "standard");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), ESCPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), ESCPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_software"), ESCPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_input"), ESCPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: Scalar? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Frame"), ESCPvarNum, CDF_CHAR, strlen("vector>xyz"), "vector>xyz");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SI_conversion"), ESCPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SC_id"), ESCPvarNum, CDF_CHAR, strlen(scid), scid);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Sig_digits"), ESCPvarNum, CDF_CHAR, strlen("3"), "3");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DISPLAY_TYPE"), ESCPvarNum, CDF_CHAR, strlen("time_series"), "time_series");
  if (status != CDF_OK) UserStatusHandler(status);

/////
// Sixth (actually forth)  write variable attributes to mmsX_psp_dcv (named spd_psp_dcv
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), PSPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_psp_dcv")), strcat(strcpy(tmp_string,sc),"sdp_psp_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMin_6[1] = { -500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), PSPvarNum, CDF_REAL4, 1, ValidMin_6);
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMax_6[1] = { 500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), PSPvarNum, CDF_REAL4, 1, ValidMax_6);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), PSPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_psp_dcv")), strcat(strcpy(tmp_string,sc),"sdp_psp_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), PSPvarNum, CDF_CHAR, strlen("V"), "V");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), PSPvarNum, CDF_CHAR, strlen("F8.3"), "F8.3");
  if (status != CDF_OK) UserStatusHandler(status);

float Fillval_6[1] = { -1.0e+31 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), PSPvarNum, CDF_REAL4, 1, Fillval_6);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), PSPvarNum, CDF_CHAR, strlen("data"), "data");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: SCALAR? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DICT_KEY"), PSPvarNum, CDF_CHAR, strlen("dc_electric_field>vector"), "dc_electric_field>vector");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), PSPvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"AVG_TYPE"), PSPvarNum, CDF_CHAR, strlen("standard"), "standard");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), PSPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), PSPvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_software"), PSPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_input"), PSPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: Scalar? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Frame"), PSPvarNum, CDF_CHAR, strlen("vector>xyz"), "vector>xyz");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SI_conversion"), PSPvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SC_id"), PSPvarNum, CDF_CHAR, strlen(scid), scid);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Sig_digits"), PSPvarNum, CDF_CHAR, strlen("3"), "3");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DISPLAY_TYPE"), PSPvarNum, CDF_CHAR, strlen("time_series"), "time_series");
  if (status != CDF_OK) UserStatusHandler(status);


/////
// Seventh (actually fifth)  write variable attributes to mmsX_delta_dcv (named spd_delta_dcv
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), SENSORvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_delta_dcv")), strcat(strcpy(tmp_string,sc),"sdp_delta_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMin_7[1] = { -500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), SENSORvarNum, CDF_REAL4, 1, ValidMin_7);
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMax_7[1] = { 500000.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), SENSORvarNum, CDF_REAL4, 1, ValidMax_7);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), SENSORvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_delta_dcv")), strcat(strcpy(tmp_string,sc),"sdp_delta_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), SENSORvarNum, CDF_CHAR, strlen("V"), "V");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), SENSORvarNum, CDF_CHAR, strlen("F8.3"), "F8.3");
  if (status != CDF_OK) UserStatusHandler(status);

float Fillval_7[1] = { -1.0e+31 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), SENSORvarNum, CDF_REAL4, 1, Fillval_7);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), SENSORvarNum, CDF_CHAR, strlen("data"), "data");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: SCALAR? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DICT_KEY"), SENSORvarNum, CDF_CHAR, strlen("dc_electric_field>vector"), "dc_electric_field>vector");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), SENSORvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"AVG_TYPE"), SENSORvarNum, CDF_CHAR, strlen("standard"), "standard");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), SENSORvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), SENSORvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_software"), SENSORvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_input"), SENSORvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);
//FIXME: FIXME: Scalar? At least not vector
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Frame"), SENSORvarNum, CDF_CHAR, strlen("vector>xyz"), "vector>xyz");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SI_conversion"), SENSORvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SC_id"), SENSORvarNum, CDF_CHAR, strlen(scid), scid);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Sig_digits"), SENSORvarNum, CDF_CHAR, strlen("3"), "3");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DISPLAY_TYPE"), SENSORvarNum, CDF_CHAR, strlen("time_series"), "time_series");
  if (status != CDF_OK) UserStatusHandler(status);




/////
// Eighth (actually sixth) write variable attributes to mmsX_sdp_psp_probes_dcv_variable
/////
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), SENSORINDvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_psp_probes_dcv")), strcat(strcpy(tmp_string,sc),"sdp_psp_probes_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMin_8[1] = { -60.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), SENSORINDvarNum, CDF_REAL4, 1, ValidMin_8);
  if (status != CDF_OK) UserStatusHandler(status);

float ValidMax_8[1] = { 60.14 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), SENSORINDvarNum, CDF_REAL4, 1, ValidMax_8);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABL_PTR_1"), SENSORINDvarNum, CDF_CHAR, strlen("DCE_LABL_1"), "DCE_LABL_1");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), SENSORINDvarNum, CDF_CHAR, strlen("V"), "V");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), SENSORINDvarNum, CDF_CHAR, strlen("F8.3"), "F8.3");
  if (status != CDF_OK) UserStatusHandler(status);

float Fillval_8[1] = { -1.0e+31 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), SENSORINDvarNum, CDF_REAL4, 1, Fillval_8);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), SENSORINDvarNum, CDF_CHAR, strlen("data"), "data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DICT_KEY"), SENSORINDvarNum, CDF_CHAR, strlen("dc_electric_field>vector"), "dc_electric_field>vector");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), SENSORINDvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"AVG_TYPE"), SENSORINDvarNum, CDF_CHAR, strlen("standard"), "standard");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), SENSORINDvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), SENSORINDvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_software"), SENSORINDvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Calib_input"), SENSORINDvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Frame"), SENSORINDvarNum, CDF_CHAR, strlen("vector>xyz"), "vector>xyz");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SI_conversion"), SENSORINDvarNum, CDF_CHAR, strlen(" "), " ");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SC_id"), SENSORINDvarNum, CDF_CHAR, strlen(scid), scid);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"Sig_digits"), SENSORINDvarNum, CDF_CHAR, strlen("3"), "3");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DISPLAY_TYPE"), SENSORINDvarNum, CDF_CHAR, strlen("time_series"), "time_series");
  if (status != CDF_OK) UserStatusHandler(status);


/////
// Nineth (actually seventh) write variable attributes to mmsX_sdp_dcv_bitmask_variable (named spd_dcv_bitmask)
/////

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FIELDNAM"), BITMASKvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_bitmask_dcv")), strcat(strcpy(tmp_string,sc),"sdp_bitmask_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int ValidMin_9[1] = { 1 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMIN"), BITMASKvarNum, CDF_UINT4, 1, ValidMin_9);
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int ValidMax_9[1] = { 262144 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VALIDMAX"), BITMASKvarNum, CDF_UINT4, 1, ValidMax_9);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"LABLAXIS"), BITMASKvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_bitmask_dcv")), strcat(strcpy(tmp_string,sc),"sdp_bitmask_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"UNITS"), BITMASKvarNum, CDF_CHAR, strlen("Bitmask"), "Bitmask");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FORMAT"), BITMASKvarNum, CDF_CHAR, strlen("I7"), "I7");
  if (status != CDF_OK) UserStatusHandler(status);

unsigned int Fillval_9[1] = { 4294967295 };
status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"FILLVAL"), BITMASKvarNum, CDF_UINT4, 1, Fillval_9);
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"VAR_TYPE"), BITMASKvarNum, CDF_CHAR, strlen("support_data"), "support_data");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"SCALETYP"), BITMASKvarNum, CDF_CHAR, strlen("linear"), "linear");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"CATDESC"), BITMASKvarNum, CDF_CHAR, strlen("Bitmask"), "Bitmask");
  if (status != CDF_OK) UserStatusHandler(status);

status = CDFputAttrzEntry (id, CDFgetAttrNum(id,"DEPEND_0"), BITMASKvarNum, CDF_CHAR, strlen(strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv")), strcat(strcpy(tmp_string,sc),"sdp_epoch_dcv"));
  if (status != CDF_OK) UserStatusHandler(status);


/////////////////////////////
// Write acutal data to file.
/////////////////////////////

char *buffer[] = {"PSP1PSP2PSP3PSP4PSP5PSP6"};

status = CDFputzVarAllRecordsByVarID (id, EPOCHvarNum, mxGetM(prhs[2]), (int64_t *)mxGetData(prhs[2]));
  if (status != CDF_OK) UserStatusHandler (status);

// NOTE: It appears as if the CDFputzVarAllRecordsByVarID for one record reads the size of one Label from buffer,
//       the size of each record is 3x4 char (3 row, 4 column) as when it was decleared above.
status = CDFputzVarAllRecordsByVarID (id, LABELvarNum, 1, buffer[0]);
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputzVarAllRecordsByVarID (id, ESCPvarNum, mxGetM(prhs[3]), (float *)mxGetData(prhs[3]));
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputzVarAllRecordsByVarID (id, PSPvarNum, mxGetM(prhs[4]), (float *)mxGetData(prhs[4]));
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputzVarAllRecordsByVarID (id, SENSORvarNum, mxGetM(prhs[5]), (float *)mxGetData(prhs[5]));
  if (status != CDF_OK) UserStatusHandler (status);

status = CDFputzVarAllRecordsByVarID (id, SENSORINDvarNum, mxGetM(prhs[6]), (float *)buffer7);
  if (status != CDF_OK) UserStatusHandler (status);


status = CDFputzVarAllRecordsByVarID (id, BITMASKvarNum, mxGetM(prhs[7]), (uint32_t *)mxGetData(prhs[7]));
  if (status != CDF_OK) UserStatusHandler (status);

//////////////////////
// Close the CDF file.
//////////////////////

status = CDFcloseCDF( id);
  if (status != CDF_OK) UserStatusHandler (status);

// DONE
 // return 0;
}
