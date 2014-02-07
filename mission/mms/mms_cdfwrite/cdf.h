/******************************************************************************
* Copyright 1996-2013 United States Government as represented by the
* Administrator of the National Aeronautics and Space Administration.
* All Rights Reserved.
******************************************************************************/
/******************************************************************************
*
*  NSSDC/CDF				CDF Header file for C/C++ applications.
*
*  Version 3.5d, 14-Dec-97, Hughes STX.
*
*  Modification history:
*
*   V1.0  22-Jan-91, R Kulkarni	Original version (for CDF V2.0).
*		     J Love
*   V2.0   3-Jun-91, J Love     Modified for CDF V2.1 enhancements,
*				namely the INTERNAL interface and the
*				MULTI/SINGLE file option.  Added
*				macros to replace C i/f functions.
*   V2.1  20-Jun-91, J Love	Added function prototypes.
*   V2.2   8-Aug-91, J Love	Increment for CDF V2.1.2.  Use
*				'CDFlib'.  Renamed a few items.
*   V3.0  19-May-92, J Love	IBM PC & HP-UX port.  CDF V2.2.
*   V3.1  23-Sep-92, J Love	CDF V2.3 (shareable/NeXT/zVar).
*   V3.1a  5-Oct-92, J Love	CDF V2.3.0a (NeXT/encoding).
*   V3.1b  6-Oct-92, J Love	CDF V2.3.0b (CDFcompare).
*   V3.1c 27-Oct-92, J Love	CDF V2.3.0c (pad values).
*   V3.2  12-Jan-94, J Love	CDF V2.4.
*   V3.2a  4-Feb-94, J Love	DEC Alpha/OpenVMS port.
*   V3.2b 22-Feb-94, J Love	Spelling lesson.
*   V3.3   8-Dec-94, J Love	CDF V2.5.
*   V3.3a  3-Mar-95, J Love	Solaris 2.3 IDL i/f.
*   V3.4  28-Mar-95, J Love	POSIX.
*   V3.4a  8-May-95, J Love	ILLEGAL_EPOCH_VALUE.
*   V3.4b  9-Jun-95, J Love	EPOCH custom format.
*   V3.4c 20-Jul-95, J Love	CDFexport-related changes.
*   V3.5  12-Sep-96, J Love	CDF V2.6.
*   V3.5a 21-Feb-97, J Love	Removed RICE.
*   V3.5b  9-Mar-97, J Love	Windows NT for MS Visual C++ 4.0 on an IBM PC.
*   V3.5c  2-Sep-97, J Love	`__STDC__' not defined for all AIX compilers.
*   V3.5d 14-Dec-97, J Love	Added ALPHAVMSi encoding.
*   V3.6  08-Apr-04, M Liu      Added  new data type CDF_EPOCH16 and some
*                               epoch functions related to the new type.
*   V3.7  28-Apr-09, M Liu      Modified MAC_ENCODING/DECODEING to PPC_ENCODING
*                               /DECODING as Mac and Linux can run on PPC box.
*   V3.8  10-Dec-10, M Liu      Added encodeEPOCH4, encodeEPOCH16_4,
*                               parseEPOCH4, parseEPOCH16_4 to handle epochs
*                               conforming to ISO 8601.
*   V3.9  04-Apr-11, M Liu      Added a few new functions for TT2000 epoch.
*   V3.10 03-Jan-12, M Liu      Added CDFgetzVarAllRecordsByVarID,
*                               CDFgetzVarRangeRecordsByVarID,
*                               CDFgetVarAllRecordsByVarName,
*                               CDFgetVarRangeRecordsByVarName functions, and
*                               a set of similar functions for put operations.
*                               Added a new error message.
*   V3.11 29-May-12, M Liu      Added CDFinsertVarRecordsByVarID and
*                               CDFinsertVarRecorsByVarName.
******************************************************************************/

#if !defined(CDFh_INCLUDEd__)
#define CDFh_INCLUDEd__

/******************************************************************************
* CDF defined types
******************************************************************************/

typedef void *CDFid;
typedef long CDFstatus;

/******************************************************************************
* CDF defined variables
******************************************************************************/

static double *TT2000NULL = 0;

/******************************************************************************
* Limits
******************************************************************************/

#define CDF_MIN_DIMS    0               /* Min number of dimensions a CDF
					   variable may have */
#define CDF_MAX_DIMS    10              /* Max number of dimensions a CDF
					   variable may have */

/******************************************************************************
* Lengths
******************************************************************************/

#define CDF_VAR_NAME_LEN        64
#define CDF_ATTR_NAME_LEN       64

#define CDF_VAR_NAME_LEN256     256
#define CDF_ATTR_NAME_LEN256    256

#define CDF_COPYRIGHT_LEN       256
#define CDF_STATUSTEXT_LEN      120
#define CDF_PATHNAME_LEN        512

#define EPOCH_STRING_LEN	24
#define EPOCH1_STRING_LEN	16
#define EPOCH2_STRING_LEN	14
#define EPOCH3_STRING_LEN	24
#define EPOCH4_STRING_LEN	23

#define EPOCH16_STRING_LEN      36
#define EPOCH16_1_STRING_LEN    24
#define EPOCH16_2_STRING_LEN    14
#define EPOCH16_3_STRING_LEN    36
#define EPOCH16_4_STRING_LEN    32

#define TT2000_0_STRING_LEN     30
#define TT2000_1_STRING_LEN     19
#define TT2000_2_STRING_LEN     14
#define TT2000_3_STRING_LEN     29

#define EPOCHx_STRING_MAX	50
#define EPOCHx_FORMAT_MAX	68

/******************************************************************************
* Data types.
******************************************************************************/

#define CDF_INT1		1L
#define CDF_INT2		2L
#define CDF_INT4		4L
#define CDF_INT8		8L
#define CDF_UINT1		11L
#define CDF_UINT2		12L
#define CDF_UINT4		14L
#define CDF_REAL4		21L
#define CDF_REAL8		22L
#define CDF_EPOCH		31L	/* Standard style. */
#define CDF_EPOCH16		32L	/* Extended style. */
#define CDF_TIME_TT2000		33L	/* One more style with leap seconds
					   and J2000 base time. */
#define CDF_BYTE		41L     /* same as CDF_INT1 (signed) */
#define CDF_FLOAT		44L     /* same as CDF_REAL4 */
#define CDF_DOUBLE		45L     /* same as CDF_REAL8 */
#define CDF_CHAR		51L     /* a "string" data type */
#define CDF_UCHAR		52L     /* a "string" data type */

/******************************************************************************
* Encoding (for data only, everything else is network encoding).
******************************************************************************/

#define NETWORK_ENCODING        1L
#define SUN_ENCODING            2L
#define VAX_ENCODING            3L
#define DECSTATION_ENCODING     4L
#define SGi_ENCODING            5L
#define IBMPC_ENCODING          6L
#define IBMRS_ENCODING          7L
#define HOST_ENCODING           8L
#define PPC_ENCODING            9L
#define HP_ENCODING             11L
#define NeXT_ENCODING           12L
#define ALPHAOSF1_ENCODING      13L
#define ALPHAVMSd_ENCODING      14L
#define ALPHAVMSg_ENCODING      15L
#define ALPHAVMSi_ENCODING	16L

/******************************************************************************
* Decodings.
******************************************************************************/

#define NETWORK_DECODING        NETWORK_ENCODING
#define SUN_DECODING            SUN_ENCODING
#define VAX_DECODING            VAX_ENCODING
#define DECSTATION_DECODING     DECSTATION_ENCODING
#define SGi_DECODING            SGi_ENCODING
#define IBMPC_DECODING          IBMPC_ENCODING
#define IBMRS_DECODING          IBMRS_ENCODING
#define HOST_DECODING           HOST_ENCODING
#define PPC_DECODING            PPC_ENCODING
#define MAC_ENCODING            PPC_ENCODING
#define MAC_DECODING            PPC_ENCODING
#define HP_DECODING             HP_ENCODING
#define NeXT_DECODING           NeXT_ENCODING
#define ALPHAOSF1_DECODING      ALPHAOSF1_ENCODING
#define ALPHAVMSd_DECODING      ALPHAVMSd_ENCODING
#define ALPHAVMSg_DECODING      ALPHAVMSg_ENCODING
#define ALPHAVMSi_DECODING	ALPHAVMSi_ENCODING

/******************************************************************************
* Variance flags
******************************************************************************/

#define VARY   (-1L)        /* TRUE record or dimension variance flag */
#define NOVARY 0L           /* FALSE record or dimension variance flag */

/******************************************************************************
* Majorities
******************************************************************************/

#define ROW_MAJOR       1L
#define COLUMN_MAJOR    2L

/******************************************************************************
* Formats.
******************************************************************************/

#define SINGLE_FILE     1L
#define MULTI_FILE      2L

/******************************************************************************
* Checksum
******************************************************************************/

#define NO_CHECKSUM     0L
#define MD5_CHECKSUM    1L
#define OTHER_CHECKSUM  2L

/******************************************************************************
* Attribute scopes
******************************************************************************/

#define GLOBAL_SCOPE            1L
#define VARIABLE_SCOPE          2L

/******************************************************************************
* Readonly modes.
******************************************************************************/

#define READONLYon              (-1L)
#define READONLYoff             0L

/******************************************************************************
* Validate data modes.
******************************************************************************/

#define VALIDATEFILEon          (-1L)
#define VALIDATEFILEoff         0L

/******************************************************************************
* zModes.
******************************************************************************/

#define zMODEoff                0L
#define zMODEon1                1L
#define zMODEon2                2L

/******************************************************************************
* Negative to positive floating point zero modes.
******************************************************************************/

#define NEGtoPOSfp0on           (-1L)
#define NEGtoPOSfp0off          0L

/******************************************************************************
* Backward file mode. 
******************************************************************************/

#define BACKWARDFILEon          1
#define BACKWARDFILEoff         0

/******************************************************************************
* Compression/sparseness constants.
******************************************************************************/

#define CDF_MAX_PARMS			5
#define NO_COMPRESSION			0L
#define RLE_COMPRESSION			1L
#define HUFF_COMPRESSION		2L
#define AHUFF_COMPRESSION		3L
/**************************************************
* Compression `4' used to be RICE.  Do not reuse! *
**************************************************/
#define GZIP_COMPRESSION		5L

#define RLE_OF_ZEROs			0L
#define OPTIMAL_ENCODING_TREES		0L
#define NO_SPARSEARRAYS			0L
#define NO_SPARSERECORDS		0L
#define PAD_SPARSERECORDS		1L
#define PREV_SPARSERECORDS		2L

/*****************************************************************************
* Invalid/reserved constants.
*****************************************************************************/

#define RESERVED_CDFID      ((CDFid) NULL)      /* Indicates that a CDF hasn't
						   been selected yet. */
#define RESERVED_CDFSTATUS  ((CDFstatus) (-1))  /* Indicates that a CDFstatus
						   hasn't been selected yet. */

#define ILLEGAL_EPOCH_VALUE	(-1.0)
#define ILLEGAL_TT2000_VALUE    (-9223372036854775805LL)
#define FILLED_TT2000_VALUE	(-9223372036854775807LL-1)

/******************************************************************************
* Status codes (CDFstatus)
*  - informatory codes are greater than CDF_OK
******************************************************************************/

#define VIRTUAL_RECORD_DATA             ((CDFstatus) 1001)
#define DID_NOT_COMPRESS		((CDFstatus) 1002)
#define VAR_ALREADY_CLOSED              ((CDFstatus) 1003)
#define SINGLE_FILE_FORMAT              ((CDFstatus) 1004)
#define NO_PADVALUE_SPECIFIED           ((CDFstatus) 1005)
#define NO_VARS_IN_CDF                  ((CDFstatus) 1006)
#define MULTI_FILE_FORMAT		((CDFstatus) 1007)
#define SOME_ALREADY_ALLOCATED		((CDFstatus) 1008)
#define PRECEEDING_RECORDS_ALLOCATED	((CDFstatus) 1009)

#define CDF_OK                          ((CDFstatus) 0)

#define ATTR_NAME_TRUNC                 ((CDFstatus) (-1001))
#define CDF_NAME_TRUNC                  ((CDFstatus) (-1002))
#define VAR_NAME_TRUNC                  ((CDFstatus) (-1003))
#define NEGATIVE_FP_ZERO		((CDFstatus) (-1004))
					/* -1005 unused. */
#define FORCED_PARAMETER		((CDFstatus) (-1006))
#define NA_FOR_VARIABLE			((CDFstatus) (-1007))

#define CDF_WARN                        ((CDFstatus) (-2000))

#define ATTR_EXISTS                     ((CDFstatus) (-2001))
#define BAD_CDF_ID                      ((CDFstatus) (-2002))
#define BAD_DATA_TYPE                   ((CDFstatus) (-2003))
#define BAD_DIM_SIZE                    ((CDFstatus) (-2004))
#define BAD_DIM_INDEX                   ((CDFstatus) (-2005))
#define BAD_ENCODING                    ((CDFstatus) (-2006))
#define BAD_MAJORITY                    ((CDFstatus) (-2007))
#define BAD_NUM_DIMS                    ((CDFstatus) (-2008))
#define BAD_REC_NUM                     ((CDFstatus) (-2009))
#define BAD_SCOPE                       ((CDFstatus) (-2010))
#define BAD_NUM_ELEMS                   ((CDFstatus) (-2011))
#define CDF_OPEN_ERROR                  ((CDFstatus) (-2012))
#define CDF_EXISTS                      ((CDFstatus) (-2013))
#define BAD_FORMAT                      ((CDFstatus) (-2014))
#define BAD_ALLOCATE_RECS		((CDFstatus) (-2015))
#define BAD_CDF_EXTENSION		((CDFstatus) (-2016))
#define NO_SUCH_ATTR                    ((CDFstatus) (-2017))
#define NO_SUCH_ENTRY                   ((CDFstatus) (-2018))
#define NO_SUCH_VAR                     ((CDFstatus) (-2019))
#define VAR_READ_ERROR                  ((CDFstatus) (-2020))
#define VAR_WRITE_ERROR                 ((CDFstatus) (-2021))
#define BAD_ARGUMENT                    ((CDFstatus) (-2022))
#define IBM_PC_OVERFLOW                 ((CDFstatus) (-2023))
#define TOO_MANY_VARS                   ((CDFstatus) (-2024))
#define VAR_EXISTS                      ((CDFstatus) (-2025))
#define BAD_MALLOC                      ((CDFstatus) (-2026))
#define NOT_A_CDF                       ((CDFstatus) (-2027))
#define CORRUPTED_V2_CDF                ((CDFstatus) (-2028))
#define VAR_OPEN_ERROR                  ((CDFstatus) (-2029))
#define BAD_INITIAL_RECS                ((CDFstatus) (-2030))
#define BAD_BLOCKING_FACTOR             ((CDFstatus) (-2031))
#define END_OF_VAR                      ((CDFstatus) (-2032))
					/* -2033 unused. */
#define BAD_CDFSTATUS                   ((CDFstatus) (-2034))
#define CDF_INTERNAL_ERROR		((CDFstatus) (-2035))
#define BAD_NUM_VARS			((CDFstatus) (-2036))
#define BAD_REC_COUNT                   ((CDFstatus) (-2037))
#define BAD_REC_INTERVAL                ((CDFstatus) (-2038))
#define BAD_DIM_COUNT                   ((CDFstatus) (-2039))
#define BAD_DIM_INTERVAL                ((CDFstatus) (-2040))
#define BAD_VAR_NUM                     ((CDFstatus) (-2041))
#define BAD_ATTR_NUM                    ((CDFstatus) (-2042))
#define BAD_ENTRY_NUM                   ((CDFstatus) (-2043))
#define BAD_ATTR_NAME                   ((CDFstatus) (-2044))
#define BAD_VAR_NAME                    ((CDFstatus) (-2045))
#define NO_ATTR_SELECTED                ((CDFstatus) (-2046))
#define NO_ENTRY_SELECTED               ((CDFstatus) (-2047))
#define NO_VAR_SELECTED                 ((CDFstatus) (-2048))
#define BAD_CDF_NAME                    ((CDFstatus) (-2049))
					/* -2050 unused. */
#define CANNOT_CHANGE                   ((CDFstatus) (-2051))
#define NO_STATUS_SELECTED              ((CDFstatus) (-2052))
#define NO_CDF_SELECTED                 ((CDFstatus) (-2053))
#define READ_ONLY_DISTRIBUTION          ((CDFstatus) (-2054))
#define CDF_CLOSE_ERROR                 ((CDFstatus) (-2055))
#define VAR_CLOSE_ERROR                 ((CDFstatus) (-2056))
					/* -2057 unused. */
#define BAD_FNC_OR_ITEM                 ((CDFstatus) (-2058))
					/* -2059 unused. */
#define ILLEGAL_ON_V1_CDF               ((CDFstatus) (-2060))
					/* -2061 unused. */
					/* -2062 unused. */
#define BAD_CACHE_SIZE                  ((CDFstatus) (-2063))
					/* -2064 unused. */
					/* -2065 unused. */
#define CDF_CREATE_ERROR                ((CDFstatus) (-2066))
#define NO_SUCH_CDF                     ((CDFstatus) (-2067))
#define VAR_CREATE_ERROR                ((CDFstatus) (-2068))
					/* -2069 unused. */
#define READ_ONLY_MODE                  ((CDFstatus) (-2070))
#define ILLEGAL_IN_zMODE                ((CDFstatus) (-2071))
#define BAD_zMODE                       ((CDFstatus) (-2072))
#define BAD_READONLY_MODE               ((CDFstatus) (-2073))
#define CDF_READ_ERROR                  ((CDFstatus) (-2074))
#define CDF_WRITE_ERROR                 ((CDFstatus) (-2075))
#define ILLEGAL_FOR_SCOPE               ((CDFstatus) (-2076))
#define NO_MORE_ACCESS                  ((CDFstatus) (-2077))
					/* -2078 unused. */
#define BAD_DECODING		        ((CDFstatus) (-2079))
					/* -2080 unused. */
#define BAD_NEGtoPOSfp0_MODE		((CDFstatus) (-2081))
#define UNSUPPORTED_OPERATION		((CDFstatus) (-2082))
#define CDF_SAVE_ERROR			((CDFstatus) (-2083))
#define VAR_SAVE_ERROR			((CDFstatus) (-2084))
					/* -2085 unused. */
#define NO_WRITE_ACCESS                 ((CDFstatus) (-2086))
#define NO_DELETE_ACCESS                ((CDFstatus) (-2087))
#define CDF_DELETE_ERROR		((CDFstatus) (-2088))
#define VAR_DELETE_ERROR		((CDFstatus) (-2089))
#define UNKNOWN_COMPRESSION		((CDFstatus) (-2090))
#define CANNOT_COMPRESS			((CDFstatus) (-2091))
#define DECOMPRESSION_ERROR		((CDFstatus) (-2092))
#define COMPRESSION_ERROR		((CDFstatus) (-2093))
					/* -2094 unused. */
					/* -2095 unused. */
#define EMPTY_COMPRESSED_CDF		((CDFstatus) (-2096))
#define BAD_COMPRESSION_PARM		((CDFstatus) (-2097))
#define UNKNOWN_SPARSENESS		((CDFstatus) (-2098))
#define CANNOT_SPARSERECORDS		((CDFstatus) (-2099))
#define CANNOT_SPARSEARRAYS		((CDFstatus) (-2100))
#define TOO_MANY_PARMS			((CDFstatus) (-2101))
#define NO_SUCH_RECORD			((CDFstatus) (-2102))
#define CANNOT_ALLOCATE_RECORDS		((CDFstatus) (-2103))
					/* -2104 unused. */
					/* -2105 unused. */
#define SCRATCH_DELETE_ERROR		((CDFstatus) (-2106))
#define SCRATCH_CREATE_ERROR		((CDFstatus) (-2107))
#define SCRATCH_READ_ERROR		((CDFstatus) (-2108))
#define SCRATCH_WRITE_ERROR		((CDFstatus) (-2109))
#define BAD_SPARSEARRAYS_PARM		((CDFstatus) (-2110))
#define BAD_SCRATCH_DIR			((CDFstatus) (-2111))
#define NOT_A_CDF_OR_NOT_SUPPORTED      ((CDFstatus) (-2113))
#define CORRUPTED_V3_CDF                ((CDFstatus) (-2223))
#define ILLEGAL_EPOCH_FIELD             ((CDFstatus) (-2224))
#define BAD_CHECKSUM                    ((CDFstatus) (-2225)) 
#define CHECKSUM_ERROR                  ((CDFstatus) (-2226))
#define CHECKSUM_NOT_ALLOWED            ((CDFstatus) (-2227))
#define IS_A_NETCDF                     ((CDFstatus) (-2228))
#define TT2000_TIME_ERROR               ((CDFstatus) (-2229))
#define UNABLE_TO_PROCESS_CDF           ((CDFstatus) (-2230))
#define ZLIB_COMPRESS_ERROR             ((CDFstatus) (-2231))
#define ZLIB_UNCOMPRESS_ERROR           ((CDFstatus) (-2232))
#define CANNOT_INSERT_RECORDS           ((CDFstatus) (-2233))

/******************************************************************************
* Functions (for INTERNAL interface).
* NOTE: These values must be different from those of the items.
******************************************************************************/

#define CREATE_			1001L
#define OPEN_			1002L
#define DELETE_			1003L
#define CLOSE_			1004L
#define SELECT_			1005L
#define CONFIRM_		1006L
#define GET_			1007L
#define PUT_			1008L

#define SAVE_                   1009L
#define BACKWARD_               1010L
#define GETCDFFILEBACKWARD_     1011L
#define CHECKSUM_               1012L
#define GETCDFCHECKSUM_         1013L
#define VALIDATE_               1014L
#define GETCDFVALIDATE_         1015L
#define GETLEAPSECONDSENVVAR_   1016L

#define NULL_			1000L

/******************************************************************************
* Items on which functions are performed (for INTERNAL interface).
* NOTE: These values must be different from those of the functions.
******************************************************************************/

#define CDF_                    1L
#define CDF_NAME_               2L
#define CDF_ENCODING_           3L
#define CDF_DECODING_		4L
#define CDF_MAJORITY_           5L
#define CDF_FORMAT_             6L
#define CDF_COPYRIGHT_          7L
#define CDF_NUMrVARS_           8L
#define CDF_NUMzVARS_           9L
#define CDF_NUMATTRS_           10L
#define CDF_NUMgATTRS_          11L
#define CDF_NUMvATTRS_          12L
#define CDF_VERSION_            13L
#define CDF_RELEASE_            14L
#define CDF_INCREMENT_          15L
#define CDF_STATUS_             16L
#define CDF_READONLY_MODE_      17L
#define CDF_zMODE_              18L
#define CDF_NEGtoPOSfp0_MODE_	19L
#define LIB_COPYRIGHT_          20L
#define LIB_VERSION_            21L
#define LIB_RELEASE_            22L
#define LIB_INCREMENT_          23L
#define LIB_subINCREMENT_       24L
#define rVARs_NUMDIMS_          25L
#define rVARs_DIMSIZES_         26L
#define rVARs_MAXREC_           27L
#define rVARs_RECDATA_		28L
#define rVARs_RECNUMBER_        29L
#define rVARs_RECCOUNT_         30L
#define rVARs_RECINTERVAL_      31L
#define rVARs_DIMINDICES_       32L
#define rVARs_DIMCOUNTS_        33L
#define rVARs_DIMINTERVALS_     34L
#define rVAR_                   35L
#define rVAR_NAME_              36L
#define rVAR_DATATYPE_          37L
#define rVAR_NUMELEMS_          38L
#define rVAR_RECVARY_           39L
#define rVAR_DIMVARYS_          40L
#define rVAR_NUMBER_            41L
#define rVAR_DATA_              42L
#define rVAR_HYPERDATA_         43L
#define rVAR_SEQDATA_           44L
#define rVAR_SEQPOS_            45L
#define rVAR_MAXREC_            46L
#define rVAR_MAXallocREC_       47L
#define rVAR_DATASPEC_          48L
#define rVAR_PADVALUE_          49L
#define rVAR_INITIALRECS_       50L
#define rVAR_BLOCKINGFACTOR_    51L
#define rVAR_nINDEXRECORDS_	52L
#define rVAR_nINDEXENTRIES_	53L
#define rVAR_EXISTENCE_		54L
#define zVARs_MAXREC_		55L
#define zVARs_RECDATA_		56L
#define zVAR_                   57L
#define zVAR_NAME_              58L
#define zVAR_DATATYPE_          59L
#define zVAR_NUMELEMS_          60L
#define zVAR_NUMDIMS_           61L
#define zVAR_DIMSIZES_          62L
#define zVAR_RECVARY_           63L
#define zVAR_DIMVARYS_          64L
#define zVAR_NUMBER_            65L
#define zVAR_DATA_              66L
#define zVAR_HYPERDATA_         67L
#define zVAR_SEQDATA_           68L
#define zVAR_SEQPOS_            69L
#define zVAR_MAXREC_            70L
#define zVAR_MAXallocREC_       71L
#define zVAR_DATASPEC_          72L
#define zVAR_PADVALUE_          73L
#define zVAR_INITIALRECS_       74L
#define zVAR_BLOCKINGFACTOR_    75L
#define zVAR_nINDEXRECORDS_	76L
#define zVAR_nINDEXENTRIES_	77L
#define zVAR_EXISTENCE_		78L
#define zVAR_RECNUMBER_         79L
#define zVAR_RECCOUNT_          80L
#define zVAR_RECINTERVAL_       81L
#define zVAR_DIMINDICES_        82L
#define zVAR_DIMCOUNTS_         83L
#define zVAR_DIMINTERVALS_      84L
#define ATTR_                   85L
#define ATTR_SCOPE_             86L
#define ATTR_NAME_              87L
#define ATTR_NUMBER_            88L
#define ATTR_MAXgENTRY_         89L
#define ATTR_NUMgENTRIES_       90L
#define ATTR_MAXrENTRY_         91L
#define ATTR_NUMrENTRIES_       92L
#define ATTR_MAXzENTRY_         93L
#define ATTR_NUMzENTRIES_       94L
#define ATTR_EXISTENCE_		95L
#define gENTRY_                 96L
#define gENTRY_EXISTENCE_       97L
#define gENTRY_DATATYPE_        98L
#define gENTRY_NUMELEMS_        99L
#define gENTRY_DATASPEC_        100L
#define gENTRY_DATA_            101L
#define rENTRY_                 102L
#define rENTRY_NAME_		103L
#define rENTRY_EXISTENCE_       104L
#define rENTRY_DATATYPE_        105L
#define rENTRY_NUMELEMS_        106L
#define rENTRY_DATASPEC_        107L
#define rENTRY_DATA_            108L
#define zENTRY_                 109L
#define zENTRY_NAME_		110L
#define zENTRY_EXISTENCE_       111L
#define zENTRY_DATATYPE_        112L
#define zENTRY_NUMELEMS_        113L
#define zENTRY_DATASPEC_        114L
#define zENTRY_DATA_            115L
#define STATUS_TEXT_            116L
#define CDF_CACHESIZE_		117L
#define rVARs_CACHESIZE_	118L
#define zVARs_CACHESIZE_	119L
#define rVAR_CACHESIZE_		120L
#define zVAR_CACHESIZE_		121L
#define zVARs_RECNUMBER_	122L
#define rVAR_ALLOCATERECS_	123L
#define zVAR_ALLOCATERECS_	124L
#define DATATYPE_SIZE_		125L
#define CURgENTRY_EXISTENCE_	126L
#define CURrENTRY_EXISTENCE_	127L
#define CURzENTRY_EXISTENCE_	128L
#define CDF_INFO_		129L
#define CDF_COMPRESSION_	130L
#define zVAR_COMPRESSION_	131L
#define zVAR_SPARSERECORDS_	132L
#define zVAR_SPARSEARRAYS_	133L
#define zVAR_ALLOCATEBLOCK_	134L
#define zVAR_NUMRECS_		135L
#define zVAR_NUMallocRECS_	136L
#define rVAR_COMPRESSION_	137L
#define rVAR_SPARSERECORDS_	138L
#define rVAR_SPARSEARRAYS_	139L
#define rVAR_ALLOCATEBLOCK_	140L
#define rVAR_NUMRECS_		141L
#define rVAR_NUMallocRECS_	142L
#define rVAR_ALLOCATEDFROM_	143L
#define rVAR_ALLOCATEDTO_	144L
#define zVAR_ALLOCATEDFROM_	145L
#define zVAR_ALLOCATEDTO_	146L
#define zVAR_nINDEXLEVELS_	147L
#define rVAR_nINDEXLEVELS_	148L
#define CDF_SCRATCHDIR_		149L
#define rVAR_RESERVEPERCENT_	150L
#define zVAR_RESERVEPERCENT_	151L
#define rVAR_RECORDS_		152L
#define zVAR_RECORDS_		153L
#define STAGE_CACHESIZE_	154L
#define COMPRESS_CACHESIZE_	155L
#define CDF_CHECKSUM_           156L

#define CDFwithSTATS_		200L	/* For CDF internal use only! */
#define CDF_ACCESS_		201L	/* For CDF internal use only! */

#define TT2000END 		-99999.999

/******************************************************************************
* C interface macros.
******************************************************************************/

#define CDFattrCreate CDFcreateAttr
#define CDFattrNum CDFgetAttrNum
#define CDFvarCreate CDFcreaterVar
#define CDFvarNum CDFgetVarNum
#define CDFerror CDFgetStatusText
#define CDFattrRename CDFrenameAttr
#define CDFopenCDF CDFopen
#define CDFdeleteCDF CDFdelete
#define CDFcloseCDF CDFclose
#define CDFselectCDF CDFselect

#define CDFattrEntryInquire(id,attrNum,entryNum,dataType,numElems) \
CDFinquireAttrEntry(id,0,attrNum,entryNum,dataType,numElems)
#define CDFinquireAttrgEntry(id,attrNum,entryNum,dataType,numElems) \
CDFinquireAttrEntry(id,1,attrNum,entryNum,dataType,numElems)
#define CDFinquireAttrrEntry(id,attrNum,entryNum,dataType,numElems) \
CDFinquireAttrEntry(id,2,attrNum,entryNum,dataType,numElems)
#define CDFinquireAttrzEntry(id,attrNum,entryNum,dataType,numElems) \
CDFinquireAttrEntry(id,3,attrNum,entryNum,dataType,numElems)

#define CDFinquireAttr1Info(id,attrNum,attrName,attrScope,maxEntry) \
CDFinquireAttrInfo(id,0,attrNum,attrName,attrScope,maxEntry)
#define CDFinquireAttr2Info(id,attrNum,attrName,attrScope,maxEntry) \
CDFinquireAttrInfo(id,1,attrNum,attrName,attrScope,maxEntry)

#define CDFattrPut(id,attrNum,entryNum,dataType,numElems,value) \
CDFputAttrEntry(id,0,attrNum,entryNum,dataType,numElems,value)
#define CDFputAttrgEntry(id,attrNum,entryNum,dataType,numElems,value) \
CDFputAttrEntry(id,1,attrNum,entryNum,dataType,numElems,value)
#define CDFputAttrrEntry(id,attrNum,entryNum,dataType,numElems,value) \
CDFputAttrEntry(id,2,attrNum,entryNum,dataType,numElems,value)
#define CDFputAttrzEntry(id,attrNum,entryNum,dataType,numElems,value) \
CDFputAttrEntry(id,3,attrNum,entryNum,dataType,numElems,value)

#define CDFattrGet(id,attrNum,entryNum,value) \
CDFgetAttrEntry(id,0,attrNum,entryNum,value)
#define CDFgetAttrgEntry(id,attrNum,entryNum,value) \
CDFgetAttrEntry(id,1,attrNum,entryNum,value)
#define CDFgetAttrrEntry(id,attrNum,entryNum,value) \
CDFgetAttrEntry(id,2,attrNum,entryNum,value)
#define CDFgetAttrzEntry(id,attrNum,entryNum,value) \
CDFgetAttrEntry(id,3,attrNum,entryNum,value)

#define CDFgetAttrgEntryDataType(id,attrNum,entryNum,dataType) \
CDFgetAttrEntryDataType(id,1,attrNum,entryNum,dataType)
#define CDFgetAttrrEntryDataType(id,attrNum,entryNum,dataType) \
CDFgetAttrEntryDataType(id,2,attrNum,entryNum,dataType)
#define CDFgetAttrzEntryDataType(id,attrNum,entryNum,dataType) \
CDFgetAttrEntryDataType(id,3,attrNum,entryNum,dataType)

#define CDFsetAttrgEntryDataSpec(id,attrNum,entryNum,dataType) \
CDFsetAttrEntryDataSpec(id,1,attrNum,entryNum,dataType,(long)-99)
#define CDFsetAttrrEntryDataSpec(id,attrNum,entryNum,dataType) \
CDFsetAttrEntryDataSpec(id,2,attrNum,entryNum,dataType,(long)-99)
#define CDFsetAttrzEntryDataSpec(id,attrNum,entryNum,dataType) \
CDFsetAttrEntryDataSpec(id,3,attrNum,entryNum,dataType,(long)-99)

#define CDFvarRename CDFrenamerVar
#define CDFrenamerVar(id,varNum,varName) CDFrenameVar(id,0,varNum,varName)
#define CDFrenamezVar(id,varNum,varName) CDFrenameVar(id,1,varNum,varName)

#define CDFinquirerVar(id,varN,varName,dataType,numElems,numDims,dimSizes,recVary,dimVarys) \
CDFinquireVar(id,0,varN,varName,dataType,numElems,numDims,dimSizes,recVary,dimVarys)
#define CDFinquirezVar(id,varN,varName,dataType,numElems,numDims,dimSizes,recVary,dimVarys) \
CDFinquireVar(id,1,varN,varName,dataType,numElems,numDims,dimSizes,recVary,dimVarys)

#define CDFvarPut CDFputrVarData
#define CDFputrVarData(id,varNum,recNum,indices,value) \
CDFputVarData(id,0,varNum,recNum,indices,value)
#define CDFputzVarData(id,varNum,recNum,indices,value) \
CDFputVarData(id,1,varNum,recNum,indices,value)

#define CDFvarGet CDFgetrVarData
#define CDFgetrVarData(id,varNum,recNum,indices,value) \
CDFgetVarData(id,0,varNum,recNum,indices,value)
#define CDFgetzVarData(id,varNum,recNum,indices,value) \
CDFgetVarData(id,1,varNum,recNum,indices,value)

#define CDFvarHyperPut CDFhyperPutrVarData
#define CDFhyperPutrVarData(id,varNum,recS,recC,recI,indices,counts,intervals,buff) \
CDFhyperPutVarData(id,0,varNum,recS,recC,recI,indices,counts,intervals,buff)
#define CDFhyperPutzVarData(id,varNum,recS,recC,recI,indices,counts,intervals,buff) \
CDFhyperPutVarData(id,1,varNum,recS,recC,recI,indices,counts,intervals,buff)

#define CDFvarHyperGet CDFhyperGetrVarData
#define CDFhyperGetrVarData(id,varNum,recS,recC,recI,indices,counts,intervals,buff) \
CDFhyperGetVarData(id,0,varNum,recS,recC,recI,indices,counts,intervals,buff)
#define CDFhyperGetzVarData(id,varNum,recS,recC,recI,indices,counts,intervals,buff) \
CDFhyperGetVarData(id,1,varNum,recS,recC,recI,indices,counts,intervals,buff)

#define CDFvarClose CDFcloserVar
#define CDFcloserVar(id,varNum) CDFcloseVar(id,0,varNum)
#define CDFclosezVar(id,varNum) CDFcloseVar(id,1,varNum)

#define CDFdeleteAttrgEntry(id,attrNum,entryNum) \
CDFdeleteAttrEntry(id,1,attrNum,entryNum)
#define CDFdeleteAttrrEntry(id,attrNum,entryNum) \
CDFdeleteAttrEntry(id,2,attrNum,entryNum)
#define CDFdeleteAttrzEntry(id,attrNum,entryNum) \
CDFdeleteAttrEntry(id,3,attrNum,entryNum)

#define CDFgetNumAttrgEntries(id,attrNum,numEntries) \
CDFgetNumAttrEntries(id,1,attrNum,numEntries)
#define CDFgetNumAttrrEntries(id,attrNum,numEntries) \
CDFgetNumAttrEntries(id,2,attrNum,numEntries)
#define CDFgetNumAttrzEntries(id,attrNum,numEntries) \
CDFgetNumAttrEntries(id,3,attrNum,numEntries)

#define CDFgetAttrMaxgEntry(id,attrNum,maxEntry) \
CDFgetAttrMaxEntry(id,1,attrNum,maxEntry)
#define CDFgetAttrMaxrEntry(id,attrNum,maxEntry) \
CDFgetAttrMaxEntry(id,2,attrNum,maxEntry)
#define CDFgetAttrMaxzEntry(id,attrNum,maxEntry) \
CDFgetAttrMaxEntry(id,3,attrNum,maxEntry)

#define CDFgetAttrgEntryNumElements(id,attrNum,entryNum,numElems) \
CDFgetAttrEntryNumElements(id,1,attrNum,entryNum,numElems)
#define CDFgetAttrrEntryNumElements(id,attrNum,entryNum,numElems) \
CDFgetAttrEntryNumElements(id,2,attrNum,entryNum,numElems)
#define CDFgetAttrzEntryNumElements(id,attrNum,entryNum,numElems) \
CDFgetAttrEntryNumElements(id,3,attrNum,entryNum,numElems)

#define CDFgetNumrVars(id,numVars) CDFgetNumVars(id,0,numVars)
#define CDFgetNumzVars(id,numVars) CDFgetNumVars(id,1,numVars)

#define CDFdeletezVar(id,varNum) CDFdeleteVar(id,1,varNum)

#define CDFdeletezVarRecords(id,varNum,sRec,eRec) \
CDFdeleteVarRecords(id,1,varNum,sRec,eRec)

#define CDFgetzVarName(id,varNum,varName) \
CDFgetVarName(id,1,varNum,varName)

#define CDFgetzVarMaxWrittenRecNum(id,varNum,maxRec) \
CDFgetVarMaxWrittenRecNum(id,1,varNum,maxRec)

#define CDFgetzVarsMaxWrittenRecNum(id,maxRec) \
CDFgetVarsMaxWrittenRecNum(id,1,maxRec)

#define CDFgetzVarMaxAllocRecNum(id,varNum,maxAllocRec) \
CDFgetVarMaxAllocRecNum(id,1,varNum,maxAllocRec)

#define CDFgetzVarDataType(id,varNum,dataType) \
CDFgetVarDataType(id,1,varNum,dataType)

#define CDFgetzVarAllocRecords(id,varNum,allocRecs) \
CDFgetVarAllocRecords(id,1,varNum,allocRecs)
#define CDFsetzVarAllocRecords(id,varNum,allocRecs) \
CDFsetVarAllocRecords(id,1,varNum,allocRecs)

#define CDFsetzVarAllocBlockRecords(id,varNum,firstRec,lastRec) \
CDFsetVarAllocBlockRecords(id,1,varNum,firstRec,lastRec)

#define CDFgetzVarBlockingFactor(id,varNum,bf) \
CDFgetVarBlockingFactor(id,1,varNum,bf)
#define CDFsetzVarBlockingFactor(id,varNum,bf) \
CDFsetVarBlockingFactor(id,1,varNum,bf)

#define CDFgetzVarCompression(id,varNum,cType,cParms,cPct) \
CDFgetVarCompression(id,1,varNum,cType,cParms,cPct)
#define CDFsetzVarCompression(id,varNum,cType,cParms) \
CDFsetVarCompression(id,1,varNum,cType,cParms)

#define CDFsetzVarDataSpec(id,varNum,dataType) \
CDFsetVarDataSpec(id,1,varNum,dataType,(long)-99)

#define CDFsetzVarDimVariances(id,varNum,dimVarys) \
CDFsetVarDimVariances(id,1,varNum,dimVarys)

#define CDFgetzVarDimVariances(id,varNum,dimVarys) \
CDFgetVarDimVariances(id,1,varNum,dimVarys)

#define CDFgetzVarNumElements(id,varNum,numEles) \
CDFgetVarNumElements(id,1,varNum,numEles)

#define CDFgetzVarNumRecsWritten(id,varNum,numRecs) \
CDFgetVarNumRecsWritten(id,1,varNum,numRecs)

#define CDFsetzVarInitialRecs(id,varNum,initRecs) \
CDFsetVarInitialRecs(id,1,varNum,initRecs)

#define CDFgetzVarPadValue(id,varNum,pad) \
CDFgetVarPadValue(id,1,varNum,pad)
#define CDFsetzVarPadValue(id,varNum,pad) \
CDFsetVarPadValue(id,1,varNum,pad)

#define CDFgetzVarRecVariance(id,varNum,recVary) \
CDFgetVarRecVariance(id,1,varNum,recVary)
#define CDFsetzVarRecVariance(id,varNum,recVary) \
CDFsetVarRecVariance(id,1,varNum,recVary)

#define CDFgetzVarSeqData(id,varNum,data) \
CDFgetVarSeqData(id,1,varNum,data)
#define CDFputzVarSeqData(id,varNum,data) \
CDFputVarSeqData(id,1,varNum,data)

#define CDFgetzVarSparseRecords(id,varNum,sprecs) \
CDFgetVarSparseRecords(id,1,varNum,sprecs)
#define CDFsetzVarSparseRecords(id,varNum,sprecs) \
CDFsetVarSparseRecords(id,1,varNum,sprecs)

#define CDFgetzVarsRecordData(id,numVars,varNames,recNum,bufferPtr) \
CDFgetVarsRecordDatabyNames(id,1,numVars,varNames,recNum,bufferPtr)
#define CDFputzVarsRecordData(id,numVars,varNames,recNum,bufferPtr) \
CDFputVarsRecordDatabyNames(id,1,numVars,varNames,recNum,bufferPtr)

#define CDFgetzVarsRecordDatabyNumbers(id,numVars,varNumbers,recNum,buffer) \
CDFgetVarsRecordDatabyNumbers(id,1,numVars,varNumbers,recNum,buffer)
#define CDFputzVarsRecordDatabyNumbers(id,numVars,varNumbers,recNum,buffer) \
CDFputVarsRecordDatabyNumbers(id,1,numVars,varNumbers,recNum,buffer)

#define CDFgetzVarRecordData(id,varNum,recNum,buffer) \
CDFgetVarRecordData(id,1,varNum,recNum,buffer)
#define CDFputzVarRecordData(id,varNum,recNum,buffer) \
CDFputVarRecordData(id,1,varNum,recNum,buffer)

#define CDFsetzVarCacheSize(id,varNum,numBuffers) \
CDFsetVarCacheSize(id,1,varNum,numBuffers)

#define CDFsetzVarsCacheSize(id,numBuffers) \
CDFsetVarsCacheSize(id,1,numBuffers)

#define CDFgetzVarCacheSize(id,varNum,numBuffers) \
CDFgetVarCacheSize(id,1,varNum,numBuffers)

#define CDFconfirmzVarExistence(id,varName) \
CDFconfirmVarExistence(id,1,varName)

#define CDFconfirmzVarPadValueExistence(id,varNum) \
CDFconfirmVarPadValueExistence(id,1,varNum)

#define CDFgetzVarReservePercent(id,varNum,percent) \
CDFgetVarReservePercent(id,1,varNum,percent)

#define CDFsetzVarReservePercent(id,varNum,percent) \
CDFsetVarReservePercent(id,1,varNum,percent)

#define CDFgetzVarSeqPos(id,varNum,recNum,indices) \
CDFgetVarSeqPos(id,1,varNum,recNum,indices)
#define CDFsetzVarSeqPos(id,varNum,recNum,indices) \
CDFsetVarSeqPos(id,1,varNum,recNum,indices)

#define CDFgetzVarAllRecordsByVarID(id,varNum,buffer) \
CDFgetVarAllRecordsByVarID(id,1,varNum,buffer)
#define CDFputzVarAllRecordsByVarID(id,varNum,numRecs,buffer) \
CDFputVarAllRecordsByVarID(id,1,varNum,numRecs,buffer)

#define CDFgetzVarRangeRecordsByVarID(id,varNum,startRec,stopRec,buffer) \
CDFgetVarRangeRecordsByVarID(id,1,varNum,startRec,stopRec,buffer)
#define CDFputzVarRangeRecordsByVarID(id,varNum,startRec,stopRec,buffer) \
CDFputVarRangeRecordsByVarID(id,1,varNum,startRec,stopRec,buffer)

#define CDFinsertzVarRecordsByVarID(id,varNum,startRec,numRecs,buffer) \
CDFinsertVarRecordsByVarID(id,1,varNum,startRec,numRecs,buffer)
#define CDFinsertrVarRecordsByVarID(id,varNum,startRec,numRecs,buffer) \
CDFinsertVarRecordsByVarID(id,0,varNum,startRec,numRecs,buffer)

/*
 * CLOSE_  *
 *         */

#define CDFclose(id) \
CDFlib (SELECT_, CDF_, id, \
        CLOSE_, CDF_, \
        NULL_)

#define CDFcloseVar(id,zVar,varNum) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        CLOSE_, (zVar? zVAR_: rVAR_), \
        NULL_)

/*
 * CONFIRM_  *
 *           */

#define CDFconfirmAttrExistence(id,attrName) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, ATTR_EXISTENCE_, attrName, \
        NULL_)

#define CDFgetCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFgetVarCacheSize(id,zVar,varNum,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        CONFIRM_, (zVar?zVAR_CACHESIZE_:rVAR_CACHESIZE_), numBuffers, \
        NULL_)

#define CDFgetDecoding(id,decoding) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_DECODING_, decoding, \
        NULL_)

#define CDFgetName(id,cdfName) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_NAME_, cdfName, \
        NULL_)

#define CDFgetNegtoPosfp0Mode(id,negtoPosfp0) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_NEGtoPOSfp0_MODE_, negtoPosfp0, \
        NULL_)

#define CDFgetReadOnlyMode(id,readOnlyMode) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_READONLY_MODE_, readOnlyMode, \
        NULL_)

#define CDFgetzMode(id,zMode) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_zMODE_, zMode, \
        NULL_)

#define CDFgetCompressionCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, COMPRESS_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFconfirmgEntryExistence(id,attrNum,entryNum) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        CONFIRM_, gENTRY_EXISTENCE_, entryNum, \
        NULL_)

#define CDFconfirmrEntryExistence(id,attrNum,entryNum) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        CONFIRM_, rENTRY_EXISTENCE_, entryNum, \
        NULL_)

#define CDFgetStageCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, STAGE_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFconfirmVarExistence(id,zVar,varName) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, (zVar?zVAR_EXISTENCE_:rVAR_EXISTENCE_), varName, \
        NULL_)

#define CDFconfirmVarPadValueExistence(id,zVar,varNum) \
CDFlib (SELECT_, CDF_, id, \
		 (zVar?zVAR_:rVAR_), varNum, \
        CONFIRM_, (zVar?zVAR_PADVALUE_:rVAR_PADVALUE_), \
        NULL_)

#define CDFgetVarReservePercent(id,zVar,varNum,percent) \
CDFlib (SELECT_, CDF_, id, \
		 (zVar?zVAR_:rVAR_), varNum, \
        CONFIRM_, (zVar?zVAR_RESERVEPERCENT_:rVAR_RESERVEPERCENT_), percent, \
        NULL_)

#define CDFgetVarSeqPos(id,zVar,varNum,recNum,indices) \
CDFlib (SELECT_, CDF_, id, \
		 (zVar?zVAR_:rVAR_), varNum, \
        CONFIRM_, (zVar?zVAR_SEQPOS_:rVAR_SEQPOS_), recNum, indices, \
        NULL_)

#define CDFconfirmzEntryExistence(id,attrNum,entryNum) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        CONFIRM_, zENTRY_EXISTENCE_, entryNum, \
        NULL_)

#define CDFconfirmChecksum(id) \
CDFlib (SELECT_, CDF_, id, \
        CONFIRM_, CDF_CHECKSUM_, \
        NULL_)

/*
 * CREATE_ *
 *         */

#define CDFcreate(CDFname,numDims,dimSizes,encoding,majority,id) \
CDFlib (CREATE_, CDF_, CDFname, numDims, dimSizes, id, \
        PUT_, CDF_ENCODING_, encoding, \
              CDF_MAJORITY_, majority, \
        NULL_)

#define CDFcreateAttr(id,attrName,attrScope,attrNum) \
CDFlib (SELECT_, CDF_, id, \
        CREATE_, ATTR_, attrName, attrScope, attrNum, \
        NULL_)

#define CDFcreaterVar(id,varName,dataType,numElements,recVary,dimVarys,varNum) \
CDFlib (SELECT_, CDF_, id, \
        CREATE_, rVAR_, varName, dataType, numElements, \
                        recVary, dimVarys, varNum, \
        NULL_)

#define CDFcreatezVar(id,varName,dataType,numElements,numDims,dimSizes,recVary,dimVarys,varNum) \
CDFlib (SELECT_, CDF_, id, \
        CREATE_, zVAR_, varName, dataType, numElements, \
                        numDims,dimSizes, recVary, dimVarys, varNum, \
        NULL_)

/*
 * DELETE_ *
 *         */

#define CDFdelete(id) \
CDFlib (SELECT_, CDF_, id, \
        DELETE_, CDF_, \
        NULL_)

#define CDFdeleteAttr(id,attrNum) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        DELETE_, ATTR_, \
        NULL_)

#define CDFdeleteVar(id,zVar,varNum) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar? zVAR_: rVAR_), varNum, \
        DELETE_, (zVar? zVAR_: rVAR_), \
        NULL_)

#define CDFdeleteVarRecords(id,zVar,varNum,firstRec,lastRec) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar? zVAR_: rVAR_), varNum, \
        DELETE_, (zVar? zVAR_RECORDS_: rVAR_RECORDS_), firstRec, lastRec, \
        NULL_)

/*
 * GET_ *
 *      */

#define CDFgetAttrName(id,attrNum,attrName) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        GET_, ATTR_NAME_, attrName, \
        NULL_)

#define CDFgetAttrScope(id,attrNum,attrScope) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        GET_, ATTR_SCOPE_, attrScope, \
        NULL_)

#define CDFgetCompression(id,cType, cParms, cPercent) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_COMPRESSION_, cType, cParms, cPercent, \
        NULL_)

#define CDFgetCopyright(id,copyright) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_COPYRIGHT_, copyright, \
        NULL_)

#define CDFgetEncoding(id,encoding) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_ENCODING_, encoding, \
        NULL_)

#define CDFgetFormat(id,format) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_FORMAT_, format, \
        NULL_)

#define CDFgetCompressionInfo(name,cType,cParms,cSize,uSize) \
CDFlib (GET_, CDF_INFO_, name, cType, cParms, cSize, uSize, \
        NULL_)

#define CDFgetMajority(id,majority) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_MAJORITY_, majority, \
        NULL_)

#define CDFgetNumAttributes(id,numAttrs) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_NUMATTRS_, numAttrs, \
        NULL_)

#define CDFgetNumgAttributes(id,numgAttrs) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_NUMgATTRS_, numgAttrs, \
        NULL_)

#define CDFgetNumVars(id,zVar,numVars) \
CDFlib (SELECT_, CDF_, id, \
        GET_, (zVar?CDF_NUMzVARS_:CDF_NUMrVARS_), numVars, \
        NULL_)

#define CDFgetNumvAttributes(id,numvAttrs) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_NUMvATTRS_, numvAttrs, \
        NULL_)

#define CDFdoc(id,version,release,copyright) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_VERSION_, version, \
              CDF_RELEASE_, release, \
              CDF_COPYRIGHT_, copyright, \
        NULL_)

#define CDFgetDataTypeSize(dataType,numBytes) \
CDFlib (GET_, DATATYPE_SIZE_, dataType, numBytes, \
        NULL_)

#define CDFgetLibraryCopyright(copyright) \
CDFlib (GET_, LIB_COPYRIGHT_, copyright, \
        NULL_)

#define CDFgetLibraryVersion(version,release,increment,subincrement) \
CDFlib (GET_, LIB_VERSION_, version, \
              LIB_RELEASE_, release, \
              LIB_INCREMENT_, increment, \
              LIB_subINCREMENT_, subincrement, \
        NULL_)

#define CDFgetVersion(id,version,release,increment) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_VERSION_, version, \
              CDF_RELEASE_, release, \
              CDF_INCREMENT_, increment, \
        NULL_)

#define CDFgetVarBlockingFactor(id,zVar,varNum,bf) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_BLOCKINGFACTOR_:rVAR_BLOCKINGFACTOR_), bf, \
        NULL_)

#define CDFgetVarCompression(id,zVar,varNum,cType,cParms,cPct) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_COMPRESSION_:rVAR_COMPRESSION_), cType, cParms, cPct, \
        NULL_)

#define CDFgetVarData(id,zVar,varNum,recNum,indices,value) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
                 (zVar?zVAR_RECNUMBER_:rVARs_RECNUMBER_), recNum, \
                 (zVar?zVAR_DIMINDICES_:rVARs_DIMINDICES_), indices, \
        GET_, (zVar?zVAR_DATA_:rVAR_DATA_), value, \
        NULL_)

#define CDFgetVarDataType(id,zVar,varNum,dataType) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_DATATYPE_:rVAR_DATATYPE_), dataType, \
        NULL_)

#define CDFgetVarDimVariances(id,zVar,varNum,dimVarys) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_DIMVARYS_:rVAR_DIMVARYS_), dimVarys, \
        NULL_)

#define CDFgetVarMaxAllocRecNum(id,zVar,varNum,maxAllocRec) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_MAXallocREC_:rVAR_MAXallocREC_), maxAllocRec, \
        NULL_)

#define CDFgetVarMaxWrittenRecNum(id,zVar,varNum,maxRec) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_MAXREC_:rVAR_MAXREC_), maxRec, \
        NULL_)

#define CDFgetVarAllocRecords(id,zVar,varNum,numAllocRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_NUMallocRECS_:rVAR_NUMallocRECS_), numAllocRecs, \
        NULL_)

#define CDFgetVarNumElements(id,zVar,varNum,numElements) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_NUMELEMS_:rVAR_NUMELEMS_), numElements, \
        NULL_)

#define CDFgetVarNumRecsWritten(id,zVar,varNum,numRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_NUMRECS_:rVAR_NUMRECS_), numRecs, \
        NULL_)

#define CDFgetVarPadValue(id,zVar,varNum,padValue) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_PADVALUE_:rVAR_PADVALUE_), padValue, \
        NULL_)

#define CDFgetVarRecVariance(id,zVar,varNum,recVary) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_RECVARY_:rVAR_RECVARY_), recVary, \
        NULL_)

#define CDFgetVarSeqData(id,zVar,varNum,seqData) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_SEQDATA_:rVAR_SEQDATA_), seqData, \
        NULL_)

#define CDFgetVarsRecordDatabyNumbers(id,zVar,numVars,varNums,recNum,buffer) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVARs_RECNUMBER_:rVARs_RECNUMBER_), recNum, \
        GET_, (zVar?zVARs_RECDATA_:rVARs_RECDATA_), numVars, varNums, buffer, \
        NULL_)

#define CDFgetVarSparseRecords(id,zVar,varNum,sparseRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_SPARSERECORDS_:rVAR_SPARSERECORDS_), sparseRecs, \
        NULL_)

#define CDFgetrVarsDimSizes(id,dimSizes) \
CDFlib (SELECT_, CDF_, id, \
        GET_, rVARs_DIMSIZES_, dimSizes, \
        NULL_)

#define CDFgetzVarDimSizes(id,varNum,dimSizes) \
CDFlib (SELECT_, CDF_, id, \
                 zVAR_, varNum, \
        GET_, zVAR_DIMSIZES_, dimSizes, \
        NULL_)

#define CDFgetVarName(id,zVar,varNum,varName) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        GET_, (zVar?zVAR_NAME_:rVAR_NAME_), varName, \
        NULL_)

#define CDFgetzVarNumDims(id,varNum,numDims) \
CDFlib (SELECT_, CDF_, id, \
                 zVAR_, varNum, \
        GET_, zVAR_NUMDIMS_, numDims, \
        NULL_)

#define CDFgetrVarsNumDims(id,numDims) \
CDFlib (SELECT_, CDF_, id, \
        GET_, rVARs_NUMDIMS_, numDims, \
        NULL_)

#define CDFgetStatusText(status,text) \
CDFlib (SELECT_, CDF_STATUS_, status, \
        GET_, STATUS_TEXT_, text, \
        NULL_)

#define CDFhyperGetVarData(id,zVar,varN,recS,recC,recI,indices,counts,intervals,buff) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varN, \
                 (zVar?zVAR_RECNUMBER_:rVARs_RECNUMBER_), recS, \
                 (zVar?zVAR_RECCOUNT_:rVARs_RECCOUNT_), recC, \
                 (zVar?zVAR_RECINTERVAL_:rVARs_RECINTERVAL_), recI, \
                 (zVar?zVAR_DIMINDICES_:rVARs_DIMINDICES_), indices, \
                 (zVar?zVAR_DIMCOUNTS_:rVARs_DIMCOUNTS_), counts, \
                 (zVar?zVAR_DIMINTERVALS_:rVARs_DIMINTERVALS_), intervals, \
        GET_, (zVar?zVAR_HYPERDATA_:rVAR_HYPERDATA_), buff, \
        NULL_)

#define CDFgetMaxWrittenRecNums(id,maxRecrVars,maxReczVars) \
CDFlib (SELECT_, CDF_, id, \
        GET_, rVARs_MAXREC_, maxRecrVars, \
              zVARs_MAXREC_, maxReczVars, \
        NULL_)

#define CDFgetVarsMaxWrittenRecNum(id,zVar,maxRecVar) \
CDFlib (SELECT_, CDF_, id, \
        GET_, (zVar?zVARs_MAXREC_:rVARs_MAXREC_), maxRecVar, \
        NULL_)

#define CDFinquireVar(id,zVar,varN,varName,dataType,numElements,numDims,dimSizes,recVary,dimVarys) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varN, \
        GET_, (zVar?zVAR_NAME_:rVAR_NAME_), varName, \
              (zVar?zVAR_DATATYPE_:rVAR_DATATYPE_), dataType, \
              (zVar?zVAR_NUMELEMS_:rVAR_NUMELEMS_), numElements, \
              (zVar?zVAR_NUMDIMS_:rVARs_NUMDIMS_), numDims, \
              (zVar?zVAR_DIMSIZES_:rVARs_DIMSIZES_), dimSizes, \
              (zVar?zVAR_RECVARY_:rVAR_RECVARY_), recVary, \
              (zVar?zVAR_DIMVARYS_:rVAR_DIMVARYS_), dimVarys, \
        NULL_)

#define CDFvarInquire(id,varN,varName,dataType,numElements,recVary,dimVarys) \
CDFlib (SELECT_, CDF_, id, \
                 rVAR_, varN, \
        GET_, rVAR_NAME_, varName, \
              rVAR_DATATYPE_, dataType, \
              rVAR_NUMELEMS_, numElements, \
              rVAR_RECVARY_, recVary, \
              rVAR_DIMVARYS_, dimVarys, \
        NULL_)
#define CDFinquire(id,numDims,dimSizes,encoding,majority,maxRec,nVars,nAttrs) \
CDFlib (SELECT_, CDF_, id, \
        GET_, rVARs_NUMDIMS_, numDims, \
              rVARs_DIMSIZES_, dimSizes, \
              CDF_ENCODING_, encoding, \
              CDF_MAJORITY_, majority, \
              rVARs_MAXREC_, maxRec, \
              CDF_NUMrVARS_, nVars, \
              CDF_NUMATTRS_, nAttrs, \
        NULL_)
#define CDFinquireCDF(id,numDims,dimSizes,encoding,majority,maxrRec,nrVars,maxzRec,nzVars,nAttrs) \
CDFlib (SELECT_, CDF_, id, \
        GET_, rVARs_NUMDIMS_, numDims, \
              rVARs_DIMSIZES_, dimSizes, \
              CDF_ENCODING_, encoding, \
              CDF_MAJORITY_, majority, \
              rVARs_MAXREC_, maxrRec, \
              CDF_NUMrVARS_, nrVars, \
              zVARs_MAXREC_, maxzRec, \
              CDF_NUMzVARS_, nzVars, \
              CDF_NUMATTRS_, nAttrs, \
        NULL_)

#define CDFgetChecksum(id,checksum) \
CDFlib (SELECT_, CDF_, id, \
        GET_, CDF_CHECKSUM_, checksum, \
        NULL_)

/*
 * OPEN_ *
 *       */

#define CDFopen(CDFname,id) \
CDFlib (OPEN_, CDF_, CDFname, id, \
        NULL_)

/*
 * PUT_ *
 *      */

#define CDFsetAttrScope(id,attrNum,attrScope) \
CDFlib (SELECT_, CDF_, id, \
                 ATTR_, attrNum, \
        PUT_, ATTR_SCOPE_, attrScope, \
        NULL_)

#define CDFsetCompression(id,cType, cParms) \
CDFlib (SELECT_, CDF_, id, \
        PUT_, CDF_COMPRESSION_, cType, cParms, \
        NULL_)

#define CDFsetEncoding(id,encoding) \
CDFlib (SELECT_, CDF_, id, \
        PUT_, CDF_ENCODING_, encoding, \
        NULL_)

#define CDFsetFormat(id,format) \
CDFlib (SELECT_, CDF_, id, \
        PUT_, CDF_FORMAT_, format, \
        NULL_)

#define CDFsetMajority(id,majority) \
CDFlib (SELECT_, CDF_, id, \
        PUT_, CDF_MAJORITY_, majority, \
        NULL_)

#define CDFrenameAttr(id,attrNum,attrName) \
CDFlib (SELECT_, CDF_, id, \
		 ATTR_, attrNum, \
	PUT_, ATTR_NAME_, attrName, \
	NULL_)

#define CDFrenameVar(id,zVar,varNum,varName) \
CDFlib (SELECT_, CDF_, id, \
		 (zVar?zVAR_:rVAR_), varNum, \
	PUT_, (zVar?zVAR_NAME_:rVAR_NAME_), varName, \
	NULL_)

#define CDFsetVarAllocRecords(id,zVar,varNum,allocRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_ALLOCATERECS_:rVAR_ALLOCATERECS_), allocRecs, \
        NULL_)

#define CDFsetVarAllocBlockRecords(id,zVar,varNum,firstRec, lastRec) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_ALLOCATEBLOCK_:rVAR_ALLOCATEBLOCK_), firstRec, \
                                                              lastRec, \
        NULL_)

#define CDFsetVarBlockingFactor(id,zVar,varNum,bf) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_BLOCKINGFACTOR_:rVAR_BLOCKINGFACTOR_), bf, \
        NULL_)

#define CDFsetVarCompression(id,zVar,varNum,cType,cParms) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_COMPRESSION_:rVAR_COMPRESSION_), cType, cParms, \
        NULL_)

#define CDFsetVarDataSpec(id,zVar,varNum,dataType,numElems) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_DATASPEC_:rVAR_DATASPEC_), dataType, numElems, \
        NULL_)

#define CDFsetVarDimVariances(id,zVar,varNum,dimVarys) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_DIMVARYS_:rVAR_DIMVARYS_), dimVarys, \
        NULL_)

#define CDFsetVarInitialRecs(id,zVar,varNum,numRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_INITIALRECS_:rVAR_INITIALRECS_), numRecs, \
        NULL_)

#define CDFsetVarPadValue(id,zVar,varNum,padValue) \
CDFlib (SELECT_, CDF_, id, \
                (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_PADVALUE_:rVAR_PADVALUE_), padValue, \
        NULL_)

#define CDFsetVarRecVariance(id,zVar,varNum,recVary) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_RECVARY_:rVAR_RECVARY_), recVary, \
        NULL_)

#define CDFputVarSeqData(id,zVar,varNum,seqData) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_SEQDATA_:rVAR_SEQDATA_), seqData, \
        NULL_)

#define CDFsetVarSparseRecords(id,zVar,varNum,sparseRecs) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
        PUT_, (zVar?zVAR_SPARSERECORDS_:rVAR_SPARSERECORDS_), sparseRecs, \
        NULL_)

#define CDFputVarData(id,zVar,varNum,recNum,indices,value) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
                 (zVar?zVAR_RECNUMBER_:rVARs_RECNUMBER_), recNum, \
                 (zVar?zVAR_DIMINDICES_:rVARs_DIMINDICES_), indices, \
        PUT_, (zVar?zVAR_DATA_:rVAR_DATA_), value, \
        NULL_)

#define CDFputVarsRecordDatabyNumbers(id,zVar,numVars,varNums,recNum,buffer) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVARs_RECNUMBER_:rVARs_RECNUMBER_), recNum, \
        PUT_, (zVar?zVARs_RECDATA_:rVARs_RECDATA_), numVars, varNums, buffer, \
        NULL_)

#define CDFhyperPutVarData(id,zVar,varN,recS,recC,recI,indices,counts,intervals,buff) \
CDFlib (SELECT_, CDF_, id, \
		 (zVar?zVAR_:rVAR_), varN, \
		 (zVar?zVAR_RECNUMBER_:rVARs_RECNUMBER_), recS, \
		 (zVar?zVAR_RECCOUNT_:rVARs_RECCOUNT_), recC, \
		 (zVar?zVAR_RECINTERVAL_:rVARs_RECINTERVAL_), recI, \
		 (zVar?zVAR_DIMINDICES_:rVARs_DIMINDICES_), indices, \
		 (zVar?zVAR_DIMCOUNTS_:rVARs_DIMCOUNTS_), counts, \
		 (zVar?zVAR_DIMINTERVALS_:rVARs_DIMINTERVALS_), intervals, \
	PUT_, (zVar?zVAR_HYPERDATA_:rVAR_HYPERDATA_), buff, \
	NULL_)

#define CDFsetChecksum(id,checksum) \
CDFlib (SELECT_, CDF_, id, \
        PUT_, CDF_CHECKSUM_, checksum, \
        NULL_)

/*
 * SELECT_ *
 *         */

#define CDFselect(id) \
CDFlib (SELECT_, CDF_, id, \
        NULL_)

#define CDFsetDecoding(id,decoding) \
CDFlib (SELECT_, CDF_, id, \
                 CDF_DECODING_, decoding, \
        NULL_)

#define CDFsetCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 CDF_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFsetVarCacheSize(id,zVar,varNum,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
                 (zVar?zVAR_CACHESIZE_:rVAR_CACHESIZE_), numBuffers, \
        NULL_)

#define CDFsetVarsCacheSize(id,zVar,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVARs_CACHESIZE_:rVARs_CACHESIZE_), numBuffers, \
        NULL_)

#define CDFsetVarSeqPos(id,zVar,varNum,recNum,indices) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
                 (zVar?zVAR_SEQPOS_:rVAR_SEQPOS_), recNum, indices, \
        NULL_)

#define CDFsetNegtoPosfp0Mode(id,negtoPosfp0) \
CDFlib (SELECT_, CDF_, id, \
                 CDF_NEGtoPOSfp0_MODE_, negtoPosfp0, \
        NULL_)

#define CDFsetReadOnlyMode(id,readOnlyMode) \
CDFlib (SELECT_, CDF_, id, \
                 CDF_READONLY_MODE_, readOnlyMode, \
        NULL_)

#define CDFsetVarReservePercent(id,zVar,varNum,percent) \
CDFlib (SELECT_, CDF_, id, \
                 (zVar?zVAR_:rVAR_), varNum, \
                 (zVar?zVAR_RESERVEPERCENT_:rVAR_RESERVEPERCENT_), percent, \
        NULL_)

#define CDFsetCompressionCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 COMPRESS_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFsetStageCacheSize(id,numBuffers) \
CDFlib (SELECT_, CDF_, id, \
                 STAGE_CACHESIZE_, numBuffers, \
        NULL_)

#define CDFsetzMode(id,zMode) \
CDFlib (SELECT_, CDF_, id, \
                 CDF_zMODE_, zMode, \
        NULL_)

/******************************************************************************
* TT2000 macros define'd
******************************************************************************/

#define parseTT2000 CDF_TT2000_from_UTC_string
#define encodeTT2000 CDF_TT2000_to_UTC_string
#define computeTT2000 CDF_TT2000_from_UTC_parts
#define breakdownTT2000 CDF_TT2000_to_UTC_parts
#define TT2000breakdown CDF_TT2000_to_UTC_parts

/******************************************************************************
* Function prototypes.
*     It is assumed that `__cplusplus' is defined for ALL C++ compilers.  If
* ANSI function prototypes are not desired (for whatever reason), define
* noPROTOs on the compile command line.  Otherwise, ANSI function prototypes
* will be used where appropriate.
******************************************************************************/

#if !defined(noPROTOs)
#  if defined(__STDC__)
#    define PROTOs_
#  else
#    if defined(vms)
#      define PROTOs_
#    endif
#    if defined(__MSDOS__) || defined(MSDOS)
#      define PROTOs_
#    endif
#    if defined(macintosh) || defined(THINK_C)
#      define PROTOs_
#    endif
#    if defined(WIN32)
#      define PROTOs_
#    endif
#    if defined(AIX)
#      define PROTOs_
#    endif
#  endif
#endif

#if defined(PROTOs_)
#  define PROTOARGs(args) args
#else
#  define PROTOARGs(args) ()
#endif

#if defined(BUILDINGforIDL)
#  define STATICforIDL static
#  define VISIBLE_PREFIX static
#else
#  if defined(WIN32) && defined(BUILDINGforDLL)
#    if defined(LIBCDF_SOURCE_)
#      define VISIBLE_PREFIX _declspec(dllexport)
#    else
#      define VISIBLE_PREFIX _declspec(dllimport)
#    endif
#  else
#    define VISIBLE_PREFIX \

#  endif
#  define STATICforIDL \

#endif

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(BUILDINGforIDL)
/* Isn't a prototype needed? */
#else
#if !defined(__CFM68K__) || defined(__USING_STATIC_LIBS__) || !defined(CFM68KDLL)
VISIBLE_PREFIX CDFstatus CDFlib PROTOARGs((long op1, ...));
#endif
#endif
VISIBLE_PREFIX CDFstatus CDFcreateCDF PROTOARGs((
  char *name, CDFid *id
));
VISIBLE_PREFIX CDFstatus CDFattrInquire PROTOARGs((
  CDFid id, long attrNum, char *attrName, long *attrScope,
  long *maxgrEntry
)); 
VISIBLE_PREFIX CDFstatus CDFinquireAttr PROTOARGs((
  CDFid id, long attrNum, char *attrName, long *attrScope,
  long *maxgEntry, long *maxrEntry, long *maxzEntry
));
VISIBLE_PREFIX CDFstatus CDFinquireAttrEntry PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, long *dataType,
  long *numElems
));
VISIBLE_PREFIX CDFstatus CDFinquireAttrInfo PROTOARGs((
  CDFid id, int zEntry, long attrNum, char *attrName, long *attrScope,
  long *maxEntry
));
VISIBLE_PREFIX CDFstatus CDFputAttrEntry PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, long dataType,
  long numElems, void *value
));
VISIBLE_PREFIX CDFstatus CDFgetAttrEntry PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, void *value
));
VISIBLE_PREFIX CDFstatus CDFdeleteAttrEntry PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum
));
VISIBLE_PREFIX CDFstatus CDFsetAttrEntryDataSpec PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, long dataType,
  long numElems
));
VISIBLE_PREFIX long CDFgetAttrNum PROTOARGs((CDFid id, char *attrName));
VISIBLE_PREFIX long CDFgetVarNum PROTOARGs((CDFid id, char *varName));
VISIBLE_PREFIX CDFstatus CDFgetNumAttrEntries PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long *numEntries
));
VISIBLE_PREFIX CDFstatus CDFgetAttrMaxEntry PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long *maxEntry
));
VISIBLE_PREFIX CDFstatus CDFgetAttrEntryDataType PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, long *dataType
));
VISIBLE_PREFIX CDFstatus CDFgetAttrEntryNumElements PROTOARGs((
  CDFid id, int grzEntry, long attrNum, long entryNum, long *numElements
));
VISIBLE_PREFIX CDFstatus CDFgetVarRecordData PROTOARGs((
  CDFid id, int zVar, long varNum, long recNum, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFputVarRecordData PROTOARGs((
  CDFid id, int zVar, long varNum, long recNum, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFgetVarsRecordDatabyNames PROTOARGs((
  CDFid id, int zVar, long numVars, char *varNames[], long recNum,
  void *buffer[]
));
VISIBLE_PREFIX CDFstatus CDFputVarsRecordDatabyNames PROTOARGs((
  CDFid id, int zVar, long numVars, char *varNames[], long recNum,
  void *buffer[]
));
VISIBLE_PREFIX CDFstatus CDFgetVarAllRecordsByVarID PROTOARGs((
  CDFid id, int zVar, long varNum, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFputVarAllRecordsByVarID PROTOARGs((
  CDFid id, int zVar, long varNum, long numRec, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFgetVarRangeRecordsByVarID PROTOARGs((
  CDFid id, int zVar, long varNum, long startRec, long stopRec,
  void *buffer
));
VISIBLE_PREFIX CDFstatus CDFputVarRangeRecordsByVarID PROTOARGs((
  CDFid id, int zVar, long varNum, long startRec, long stopRec,
  void *buffer
));
VISIBLE_PREFIX CDFstatus CDFgetVarAllRecordsByVarName PROTOARGs((
  CDFid id, char *varName, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFputVarAllRecordsByVarName PROTOARGs((
  CDFid id, char *varName, long numRecs, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFgetVarRangeRecordsByVarName PROTOARGs((
  CDFid id, char *varName, long startRec, long stopRec, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFputVarRangeRecordsByVarName PROTOARGs((
  CDFid id, char *varName, long startRec,long stopRec, void *buffer
));
VISIBLE_PREFIX CDFstatus CDFinsertVarRecordsByVarID PROTOARGs((
  CDFid id, int zVar, long varNum, long startRec,long numRecs,
  void *buffer
));
VISIBLE_PREFIX CDFstatus CDFinsertVarRecordsByVarName PROTOARGs((
  CDFid id, char *varName, long startRec,long numRecs, void *buffer
));
VISIBLE_PREFIX void CDFsetFileBackward PROTOARGs((
  int flag
));
VISIBLE_PREFIX void CDFsetFileBackwardFlag PROTOARGs((
  int flag
));
VISIBLE_PREFIX int CDFgetFileBackward PROTOARGs(());
VISIBLE_PREFIX int CDFgetFileBackwardFlag PROTOARGs(());
VISIBLE_PREFIX void CDFsetChecksumMode PROTOARGs((
  long flag
));
VISIBLE_PREFIX long CDFgetChecksumMode PROTOARGs(());
VISIBLE_PREFIX int CDFgetFileBackwardEnvVar PROTOARGs(());
VISIBLE_PREFIX void CDFsetValidate PROTOARGs((long mode));
VISIBLE_PREFIX int CDFgetValidate PROTOARGs(());
VISIBLE_PREFIX int CDFgetValidateDebug PROTOARGs(());
#if !defined(__CFM68K__) || defined(__USING_STATIC_LIBS__) || !defined(CFM68KDLL)
VISIBLE_PREFIX void EPOCHbreakdown PROTOARGs((
  double epoch, long *year, long *month, long *day, long *hour, long *minute,
  long *second, long *msec
));
VISIBLE_PREFIX double computeEPOCH PROTOARGs((
  long year, long month, long day, long hour, long minute, long second,
  long msec
));
VISIBLE_PREFIX double parseEPOCH PROTOARGs((char *inString));
VISIBLE_PREFIX double parseEPOCH1 PROTOARGs((char *inString));
VISIBLE_PREFIX double parseEPOCH2 PROTOARGs((char *inString));
VISIBLE_PREFIX double parseEPOCH3 PROTOARGs((char *inString));
VISIBLE_PREFIX double parseEPOCH4 PROTOARGs((char *inString));
VISIBLE_PREFIX void encodeEPOCH PROTOARGs((
  double epoch, char epString[EPOCH_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH1 PROTOARGs((
  double epoch, char epString[EPOCH1_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH2 PROTOARGs((
  double epoch, char epString[EPOCH2_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH3 PROTOARGs((
  double epoch, char epString[EPOCH3_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH4 PROTOARGs((
  double epoch, char epString[EPOCH4_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCHx PROTOARGs((
  double epoch, char format[EPOCHx_FORMAT_MAX],
  char encoded[EPOCHx_STRING_MAX]
));
VISIBLE_PREFIX void EPOCH16breakdown PROTOARGs((
  double *epoch, long *year, long *month, long *day, long *hour,
  long *minute, long *second, long *msec, long *usec, long *nsec, long *psec
));
VISIBLE_PREFIX double computeEPOCH16 PROTOARGs((
  long year, long month, long day, long hour, long minute, long second,
  long msec, long usec, long nsec, long psec, double *epoch
));
VISIBLE_PREFIX double parseEPOCH16 PROTOARGs((char *inString,
  double *epoch
));
VISIBLE_PREFIX double parseEPOCH16_1 PROTOARGs((char *inStringch,
  double *epoch
));
VISIBLE_PREFIX double parseEPOCH16_2 PROTOARGs((char *inStringch,
  double *epoch
));
VISIBLE_PREFIX double parseEPOCH16_3 PROTOARGs((char *inStringch,
  double *epoch
));
VISIBLE_PREFIX double parseEPOCH16_4 PROTOARGs((char *inStringch,
  double *epoch
));
VISIBLE_PREFIX void encodeEPOCH16 PROTOARGs((
  double *epoch, char epString[EPOCH16_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH16_1 PROTOARGs((
  double *epoch, char epString[EPOCH16_1_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH16_2 PROTOARGs((
  double *epoch, char epString[EPOCH16_2_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH16_3 PROTOARGs((
  double *epoch, char epString[EPOCH16_3_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH16_4 PROTOARGs((
  double *epoch, char epString[EPOCH16_4_STRING_LEN+1]
));
VISIBLE_PREFIX void encodeEPOCH16_x PROTOARGs((
  double *epoch, char format[EPOCHx_FORMAT_MAX],
  char encoded[EPOCHx_STRING_MAX]
));
/******************************************************************************
* A new set of functions to handle CDF_TIME_TT2000 time type.
******************************************************************************/

VISIBLE_PREFIX void CDF_TT2000_to_UTC_parts PROTOARGs((
  long long tt2000, double *year, double *month, double *day, ...
));
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_parts PROTOARGs((
  double year, double month, double day, ...
));
VISIBLE_PREFIX double CDF_TT2000_to_UTC_EPOCH PROTOARGs((
  long long time
));
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_EPOCH PROTOARGs((
  double epoch
));
VISIBLE_PREFIX double CDF_TT2000_to_UTC_EPOCH16 PROTOARGs((
  long long time, double *epoch16
));
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_EPOCH16 PROTOARGs((
  double *epoch16
));
VISIBLE_PREFIX void CDF_TT2000_to_UTC_string PROTOARGs((
  long long time, char *string, ...
));
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_string PROTOARGs((
  char *string
));
VISIBLE_PREFIX void CDFClearLeapSecondsTable PROTOARGs(());
VISIBLE_PREFIX void CDFfillLeapSecondsTable PROTOARGs((
  double **table
));
VISIBLE_PREFIX int CDFgetRowsinLeapSecondsTable PROTOARGs(());
#if defined(vms)
VISIBLE_PREFIX void CDFgetLastDateinLeapSecondsTBL PROTOARGs((
  long *year, long *month, long *day
));
#else
VISIBLE_PREFIX void CDFgetLastDateinLeapSecondsTable PROTOARGs((
  long *year, long *month, long *day
));
#endif
VISIBLE_PREFIX char *CDFgetLeapSecondsTableEnvVar PROTOARGs(());
VISIBLE_PREFIX int CDFgetLeapSecondsTableStatus PROTOARGs(());
#endif
#if defined(__cplusplus)
}
#endif

/******************************************************************************
* Synonyms for compatibility with older releases.
******************************************************************************/

#define CDF_DOCUMENT_LEN	        CDF_COPYRIGHT_LEN
#define CDF_ERRTEXT_LEN         	CDF_STATUSTEXT_LEN
#define CDF_NUMDIMS_            	rVARs_NUMDIMS_
#define CDF_DIMSIZES_           	rVARs_DIMSIZES_
#define CDF_MAXREC_             	rVARs_MAXREC_
#define CDF_RECNUMBER_          	rVARs_RECNUMBER_
#define CDF_RECCOUNT_           	rVARs_RECCOUNT_
#define CDF_RECINTERVAL_        	rVARs_RECINTERVAL_
#define CDF_DIMINDICES_         	rVARs_DIMINDICES_
#define CDF_DIMCOUNTS_          	rVARs_DIMCOUNTS_
#define CDF_DIMINTERVALS_       	rVARs_DIMINTERVALS_
#define CDF_NUMVARS_            	CDF_NUMrVARS_
#define VAR_                    	rVAR_
#define VAR_NAME_               	rVAR_NAME_
#define VAR_DATATYPE_           	rVAR_DATATYPE_
#define VAR_NUMELEMS_           	rVAR_NUMELEMS_
#define VAR_RECVARY_            	rVAR_RECVARY_
#define VAR_DIMVARYS_           	rVAR_DIMVARYS_
#define VAR_NUMBER_             	rVAR_NUMBER_
#define VAR_DATA_               	rVAR_DATA_
#define VAR_HYPERDATA_          	rVAR_HYPERDATA_
#define VAR_SEQDATA_            	rVAR_SEQDATA_
#define VAR_SEQPOS_             	rVAR_SEQPOS_
#define VAR_MAXREC_             	rVAR_MAXREC_
#define VAR_DATASPEC_           	rVAR_DATASPEC_
#define VAR_FILLVALUE_          	rVAR_PADVALUE_
#define VAR_INITIALRECS_        	rVAR_INITIALRECS_
#define VAR_EXTENDRECS_         	rVAR_BLOCKINGFACTOR_
#define ATTR_MAXENTRY_          	ATTR_MAXrENTRY_
#define ATTR_NUMENTRIES_        	ATTR_NUMrENTRIES_
#define ENTRY_                  	rENTRY_
#define ENTRY_DATATYPE_         	rENTRY_DATATYPE_
#define ENTRY_NUMELEMS_         	rENTRY_NUMELEMS_
#define ENTRY_DATA_             	rENTRY_DATA_
#define MIPSEL_ENCODING			DECSTATION_ENCODING
#define MIPSEB_ENCODING			SGi_ENCODING
#define rVAR_EXISTANCE_			rVAR_EXISTENCE_
#define zVAR_EXISTANCE_			zVAR_EXISTENCE_
#define ATTR_EXISTANCE_			ATTR_EXISTENCE_
#define gENTRY_EXISTANCE_		gENTRY_EXISTENCE_
#define rENTRY_EXISTANCE_		rENTRY_EXISTENCE_
#define zENTRY_EXISTANCE_		zENTRY_EXISTENCE_
#define GLOBAL_SCOPE_ASSUMED		GLOBAL_SCOPE
#define VARIABLE_SCOPE_ASSUMED		VARIABLE_SCOPE
#define BAD_EXTEND_RECS			BAD_BLOCKING_FACTOR
#define rVAR_EXTENDRECS_		rVAR_BLOCKINGFACTOR_
#define zVAR_EXTENDRECS_		zVAR_BLOCKINGFACTOR_
#define COL_MAJOR			COLUMN_MAJOR
#define NONE_CHECKSUM			NO_CHECKSUM

#define StrlaststrIgCase		StrLaststrIgCase     
#define Strlaststr                      StrLaststr

/*****************************************************************************/

#endif
