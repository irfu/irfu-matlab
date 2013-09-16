
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbDef_h Deffiniton of Isdat/Matlab interface
% maximum number of dimensions
% changing the constant is not enough, a lot of changes needs to be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbDef_h Deffiniton of Isdat/Matlab interface
% maximum number of dimensions
% changing the constant is not enough, a lot of changes needs to be
% done in Dbproto.h as well

NULL			= 0;
DbMAX_DIMS		= 32;
DbSPEC_NAME_DIM		= 16;

% Used in DbGetData() and DbQuery()
DbUNKNOWN		= -1;

% Valid for type & filter
DbUNDEF			= 0;
DbUNUSED		= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project

 DbCLUSTER		= 1;
 DbCRRES		= 2;
 DbEISCAT		= 3;
 DbFAST			= 4;
 DbFREJA		= 5;
 DbGALILEO		= 6;
 DbGEOTAIL		= 7;
 DbISEE			= 8;
 DbMAGNETOMETER		= 9;
 DbPOLAR		= 10;
 DbPROTO		= 11;
 DbRIOMETER		= 12;
 DbSIMCLU		= 13;
 DbVIKING		= 14;
 DbTEST			= 15;		% to get test data 
 DbCSDS_PP		= 16;
 DbCSDS_SP		= 17;
 DbCSDS_LOC		= 18;
 DbCASSINI		= 19;
 DbMARS96		= 20;
 DbASTRID_2             = 21;
 DbFLATFILE		= 22;
 DbCLUSTER_DSN		= 23;
 DbNUM_PROJECTS		= 24;		% count zero even if it's not used 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project groups 
 DbISTP			= 1000;
 DbGROUP_SIZE		= 100;

% sensor 
 DbSTATUS		= 100;		% instrument specific status 
 DbATTITUDE		= 101;		% instrument specific attitude 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% units 

% general units 
 DbUN_TM		= 1;
 DbUN_CORR		= 2;
 DbUN_PHYS		= 3;		% must be last of the general units 

% variable meta unit 
 DbUN_VAR_META		= 4;

% global meta unit 
 DbUN_GLOBAL_META	= 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific units (must have higher numbers than the general units) 
 DbUN_METER		= 10;		% (m) 
 DbUN_KILOMETER		= 11;		% (km) 
 DbUN_M_PER_S		= 12;		% (m/s) 
 DbUN_VOLT		= 13;		% (V) 
 DbUN_V_PER_M		= 14;		% (V/m) 
 DbUN_MV_PER_M		= 15;		% (mV/m) 
DbUN_V_PER_M_SQR_PER_HZ	= 16;		% (V/m)^2/Hz 
DbUN_MV_PER_M_SQR_PER_HZ= 17;		% (mV/m)^2/Hz 
 DbUN_AMPERE		= 18;		% (A) 
 DbUN_MICRO_AMP		= 19;		% (uA) 
 DbUN_TESLA		= 20;		% (T) 
 DbUN_NANO_TESLA	= 21;		% (nT) 
 DbUN_KELVIN		= 22;		% (K) 
 DbUN_HZ		= 23;		% (Hz) 
 DbUN_PROCENT		= 24;		% (%) 
 DbUN_DECIBELL		= 25;		% (dB) 
 DbUN_CELCIUS		= 26;		% (degC) 
 DbUN_EV		= 27;		% (eV) 
 DbUN_KILO_EV		= 28;		% (keV) 
 DbUN_DEGREES		= 29;		% (deg) 
 DbUN_KM_PER_S		= 30;		% (km/s) 
 DbUN_NT_SQR_PER_HZ	= 31;		% (nT)^2/Hz 
 DbUN_SECONDS		= 32;		% (sec) 

% quantity 
 DbQTY_AMPLITUDE	= 1;
 DbQTY_FREQUENCY	= 2;
 DbQTY_POWER		= 3;
 DbQTY_COUNTS		= 4;
 DbQTY_ENERGY		= 5;
 DbQTY_ANGLE		= 6;
 DbQTY_MASS		= 7;
 DbQTY_CHARGE		= 8;
 DbQTY_COUNTSSQUARED	= 9;
 DbQTY_TIME		= 10;
 DbQTY_VECTOR 		=77;		% koko, hack to make Staff happy 

% reduction 
 DbRED_NONE		= 1;
 DbRED_AVERAGE		= 2;
 DbRED_RESAMPLE		= 3;
 DbRED_SKIP		= 4;
 DbRED_MIN		= 5;
 DbRED_MAX		= 6;

% gap fill strategy (DbGetData() only) 
 DbGAP_NAN		= 1;
 DbGAP_INTERPOL		= 2;
 DbGAP_ZERO		= 3;
 DbGAP_MERGE		= 4;

% scale values 
 DbSCALE_LIN	 	= 1;
 DbSCALE_LOG		= 2;
 DbSCALE_IRREGULAR	= 3;

% pack values 
 DbPACK_TIMETAG		= 1;
 DbPACK_SEGMENT		= 2;
 DbPACK_FILL		= 3;

% DbGetContent() event values 
 DbEVENT_SWEEP		= 1;
 DbEVENT_CALIBRATION	= 2;
 DbEVENT_SOUNDER	= 3;

% DbQuery() mode values 
 DbQUERY_ALL		= 1;
 DbQUERY_ONLINE		= 2;

% DbGetContent(), DbGetInfo() and DbQuery() level values 
 DbLEVEL_PROJECT	= 0;
 DbLEVEL_MEMBER		= 1;
 DbLEVEL_INSTRUMENT	= 2;
 DbLEVEL_SENSOR		= 3;
 DbLEVEL_SIGNAL		= 4;
 DbLEVEL_CHANNEL	= 5;
 DbLEVEL_PARAMETER	= 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbSearch() definitions 

% quantity types 
 DbSEARCH_SPEC		= 1;
 DbSEARCH_CONST_INT	= 2;
 DbSEARCH_CONST_DOUBLE	= 3;

% Available operators: ( ) ! / * + - << >> < > == != & | && || 
% Single character operators are defined using ascii representation 
 DbSEARCH_PAREN_LEFT	= '(';		% ( 
 DbSEARCH_PAREN_RIGHT	= ')';		% ) 
 DbSEARCH_LOGIC_NOT	= '!';		% ! 
 DbSEARCH_DIVIDE	= '/';		% / 
 DbSEARCH_MULTIPLY	= '*';		% * 
 DbSEARCH_PLUS		= '+';		% + 
 DbSEARCH_MINUS		= '-';		% - 
 DbSEARCH_LESS_THAN	= '<';		% < 
 DbSEARCH_GREATER_THAN	= '>';		% > 
 DbSEARCH_BIT_AND	= '&';		% & 
 DbSEARCH_BIT_OR	= '|';		% | 

 DbSEARCH_SHIFT_LEFT	= 1;		% << 
 DbSEARCH_SHIFT_RIGHT	= 2;		% >> 
 DbSEARCH_EQUAL		= 3;		% == 
 DbSEARCH_NOT_EQUAL	= 4;		% != 
 DbSEARCH_LOGIC_AND	= 5;		% && 
 DbSEARCH_LOGIC_OR	= 6;		% || 


% Few Database ERROR CODES provided with Matlab_Isdat.Interface
 DbSUCCESS	   = 0;	 	% everything's okay
 DbBAD_DROP	   = 1;  % must be same as DbWARN_DROP
 DbBAD_TIME        = 2;  % requested time not on disc
 DbBAD_PROJECT     = 3;
 DbBAD_MEMBER      = 4;   
 DbBAD_EOF         = 5;  % must be same as DbWARN_EOF
 DbBAD_INSTRUMENT  = 6;   
 DbBAD_SENSOR      = 7;  %  mux1 or mux2 doesn't correspond to requested sensor
 DbBAD_SIGNAL      = 8;   
 DbBAD_CHANNEL     = 9;   
 DbBAD_PARAMETER   = 10;   
 DbBAD_UNITS       = 11;
 DbBAD_REDUCTION   = 12;
 DbBAD_GAPFILL     = 13;
 DbBAD_ALLOC       = 14;
 DbBAD_INTERNAL    = 15;
 DbBAD_GAP         = 16;  % must be same as DbWARN_GAP
 DbBAD_ZONE        = 17;  % requested interval is between two samples
 DbSUSPEND         = 20;  % for internal server use only
 DbBAD_REQUEST     = 21;  % bad request code
 DbBAD_VALUE       = 22;  % int parameter out of range
 DbBAD_ACCESS      = 23;  %  depending on context:
 			  %	- attempt to modify the access control
			  %	   list from other than the local hosts.
 DbBAD_LENGTH      = 24;  % Request length incorrect 
 DbNOT_IMPLEMENTED = 25;  % request not implemented for specified project 
 DbBAD_INTERVAL    = 26;  % invalid interval, eg. negative or too large
 DbBAD_MESSAGE     = 27;  % look in errorMessage string for details
 DbBAD_MEMLIMIT    = 28;  % memory limit exceeded

















