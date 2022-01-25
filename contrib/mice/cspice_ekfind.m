%-Abstract
%
%   CSPICE_EKFIND finds E-kernel data that satisfy a set of constraints.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      query    a string specifying the data to locate from data
%               available in all loaded EK files.
%
%               [1,c1] = size(query); char = class(query)
%
%                  or
%
%               [1,1] = size(query); cell = class(query)
%
%               The general form of a query general form:
%
%                  SELECT   <column list>
%                  FROM     <table list>
%                  [WHERE    <constraint list>]
%                  [ORDER BY <ORDER BY column list>]
%
%               (WHERE and ORDER BY are optional parameters)
%
%   the call:
%
%      [nmrows, error, errmsg] = cspice_ekfind( query )
%
%   returns:
%
%      query    the number of rows matching the query.
%
%               [1,1] = size(query); int32 = class(query)
%
%      error    a boolean indicating whether the query parsed correctly
%               (false) or not (true).
%
%               [1,1] = size(error); logical = class(error)
%
%      errmsg   the description of the parse error, should one occur,
%               otherwise the string returns as blank.
%
%               [1,c2] = size(errmsg); char = class(errmsg)
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Perform a query on an EK file that contains a database with
%      the different commands of the Deep Impact spacecraft subsystem,
%      and a table with the subsystem id, the parameter name, and the
%      description of that parameters, ordered by subsystem name. Print
%      the number of records of that table.
%
%      Use the EK kernel below to load the Deep Impact spacecraft
%      subsystem commands dictionary.
%
%         dif_cmdict_128_20050620.bdb
%
%
%      Example code begins here.
%
%
%      function ekfind_ex1()
%
%         %
%         % Assign an EK file to load.
%         %
%         EK = 'dif_cmdict_128_20050620.bdb';
%
%         %
%         % Load the EK.
%         %
%         cspice_furnsh( EK )
%
%         %
%         % The EK file contains the table 'DIF_COMMANDS',
%         % and that 'DIF_COMMANDS' contains columns named:
%         %
%         %   SUBSYSTEM, COMMAND, PARAMETER_NAME, DESCRIPTION
%         %
%         % Define a set of constraints to perform a query on all
%         % loaded EK files (the SELECT clause).
%         %
%         query = [ 'Select SUBSYSTEM, COMMAND, PARAMETER_NAME, '    ...
%                   'DESCRIPTION from DIF_COMMANDS order by SUBSYSTEM' ];
%
%         %
%         % Query the EK system for data rows matching the
%         % SELECT constraints.
%         %
%         [nmrows, error, errmsg] = cspice_ekfind( query );
%
%         %
%         % Check whether an error occurred while processing the
%         % SELECT clause. If so, output the error message.
%         %
%         if ( error )
%
%            printf( 'SELECT clause error: %s\n', errmsg );
%
%         else
%
%            %
%            % If no error occurred, `nmrows' contains the number of rows
%            % matching the constraints specified in the query string.
%            %
%            fprintf( 'Number of matching rows: %d\n', nmrows )
%
%         end
%
%         %
%         % Clear the kernel pool and database. Note, you don't normally
%         % unload an EK after a query, rather at the end of a program.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Number of matching rows: 5798
%
%
%   2) Examples of strings containing syntactically valid queries:
%
%          SELECT COL1 FROM TAB1
%
%          select col1 from tab1 where col1 gt 5
%
%          SELECT COL2 FROM TAB1 WHERE COL2 > 5.7D0 ORDER BY COL2
%
%          SELECT COL2 FROM TAB1 WHERE COL1 != 5
%
%          SELECT COL2 FROM TAB1 WHERE COL1 GE COL2
%
%          SELECT COL1, COL2, COL3 FROM TAB1 ORDER BY COL1
%
%          SELECT COL3 FROM TAB1 WHERE COL5 EQ "ABC"
%
%          SELECT COL3 FROM TAB1 WHERE COL5 = 'ABC'
%
%          SELECT COL3 FROM TAB1 WHERE COL5 LIKE 'A*'
%
%          SELECT COL3 FROM TAB1 WHERE COL5 LIKE 'A%%'
%
%          SELECT COL4 FROM TAB1 WHERE COL4 = '1995 JAN 1 12:38:09.7'
%
%          SELECT COL4 FROM TAB1 WHERE COL4 = "1995 JAN 1 12:38:09.7"
%
%          SELECT COL4 FROM TAB1 WHERE
%          COL4 NE 'GLL SCLK 02724646:67:7:2'
%
%          SELECT COL1 FROM TAB1 WHERE COL1 != NULL
%
%          SELECT COL1 FROM TAB1 WHERE COL1 IS NULL
%
%          SELECT COL1 FROM TAB1 WHERE COL1 IS NOT NULL
%
%          SELECT COL1, COL2, COL3 FROM TAB1
%          WHERE (COL1 BETWEEN 4 AND 6) AND (COL3 NOT LIKE "A%%")
%          ORDER BY COL1, COL3
%
%          SELECT COL4 FROM TAB1
%          WHERE COL4 BETWEEN "1995 JAN 1 12:38" AND
%          "October 23, 1995"
%
%          SELECT COL1, COL2 FROM TAB1 WHERE
%          NOT (    ( ( COL1 <  COL2 ) AND ( COL1 > 5   ) )  OR
%                   ( ( COL1 >= COL2 ) AND ( COL2 <= 10 ) )      )
%
%
%          SELECT T1.COL1, T1.COL2, T2.COL2, T2.COL3
%          FROM TABLE1 T1, TABLE2 T2
%          WHERE T1.COL1 = T2.COL1
%          AND T1.COL2 > 5
%          ORDER BY T1.COL1, T2.COL2
%
%
%   3) Examples of syntactically invalid queries:
%
%          SELECT TIME WHERE TIME
%          LT 1991 JAN 1                      {FROM clause is absent}
%
%          select time from table1 where
%          time lt 1991 jan 1                 {time string is not
%                                              quoted}
%
%          select time from table1
%          where time .lt. '1991 jan 1'       {operator should be lt}
%
%          select cmd from table1
%          where "cmd,6tmchg" != cmd          {value is on left side
%                                              of operator}
%
%          select event_type from table1
%          where event_type eq ""             {quoted string is empty
%                                              ---use " " to indicate
%                                              a blank string}
%
%          select event_type from table1
%          where event_type = "COMMENT"
%          order TIME                         {ORDER BY phrase is
%                                              lacking BY keyword}
%
%          select COL1 from table where
%          where COL1 eq MOC_EVENT            {literal string on
%                                              right-hand-side of
%                                              operator is not quoted}
%
%
%
%       In the following examples, we'll assume that the program
%       calling cspice_ekfind has loaded an EK containing two segments
%       having columns having the following names and attributes:
%
%
%        TABLE1:
%        ==========
%
%          Column name        Data type         Size       Indexed?
%          -----------        ---------         ----       --------
%          EVENT_TYPE         CHARACTER*32      1          YES
%          EVENT_PARAMETERS   CHARACTER*(*)     1          NO
%          COMMENT            CHARACTER*80      VARIABLE   NO
%
%
%        TABLE2:
%        ==========
%
%          Column name        Data type         Size       Indexed?
%          -----------        ---------         ----       --------
%          EVENT_TYPE         CHARACTER*32      1          YES
%          EVENT_PARAMETERS   CHARACTER*80      1          NO
%          COMMENT            CHARACTER*80      VARIABLE   NO
%          COMMAND            CHARACTER*80      1          YES
%
%
%       Then the following queries are semantically invalid:
%
%          SELECT EVENT_PARAMETERS
%          FROM TABLE1
%          WHERE EVENT_DURATION = 7.0         {No column called
%                                              EVENT_DURATION
%                                             is present in a loaded
%                                             EK}
%
%          SELECT COMMENT FROM TABLE2
%          WHERE COMMENT EQ "N/A"             {The COMMENT column does
%                                              not have size 1 and
%                                              therefore cannot be
%                                              referenced in a query}
%
%-Particulars
%
%   This routine operates almost entirely by side effects: it
%   prepares the EK fetch routines to return event data that
%   satisfy the input query. See the EK Required Reading for examples
%   of use of this routine in conjunction with the EK fetch routines.
%
%-Exceptions
%
%   1)  Most of the exceptions that can occur on a call to
%       cspice_ekfind are caused by errors in the input query. cspice_ekfind
%       attempts to diagnose these via the output error flag and
%       error message, instead of signaling errors. The following
%       classes of errors are detected:
%
%          Scanning errors---these result from badly formed query
%          in which cspice_ekfind could not identify all of the tokens.
%          When these errors occur, cspice_ekfind may be too confused to
%          give a helpful diagnostic message.
%
%          Parsing errors---these result from a badly formed
%          query that cspice_ekfind was able to separate into tokens
%          but that cspice_ekfind determined to be syntactically invalid:
%
%          Name resolution errors---these result from referencing
%          invalid or ambiguous column or table names in a query.
%
%          Time resolution errors---these result from use of time
%          strings that cannot be parsed.
%
%          Semantic errors---these result from a syntactically
%          valid query that violates a limit or a restriction on
%          values used in a query.
%
%
%   Some problems with queries are not trapped by cspice_ekfind but
%   instead cause errors to be signaled. These are listed below.
%
%   2)  If no E-kernels are loaded at the time this routine is called,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  If a leapseconds kernel is is not loaded before this routine
%       is called, UTC time values may not be used in queries. If they
%       are, an error is signaled by a routine in the call tree of
%       this routine.
%
%   4)  If an SCLK kernel for the appropriate spacecraft clock has not
%       been loaded before this routine is called, SCLK values for
%       that clock may not be used in queries. If they are, an error
%       is signaled by a routine in the call tree of this routine.
%
%   5)  If the input argument `query' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   6)  If the input argument `query' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  A leapseconds kernel must be loaded before this routine may
%       be called, if UTC time values are used in input queries.
%
%   2)  An appropriate SCLK kernel must be loaded before this routine
%       may be called, if SCLK values are used in input queries.
%
%-Required_Reading
%
%   MICE.REQ
%   EK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.3.0, 13-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "ok" to "error".
%
%       Edited the header to comply with NAIF standard. Fixed description
%       of output argument "error". Added example's problem statement and
%       and example's EK. Updated example code to work with provided EK.
%       Added examples of syntactically valid and invalid queries.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.2.1, 03-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.2.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.0, 10-APR-2010 (EDW)
%
%-Index_Entries
%
%   find EK data
%   issue EK query
%
%-&

function [nmrows, error, errmsg] = cspice_ekfind(query)

   switch nargin
      case 1

         sample = zzmice_str(query);

      otherwise

         error ( [ 'Usage: [nmrows, error, `errmsg`] = ' ...
                   'cspice_ekfind( `query` )' ])

   end

   %
   % Call the MEX library.
   %
   try
      [nmrows, error, errmsg] = mice('ekfind_c', query );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      error = zzmice_logical(error);
   catch spiceerr
      rethrow(spiceerr)
   end
