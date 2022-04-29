%-Abstract
%
%   CSPICE_EKNELT returns the number of elements in a specified column entry
%   in the current row.
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
%      selidx   the index of the column in the SELECT from which to retrieve
%               data. The range of 'selidx' is 1:nsel inclusive, where
%               'nsel' is the number of items in the SELECT clause of the
%               current query.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%      row      the index of the row containing the element.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%               This number refers to a member of the set of rows matching
%               a query. 'row' must be in the range
%
%                  1:nmrows
%
%               where 'nmrows' is the matching row count returned
%               by cspice_ekfind.
%
%   the call:
%
%       [nelt] = cspice_eknelt( selidx, row )
%
%   returns:
%
%      nelt    the number of elements in the column entry belonging to the
%              specified column in the specified row.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%
%              Null entries in variable-size columns are considered to have
%              size 1.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Perform a query on an EK file that contains a database with
%      the different commands of the Deep Impact spacecraft subsystem,
%      and a table with the commands, their parameter names and the
%      parameters' ranges. The parameter's ranges are provided in a
%      a character column with variable number of items. Obtain the
%      parameter ranges for given parameter name for all commands that
%      use that parameter.
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
%      function eknelt_ex1()
%
%         %
%         % Assign an EK file to load.
%         %
%         EK = 'dif_cmdict_128_20050620.bdb';
%         MAXSTR = 1025;
%
%         %
%         % Load the EK.
%         %
%         cspice_furnsh( EK )
%
%         %
%         % The file "test_file.EK" contains the table 'vector_1', and
%         % 'vector_1' has the column named 'd_col_1', a vector of double
%         % precision values.
%         %
%
%         %
%         % Define a set of constraints to perform a query on all
%         % loaded EK files (the SELECT clause). In this case select
%         % the column "d_col_1" from table "vector_1."
%         %
%         query = [ 'Select COMMAND, PARAMETER_RANGE from DIF_COMMANDS ',  ...
%                   'where PARAMETER_NAME="DEVICE_SELECT"' ];
%
%         %
%         % Query the EK system for data rows matching the
%         % SELECT restraints.
%         %
%         [nmrows, ok, errmsg] = cspice_ekfind( query );
%
%         %
%         % Check whether an error occurred while processing the
%         % SELECT clause. If so, output the error message.
%         %
%         if ( ok )
%            printf( 'SELECT clause error: %s\n', errmsg );
%         end
%
%         %
%         % Loop over each row found matching the query.
%         %
%         for rowno = 1:nmrows
%
%            %
%            % Fetch the command name. We know it's only one element.
%            %
%            selidx = 1;
%            [cdata, isnull, found] = cspice_ekgc( selidx, rowno, 1, MAXSTR );
%
%            if  ~isnull
%               fprintf( 'Row: %d\nCommand: %s\n', rowno, cdata );
%            end
%
%
%            %
%            % Fetch now the parameter range data. We know the query
%            % returned these data in the second column, which may have
%            % multiple elements.
%            %
%            selidx = 2;
%            nelt   = cspice_eknelt( selidx, rowno);
%
%            %
%            % Use cspice_ekgc to retrieve the value from
%            % row/column position.
%            %
%            for eltidx = 1:nelt
%
%               [cdata, isnull, found] = cspice_ekgc( selidx, rowno,       ...
%                                                     eltidx, MAXSTR );
%
%               %
%               % Output the value, if non-null data exist at the
%               % requested position.
%               %
%               if  ~isnull
%                  fprintf( '  Range (elm %d): %s\n', eltidx, cdata );
%               end
%
%            end
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
%      Row: 1
%      Command: PWFPSW_PRI_SEL
%        Range (elm 1): 'RW_SET'='7', 'INST_SET'='21', 'SBAND_RX'='18', 'SBA***
%        Range (elm 2): 'IMPACTOR_SOH'='23'
%      Row: 2
%      Command: PWFPSW_AVAIL
%        Range (elm 1): 'PI'='0', 'CTB'='1', 'RIU'='2', 'XBT'='3', 'HTRS_ESB***
%        Range (elm 2): 'HTRS_NEB'='6', 'RW1'='7', 'RW2'='8', 'RW3'='9', 'RW***
%        Range (elm 3): 'STRTKR'='12', 'THRSTRS'='13', 'LVALVES'='14', 'SPAR***
%        Range (elm 4): 'TWTA'='17', 'SBANDRX'='18', 'SBANDTX'='19', 'SPARE2***
%        Range (elm 5): 'HRI'='22', 'IMP_SOH'='23'
%      Row: 3
%      Command: PWFCTL_CTB
%        Range (elm 1): 'CTB_A'='2', 'CTB_B'='3'
%      Row: 4
%      Command: PWFCTL_HTRS_ESB
%        Range (elm 1): 'SRVHTR_A'='2', 'SRVHTR_B'='3'
%      Row: 5
%      Command: PWFCTL_SCUCPU
%        Range (elm 1): 'SCUCPU_A'='2', 'SCUCPU_B'='3'
%      Row: 6
%      Command: PWFCTL_HTRS_NEB
%        Range (elm 1): 'HTRS_NEB_A'='2', 'HTRS_NEB_B'='3'
%
%
%      Warning: incomplete output. 5 lines extended past the right
%      margin of the header and have been truncated. These lines are
%      marked by "***" at the end of each line.
%
%
%-Particulars
%
%   This routine is meant to be used in conjunction with the EK fetch
%   entry points cspice_ekgc, cspice_ekgd, and cspice_ekgi. This routine
%   allows the caller of those routines to determine appropriate
%   loop bounds to use to fetch each column entry in the current row.
%
%-Exceptions
%
%   1)  If this routine is called when no E-kernels have been loaded,
%       the error SPICE(NOLOADEDFILES) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If `selidx' is outside of the range established by the last query
%       passed to the EK search engine, the error SPICE(INVALIDINDEX) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If `row' is outside of the range established by the last query passed
%       to the EK search engine, the error SPICE(INVALIDINDEX) is signaled by
%       a routine in the call tree of this routine.
%
%   4)  If any of the input arguments, `selidx' or `row', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `selidx' or `row', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   This routine reads binary "sequence component" EK files.
%   In order for a binary EK file to be accessible to this routine,
%   the file must be "loaded" via a call to the routine cspice_furnsh.
%
%   Text format EK files cannot be used by this routine; they must
%   first be converted by binary format by the NAIF Toolkit utility
%   SPACIT.
%
%-Restrictions
%
%   None.
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Changed the example
%       code to work with Deep Impact PDS archived data, and added the
%       corresponding problem statement.
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
%   -Mice Version 1.0.1, 03-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 10-APR-2010 (EDW)
%
%-Index_Entries
%
%   return the number of elements in a column entry
%
%-&

function [nelt] = cspice_eknelt( selidx, row )

   switch nargin

      case 2

         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);

      otherwise

         error ( 'Usage: [nelt] = cspice_eknelt( selidx, row )' )

   end

   try
      [nelt] = mice('eknelt_c', selidx, row );
   catch spiceerr
      rethrow(spiceerr)
   end
