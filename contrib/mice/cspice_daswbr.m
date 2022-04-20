%-Abstract
%
%   CSPICE_DASWBR writes out all buffered records of a specified DAS file.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      handle   the handle of a DAS file opened for writing.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_daswbr( handle )
%
%   returns:
%
%      None.
%
%      See -Particulars for a description of the action of this routine.
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
%   1) Write a DAS file by adding data to it over multiple passes.
%      Avoid spending time on file segregation between writes.
%
%      Each pass opens the file, adds character, double precision,
%      and integer data to the file, writes out buffered data by
%      calling cspice_daswbr, and closes the file without segregating the
%      data by calling cspice_dasllc.
%
%      The program also checks the file: after the final write,
%      the program reads the data and compares it to expected values.
%
%      Note that most user-oriented applications should segregate a
%      DAS file after writing it, since this greatly enhances file
%      reading efficiency. The technique demonstrated here may be
%      useful for cases in which a file will be written via many
%      small data additions, and in which the file is read between
%      write operations.
%
%
%      Example code begins here.
%
%
%      function daswbr_ex1()
%
%         %
%         % Local parameters
%         %
%         CHRLEN =   50;
%         IBUFSZ =   20;
%         DBUFSZ =   30;
%
%         %
%         % Local variables
%         %
%         dpbuf  = zeros(DBUFSZ,1);
%         xdpbuf = zeros(DBUFSZ,1);
%
%         intbuf = zeros(IBUFSZ,1, 'int32');
%         xintbf = zeros(IBUFSZ,1, 'int32');
%
%         bytbuf = zeros( 1, CHRLEN, 'uint8' );
%         xbytbf = zeros( 1, CHRLEN, 'uint8' );
%
%         %
%         % Initial values
%         %
%         fname = 'daswbr_ex1.das';
%         ftype = 'ANG';
%         ncall = 1000;
%         ncomr = 10;
%         npass = 3;
%
%         %
%         % Open a new DAS file. We'll allocate `ncomr' records
%         % for comments. The file type is not one of the standard
%         % types recognized by SPICE; however it can be used to
%         % ensure the database file is of the correct type.
%         %
%         % We'll use the file name as the internal file name.
%         %
%         [handle] = cspice_dasonw( fname, ftype, fname, ncomr );
%
%         %
%         % Add data of character, integer, and double precision
%         % types to the file in interleaved fashion. We'll add to
%         % the file over `npass' "passes," in each of which we close
%         % the file after writing.
%         %
%         for passno=1:npass
%
%            if ( passno > 1 )
%
%               fprintf( 'Opening file for write access...\n' )
%
%               [handle] = cspice_dasopw( fname );
%
%            end
%
%            for i=1:ncall
%
%               %
%               % Add string data to the file.
%               %
%               chrbuf   = 'Character value #';
%               [chrbuf] = cspice_repmi( chrbuf, '#', i );
%
%               bytbuf(1:length(chrbuf)) = uint8(chrbuf);
%               cspice_dasadc( handle, CHRLEN, 1, CHRLEN, bytbuf );
%
%               %
%               % Add double precision data to the file.
%               %
%               for j=1:DBUFSZ
%
%                  dpbuf(j) = double( 100000000*passno + 100*i + j );
%
%               end
%
%               cspice_dasadd( handle, dpbuf );
%
%               %
%               % Add integer data to the file.
%               %
%               for j=1:IBUFSZ
%
%                  intbuf(j) = 100000000*passno  +  100 * i  +  j;
%
%               end
%
%               cspice_dasadi( handle, intbuf );
%
%            end
%
%            %
%            % Write buffered data to the file.
%            %
%            fprintf( 'Writing buffered data...\n' )
%            cspice_daswbr( handle );
%
%            %
%            % Close the file without segregating it.
%            %
%            fprintf( 'Closing DAS file...\n' )
%            cspice_dasllc( handle );
%
%         end
%
%         fprintf( 'File write is done.\n' )
%
%         %
%         % Check file contents.
%         %
%         [handle] = cspice_dasopr( fname );
%
%         %
%         % Read data from the file; compare to expected values.
%         %
%         % Initialize end addresses.
%         %
%         lastc = 0;
%         lastd = 0;
%         lasti = 0;
%
%         for passno=1:npass
%
%            for i=1:ncall
%
%               %
%               % Check string data.
%               %
%               xchrbf   = 'Character value #';
%               [xchrbf] = cspice_repmi( xchrbf, '#', i );
%
%               firstc   = lastc + 1;
%               lastc    = lastc + CHRLEN;
%
%               xbytbf(1:length(xchrbf)) = uint8(xchrbf);
%               [chrbuf] = cspice_dasrdc( handle, firstc, lastc,           ...
%                                         1,      CHRLEN, bytbuf );
%
%               if ( bytbuf != xbytbf )
%
%                  fprintf( 'Character data mismatch:\n' )
%                  fprintf( 'PASS     = %d\n', passno )
%                  fprintf( 'I        = %d\n', i )
%                  fprintf( 'Expected = %s\n', char(xbytbf) )
%                  fprintf( 'Actual   = %s\n', char(bytbuf) )
%                  exit;
%
%               end
%
%               %
%               % Check double precision data.
%               %
%               for j=1:DBUFSZ
%
%                  xdpbuf(j) = double( 100000000*passno + 100*i + j );
%
%               end
%
%               firstd  = lastd + 1;
%               lastd   = lastd + DBUFSZ;
%
%               [dpbuf] = cspice_dasrdd( handle, firstd, lastd );
%
%               for j=1:DBUFSZ
%
%                  if ( dpbuf(j) != xdpbuf(j) )
%
%                     fprintf( 'Double precision data mismatch:\n' )
%                     fprintf( 'PASS     = %d\n', passno )
%                     fprintf( 'I        = %d\n', i )
%                     fprintf( 'J        = %d\n', j )
%                     fprintf( 'Expected = %f\n', xdpbuf(j) )
%                     fprintf( 'Actual   = %f\n', dpbuf(j) )
%                     exit;
%
%                  end
%
%               end
%
%               %
%               % Check integer data.
%               %
%               for j=1:IBUFSZ
%
%                  xintbf(j) = 100000000*passno  +  100 * i  +  j;
%
%               end
%
%               firsti   = lasti + 1;
%               lasti    = lasti + IBUFSZ;
%
%               [intbuf] = cspice_dasrdi( handle, firsti, lasti );
%
%               for j=1:IBUFSZ
%
%                  if ( intbuf(j) != xintbf(j) )
%
%                     fprintf( 'Integer data mismatch:\n' )
%                     fprintf( 'PASS     = %d\n', passno )
%                     fprintf( 'I        = %d\n', i )
%                     fprintf( 'J        = %d\n', j )
%                     fprintf( 'Expected = %d\n', xintbf(j) )
%                     fprintf( 'Actual   = %d\n', intbuf(j) )
%                     exit;
%
%                  end
%
%               end
%
%            end
%
%         end
%
%         fprintf( 'File check is done.\n' )
%
%         %
%         % Close the file.
%         %
%         cspice_dascls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Writing buffered data...
%      Closing DAS file...
%      Opening file for write access...
%      Writing buffered data...
%      Closing DAS file...
%      Opening file for write access...
%      Writing buffered data...
%      Closing DAS file...
%      File write is done.
%      File check is done.
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   This routine writes buffered records out to the DAS file to which
%   they correspond.
%
%   Because the DAS system buffers records that are written as well as those
%   that are read, data supplied to the DAS add data (cspice_dasadc,
%   cspice_dasadd, cspice_dasadi) and DAS update data (cspice_dasudc,
%   cspice_dasudd, cspice_dasudi) routines on input has not necessarily been
%   physically written to the DAS file specified by the caller of those
%   routines, at the time those routines return. Before closing a DAS file
%   that has been opened for writing, the DAS system must write out to the
%   file any updated records present in the DAS buffers. The Mice routine
%   cspice_dascls uses this routine to perform this function. The routines
%   DASAC and DASDC, through the use of the SPICELIB routines DASACR and
%   DASRCR, which respectively add comment records to or delete comment
%   records from a DAS file, use this routine to ensure that the SPICELIB
%   routine DASRWR record buffers don't become out of sync with the file they
%   operate upon.
%
%   In addition, this routine can be used by application programs
%   that create or update DAS files. The reason for calling this
%   routine directly would be to provide a measure of safety when
%   writing a very large file: if the file creation or update were
%   interrupted, the amount of work lost due to the loss of buffered,
%   unwritten records could be reduced.
%
%   However, routines outside of Mice will generally not need to
%   call this routine directly.
%
%-Exceptions
%
%   1)  If the input file handle is invalid, an error is signaled
%       by a routine in the call tree of this routine. The indicated
%       file will not be modified.
%
%   2)  If a write operation attempted by this routine fails, an
%       error is signaled by a routine in the call tree of this
%       routine. The status of the DAS file written to is uncertain
%       in this case.
%
%   3)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   4)  If the input argument `handle' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See the description of the argument `handle' in -I/O.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   DAS.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 29-JUN-2021 (JDR)
%
%-Index_Entries
%
%   write buffered records to a DAS file
%
%-&
function cspice_daswbr( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_daswbr( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('daswbr_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
