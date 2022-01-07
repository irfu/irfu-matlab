%-Abstract
%
%   CSPICE_DASLLC closes the DAS file associated with a given handle, without
%   flushing buffered data or segregating the file.
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
%      handle   the handle of a previously opened DAS file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dasllc( handle )
%
%   returns:
%
%      None.
%
%      See -Particulars for a description of the effect of this routine.
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
%      function dasllc_ex1()
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
%         fname = 'dasllc_ex1.das';
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
%   Normally, routines outside of Mice will not need to call this
%   routine. Application programs should close DAS files by calling
%   the Mice routine cspice_dascls. This routine is a lower-level
%   routine that is called by cspice_dascls, but (obviously) does not have
%   the full functionality of cspice_dascls.
%
%   This routine closes a DAS file and updates the DAS file manager's
%   bookkeeping information on open DAS files. Because the DAS file
%   manager must keep track of which files are open at any given time,
%   it is important that DAS files be closed only with cspice_dascls or
%   cspice_dasllc, to prevent the remaining DAS routines from failing,
%   sometimes mysteriously.
%
%   Note that when a file is opened more than once for read or write
%   access, cspice_dasopr returns the same handle each time it is re-opened.
%   Each time the file is closed, cspice_dasllc checks to see if any other
%   claims on the file are still active before physically closing
%   the file.
%
%   Unlike cspice_dascls, this routine does not force a write of updated,
%   buffered records to the indicated file, nor does it segregate the
%   data records in the file.
%
%-Exceptions
%
%   1)  If the specified handle does not belong to a DAS file that is
%       currently open, this routine returns without signaling an
%       error.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
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
%   close a DAS file
%
%-&
function cspice_dasllc( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dasllc( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dasllc_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
