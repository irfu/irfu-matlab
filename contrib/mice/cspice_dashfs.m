%-Abstract
%
%   CSPICE_DASHFS returns a file summary for a specified DAS file.
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
%               The file may be open for read or write access.
%
%   the call:
%
%      [nresvr, nresvc, ncomr,  ncomc,                                     ...
%       free,   lastla, lastrc, lastwd] = cspice_dashfs( handle )
%
%   returns:
%
%      nresvr   the number of reserved records in a specified DAS file.
%
%               [1,1] = size(nresvr); int32 = class(nresvr)
%
%      nresvc   the number of characters in use in the reserved record area
%               of a specified DAS file.
%
%               [1,1] = size(nresvc); int32 = class(nresvc)
%
%      ncomr    the number of comment records in a specified DAS file.
%
%               [1,1] = size(ncomr); int32 = class(ncomr)
%
%      ncomc    the number of characters in use in the comment area of a
%               specified DAS file.
%
%               [1,1] = size(ncomc); int32 = class(ncomc)
%
%      free     the 1-based record number of the first free record in a
%               specified DAS file.
%
%               [1,1] = size(free); int32 = class(free)
%
%      lastla   an array containing the highest current 1-based logical
%               addresses, in the specified DAS file, of data of character,
%               double precision, and integer types, in that order.
%
%               [3,1] = size(lastla); int32 = class(lastla)
%
%      lastrc   an array containing the 1-based record numbers, in the
%               specified DAS file, of the directory records containing the
%               current last descriptors of clusters of character, double
%               precision, and integer data records, in that order.
%
%               [3,1] = size(lastrc); int32 = class(lastrc)
%
%      lastwd   an array containing the 1-based word indices, within the
%               respective descriptor records identified by the elements of
%               `lastrc', of the current last descriptors of clusters of
%               character, double precision, and integer data records, in
%               that order.
%
%               [3,1] = size(lastwd); int32 = class(lastwd)
%
%-Parameters
%
%   See parameter definitions file MiceDAS.m for declarations and
%   descriptions of parameters used throughout the DAS system.
%
%   SPICE_DAS_CHARDT,
%   SPICE_DAS_DPDT,
%   SPICE_DAS_INTDT
%
%               are data type specifiers which indicate SpiceChar
%               (uint8), SpiceDouble (double), and SpiceInt (int32)
%               respectively. These parameters are used in all DAS
%               routines that require a data type specifier.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Create a DAS file containing 10 integers, 5 double precision
%      numbers, and 4 characters. Print the summary of the file and
%      dump its contents.
%
%
%      Example code begins here.
%
%
%      function dashfs_ex1()
%
%         %
%         % MiceUser is a file that makes certain variables global.
%         % You must call MiceUser to have access to the parameters used
%         % in this example.
%         %
%         MiceUser;
%
%         %
%         % Local parameters.
%         %
%         FNAME  =   'dashfs_ex1.das';
%         LINLEN =   2;
%
%         %
%         % Local variables.
%         %
%         line = zeros( 1, LINLEN, 'uint8' );
%
%         %
%         % Open a new DAS file. Reserve no records for comments.
%         %
%         type     = 'TEST';
%         ifname   = 'TEST.DAS/NAIF/NJB/11-NOV-1992-20:12:20';
%
%         [handle] = cspice_dasonw( FNAME, type, ifname, 0 );
%
%         %
%         % Obtain the file summary.
%         %
%         [nresvr, nresvc,                                                 ...
%          ncomr,  ncomc,                                                  ...
%          free,   lastla,                                                 ...
%          lastrc, lastwd] = cspice_dashfs( handle );
%
%         %
%         % Print the summary of the new file.
%         %
%         fprintf( 'Summary before adding data:\n' )
%         fprintf( '   Number of reserved records     : %d\n', nresvr )
%         fprintf( '   Characters in reserved records : %d\n', nresvc )
%         fprintf( '   Number of comment records      : %d\n', ncomr )
%         fprintf( '   Characters in comment area     : %d\n', ncomc )
%         fprintf( '   Number of first free record    : %d\n', free )
%         fprintf( '   Last logical character address : %d\n',             ...
%                                              lastla(SPICE_DAS_CHARDT) )
%         fprintf( '   Last logical d.p. address      : %d\n',             ...
%                                              lastla(SPICE_DAS_DPDT) )
%         fprintf( '   Last logical integer address   : %d\n',             ...
%                                              lastla(SPICE_DAS_INTDT) )
%         fprintf( '   Last character descriptor      : %d\n',             ...
%                                              lastrc(SPICE_DAS_CHARDT) )
%         fprintf( '   Last d.p descriptor            : %d\n',             ...
%                                              lastrc(SPICE_DAS_DPDT) )
%         fprintf( '   Last integer descriptor        : %d\n',             ...
%                                              lastrc(SPICE_DAS_INTDT) )
%         fprintf( '   Character word position in desc: %d\n',             ...
%                                              lastwd(SPICE_DAS_CHARDT) )
%         fprintf( '   d.p. word position in desc     : %d\n',             ...
%                                              lastwd(SPICE_DAS_DPDT) )
%         fprintf( '   Integer word position in desc  : %d\n',             ...
%                                              lastwd(SPICE_DAS_INTDT) )
%         fprintf( '\n' )
%
%         %
%         % Add the data.
%         %
%         for i=1:10
%
%            cspice_dasadi( handle, i );
%
%         end
%
%         for i=1:5
%
%            cspice_dasadd( handle, double(i) );
%
%         end
%
%         %
%         % Add character data to the file. DAS character data are
%         % treated as a character array, not as a string. The
%         % following call adds only the first 4 characters to the
%         % DAS file.
%         %
%         cspice_dasadc( handle, 4, 1, 4, uint8('SPUDWXY') );
%
%         %
%         % Close the file and open it for reading.
%         %
%         cspice_dascls( handle );
%         [handle] = cspice_dasopr( FNAME );
%
%         %
%         % Obtain again the file summary.
%         %
%         [nresvr, nresvc, ncomr,  ncomc,                                  ...
%          free,   lastla, lastrc, lastwd] = cspice_dashfs( handle );
%
%         fprintf( 'Summary after adding data:\n' )
%         fprintf( '   Number of reserved records     : %d\n', nresvr )
%         fprintf( '   Characters in reserved records : %d\n', nresvc )
%         fprintf( '   Number of comment records      : %d\n', ncomr )
%         fprintf( '   Characters in comment area     : %d\n', ncomc )
%         fprintf( '   Number of first free record    : %d\n', free )
%         fprintf( '   Last logical character address : %d\n',             ...
%                                              lastla(SPICE_DAS_CHARDT) )
%         fprintf( '   Last logical d.p. address      : %d\n',             ...
%                                              lastla(SPICE_DAS_DPDT) )
%         fprintf( '   Last logical integer address   : %d\n',             ...
%                                              lastla(SPICE_DAS_INTDT) )
%         fprintf( '   Last character descriptor      : %d\n',             ...
%                                              lastrc(SPICE_DAS_CHARDT) )
%         fprintf( '   Last d.p descriptor            : %d\n',             ...
%                                              lastrc(SPICE_DAS_DPDT) )
%         fprintf( '   Last integer descriptor        : %d\n',             ...
%                                              lastrc(SPICE_DAS_INTDT) )
%         fprintf( '   Character word position in desc: %d\n',             ...
%                                              lastwd(SPICE_DAS_CHARDT) )
%         fprintf( '   d.p. word position in desc     : %d\n',             ...
%                                              lastwd(SPICE_DAS_DPDT) )
%         fprintf( '   Integer word position in desc  : %d\n',             ...
%                                              lastwd(SPICE_DAS_INTDT) )
%         fprintf( '\n' )
%
%         %
%         % Read the integers and dump them.
%         %
%         fprintf( 'Integer data in the DAS file:\n' )
%         for i=1:lastla(SPICE_DAS_INTDT)
%
%            [n] = cspice_dasrdi( handle, i, i );
%            fprintf( '   %d\n', n )
%
%         end
%
%         %
%         % Now the d.p. numbers:
%         %
%         fprintf( '\n' )
%         fprintf( 'Double precision data in the DAS file:\n' )
%         for i=1:lastla(SPICE_DAS_DPDT)
%
%            [x] = cspice_dasrdd( handle, i, i );
%            fprintf( '   %f\n', x )
%
%         end
%
%         %
%         % Now the characters. In this case, we read the
%         % data one line at a time.
%         %
%         first  =  0;
%         last   =  0;
%         remain =  lastla(SPICE_DAS_CHARDT);
%
%         fprintf( '\n' )
%         fprintf( 'Character data in the DAS file:\n' )
%         while remain > 0
%
%            nread  = min ( [LINLEN, remain] );
%            first  = last + 1;
%            last   = last + nread;
%
%            [line] = cspice_dasrdc( handle, first, last, 1, LINLEN, line );
%
%            fprintf( '   %s\n', char(line) )
%
%            remain = remain - nread;
%
%         end
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
%      Summary before adding data:
%         Number of reserved records     : 0
%         Characters in reserved records : 0
%         Number of comment records      : 0
%         Characters in comment area     : 0
%         Number of first free record    : 3
%         Last logical character address : 0
%         Last logical d.p. address      : 0
%         Last logical integer address   : 0
%         Last character descriptor      : 0
%         Last d.p descriptor            : 0
%         Last integer descriptor        : 0
%         Character word position in desc: 0
%         d.p. word position in desc     : 0
%         Integer word position in desc  : 0
%
%      Summary after adding data:
%         Number of reserved records     : 0
%         Characters in reserved records : 0
%         Number of comment records      : 0
%         Characters in comment area     : 0
%         Number of first free record    : 6
%         Last logical character address : 4
%         Last logical d.p. address      : 5
%         Last logical integer address   : 10
%         Last character descriptor      : 2
%         Last d.p descriptor            : 2
%         Last integer descriptor        : 2
%         Character word position in desc: 10
%         d.p. word position in desc     : 11
%         Integer word position in desc  : 12
%
%      Integer data in the DAS file:
%         1
%         2
%         3
%         4
%         5
%         6
%         7
%         8
%         9
%         10
%
%      Double precision data in the DAS file:
%         1.000000
%         2.000000
%         3.000000
%         4.000000
%         5.000000
%
%      Character data in the DAS file:
%         SP
%         UD
%
%
%      Note that after run completion, a new DAS file exists in the
%      output directory.
%
%-Particulars
%
%   The quantities
%
%      nresvr
%      nresvc
%      ncomr
%      ncomc
%      free
%      lastla
%      lastrc
%      lastwd
%
%   define the "state" of a DAS file, and in particular the state of
%   the directory structure of the file. This information is needed by
%   other DAS routines, but application programs will usually have no
%   need for it. The one exception is the array of "last" logical
%   addresses `lastla': these addresses indicate how many words of data
%   of each type are contained in the specified DAS file. The elements
%   of `lastla' can be conveniently retrieved by calling cspice_daslla.
%
%-Exceptions
%
%   1)  If the specified handle does not belong to any file that is
%       currently known to be open, the error SPICE(DASNOSUCHHANDLE)
%       is signaled by a routine in the call tree of this routine.
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
%   None.
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
%   -Mice Version 1.0.0, 23-JUL-2021 (JDR)
%
%-Index_Entries
%
%   return the file summary of a DAS file
%   find the amount of data in a DAS file
%
%-&
function [nresvr, nresvc, ncomr,  ncomc,                                   ...
          free,   lastla, lastrc, lastwd] = cspice_dashfs( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( [ 'Usage: [nresvr, nresvc, ncomr, ncomc, free, '          ...
                   'lastla(3), lastrc(3), lastwd(3)] = '                   ...
                   'cspice_dashfs( handle )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [nresvr, nresvc, ncomr,  ncomc,                                      ...
       free,   lastla, lastrc, lastwd] = mice('dashfs_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
