%-Abstract
%
%   CSPICE_DAFBBS initiates a backwards search for arrays in a DAF.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dafbbs( handle )
%
%   starts a backwards search, i.e. end of file to start of file,
%   on a DAF.
%
%   returns:
%
%      None.
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
%   1) Create a simple program to output the double precision and
%      integer values stored in an SPK's segments' descriptors. This
%      program opens a DAF for read, performs a backward search for
%      the DAF arrays, prints the segment descriptor for each array
%      found, then closes the DAF.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         de421.bsp
%
%
%      Example code begins here.
%
%
%      function dafbbs_ex1()
%
%         %
%         % Open a DAF for read. Return a 'handle' referring
%         % to the file.
%         %
%         kernel = 'de421.bsp';
%         handle = cspice_dafopr( kernel );
%
%         %
%         % Define the summary parameters appropriate
%         % for an SPK file.
%         %
%         ND = 2;
%         NI = 6;
%
%         %
%         % Begin a backwards search on the file.
%         %
%         cspice_dafbbs( handle );
%
%         %
%         % Search until a DAF array is found.
%         %
%         found = cspice_daffpa;
%
%         %
%         % Loop while the search finds previous DAF arrays.
%         %
%         while found
%
%            [dc, ic ] = cspice_dafgs( ND, NI );
%
%            fprintf( 'Doubles:  ' )
%            fprintf( '%f   ', dc )
%            fprintf( '\n' )
%
%            fprintf( 'Integers: ' )
%            fprintf( '%d   ', ic )
%            fprintf( '\n\n' )
%
%
%            %
%            % Check for another segment.
%            %
%            found = cspice_daffpa;
%
%         end
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 499   4   1   2   2098633   2098644
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 299   2   1   2   2098621   2098632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 199   1   1   2   2098609   2098620
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 399   3   1   2   1521325   2098608
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 301   3   1   2   944041   1521324
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 10   0   1   2   820837   944040
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 9   0   1   2   785633   820836
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 8   0   1   2   750429   785632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 7   0   1   2   715225   750428
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 6   0   1   2   674741   715224
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 5   0   1   2   628977   674740
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 4   0   1   2   567373   628976
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 3   0   1   2   423049   567372
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 2   0   1   2   310405   423048
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 1   0   1   2   641   310404
%
%
%      Note, the specific contents of 'ic' and 'dc' depend on the
%      type of DAF.
%
%      Note, the final entries in the integer array contain the segment
%      start/end indexes. The output indicates the search proceeded
%      from the end of the file (high value index) towards the beginning
%      (low value index).
%
%-Particulars
%
%   The DAF search routines are:
%
%      cspice_dafbfs       Begin forward search.
%      cspice_daffna       Find next array.
%
%      cspice_dafbbs       Begin backward search.
%      cspice_daffpa       Find previous array.
%
%      cspice_dafgs        Get summary.
%      cspice_dafgn        Get name.
%
%      cspice_dafcs        Continue search.
%
%   The main function of these entry points is to allow the
%   contents of any DAF to be examined on an array-by-array
%   basis.
%
%   Conceptually, the arrays in a DAF form a doubly linked list,
%   which can be searched in either of two directions: forward or
%   backward. It is possible to search multiple DAFs simultaneously.
%
%   cspice_dafbfs (begin forward search) and daffna are used to search the
%   arrays in a DAF in forward order. In applications that search a
%   single DAF at a time, the normal usage is
%
%      cspice_dafbfs( handle );
%      [found] = cspice_daffna;
%
%      while found
%
%         [dc, ic] = cspice_dafgs( ND, NI );
%         [sum]    = cspice_dafps( dc, ic );
%         [name]   = cspice_dafgn;
%                  .
%                  .
%         [found]  = cspice_daffna;
%
%      end
%
%
%   cspice_dafbbs (begin backward search) and cspice_daffpa are used to
%   search the arrays in a DAF in backward order. In applications that search
%   a single DAF at a time, the normal usage is
%
%      cspice_dafbbs( handle );
%      [found] = cspice_daffpa;
%
%      while found
%
%         [dc, ic] = cspice_dafgs( ND, NI );
%         [sum]    = cspice_dafps( dc, ic );
%         [name]   = cspice_dafgn;
%                  .
%                  .
%         [found]  = cspice_daffpa;
%
%       end
%
%
%   In applications that conduct multiple searches simultaneously, the above
%   usage must be modified to specify the handle of the file to operate on,
%   in any case where the file may not be the last one specified by
%   cspice_dafbfs or cspice_dafbbs. The routine cspice_dafcs (DAF, continue
%   search) is used for this purpose. Below, we give an example of an
%   interleaved search of two files specified by the handles handl1 and
%   handl2. The directions of searches in different DAFs are independent;
%   here we conduct a forward search on one file and a backward search on the
%   other. Throughout, we use dafcs to specify which file to operate on,
%   before calling cspice_daffna, cspice_daffpa, cspice_dafgs, or
%   cspice_dafgn.
%
%
%      cspice_dafbfs( handl1 );
%      cspice_dafbbs( handl2 );
%
%      cspice_dafcs( handl1 );
%      [found1] = cspice_daffna;
%
%      cspice_dafcs( handl2 );
%      [found2] = cspice_daffpa;
%
%      while ( found1 | found2 )
%
%         if found1
%
%            cspice_dafcs( handl1 );
%            [dc, ic] = cspice_dafgs( ND, NI );
%            [sum]    = cspice_dafps( dc, ic );
%            [name]   = cspice_dafgn;
%                     .
%                     .
%            cspice_dafcs( handl1 );
%            [found1] = cspice_daffna;
%
%         end
%
%         if found2
%
%            cspice_dafcs( handl2 );
%            [dc, ic] = cspice_dafgs( ND, NI );
%            [sum]    = cspice_dafps( dc, ic );
%            [name]   = cspice_dafgn;
%                     .
%                     .
%            cspice_dafcs( handl2 );
%            [found2] = cspice_daffpa;
%
%         end
%
%      end
%
%
%   At any time, the latest array found (whether by cspice_daffna or
%   cspice_daffpa) is regarded as the 'current' array for the file in which
%   the array was found. The last DAF in which a search was started,
%   executed, or continued by any of cspice_dafbfs, cspice_dafbbs,
%   cspice_daffna, cspice_daffpa or cspice_dafcs is regarded as the 'current'
%   DAF. The summary and name for the current array in the current DAF can be
%   obtained separately, as shown above, by calls to cspice_dafgs (get
%   summary) and cspice_dafgn (get name).
%
%   Once a search has been begun, it may be continued in either
%   direction. That is, cspice_daffpa may be used to back up during a
%   forward search, and cspice_daffna may be used to advance during a
%   backward search.
%
%-Exceptions
%
%   1)  If the input handle is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the summary record of the last record in the DAF file
%       cannot be read, the error SPICE(RECORDNOTFOUND) is signaled by
%       a routine in the call tree of this routine.
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
%   See argument `handle'.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DAF.REQ
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
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard.
%       Modified code example to hardcode the input DAF file.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 10-JUL-2012 (EDW)
%
%-Index_Entries
%
%   begin DAF backward search
%
%-&

function cspice_dafbbs( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafbbs(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafbbs_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end




