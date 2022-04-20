%-Abstract
%
%   CSPICE_DAFGN returns the name for current array in the current
%   DAF being searched
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
%      None.
%
%   the call:
%
%      name = cspice_dafgn
%
%   returns:
%
%      name     the name of the current DAF array - that array
%               found by a previous call to cspice_daffna or cspice_daffpa.
%
%               [1,c1] = size(name); char = class(name)
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
%   1) Read the name for each array in a DAF file.
%
%      Use the SPK kernel below as input DAF file for the example.
%
%         OUTERPLANETS_V0002.BSP
%
%
%      Example code begins here.
%
%
%      function dafgn_ex1()
%
%         %
%         % Define a DAF from which to read the name for each array in
%         % the DAF.
%         %
%         DAF = 'OUTERPLANETS_V0002.BSP';
%         NI  = 6;
%         ND  = 2;
%
%         %
%         % Open the DAF for read
%         %
%         handle = cspice_dafopr( DAF );
%
%         %
%         % Begin a forward search on DAF.
%         %
%         cspice_dafbfs( handle )
%         found = cspice_daffna;
%
%         %
%         % Loop while found
%         %
%         while ( found )
%
%            [dc, ic] = cspice_dafgs( ND, NI );
%            name = cspice_dafgn;
%
%            %
%            % Output each array name.
%            %
%            fprintf( '%s\n', name)
%
%            %
%            % Check for a next segment.
%            %
%            found = cspice_daffna;
%
%         end
%
%         %
%         % SAFELY close the file.
%         %
%         cspice_dafcls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      JUP230
%      SAT261xl
%      URA083
%      NEP016.6
%
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
%   The main function of these routines is to allow the
%   contents of any DAF to be examined on an array-by-array
%   basis.
%
%   Conceptually, the arrays in a DAF form a doubly linked list,
%   which can be searched in either of two directions: forward or
%   backward. It is possible to search multiple DAFs simultaneously.
%
%   cspice_dafbfs (begin forward search) and cspice_daffna are used to
%   search the arrays in a DAF in forward order. In applications that
%   search a single DAF at a time, the normal usage is
%
%      cspice_dafbfs( handle );
%      found = cspice_daffna;
%
%      while found
%
%         [dc, ic] = cspice_dafgs( ND, NI );
%         name = cspice_dafgn;
%
%          .
%          .
%
%         found = cspice_daffna;
%
%      end
%
%
%   cspice_dafbbs (begin backward search) and cspice_daffpa are used to
%   search the arrays in a DAF in backward order. In applications that
%   search a single DAF at a time, the normal usage is
%
%      cspice_dafbbs( handle );
%      found = cspice_daffpa;
%
%      while found
%
%         [dc, ic] = cspice_dafgs( ND, NI );
%         name = cspice_dafgn;
%
%          .
%          .
%
%         found = cspice_daffpa;
%
%      end
%
%
%   In applications that conduct multiple searches simultaneously, the above
%   usage must be modified to specify the handle of the file to operate on,
%   in any case where the file may not be the last one specified by
%   cspice_dafbfs or cspice_dafbbs. The routine cspice_dafcs (DAF, continue
%   search) is used for this purpose. Below, we give an example of an
%   interleaved search of two files specified by the handles `handl1' and
%   `handl2'. The directions of searches in different DAFs are independent;
%   here we conduct a forward search on one file and a backward search on the
%   other. Throughout, we use cspice_dafcs to specify which file to operate
%   on, before calling cspice_daffna, cspice_daffpa, cspice_dafgs, or
%   cspice_dafgn.
%
%
%      cspice_dafbfs( handl1 );
%      cspice_dafbbs( handl2 );
%
%      cspice_dafcs( handl1 );
%      found1 = cspice_daffna;
%
%      cspice_dafcs( handl2 );
%      found2 = cspice_daffpa;
%
%      while ( found1 | found2 )
%
%         if found1
%
%            cspice_dafcs( handl1 );
%            [dc, ic] = cspice_dafgs( ND, NI );
%            name = cspice_dafgn;
%
%             .
%             .
%
%            cspice_dafcs( handl1 );
%            found1 = cspice_daffna;
%
%         end
%
%         if found2
%
%            cspice_dafcs( handl2 );
%            [dc, ic] = cspice_dafgs( ND, NI );
%            name = cspice_dafgn;
%
%             .
%             .
%
%            cspice_dafcs( handl2 );
%            found2 = cspice_daffna;
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
%   DAF. The summary and name for the current array in the current DAF can
%   be obtained separately, as shown above, by calls to cspice_dafgs
%   (get summary) and cspice_dafgn (get name).
%
%   Once a search has been begun, it may be continued in either
%   direction. That is, cspice_daffpa may be used to back up during a
%   forward search, and cspice_daffna may be used to advance during a
%   backward search.
%
%-Exceptions
%
%   1)  If this routine is called when no search is in progress in the
%       current DAF, the error SPICE(DAFNOSEARCH) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the DAF for which the "current" array's name is to be
%       returned has actually been closed, an error is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If no array is current in the current DAF, the error
%       SPICE(NOCURRENTARRAY) is signaled by a routine in the call tree of
%       this routine. There is no current array when a search is started by
%       cspice_dafbfs or cspice_dafbbs, but no calls to cspice_daffna or
%       cspice_daffpa have been made yet, or whenever cspice_daffna or
%       cspice_daffpa return the value false in the `found' argument.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement and modified input DAF file.
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
%   -Mice Version 1.0.0, 11-JUN-2013 (EDW)
%
%-Index_Entries
%
%   get DAF array name
%
%-&

function [name] = cspice_dafgn

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [`name`] = cspice_dafgn' )

   end


   %
   % Call the MEX library.
   %
   try
      [name] = mice( 'dafgn_c' );
   catch spiceerr
      rethrow(spiceerr)
   end
