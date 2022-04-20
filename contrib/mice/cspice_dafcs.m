%-Abstract
%
%   CSPICE_DAFCS sets the active DAF to search. A search must be
%   in progress for the DAF.
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
%      handle   the file handle referring to a DAF to
%               set as the "active" file for a search.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dafcs( handle )
%
%   causes all DAF search activity apply to the file
%   referred to by 'handle'.
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
%   1) Given two SPK files, begin a forward search on one and, right
%      after, a backward search on the other. Restitute the search on
%      the first SPK, and print all the segment IDs found in the file.
%      Go back to the second file, and print all its segment IDs.
%
%      Use the SPK kernel below as first input DAF file for the program.
%
%         OUTERPLANETS_V0002.BSP
%
%      Use the SPK kernel below as second input DAF file for the program.
%
%         sat382-rocks-merge.bsp
%
%
%      Example code begins here.
%
%
%      function dafcs_ex1()
%
%         %
%         % Define two SPK test files.
%         %
%         SPK1 = 'OUTERPLANETS_V0002.BSP';
%         SPK2 = 'sat382-rocks-merge.bsp';
%
%         %
%         % Open the DAFs for read
%         %
%         han1 = cspice_dafopr( SPK1 );
%         han2 = cspice_dafopr( SPK2 );
%
%         %
%         % Begin a forward search on SPK1
%         %
%         cspice_dafbfs( han1 )
%         found = cspice_daffna;
%
%         %
%         % Begin a backwards search on SPK2
%         %
%         cspice_dafbbs( han2 )
%         found2 = cspice_daffpa;
%
%         %
%         % Reinstitute the search on han1, loop
%         % so long as segment data are found.
%         %
%         cspice_dafcs( han1 );
%         fprintf( 'Segment IDs found on forward search of: %s\n', SPK1 );
%
%         while ( found )
%
%            segid    = cspice_dafgn;
%            found    = cspice_daffna;
%
%            %
%            % Output each segment ID.
%            %
%            fprintf( '%s\n', segid )
%
%         end
%
%         %
%         % Reinstitute the search on han2, loop
%         % so long as segment data are found.
%         %
%         cspice_dafcs( han2 )
%         fprintf( '\nSegment IDs found on backward search of: %s\n', SPK2 );
%
%         while ( found2 )
%
%            segid    = cspice_dafgn;
%            found2   = cspice_daffpa;
%
%            %
%            % Output each segment ID.
%            %
%            fprintf( '%s\n', segid )
%
%         end
%
%         %
%         % Close the files.
%         %
%         cspice_dafcls( han1 )
%         cspice_dafcls( han2 )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Segment IDs found on forward search of: OUTERPLANETS_V0002.BSP
%      JUP230
%      SAT261xl
%      URA083
%      NEP016.6
%
%      Segment IDs found on backward search of: sat382-rocks-merge.bsp
%      SAT375
%      DE-0431LE-0431
%      DE-0431LE-0431
%      DE-0431LE-0431
%      DE-0431LE-0431
%      PAN
%      DAPHNIS
%      PAN
%      DAPHNIS
%
%
%-Particulars
%
%   cspice_dafcs supports simultaneous searching of multiple DAFs. In
%   applications that use this capability, cspice_dafcs should be called
%   prior to each call to cspice_daffna, cspice_daffpa, cspice_dafgn, or
%   cspice_dafgs to specify which DAF is to be acted upon.
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
%         [sum]    = cspice_dafps( dc, ic );
%         name     = cspice_dafgn;
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
%         [sum]    = cspice_dafps( dc, ic );
%         name     = cspice_dafgn;
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
%            [sum]    = cspice_dafps( dc, ic );
%            name     = cspice_dafgn;
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
%            [sum]    = cspice_dafps( dc, ic );
%            name     = cspice_dafgn;
%
%             .
%             .
%
%            cspice_dafcs( handl2 );
%            found2 = cspice_daffpa;
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
%   1)  If the input handle is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If this routine is called when no search is in progress in the
%       the current DAF, the error SPICE(DAFNOSEARCH) is signaled by a
%       routine in the call tree of this routine.
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
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard.
%       Update code example to output segment IDs for both input SPK
%       files. Added example's problem statement.
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
%   -Mice Version 1.0.0, 10-JUL-2012 (EDW)
%
%-Index_Entries
%
%   select a DAF to continue searching
%
%-&

function cspice_dafcs( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafcs(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafcs_c', handle );
   catch spiceerr
      rethrow(spiceerr)
   end




