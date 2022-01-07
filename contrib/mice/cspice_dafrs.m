%-Abstract
%
%   CSPICE_DAFRS changes the summary for the current array in the current
%   DAF.
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
%      dc       the double precision components of the summary.
%
%               [1,nd] = size(dc); double = class(dc)
%
%      ic       the integer components of the summary.
%
%               [1,ni] = size(ic); int32 = class(ic)
%
%   the call:
%
%      cspice_dafrs( dc, ic )
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
%   1) Replace the body ID code 301 (Moon) with a test body ID,
%      e.g. -999, in every descriptor of an SPK file.
%
%
%      Example code begins here.
%
%
%      function dafrs_ex1()
%
%         %
%         % Local parameters.
%         %
%         MAXOBJ  = 1000;
%
%         ND      = 2;
%         NI      = 6;
%
%         NEWCODE = ( -999 );
%         OLDCODE = ( 301 );
%
%         %
%         % Get the SPK file name.
%         %
%         fname = input( 'Enter name of the SPK file > ', 's' );
%
%         %
%         % Open for writing the SPK file.
%         %
%         [handle] = cspice_dafopw( fname );
%
%         %
%         % Search the file in forward order.
%         %
%         cspice_dafbfs( handle );
%         [found] = cspice_daffna;
%
%         while found
%
%            %
%            % Fetch and unpack the descriptor (aka summary)
%            % of the current segment.
%            %
%            [dc, ic] = cspice_dafgs( ND, NI );
%
%            %
%            % Replace ID codes if necessary.
%            %
%            if ( ic(1) == OLDCODE )
%               ic(1) = NEWCODE;
%            end
%            if ( ic(2) == OLDCODE )
%               ic(2) = NEWCODE;
%            end
%
%            %
%            % Re-pack the descriptor; replace the descriptor
%            % in the file.
%            %
%            cspice_dafrs( dc, ic );
%
%            %
%            % Find the next segment.
%            %
%            [found] = cspice_daffna;
%
%         end
%
%         %
%         % Close the file.
%         %
%         cspice_dafcls( handle );
%
%         %
%         % Find the set of objects in the SPK file.
%         %
%         [ids] = cspice_spkobj( fname, MAXOBJ );
%
%         fprintf( 'Objects in the DAF file:\n' )
%         fprintf( '\n' )
%         for i=1:numel( ids )
%
%            obj =  ids(i);
%            fprintf( '%4d', obj )
%
%         end
%
%         fprintf( '\n' )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the SPK file named de430.bsp, the output was:
%
%
%      Enter name of the SPK file > de430.bsp
%      Objects in the DAF file:
%
%      -999   1   2   3   4   5   6   7   8   9  10 199 299 399
%
%
%-Particulars
%
%   This function writes the contents of `dc' and `ic' to the current
%   DAF segment summary.
%
%   A single call to cspice_dafrs equates to the C calls:
%
%      dafps_c ( nd, ni, dc, ic, sum );
%      dafrs_c ( sum );
%
%   without use of the `sum' variable.
%
%   The summary of the current array in the current DAF can
%   be updated by providing new ones through cspice_dafrs. This feature
%   should not be used except to correct errors that occurred during
%   the creation of a file. Note that changes can only be made to
%   files opened for write access.
%
%-Exceptions
%
%   1)  If `nd' is zero or negative, no double precision components are
%       stored.
%
%   2)  If `ni' is zero or negative, no integer components are stored.
%
%   3)  If the total size of the summary, computed as
%
%               (ni - 1)
%          nd + -------- + 1
%                   2
%
%       is greater than 125 double precision words, some components may
%       not be stored.
%
%   4)  If this routine is called when no search is in progress in the
%       the current DAF, the error SPICE(DAFNOSEARCH) is signaled by a
%       routine in the call tree of this routine.
%
%   5)  If the DAF containing the "current" array has actually been
%       closed, an error is signaled by a routine in the call tree of
%       this routine.
%
%   6)  If the DAF containing the "current" array is not open for
%       writing, an error is signaled by a routine in the call tree of
%       this routine.
%
%   7)  If no array is current in the current DAF, the error
%       SPICE(NOCURRENTARRAY) is signaled by a routine in the call tree
%       of this routine. There is no current array when a search is
%       started by cspice_dafbfs or cspice_dafbbs, but no calls to
%       cspice_daffna or cspice_daffpa have been made yet, or whenever
%       cspice_daffna or cspice_daffpa return the value False.
%
%   8)  If any of the input arguments, `nd', `ni', `dc' or `ic', is
%       undefined, an error is signaled by the Matlab error handling system.
%
%   9)  If any of the input arguments, `nd', `ni', `dc' or `ic', is not
%       of the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   DAF.REQ
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
%   -Mice Version 1.0.0, 26-NOV-2021 (JDR)
%
%-Index_Entries
%
%   replace DAF summary
%
%-&
function cspice_dafrs( dc, ic )

   switch nargin
      case 2

         dc  = zzmice_dp(dc);
         ic  = zzmice_int(ic);

      otherwise

         error ( 'Usage: cspice_dafrs( dc(nd), ic(ni) )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dafrs_c', dc, ic);
   catch spiceerr
      rethrow(spiceerr)
   end
