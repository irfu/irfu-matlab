%-Abstract
%
%   CSPICE_DPMAX returns the value of the largest (positive) number
%   representable in a double precision variable.
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
%      None.
%
%   the call:
%
%      [dpmax] = cspice_dpmax
%
%   returns:
%
%      dpmax    the value of the largest (positive) number that can be
%               represented in a double precision variable.
%
%               [1,1] = size(dpmax); double = class(dpmax)
%
%               This varies from machine to machine.
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
%   1) Return the value of the largest positive number representable
%      in a double precision variable.
%
%      Example code begins here.
%
%
%      function dpmax_ex1()
%
%         %
%         % Display the largest positive number representable
%         % in a double precision variable.
%         %
%         fprintf( 'Maximum Double Precision Number:  %21.14e\n',          ...
%                                                    cspice_dpmax )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Maximum Double Precision Number:  1.79769313486232e+308
%
%
%-Particulars
%
%   The function always returns a constant value, set by the user
%   prior to compilation.
%
%-Exceptions
%
%   Error free.
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
%
%-Literature_References
%
%   [1]  "Programming in VAX FORTRAN", Digital Equipment Corporation,
%        September 1984, Appendix C, FORTRAN Data Representation,
%        page C-3.
%
%   [2]  "Microsoft FORTRAN Reference", Microsoft Corporation
%        1989, Section 1.3.1, page 13.
%
%   [3]  "Sun FORTRAN Programmer's Guide", Sun Microsystems,
%        Revision A of 6 May 1988, Appendix F, Manual Pages for
%        FORTRAN, page 306 (LIBM_DOUBLE). This routine includes
%        the function D_MAX_NORMAL.
%
%   [4]  "FORTRAN Programmer's Guide for the Sun Workstation",
%        Sun Microsystems, Revision E of 17 February 1986,
%        Appendix E, Manual Pages for FORTRAN, page 184 (RANGE).
%        This routine includes the function DFLMAX.
%
%   [5]  "Language Systems FORTRAN Reference Manual", Language Systems
%        Corporation, version 1.2.1, page 3-12.
%
%   [6]  "Lahey F77L EM/32 Programmers Reference Manual", version 4.0,
%        page 95.
%
%   [7]  "FORTRAN/9000 Reference HP 9000 Series 700 Computers",
%        First Edition, June 1991, Hewlett Packard Company, page 4-5.
%
%   [8]  "SGI Fortran 77 Programmer's Guide", Document number
%        007-0711-030, page 2-2.
%
%   [9]  "Language Reference Manual", Absoft Fortran V3.2, 1993,
%        page 3-16, section 3.6.1.5. (for the NeXT)
%
%   [10] "Unix/VMS Compatibility Libraries", Absoft Fortran V3.2,
%        1993; Chapter 3, Support Libraries, page 3-5, dflmax. (for
%        the NeXT)
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   largest d.p. number
%
%-&
function [dpmax] = cspice_dpmax

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [dpmax] = cspice_dpmax' )

   end

   %
   % Call the MEX library.
   %
   try
      [dpmax] = mice('dpmax_c');
   catch spiceerr
      rethrow(spiceerr)
   end
