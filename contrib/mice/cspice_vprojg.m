%-Abstract
%
%   CSPICE_VPROJG computes the projection of one vector onto another vector.
%   All vectors are of arbitrary dimension.
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
%      a        a double precision vector of arbitrary dimension.
%
%               [n,1] = size(a); double = class(a)
%
%               This vector is to be projected onto the vector `b'.
%
%      b        a double precision vector of arbitrary dimension.
%
%               [n,1] = size(b); double = class(b)
%
%               This vector is the vector which receives the projection.
%
%   the call:
%
%      [p] = cspice_vprojg( a, b )
%
%   returns:
%
%      p        a double precision vector of arbitrary dimension containing
%               the projection of `a' onto `b'.
%
%               [n,1] = size(p); double = class(p)
%
%               (`p' is necessarily parallel to `b'.)
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
%   1) Define two sets of vectors and compute the projection of
%      each vector of the first set on the corresponding vector of
%      the second set.
%
%      Example code begins here.
%
%
%      function vprojg_ex1()
%
%         %
%         % Local parameters.
%         %
%         SETSIZ = 4;
%
%         %
%         % Define the two vector sets.
%         %
%         seta = [ [6.0,  6.0,  6.0,  0.0]',                               ...
%                  [6.0,  6.0,  6.0,  0.0]',                               ...
%                  [6.0,  6.0,  0.0,  0.0]',                               ...
%                  [6.0,  0.0,  0.0,  0.0]' ];
%
%         setb = [ [2.0,  0.0,  0.0,  0.0]',                               ...
%                  [-3.0,  0.0,  0.0,  0.0]',                              ...
%                  [0.0,  7.0,  0.0,  0.0]',                               ...
%                  [0.0,  0.0,  9.0,  0.0]' ];
%
%         %
%         % Calculate the projection
%         %
%         for i=1:SETSIZ
%            [pvec] = cspice_vprojg( seta(:,i), setb(:,i) );
%            fprintf( 'Vector A  :  %4.1f %4.1f %4.1f %4.1f\n',            ...
%                     seta(1,i), seta(2,i), seta(3,i), seta(4,i) )
%            fprintf( 'Vector B  :  %4.1f %4.1f %4.1f %4.1f\n',            ...
%                     setb(1,i), setb(2,i), setb(3,i), setb(4,i) )
%            fprintf( 'Projection:  %4.1f %4.1f %4.1f %4.1f\n',            ...
%                            pvec(1), pvec(2), pvec(3), pvec(4) )
%            fprintf( ' \n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Vector A  :   6.0  6.0  6.0  0.0
%      Vector B  :   2.0  0.0  0.0  0.0
%      Projection:   6.0  0.0  0.0  0.0
%
%      Vector A  :   6.0  6.0  6.0  0.0
%      Vector B  :  -3.0  0.0  0.0  0.0
%      Projection:   6.0 -0.0 -0.0 -0.0
%
%      Vector A  :   6.0  6.0  0.0  0.0
%      Vector B  :   0.0  7.0  0.0  0.0
%      Projection:   0.0  6.0  0.0  0.0
%
%      Vector A  :   6.0  0.0  0.0  0.0
%      Vector B  :   0.0  0.0  9.0  0.0
%      Projection:   0.0  0.0  0.0  0.0
%
%
%-Particulars
%
%   The projection of a vector `a' onto a vector `b' is, by definition,
%   that component of `a' which is parallel to `b'. To find this
%   component it is enough to find the scalar ratio of the length of
%   `b' to the projection of `a' onto `b', and then use this number to
%   scale the length of `b'. This ratio is given by
%
%       ratio = < a, b > / < b, b >
%
%   where <,> denotes the general vector dot product. This routine
%   does not attempt to divide by zero in the event that `b' is the
%   zero vector.
%
%-Exceptions
%
%   1)  If any of the input arguments, `a' or `b', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `a' or `b', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   3)  If the input vector arguments `a' and `b' do not have the same
%       dimension (n), an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  No error detection or recovery schemes are incorporated into
%       this routine except to insure that no attempt is made to
%       divide by zero. Thus, the user is required to make sure that
%       the vectors `a' and `b' are such that no floating point overflow
%       will occur when the dot products are calculated.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  G. Thomas and R. Finney, "Calculus and Analytic Geometry,"
%        7th Edition, Addison Wesley, 1988.
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
%   n-dimensional vector projection
%
%-&
function [p] = cspice_vprojg( a, b )

   switch nargin
      case 2

         a = zzmice_dp(a);
         b = zzmice_dp(b);

      otherwise

         error ( 'Usage: [p(n)] = cspice_vprojg( a(n), b(n) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [p] = mice('vprojg_c', a, b);
   catch spiceerr
      rethrow(spiceerr)
   end
