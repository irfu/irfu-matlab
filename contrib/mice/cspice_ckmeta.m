%-Abstract
%
%   CSPICE_CKMETA returns (depending upon the user's request) the ID code of
%   either the spacecraft or spacecraft clock associated with a C-Kernel ID
%   code.
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
%      ckid     the ID code for some object whose attitude and possibly
%               angular velocity are stored in some C-kernel.
%
%               [1,1] = size(ckid); int32 = class(ckid)
%
%      meta     a character string that indicates which piece of meta data to
%               fetch.
%
%               [1,c1] = size(meta); char = class(meta)
%
%                  or
%
%               [1,1] = size(meta); cell = class(meta)
%
%               Acceptable values are 'SCLK' and 'SPK'. The routine is case
%               insensitive. Leading and trailing blanks are insignificant.
%               However, blanks between characters are regarded as being
%               significant.
%
%   the call:
%
%      [idcode] = cspice_ckmeta( ckid, meta )
%
%   returns:
%
%      idcode   if `meta' is 'SCLK' then the value returned in `idcode' is
%               the ID code of the spacecraft clock used for converting ET to
%               TICKS and TICKS to ET for the C-kernel used to represent the
%               attitude of the object with ID code `ckid'.
%
%               [1,1] = size(idcode); int32 = class(idcode)
%
%               If `meta' is 'SPK' then the value returned in `idcode' is the
%               ID code of the spacecraft on which the platform indicated
%               by `ckid' is mounted.
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
%   1) Suppose you would like to look up the attitude of an object
%      in a C-kernel but have ET and seconds as your input time and
%      tolerance.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: ckmeta_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name              Contents
%            --------------------   -----------------------
%            cas00071.tsc           CASSINI SCLK
%            naif0012.tls           Leapseconds
%            04153_04182ca_ISS.bc   CASSINI image navigated
%                                   spacecraft CK
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'cas00071.tsc'
%                                '04153_04182ca_ISS.bc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function ckmeta_ex1()
%
%         %
%         % Local parameters.
%         %
%         % -- The code for CASSINI spacecraft reference frame is
%         %    -82000.
%         %
%         % -- The reference frame we want is J2000.
%         %
%         REF  =   'J2000';
%         CKID =   -82000;
%
%         %
%         % Initial values.
%         %
%         et     = [141162208.034340]';
%         sectol = [0.5]';
%
%         %
%         % First load the CK, LSK and SCLK files.
%         %
%         cspice_furnsh( 'ckmeta_ex1.tm' );
%
%         %
%         % Get the SCLK identifier of the spacecraft clock required
%         % to convert from `et' to `ticks'.
%         %
%         [idcode] = cspice_ckmeta( CKID, 'SCLK' );
%
%         %
%         % Convert `et' and et+sectol to spacecraft clock ticks.
%         %
%         [ticks] = cspice_sce2c( idcode, et );
%         [tick2] = cspice_sce2c( idcode, et+sectol );
%
%         %
%         % Compute the tolerance in spacecraft clock ticks.
%         %
%         tol = tick2 - ticks;
%
%         %
%         % Look the attitude up.
%         %
%         [cmat, av, clkout, found] = cspice_ckgpav( CKID, ticks, tol, REF );
%
%         fprintf( 'Input ET:             %19.6f\n', et )
%
%         if ( found )
%
%            [etout] = cspice_sct2e( idcode, clkout );
%            fprintf( 'Attitude found at ET: %19.6f\n', etout )
%
%         else
%
%            fprintf( 'No attitude found at ET.\n' )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Input ET:                141162208.034340
%      Attitude found at ET:    141162208.034586
%
%
%-Particulars
%
%   This is a utility routine for mapping C-kernels to associated
%   spacecraft clocks.
%
%   An association of an SCLK ID and spacecraft ID with a CK frame
%   class ID may be made by placing in a text kernel the kernel
%   variable assignments
%
%      CK_<ck_frame_class_ID_code>_SCLK = <ID code of SCLK>
%      CK_<ck_frame_class_ID_code>_SPK  = <SPK ID code>
%
%   See the Frames Required Reading section on CK frames.
%
%-Exceptions
%
%   1)  If the variable `meta' is not recognized to be one of the
%       inputs 'SPK' or 'SCLK', the error SPICE(UNKNOWNCKMETA)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If `ckid' is greater than -1000, the associated SCLK and SPK
%       IDs must be in the kernel pool. If they are not present
%       a value of zero is returned for the requested item. Zero
%       is never the valid ID of a spacecraft clock.
%
%   3)  If any of the input arguments, `ckid' or `meta', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   4)  If any of the input arguments, `ckid' or `meta', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
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
%   CK.REQ
%   FRAMES.REQ
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
%   -Mice Version 1.0.0, 21-JUN-2021 (JDR)
%
%-Index_Entries
%
%   Map C-kernel ID to SCLK and SPK ID
%
%-&
function [idcode] = cspice_ckmeta( ckid, meta )

   switch nargin
      case 2

         ckid = zzmice_int(ckid);
         meta = zzmice_str(meta);

      otherwise

         error ( 'Usage: [idcode] = cspice_ckmeta( ckid, `meta` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [idcode] = mice('ckmeta_c', ckid, meta);
   catch spiceerr
      rethrow(spiceerr)
   end
