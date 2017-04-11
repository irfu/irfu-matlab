%-Abstract
%
%   CSPICE_CCIFRM returns the frame name, frame id, and center associated with
%   a given frame class and class id.
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
%     frclss   is the class or type of the frame. This identifies which
%              subsystem will be used to perform frame transformations.
%
%              [1,n] = size(frclss); int32 = class(frclss)
%
%     clssid   is the ID code used for the frame within its class. This
%              may be different from the frame ID code.
%
%              [1,n] = size(clssid); int32 = class(clssid)
%
%   the call:
%
%      [frcode, frname, center, found] = cspice_ccifrm( frclss, clssid )
%
%   returns:
%
%      frcode  is the frame ID code for the reference frame
%              identified by `frclss' and `clssid'.
%
%              [1,n] = size(frcode); int32 = class(frcode)
%
%      frname  is the name of the frame identified by
%              `frclss' and `clssid'.
%
%              [1,n] = size(frname); char = class(frname)
%
%              If `frname' does not have enough room to hold
%              the full name of the frame, the name will be truncated
%              on the right.
%
%      center  is the body ID code for the center of the reference
%              frame identified  by `frclss' and `clssid'.
%
%              [1,n] = size(center); int32 = class(center)
%
%      found   is true if a valid frame specification
%              corresponding to the input frame class and frame class
%              ID is available, in which case the other outputs are
%              valid. Otherwise, `found' is returned with the value
%              false.
%
%              [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      function ccifrm_t
%
%         frcode1 = cspice_namfrm( 'ITRF93');
%         [center1, clss, clss_ID, found]    = cspice_frinfo( frcode1 );
%         [frcode2,  frname, center2, found] = cspice_ccifrm( clss, clss_ID );
%
%         if ( ~found )
%
%               error('No joy' )
%
%         end
%
%         fprintf( 'Class     : %d\n', clss    )
%         fprintf( 'Class ID  : %d\n', clss_ID )
%         fprintf( 'Fame name : %s\n', frname  )
%         fprintf( 'Frame Code: %d\n', frcode2 )
%         fprintf( 'Center ID : %d\n', center2 )
%
%      end
%
%   Matlab outputs:
%
%      Class     : 2
%      Class ID  : 3000
%      Fame name : ITRF93
%      Frame Code: 13000
%      Center ID : 399
%
%-Particulars
%
%   This routine allows the user to determine the frame associated
%   with a given frame class and class ID code. The kernel pool is
%   searched first for a matching frame; if no match is found, then
%   the set of built-in frames is searched.
%
%   Since the neither the frame class nor the class ID are primary
%   keys, searching for matching frames is a linear (and therefore
%   typically slow) process.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ccifrm_c.
%
%   MICE.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-03-2017, ML (JPL), EDW (JPL)
%
%-Index_Entries
%
%   Find info associated with a frame class and class id
%   Map frame class and class id to frame info
%   Map frame class and class id to frame name, id, and center
%
%-&

function [frcode, frname, center, found] = cspice_ccifrm( frclss, clssid )

   switch nargin
      case 2

         frclss = zzmice_int(frclss);
         clssid = zzmice_int(clssid);

      otherwise

         error( ['Usage: [ frcode, `frname`, center, found] = ' ...
                'cspice_ccifrm( frclss, clssid)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ccifrm, center] = mice( 'ccifrm_s', frclss, clssid );
      frcode   = reshape( [ccifrm.code],  1, [] );
      frname   = char( ccifrm.name );
      found    = reshape( [ccifrm.found], 1, [] );
   catch
      rethrow(lasterror)
   end

