%-Abstract
%
%   CSPICE_DASOPR opens a DAS file for reading.
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
%      fname    is the name of a DAS file to be opened with read
%               access.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      [handle] = cspice_dasopr( fname )
%
%   returns:
%
%      handle   is the handle that is  associated with the file. This
%               handle is used to identify the file in subsequent
%               calls to other DAS routines.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%-Examples
%
%   Refer to example in cspice_dskb02.
%
%-Particulars
%
%   Most DAS files require only read access. If you do not need to
%   change the contents of a file, you should open it using cspice_dasopr.
%
%-Required Reading
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 28-APR-2014, NJB (JPL), EDW (JPL)
%
%-Index_Entries
%
%   open a DAS file for reading
%   open a DAS file for read access
%
%-&

function [handle] = cspice_dasopr( fname )

   switch nargin
      case 1

         fname  = zzmice_str(fname);

      otherwise

         error ( 'Usage: [handle] = cspice_dafopr(`fname`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'dasopr_c', fname );
   catch
      rethrow(lasterror)
   end




