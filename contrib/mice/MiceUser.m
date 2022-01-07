%-Abstract
%
%   MiceUser.m declares global parameter definitions. 
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
%   None.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   In order to include all Mice parameter definitions, use the call:
%
%      MiceUser;
%
%   from the user's application code.
%
%-Particulars
%
%   This file is an umbrella file that includes all of the parameters
%   required to support the Mice application programming interface
%   (API). Users' application code that calls Mice only needs to include
%   this single file. 
%
%   The Mice definition files included in this file are not part of 
%   the Mice API. The set of files included may change without notice.
%   Users should not include these files directly in their own
%   application code.
%
%-Exceptions
%
%   None.
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
%   None.
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.3.0, 11-JUN-2020 (JDR)
%
%       Added MiceDAS, MiceDtl, MiceFrm, MiceGF, and MiceOsc for
%       include.
%
%       Renamed DLAMice and DSKMiceUser to MiceDLA and MiceDSK.
%
%       Updated the header for compliance with NAIF standard. Added
%       -Parameters, -Exceptions, -Files, -Restrictions, -Literature_References
%       and -Author_and_Institution sections.
%
%   -Mice Version 1.2.0, 30-MAR-2017 (EDW)
%
%       Added DLAMice for include.
%
%   -Mice Version 1.1.0, 11-MAR-2016 (EDW)
%
%       Added DSKMiceUser for include.
%
%   -Mice Version 1.0.0, 14-FEB-2012 (SCK)
%
%-Index_Entries
%
%   include global parameters
%   global parameter definitions 
%
%-&


%
% Include needed definitions.
%
   MiceDAS;
   MiceDLA;
   MiceDSK;
   MiceDtl;
   MiceFrm;
   MiceGF;
   MiceOccult;
   MiceOsc;
