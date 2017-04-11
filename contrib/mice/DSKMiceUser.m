%-Abstract
%
%   DSKMiceUser.m declares global variables for use with
%   DSKMice APIs.
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
%   None.
%
%-Examples
%
%   The call:
%
%       DSKMiceUser
%
%   makes DSKMice parameter definitions available to the calling program.
%
%   Please see the example program in cspice_dski02.
%
%-Particulars
%
%   This file is an umbrella file that includes all of the required
%   variables required to support the DSKMice application programming interface
%   (API). Users' application code that calls DSKMice only needs to include
%   this single file.
%
%-Required Reading
%
%   None.
%
%-Version
%
%   -Mice Version 1.0.0, 10-MAR-2016, NJB (JPL), SCK (JPL), EDW (JPL)
%
%-Index_Entries
%
%   Include DSKMice parameters
%
%-&


%
% Include:
%    dsk type 2 specific definitions.
%    dsk tolerence parameters
%
   DSKMice02
   DSKtol

   %
   %API-specific definitions
   %========================
   %
   %Parameters for dskxsi_c:
   %
   SPICE_DSKXSI_DCSIZE  =          1;
   SPICE_DSKXSI_ICSIZE  =          1;

%
% End of parameter declaration file DSKMiceUser.pro
%
