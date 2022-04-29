%-Abstract
%
%   MiceFrm.m performs Mice definitions for frame subsystem APIs.
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
%-Parameters
%
%   None.
%
%-Examples
%
%   Include these definitions by using the call:
%
%      MiceUser;
%
%-Particulars
%
%   This file defines constants that may be referenced in
%   application code that calls Mice frame subsystem APIs.
%
%
%      CONSTANTS
%      ==========
%
%
%         Frame counts
%         ------------
%
%         The following parameter are counts of built-in frames. These
%         parameters correspond to those defined in the SPICELIB Fortran
%         INCLUDE files
%
%            ninert.inc
%            nninrt.inc
%
%         and in the CSPICE header file
%
%            SpiceFrm.h
%
%
%         Name                   Description
%         ----                   ----------
%
%         SPICE_NFRAME_NINERT    Number of built-in inertial frames.
%         SPICE_NFRAME_NNINRT    Number of built-in non-inertial frames.
%
%
%
%         Frame classes
%         -------------
%
%         The following parameters identify SPICE frame classes. These
%         parameters correspond to those defined in the SPICELIB Fortran
%         INCLUDE file frmtyp.inc and in the CSPICE header file SpiceFrm.h
%
%         See the Frames Required Reading for a detailed discussion of
%         frame classes.
%
%
%         Name                   Description
%         ----                   ----------
%
%         SPICE_FRMTYP_INERTL    an inertial frame that is built-in into
%                                the SPICE frame subsystem and that requires
%                                no external file to compute the
%                                transformation from or to any other
%                                inertial frame.
%
%
%         SPICE_FRMTYP_PCK       a frame that is specified relative
%                                to some built-in, inertial frame (of
%                                class SPICE_FRMTYP_INERTL) and that has
%                                an IAU model that may be retrieved from
%                                the PCK system via a call to the
%                                CSPICE routine tisbod_c.
%
%
%         SPICE_FRMTYP_CK        a frame defined by a C-kernel.
%
%
%         SPICE_FRMTYP_TK        a "text kernel" frame. These frames
%                                are offset from their associated
%                                "relative" frames by a constant
%                                rotation.
%
%
%         SPICE_FRMTYP_DYN       a "dynamic" frame. These currently
%                                are limited to parameterized frames
%                                where the full frame definition depends
%                                on parameters supplied via a frame
%                                kernel.
%
%         SPICE_FRMTYP_SWTCH     is a "switch" frame. These frames have
%                                orientation defined by their alignment with
%                                base frames selected from a prioritized list.
%                                The base frames optionally have associated
%                                time intervals of applicability.
%
%         SPICE_FRMTYP_ALL       indicates any of the above classes.
%                                This parameter is used in APIs that
%                                fetch information about frames of a
%                                specified class.
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
%   CK.REQ
%   FRAMES.REQ
%   PCK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 26-AUG-2021 (JDR) (NJB)
%
%-Index_Entries
%
%   Include frame subsystem definitions parameters
%
%-&


%
% Frame counts:
%

%
% Number of built-in inertial frames. This number must be kept in
% sync with that defined in the SPICELIB include file ninert.inc.
%
   SPICE_NFRAME_NINERT  =  21;

%
% Number of built-in non-inertial frames. This number must be kept in
% sync with that defined in the SPICELIB include file nninrt.inc.
%
   SPICE_NFRAME_NNINRT  =  124;


%
% The frame class codes defined here are identical
% to those used in SPICELIB.
%

%
% Inertial, built-in frames:
%
   SPICE_FRMTYP_INERTL  =  1;

%
% PCK frames:
%
   SPICE_FRMTYP_PCK     =  2;

%
% CK frames:
%
   SPICE_FRMTYP_CK      =  3;

%
% TK frames:
%
   SPICE_FRMTYP_TK      =  4;

%
% Dynamic frames:
%
   SPICE_FRMTYP_DYN     =  5;

%
% Switch frames:
%
   SPICE_FRMTYP_SWTCH   =  6;

%
% All frame classes:
%
   SPICE_FRMTYP_ALL     = -1;

