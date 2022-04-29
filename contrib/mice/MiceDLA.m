%-Abstract
%
%   MiceDLA.m declares DLA-specific constants.
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
%
%   This header parameters that may be referenced in application code that
%   calls DLA functions.
%
%      Dimensions
%      ----------
%
%         Name                  Description
%         ----                  -----------
%         SPICE_DLA_DSCSIZ      Size of a SPICELIB DLA descriptor.
%
%
%      DLA File Offsets
%      ----------------
%
%         These parameters are provided to support CSPICE wrapper testing.
%
%         Name                  Description
%         ----                  -----------
%         SPICE_DLA_VERIDX      DAS integer address of DLA version code.
%
%         SPICE_DLA_LLBIDX      DAS integer addresses of first segment linked
%                               list pointer.
%
%         SPICE_DLA_LLEIDX      DAS integer addresses of last segment linked
%                               list pointer.
%
%
%      Structure Offsets
%      -------------------------
%
%         Name                  Description
%         ----                  ----------
%         SPICE_DLA_BWDIDX      Backward pointer index in a DLA
%                               descriptor.
%
%         SPICE_DLA_FWDIDX      Forward pointer index in a DLA
%                               descriptor.
%
%         SPICE_DLA_IBSIDX      Integer base address index in a 
%                               DLA descriptor.
%
%         SPICE_DLA_ISZIDX      Integer component size index in a
%                               DLA descriptor.
%
%         SPICE_DLA_DBSIDX      D.p. base address index in a DLA
%                               descriptor.
%
%         SPICE_DLA_DSZIDX      D.p. component size index in a 
%                               DLA descriptor.
%
%         SPICE_DLA_CBSIDX      Character base address index in a
%                               DLA descriptor.
%
%         SPICE_DLA_CSZIDX      Character component size index in a
%                               DLA descriptor.
%
%
%      Other DLA parameters
%      --------------------
%
%         Name                  Description
%         ----                  -----------
%         SPICE_DLA_NULPTR      Null pointer parameter.
%
%         SPICE_DLA_FMTVER      DLA format version.
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
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-SEP-2021, (JDR)
%
%       Added "DLA File Offsets" and "Other DLA parameters".
%
%       Updated the example section with the correct Mice call. Corrected
%       Index_Entries entry.
%
%       Updated the header for compliance with NAIF standard. Added
%       -Parameters, -Exceptions, -Files, -Restrictions, -Literature_References
%       and -Author_and_Institution sections.
%
%   -Mice Version 1.0.0, 26-MAR-2017 (NJB) (EDW)
%
%-Index_Entries
%
%   Include DLA parameters
%
%-&


%
% DAS integer address of DLA version code.
%
   SPICE_DLA_VERIDX  =              1;

%
% Linked list parameters
%
% Logical arrays (aka "segments") in a DAS linked array (DLA) file
% are organized as a doubly linked list. Each logical array may
% actually consist of character, double precision, and integer
% components. A component of a given data type occupies a
% contiguous range of DAS addresses of that type. Any or all
% array components may be empty.
%
% The segment descriptors in a SPICE DLA (DAS linked array) file
% are connected by a doubly linked list. Each node of the list is
% represented by a pair of integers acting as forward and backward
% pointers. Each pointer pair occupies the first two integers of a
% segment descriptor in DAS integer address space. The DLA file
% contains pointers to the first integers of both the first and
% last segment descriptors.
%
% At the DLA level of a file format implementation, there is
% no knowledge of the data contents. Hence segment descriptors
% provide information only about file layout (in contrast with
% the DAF system). Metadata giving specifics of segment contents
% are stored within the segments themselves in DLA-based file
% formats.
%
%
% Parameter declarations follow.
%
% DAS integer addresses of first and last segment linked list
% pointer pairs. The contents of these pointers
% are the DAS addresses of the first integers belonging
% to the first and last link pairs, respectively.
%
% The acronyms "LLB" and "LLE" denote "linked list begin"
% and "linked list end" respectively.
%
   SPICE_DLA_LLBIDX  =              SPICE_DLA_VERIDX + 1;
   SPICE_DLA_LLEIDX  =              SPICE_DLA_LLBIDX + 1;

%
% Null pointer parameter.
%
   SPICE_DLA_NULPTR  =              -1;

%
% DLA descriptor dimension:
%
   SPICE_DLA_DSCSIZ  =              8;

%
% DLA descriptor index parameters:
%
   SPICE_DLA_BWDIDX  =              1;
   SPICE_DLA_FWDIDX  =              2;
   SPICE_DLA_IBSIDX  =              3;
   SPICE_DLA_ISZIDX  =              4;
   SPICE_DLA_DBSIDX  =              5;
   SPICE_DLA_DSZIDX  =              6;
   SPICE_DLA_CBSIDX  =              7;
   SPICE_DLA_CSZIDX  =              8;

%
% Other DLA parameters:
%
%
% DLA format version. (This number is expected to occur very
% rarely at integer address SPICE_DLA_VERIDX in uninitialized DLA files.)
%
   SPICE_DLA_FMTVER  =              1000000;
