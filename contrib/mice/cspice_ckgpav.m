%-Abstract
%
%   CSPICE_CKGPAV returns pointing (attitude) and angular velocity
%   for a specified object at a user specified spacecraft clock time.
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
%      inst     NAIF ID for the instrument, spacecraft, or other structure for
%               which pointing is requested.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%               The frame fixed to this object is called the "instrument
%               frame" or "instrument-fixed" frame.
%
%      sclkdp   encoded spacecraft clock time(s) for which pointing is
%               requested.
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%      tol      time tolerance given in ticks (+/-), the units of encoded
%               spacecraft clock time, about `sclkdp'.
%
%               [1,1] = size(tol); double = class(tol)
%
%               The C-matrix returned by cspice_ckgpav, if any, is the one
%               whose time tag is closest to `sclkdp' and within `tol' units
%               of `sclkdp'.
%
%               In general, because using a non-zero tolerance affects
%               selection of the segment from which the data is obtained,
%               users are strongly discouraged from using a non-zero
%               tolerance when reading CKs with continuous data. Using
%               a non-zero tolerance should be reserved exclusively to
%               reading CKs with discrete data because in practice
%               obtaining data from such CKs using a zero tolerance is
%               often not possible due to time round off.
%
%      ref      naming the desired reference frame for the returned pointing.
%
%               [1,c1] = size(ref); char = class(ref)
%
%   the call:
%
%      [cmat, av, clkout, found] = cspice_ckgpav( inst, sclkdp, tol, ref )
%
%   returns:
%
%      cmat     rotation matrix(ces) that transform components of a vector
%               expressed in the frame specified by `ref' to components
%               expressed in the frame tied to the instrument, spacecraft, or
%               other structure at time(s) `clkout'.
%
%               If [1,1] = size(sclkdp) then [3,3]   = size(cmat)
%               If [1,n] = size(sclkdp) then [3,3,n] = size(cmat)
%                                             double = class(cmat)
%
%      av       angular velocity measured in radians per second (this is the
%               axis about which the reference frame tied to the instrument is
%               rotating in the right-handed sense at time `clkout').
%
%               [3,n] = size(av); double = class(av)
%
%      clkout   encoded spacecraft clock time(s) associated with the returned
%               C-matrix `cmat' (this value may differ from the requested
%               time, but never by more than the input tolerance `tol').
%
%               [1,n] = size(clkout); double = class(clkout)
%
%      found    the flag(s) indicating if the requested pointing is found.
%
%               [1,n] = size(found); logical = class(found)
%
%               `cmat', `av', `clkout', and `found' return with the same
%               vectorization measure (N) as `sclkdp'.
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
%   1) The following example code uses ckgpav_c to get C-matrices and
%      associated angular velocity vectors for a set of images whose
%      SCLK counts (un-encoded character string versions) are known.
%
%      For each C-matrix, a unit pointing vector is constructed and
%      printed along with the angular velocity vector.
%
%      Note: if the C-kernels of interest do not contain angular velocity
%      data, then the CSPICE routine cspice_ckgp should be used to read the
%      pointing data. An example program in the header of the Mice
%      function cspice_ckgp demonstrates this.
%
%      We need to load also an SCLK kernel to convert from clock string
%      to "ticks." Although not required for older spacecraft clocks,
%      most modern spacecraft ones require a leapseconds kernel to be
%      loaded in addition to an SCLK kernel.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: ckgpav_ex1.tm
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
%            04153_04182ca_ISS.bc   CASSINI image navigated
%                                   spacecraft CK
%
%
%         \begindata
%
%           KERNELS_TO_LOAD = ( 'cas00071.tsc'
%                               '04153_04182ca_ISS.bc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function ckgpav_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'ckgpav_ex1.tm' )
%
%         %
%         % The code for the Cassini spacecraft is -82
%         %
%         SC   = -82;
%
%         % The code for the Cassini spacecraft bus is -82000.
%         %
%         INST = -82000;
%
%         %
%         % The reference frame we want is J2000.
%         %
%         REF  = 'J2000';
%
%         %
%         % The CASSINI ISS camera boresight
%         % in the spacecraft frame is
%         % (0.0005760, -0.99999982, -0.0001710).
%         %
%         BORE   = [ 0.0005760; -0.99999982; -0.0001710];
%
%         %
%         % Spacecraft clock times for successive CASSINI
%         % navigation images always differ by more than 1.0 seconds.
%         % This is an acceptable tolerance, and must be
%         % converted to "ticks" (units of encoded `sclk') for
%         % input to cspice_ckgp.
%         %
%         TOL  = '1.0';
%
%         %
%         % Two CASSINI clock strings of interest.
%         %
%         SCLKCH =  strvcat( '1465644281.0', '1465644351.0' );
%
%         %
%         % Convert tolerance from CASSINI formatted character string
%         % SCLK to ticks, which are units of encoded SCLK.
%         %
%         toltik = cspice_sctiks( SC, TOL );
%
%         %
%         % cspice_ckgpav requires encoded spacecraft clock time.
%         %
%         sclkdp = cspice_scencd( SC, SCLKCH );
%
%         %
%         % Retrieve the 'REF' reference frame to 'INST' reference frame
%         % transformation matrix at time sclkdp with a tolerance
%         % 'toltik'.
%         %
%         %   [INST] = [cmat][ref]
%         %
%         [cmat, av, clkout, found] = cspice_ckgpav( INST,   sclkdp,  ...
%                                                    toltik, REF   );
%
%         for n=1:2
%
%            if( found(n) )
%
%               %
%               % Transform the 'BORE' vector from 'INST' reference frame to
%               % 'REF' frame.
%               %                T
%               %  [ref] = [cmat] [INST]
%               %
%               bore_ref = cmat(:,:,n)' * BORE;
%               [clkch] = cspice_scdecd( SC, clkout(n) );
%
%               fprintf( 'Requested SCLK time : %s\n', SCLKCH(n,:) );
%               fprintf( '   CASSINI SCLK time: %s\n', clkch       );
%               fprintf( ['   CASSINI ISS boresight  :',        ...
%                          '%10.7f %10.7f %10.7f\n' ], bore_ref    );
%               fprintf( ['   Angular velocity vector:',        ...
%                          '%10.7f %10.7f %10.7f\n\n'], av(:,n)    );
%
%            else
%
%               txt = sprintf( ['At CASSINI SCLK time: %s ',         ...
%                               'pointing not found'], SCLKCH(n,:) );
%               disp( txt )
%               disp( ' ' )
%
%            end
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
%      Requested SCLK time : 1465644281.0
%         CASSINI SCLK time: 1/1465644281.171
%         CASSINI ISS boresight  : 0.9376789  0.3444125  0.0462419
%         Angular velocity vector: 0.0000000  0.0000000  0.0000000
%
%      Requested SCLK time : 1465644351.0
%         CASSINI SCLK time: 1/1465644351.071
%         CASSINI ISS boresight  : 0.9376657  0.3444504  0.0462266
%         Angular velocity vector: 0.0000000  0.0000000  0.0000000
%
%
%-Particulars
%
%   How the tolerance argument is used
%   ==================================
%
%
%   Reading a type 1 CK segment (discrete pointing instances)
%   ---------------------------------------------------------
%
%   In the diagram below
%
%      - '0' is used to represent discrete pointing instances
%        (quaternions, angular velocity vectors, and associated
%        time tags).
%
%      - '( )' are used to represent the end points of the time
%        interval covered by a segment in a CK file.
%
%      - `sclkdp' is the time at which you requested pointing.
%        The location of `sclkdp' relative to the time tags of the
%        pointing instances is indicated by the '+' sign.
%
%      - `tol' is the time tolerance specified in the pointing
%        request. The square brackets '[ ]' represent the
%        endpoints of the time interval
%
%           sclkdp-tol : sclkdp+tol
%
%      - The quaternions occurring in the segment need not be
%        evenly spaced in time.
%
%
%   Case 1: pointing is available
%   ------------------------------
%
%                            sclkdp
%                                 \   tol
%                                  | /
%                                  |/\
%   Your request                [--+--]
%                               .  .  .
%   Segment      (0-----0--0--0--0--0--0---0--0------------0--0--0--0)
%                                   ^
%                                   |
%                       cspice_ckgpav returns this instance.
%
%
%   Case 2: pointing is not available
%   ----------------------------------
%
%                                                 sclkdp
%                                                    \   tol
%                                                     | /
%                                                     |/\
%   Your request                                   [--+--]
%                                                  .  .  .
%   Segment      (0-----0--0--0--0--0--0---0--0--0---------0--0--0--0)
%
%
%                       cspice_ckgpav returns no pointing; the output
%                       `found' flag is set to false.
%
%
%
%   Reading a type 2, 3, 4, or 5 CK segment (continuous pointing)
%   -------------------------------------------------------------
%
%   In the diagrams below
%
%      - '==' is used to represent periods of continuous pointing.
%
%      - '--' is used to represent gaps in the pointing coverage.
%
%      - '( )' are used to represent the end points of the time
%        interval covered by a segment in a CK file.
%
%      - `sclkdp' is the time at which you requested pointing.
%        The location of `sclkdp' relative to the time tags of the
%        pointing instances is indicated by the '+' sign.
%
%      - `tol' is the time tolerance specified in the pointing
%        request. The square brackets '[ ]' represent the
%        endpoints of the time interval
%
%           sclkdp-tol : sclkdp+tol
%
%      - The quaternions occurring in the periods of continuous
%        pointing need not be evenly spaced in time.
%
%
%   Case 1: pointing is available at the request time
%   --------------------------------------------------
%
%                           sclkdp
%                                 \   tol
%                                  | /
%                                  |/\
%   Your request                [--+--]
%                               .  .  .
%                               .  .  .
%                               .  .  .
%   Segment            (==---===========---=======----------===--)
%                                  ^
%                                  |
%
%                 The request time lies within an interval where
%                 continuous pointing is available. cspice_ckgpav returns
%                 pointing at the requested epoch.
%
%
%   Case 2: pointing is available 'near' the request time
%   ------------------------------------------------------
%
%                                  sclkdp
%                                        \   tol
%                                         | /
%                                         |/\
%   Your request                       [--+--]
%                                      .  .  .
%   Segment            (==---===========----=======---------===--)
%                                           ^
%                                           |
%
%                 The request time lies in a gap: an interval where
%                 continuous pointing is *not* available. cspice_ckgpav
%                 returns pointing for the epoch closest to the
%                 request time `sclkdp'.
%
%
%   Case 3: pointing is not available
%   ----------------------------------
%
%                                               sclkdp
%                                                     \   tol
%                                                      | /
%                                                      |/\
%   Your request                                    [--+--]
%                                                   .  .  .
%   Segment            (==---===========----=======---------===--)
%
%                       cspice_ckgpav returns no pointing; the output
%                       `found' flag is set to false.
%
%
%
%   Tolerance and segment priority
%   ==============================
%
%   cspice_ckgpav searches through loaded C-kernels to satisfy a pointing
%   request. Last-loaded files are searched first. Individual files are
%   searched in backwards order, so that between competing segments
%   (segments containing data for the same object, for overlapping time
%   ranges), the one closest to the end of the file has highest
%   priority. cspice_ckgpav considers only those segments that contain both
%   pointing and angular velocity data, as indicated by the segment
%   descriptor.
%
%   The search ends when a segment is found that can provide pointing
%   and angular velocity for the specified instrument at a time
%   falling within the specified tolerance on either side of the
%   request time. Within that segment, the instance closest to the
%   input time is located and returned.
%
%   The following four cases illustrate this search procedure. Segments
%   A and B are in the same file, with segment A located further
%   towards the end of the file than segment B. Both segments A and B
%   contain discrete pointing data, indicated by the number 0.
%
%
%   Case 1: Pointing is available in the first segment searched.
%            Because segment A has the highest priority and can
%            satisfy the request, segment B is not searched.
%
%
%                                sclkdp
%                                      \  tol
%                                       | /
%                                       |/\
%   Your request                     [--+--]
%                                    .  .  .
%   Segment A          (0-----------------0--------0--0-----0)
%                                         ^
%                                         |
%                                         |
%                             cspice_ckgpav returns this instance
%
%   Segment B     (0--0--0--0--0--0--0--0--0--0--0--0--0--0--0--0--0)
%
%
%
%   Case 2: Pointing is not available in the first segment searched.
%            Because segment A cannot satisfy the request, segment B
%            is searched.
%
%
%                           sclkdp
%                                \   tol
%                                 | /
%                                 |/\
%   Your request               [--+--]
%                              .  .  .
%   Segment A          (0-----------------0--------0--0-----0)
%                              .  .  .
%   Segment B     (0--0--0--0--0--0--0--0--0--0--0--0--0--0--0--0--0)
%                                 ^
%                                 |
%                     cspice_ckgpav returns this instance
%
%
%   Segments that contain continuous pointing data are searched in the
%   same manner as segments containing discrete pointing data. For
%   request times that fall within the bounds of continuous intervals,
%   cspice_ckgpav will return pointing at the request time. When the request
%   time does not fall within an interval, then a time at an endpoint of
%   an interval may be returned if it is the closest time in the segment
%   to the user request time and is also within the tolerance.
%
%   In the following examples, segment A is located further towards the
%   end of the file than segment C. Segment A contains discrete pointing
%   data and segment C contains continuous data, indicated by the '='
%   character.
%
%
%   Case 3: Pointing is not available in the first segment searched.
%            Because segment A cannot satisfy the request, segment C
%            is searched.
%
%                           sclkdp
%                                 \  tol
%                                  | /
%                                  |/\
%   Your request                [--+--]
%                               .  .  .
%                               .  .  .
%   Segment A          (0-----------------0--------0--0-----0)
%                               .  .  .
%                               .  .  .
%   Segment C          (---=============-----====--------==--)
%                                  ^
%                                  |
%                                  |
%                       cspice_ckgpav returns this instance
%
%
%   In the next case, assume that the order of segments A and C in the
%   file is reversed: A is now closer to the front, so data from
%   segment C are considered first.
%
%
%   Case 4: Pointing is available in the first segment searched.
%            Because segment C has the highest priority and can
%            satisfy the request, segment A is not searched.
%
%                                           sclkdp
%                                          /
%                                         |  tol
%                                         | /
%                                         |/\
%   Your request                       [--+--]
%                                      .  .  .
%                                      .  .  .
%   Segment C          (---=============-----====--------==--)
%                                           ^
%                                           |
%                              cspice_ckgpav returns this instance
%
%   Segment A          (0-----------------0--------0--0-----0)
%                                         ^
%                                         |
%                                   'Best' answer
%
%
%   The next case illustrates an unfortunate side effect of using
%   a non-zero tolerance when reading multi-segment CKs with
%   continuous data. In all cases when the look-up interval
%   formed using tolerance overlaps a segment boundary and
%   the request time falls within the coverage of the lower
%   priority segment, the data at the end of the higher priority
%   segment will be picked instead of the data from the lower
%   priority segment.
%
%
%   Case 5: Pointing is available in the first segment searched.
%            Because segment C has the highest priority and can
%            satisfy the request, segment A is not searched.
%
%                                           sclkdp
%                                          /
%                                         |  tol
%                                         | /
%                                         |/\
%   Your request                       [--+--]
%                                      .  .  .
%                                      .  .  .
%   Segment C                                (===============)
%                                            ^
%                                            |
%                              cspice_ckgpav returns this instance
%
%   Segment A          (=====================)
%                                         ^
%                                         |
%                                   'Best' answer
%
%-Exceptions
%
%   1)  If a C-kernel file has not been loaded using cspice_furnsh prior to
%       a call to this routine, an error is signaled by a routine in
%       the call tree of this routine.
%
%   2)  If `tol' is negative, found is set to false.
%
%   3)  If `ref' is not a supported reference frame, an error is
%       signaled by a routine in the call tree of this routine and
%       `found' is set to false.
%
%   4)  If any of the input arguments, `inst', `sclkdp', `tol' or
%       `ref', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   5)  If any of the input arguments, `inst', `sclkdp', `tol' or
%       `ref', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   cspice_ckgpav searches through files loaded by cspice_furnsh to locate a
%   segment that can satisfy the request for pointing and angular
%   velocity for instrument `inst' at time `sclkdp'. You must load a
%   C-kernel file using cspice_furnsh prior to calling this routine.
%
%-Restrictions
%
%   1)  Only loaded C-kernel segments containing both pointing and
%       angular velocity data will be searched by this reader.
%       Segments containing only pointing data will be skipped over.
%
%-Required_Reading
%
%   MICE.REQ
%   CK.REQ
%   SCLK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   B.V. Semenov        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.3.0, 26-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Updated example to
%       load the required kernels using meta-kernel and reformatted its
%       example output. Added problem statement. Modified input times and
%       kernel set to work with PDS archived CASSINI data. Added call to
%       cspice_kclear.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.2.0, 08-NOV-2012 (EDW) (SCK)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.1.1, 03-JUN-2010 (BVS)
%
%      Edits to header. Added warning regarding non-zero tolerance.
%
%   -Mice Version 1.1.0, 23-FEB-2009 (EDW)
%
%      Added zzmice_str call on input 'ref' to convert string cells to
%      character arrays if 'ref' has type string cells. Added proper
%      markers for usage string variable types.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   get CK pointing and angular velocity
%
%-&

function [cmat, av, clkout, found] = cspice_ckgpav(inst, sclkdp, tol, ref)

   switch nargin
      case 4

         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         tol    = zzmice_dp(tol);
         ref    = zzmice_str(ref);

      otherwise

         error ( [ 'Usage: [_cmat(3,3)_, _av(3)_, _clkout_, _found_] = ' ...
                   'cspice_ckgpav(inst, _sclkdp_, tol, `ref`)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cmat, av, clkout, found] = mice('ckgpav_c', inst, sclkdp, tol, ref);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end



