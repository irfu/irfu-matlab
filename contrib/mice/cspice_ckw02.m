%-Abstract
%
%   CSPICE_CKW02 adds a type 2 segment to a C-kernel.
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
%      handle   file handle for an open CK file, returned from cspice_ckopn.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      begtim   encoded SCLK segment begin time.
%
%               [1,1] = size(begtim); double = class(begtim)
%
%      endtim   encoded SCLK segment end time.
%
%               [1,1] = size(endtim); double = class(endtim)
%
%      inst     NAIF instrument ID code.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%      ref      name of the reference frame for the segment.
%
%               [1,c1] = size(ref); char = class(ref)
%
%                  or
%
%               [1,1] = size(ref); cell = class(ref)
%
%      segid    name to identify the segment.
%
%               [1,c2] = size(segid); char = class(segid)
%
%                  or
%
%               [1,1] = size(segid); cell = class(segid)
%
%      start    an array containing encoded SCLK interval start times.
%
%               [n,1] = size(start); double = class(start)
%
%      stop     an array containing the encoded SCLK interval stop times.
%
%               [n,1] = size(stop); double = class(stop)
%
%      quats    array of SPICE style quaternions representing instrument
%               pointing.
%
%               [4,n] = size(quats); double = class(quats)
%
%      avvs     array of angular velocity vectors in units of radians per
%               second.
%
%               [3,n] = size(avvs); double = class(avvs)
%
%      rates    the number of seconds per tick for each interval.
%
%               [n,1] = size(rates); double = class(rates)
%
%   the call:
%
%      cspice_ckw02( handle, ...
%                    begtim, ...
%                    endtim, ...
%                    inst,   ...
%                    ref,    ...
%                    segid,  ...
%                    start,  ...
%                    stop,   ...
%                    quats,  ...
%                    avvs,   ...
%                    rates )
%
%   returns:
%
%      None.
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
%   1) The following example creates a CK file with a type-2 segment,
%      with data for a simple time dependent rotation and angular
%      velocity.
%
%      Example code begins here.
%
%
%      function ckw02_ex1()
%
%         %
%         % Define needed parameters
%         %
%         CK2        = 'ckw02_ex1.bc';
%         INST       = -77702;
%         MAXREC     = 21;
%         SECPERTICK = 0.001;
%         IFNAME     = 'Test CK type 2 segment created by cspice_ckw02';
%         SEGID      = 'Test type 2 CK segment';
%
%         %
%         % NCOMCH is the number of characters to reserve for the kernel's
%         % comment area. This example doesn't write comments, so set to
%         % zero.
%         %
%         NCOMCH = 0;
%
%         %
%         % The base reference from for the rotation data.
%         %
%         REF = 'J2000';
%
%         %
%         % Time spacing in encoded ticks.
%         %
%         SPACING_TICKS = 10.;
%
%         %
%         % Time spacing in seconds
%         %
%         SPACING_SECS = SPACING_TICKS * SECPERTICK;
%
%         %
%         % Declare an angular rate in radians per sec.
%         %
%         RATE = 1.e-2;
%
%         %
%         % Create a 4xMAXREC matrix for quaternions, and a
%         % 3xMAXREC for expavs.
%         %
%         quats = zeros( [4, MAXREC] );
%         av    = zeros( [3, MAXREC] );
%
%         %
%         % Create a 3x3 double precision identity matrix.
%         %
%         work_mat = eye( 3 );
%
%         %
%         % Convert the `work_mat' to quaternion.
%         %
%         work_quat = cspice_m2q( work_mat);
%
%         %
%         % Copy the work quaternion to the first row of
%         % `quats'.
%         %
%         quats(:,1) = work_quat;
%
%         %
%         % Create an angular velocity vector. Copy to the third (Z) row
%         % of `av'. This vector is in the `REF' reference frame and
%         % indicates a constant rotation about the Z axis.
%         %
%         av(3,:) = RATE;
%
%         %
%         % Create arrays of interval start and stop times.  The interval
%         % associated with each quaternion will start at the epoch of
%         % the quaternion and will extend 0.8 * SPACING_TICKS forward in
%         % time, leaving small gaps between the intervals.
%         %
%         % Fill in the clock rates array with a constant `SECPERTICK' for
%         % all values.
%         %
%         rates  = zeros( [MAXREC,1] ) + SECPERTICK;
%
%         %
%         % Create an array of encoded tick values in increments of
%         % `SPACING_TICKS' with an initial value of 1000 ticks.
%         %
%         sclkdp = [0:MAXREC-1]' * SPACING_TICKS + 1000;
%
%         starts = sclkdp;
%         stops  = sclkdp + ( 0.8 * SPACING_TICKS );
%
%         %
%         % Fill the rest of the `av' and `quats' matrices
%         % with simple data.
%         %
%         for i = 2:MAXREC
%
%            %
%            % Create the transformation matrix for a rotation of `theta'
%            % about the Z axis. Calculate `theta' from the constant
%            % angular rate RATE at increments of SPACING_SECS.
%            %
%            theta    = (i-1) * RATE * SPACING_SECS;
%            work_mat = cspice_rotmat( work_mat, theta, 3 );
%
%            %
%            % Convert the `work_mat' matrix to SPICE type quaternion.
%            %
%            work_quat = cspice_m2q( work_mat );
%
%            %
%            % Store the quaternion in the `quats' matrix.
%            %
%            quats(:,i) = work_quat;
%
%         end
%
%         %
%         % Set the segment boundaries equal to the first and last
%         % time in the segment.
%         %
%         begtim = starts(1);
%         endtim = stops(MAXREC);
%
%         %
%         % All information ready to write. Write to a CK type 2 segment
%         % to the file indicated by `handle'.
%         %
%         try
%            handle = cspice_ckopn( CK2, IFNAME, NCOMCH );
%            cspice_ckw02(  handle, ...
%                           begtim, ...
%                           endtim, ...
%                           INST,   ...
%                           REF,    ...
%                           SEGID,  ...
%                           starts, ...
%                           stops,  ...
%                           quats,  ...
%                           av,     ...
%                           rates )
%         catch
%            error( [ 'Failure: ' lasterr] )
%         end
%
%         %
%         % SAFELY close the file.
%         %
%         cspice_ckcls( handle )
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new CK file exists in the
%      output directory.
%
%-Particulars
%
%   For a detailed description of a type 2 CK segment please see the
%   CK Required Reading.
%
%   This routine relieves the user from performing the repetitive
%   calls to the DAF routines necessary to construct a CK segment.
%
%
%   Quaternion Styles
%   -----------------
%
%   There are different "styles" of quaternions used in
%   science and engineering applications. Quaternion styles
%   are characterized by
%
%   -  The order of quaternion elements
%
%   -  The quaternion multiplication formula
%
%   -  The convention for associating quaternions
%      with rotation matrices
%
%   Two of the commonly used styles are
%
%      - "SPICE"
%
%         > Invented by Sir William Rowan Hamilton
%         > Frequently used in mathematics and physics textbooks
%
%      - "Engineering"
%
%         > Widely used in aerospace engineering applications
%
%
%   Mice routine interfaces ALWAYS use SPICE quaternions.
%   Quaternions of any other style must be converted to SPICE
%   quaternions before they are passed to Mice routines.
%
%
%   Relationship between SPICE and Engineering Quaternions
%   ------------------------------------------------------
%
%   Let `m' be a rotation matrix such that for any vector `v',
%
%      m*v
%
%   is the result of rotating `v' by theta radians in the
%   counterclockwise direction about unit rotation axis vector `a'.
%   Then the SPICE quaternions representing `m' are
%
%      (+/-) (  cos(theta/2),
%               sin(theta/2) a(1),
%               sin(theta/2) a(2),
%               sin(theta/2) a(3)  )
%
%   while the engineering quaternions representing `m' are
%
%      (+/-) ( -sin(theta/2) a(1),
%              -sin(theta/2) a(2),
%              -sin(theta/2) a(3),
%               cos(theta/2)       )
%
%   For both styles of quaternions, if a quaternion q represents
%   a rotation matrix `m', then -q represents `m' as well.
%
%   Given an engineering quaternion
%
%      qeng   = ( q0,  q1,  q2,  q3 )
%
%   the equivalent SPICE quaternion is
%
%      qspice = ( q3, -q0, -q1, -q2 )
%
%
%   Associating SPICE Quaternions with Rotation Matrices
%   ----------------------------------------------------
%
%   Let `from' and `to' be two right-handed reference frames, for
%   example, an inertial frame and a spacecraft-fixed frame. Let the
%   symbols
%
%      v    ,   v
%       from     to
%
%   denote, respectively, an arbitrary vector expressed relative to
%   the `from' and `to' frames. Let `m' denote the transformation matrix
%   that transforms vectors from frame `from' to frame `to'; then
%
%      v   =  m * v
%       to         from
%
%   where the expression on the right hand side represents left
%   multiplication of the vector by the matrix.
%
%   Then if the unit-length SPICE quaternion q represents `m', where
%
%      q = (q0, q1, q2, q3)
%
%   the elements of `m' are derived from the elements of q as follows:
%
%        .-                                                         -.
%        |           2    2                                          |
%        | 1 - 2*( q2 + q3 )   2*(q1*q2 - q0*q3)   2*(q1*q3 + q0*q2) |
%        |                                                           |
%        |                                                           |
%        |                               2    2                      |
%    m = | 2*(q1*q2 + q0*q3)   1 - 2*( q1 + q3 )   2*(q2*q3 - q0*q1) |
%        |                                                           |
%        |                                                           |
%        |                                                   2    2  |
%        | 2*(q1*q3 - q0*q2)   2*(q2*q3 + q0*q1)   1 - 2*( q1 + q2 ) |
%        |                                                           |
%        `-                                                         -'
%
%   Note that substituting the elements of -q for those of q in the
%   right hand side leaves each element of `m' unchanged; this shows
%   that if a quaternion q represents a matrix `m', then so does the
%   quaternion -q.
%
%   To map the rotation matrix `m' to a unit quaternion, we start by
%   decomposing the rotation matrix as a sum of symmetric
%   and skew-symmetric parts:
%
%                                      2
%      m = [ i  +  (1-cos(theta)) omega  ] + [ sin(theta) omega ]
%
%                   symmetric                   skew-symmetric
%
%
%   `omega' is a skew-symmetric matrix of the form
%
%                 .-             -.
%                 |  0   -n3   n2 |
%                 |               |
%       omega  =  |  n3   0   -n1 |
%                 |               |
%                 | -n2   n1   0  |
%                 `-             -'
%
%   The vector N of matrix entries (n1, n2, n3) is the rotation axis
%   of `m' and theta is M's rotation angle. Note that N and theta
%   are not unique.
%
%   Let
%
%      C = cos(theta/2)
%      s = sin(theta/2)
%
%   Then the unit quaternions `q' corresponding to `m' are
%
%      `q' = +/- ( C, s*n1, s*n2, s*n3 )
%
%   The mappings between quaternions and the corresponding rotations
%   are carried out by the Mice routines
%
%      cspice_q2m {quaternion to matrix}
%      cspice_m2q {matrix to quaternion}
%
%   cspice_m2q always returns a quaternion with scalar part greater than
%   or equal to zero.
%
%
%   SPICE Quaternion Multiplication Formula
%   ---------------------------------------
%
%   Given a SPICE quaternion
%
%      q = ( q0, q1, q2, q3 )
%
%   corresponding to rotation axis `a' and angle theta as above, we can
%   represent `q' using "scalar + vector" notation as follows:
%
%      s =   q0           = cos(theta/2)
%
%      v = ( q1, q2, q3 ) = sin(theta/2) * a
%
%      q = s + v
%
%   Let `q1' and `q2' be SPICE quaternions with respective scalar
%   and vector parts s1, s2 and v1, v2:
%
%      q1 = s1 + v1
%      q2 = s2 + v2
%
%   We represent the dot product of v1 and v2 by
%
%      <v1, v2>
%
%   and the cross product of v1 and v2 by
%
%      v1 x v2
%
%   Then the SPICE quaternion product is
%
%      q1*q2 = s1*s2 - <v1,v2>  + s1*v2 + s2*v1 + (v1 x v2)
%
%   If `q1' and `q2' represent the rotation matrices `m1' and `m2'
%   respectively, then the quaternion product
%
%      q1*q2
%
%   represents the matrix product
%
%      m1*m2
%
%-Exceptions
%
%   1)  If `handle' is not the handle of a C-kernel opened for writing,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If `segid' is more than 40 characters long, the error
%       SPICE(SEGIDTOOLONG) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If `segid' contains any nonprintable characters, the error
%       SPICE(NONPRINTABLECHARS) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the first `start' time is negative, the error
%       SPICE(INVALIDSCLKTIME) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the second or any subsequent `start' times are negative, the
%       error SPICE(TIMESOUTOFORDER) is signaled by a routine in the
%       call tree of this routine.
%
%   6)  If any of the `stop' times are negative, the error
%       SPICE(DEGENERATEINTERVAL) is signaled by a routine in the call
%       tree of this routine.
%
%   7)  If the `stop' time of any of the intervals is less than or equal
%       to the `start' time, the error SPICE(DEGENERATEINTERVAL) is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the `start' times are not strictly increasing, the error
%       SPICE(TIMESOUTOFORDER) is signaled by a routine in the call
%       tree of this routine.
%
%   9)  If the `stop' time of one interval is greater than the `start'
%       time of the next interval, the error SPICE(BADSTOPTIME)
%       is signaled by a routine in the call tree of this routine.
%
%   10) If `begtim' is greater than start(1) or `endtim' is less than
%       stop(nrec), where `nrec' is the number of pointing records,
%       the error SPICE(INVALIDDESCRTIME) is signaled by a routine in
%       the call tree of this routine.
%
%   11) If the name of the reference frame is not one of those
%       supported by the routine cspice_namfrm, the error
%       SPICE(INVALIDREFFRAME) is signaled by a routine in the call
%       tree of this routine.
%
%   12) If `nrec', the number of pointing records, is less than or
%       equal to 0, the error SPICE(INVALIDNUMRECS) is signaled by a
%       routine in the call tree of this routine.
%
%   13) If any quaternion has magnitude zero, the error
%       SPICE(ZEROQUATERNION) is signaled by a routine in the call
%       tree of this routine.
%
%   14) If any of the input arguments, `handle', `begtim', `endtim',
%       `inst', `ref', `segid', `start', `stop', `quats', `avvs' or
%       `rates', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   15) If any of the input arguments, `handle', `begtim', `endtim',
%       `inst', `ref', `segid', `start', `stop', `quats', `avvs' or
%       `rates', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   16) If the input vector arguments `start', `stop', `quats', `avvs'
%       and `rates' do not have the same dimension (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   This routine adds a type 2 segment to a C-kernel. The C-kernel
%   may be either a new one or an existing one opened for writing.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   CK.REQ
%   DAF.REQ
%   SCLK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Changed input argument names "begtime" and "endtime" to "begtim"
%       and "endtim".
%
%       Edited the header to comply with NAIF standard. Added example's
%       problem statement.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 29-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 11-JUN-2013 (EDW)
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.1, 30-DEC-2008 (EDW)
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 04-JAN-2008 (EDW)
%
%-Index_Entries
%
%   write CK type_2 pointing data segment
%
%-&

function cspice_ckw02( handle, ...
                       begtim, ...
                       endtim, ...
                       inst,   ...
                       ref,    ...
                       segid,  ...
                       start,  ...
                       stop,   ...
                       quats,  ...
                       avvs,   ...
                       rates )

   switch nargin
      case 11

         handle = zzmice_int(handle);
         begtim = zzmice_dp(begtim);
         endtim = zzmice_dp(endtim);
         inst   = zzmice_int(inst);
         ref    = zzmice_str(ref);
         segid  = zzmice_str(segid);
         start  = zzmice_dp(start);
         stop   = zzmice_dp(stop);
         quats  = zzmice_dp(quats);
         avvs   = zzmice_dp(avvs);
         rates  = zzmice_dp(rates);

      otherwise

         error ( [ 'Usage: '                    ...
                   'cspice_ckw02( handle, '     ...
                                 'begtim, '     ...
                                 'endtim, '     ...
                                 'inst, '       ...
                                 'ref, '        ...
                                 'segid, '      ...
                                 'start(N), '   ...
                                 'stop(N), '    ...
                                 'quats(4,N), ' ...
                                 'avvs(3,N)) '  ...
                                 'rates(N)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'ckw02_c', handle, ...
                       begtim, ...
                       endtim, ...
                       inst,   ...
                       ref,    ...
                       segid,  ...
                       start,  ...
                       stop,   ...
                       quats,  ...
                       avvs,   ...
                       rates )
   catch spiceerr
      rethrow(spiceerr)
   end


