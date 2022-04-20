%-Abstract
%
%   CSPICE_CKW03 adds a type 3 segment to a C-kernel.
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
%      avflag   a boolean signifying if the segment will contain
%               angular velocity.
%
%               [1,1] = size(avflag); logical = class(avflag)
%
%      segid    name to identify the segment.
%
%               [1,c2] = size(segid); char = class(segid)
%
%                  or
%
%               [1,1] = size(segid); cell = class(segid)
%
%      sclkdp   array containing the encoded SCLK times for the data.
%
%               [n,1] = size(sclkdp); double = class(sclkdp)
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
%      starts   array containing the encoded SCLK interval start times of
%               each interpolation interval, the times must be strictly
%               increasing and coincide with pointing data times.
%
%               [m,1] = size(starts); double = class(starts)
%
%   the call:
%
%      cspice_ckw03( handle, ...
%                    begtim, ...
%                    endtim, ...
%                    inst,   ...
%                    ref,    ...
%                    avflag, ...
%                    segid,  ...
%                    sclkdp, ...
%                    quats,  ...
%                    avvs,   ...
%                    starts)
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
%   1) The following example creates a CK file with a type-3 segment.
%
%      Example code begins here.
%
%
%      function ckw03_ex1()
%
%         INST3      = -77703;
%         NCOMCH     = 10;
%         REF        = 'J2000';
%         SEGID3     = 'Test type 3 test CK';
%         SECPERTICK = 0.001;
%         SPACING    = 10.0;
%         MAXREC     = 50;
%
%         %
%         % Note, `sclkdp' is a vector input, not a vectorized scalar.
%         %
%         sclkdp    = [1:MAXREC]';
%         sclkdp    = (sclkdp - 1)*SPACING;
%
%         spinrate  = [1:MAXREC]*1.e-6;
%
%         theta     = [0:MAXREC-1]*SPACING;
%         theta     = theta .* spinrate;
%
%         %
%         % Create a zero-filled array for the angular velocity
%         % vectors. This allocates the needed memory and
%         % defines a variable of the correct shape.
%         %
%         expavvs = zeros( [3 MAXREC] );
%
%         a1 = zeros( [1 MAXREC] );
%         a2 = a1;
%
%         r  = cspice_eul2m( theta, a2, a1, 3, 1 ,3 );
%         q  = cspice_m2q( r );
%
%         %
%         % Fill the z component of the expavvs vectors with the
%         % corresponding spinrate element scaled to SECPERTICK.
%         %
%         expavvs(3,:) = spinrate/SECPERTICK;
%
%         begtim = sclkdp(1);
%         endtim = sclkdp(MAXREC);
%         avflag = 1;
%
%         starts = [1:(MAXREC/2)]';
%         starts = (starts-1)*2*SPACING;
%
%         %
%         % Open a new CK, write the data, catch any errors.
%         %
%         try
%            handle = cspice_ckopn( 'ckw03_ex1.ck', 'ck', 0);
%            cspice_ckw03( handle,  ...
%                          begtim,  ...
%                          endtim,  ...
%                          INST3,   ...
%                          REF,     ...
%                          avflag,  ...
%                          SEGID3,  ...
%                          sclkdp,  ...
%                          q,       ...
%                          expavvs, ...
%                          starts )
%         catch
%
%            error( [ 'Failure: ' lasterr] )
%         end
%
%         cspice_ckcls(handle)
%
%
%      When this program is executed, no output is presented on
%      screen. After run completion, a new CK file exists in the
%      output directory.
%
%-Particulars
%
%   For a detailed description of a type 3 CK segment please see the
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
%   3)  If `segid' contains any non-printable characters, the error
%       SPICE(NONPRINTABLECHARS) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the first encoded SCLK time is negative, the error
%       SPICE(INVALIDSCLKTIME) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the second encoded SCLK or any subsequent times, or if the
%       encoded SCLK times are not strictly increasing, the error
%       SPICE(TIMESOUTOFORDER) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If `begtim' is greater than sclkdp(1) or `endtim' is less than
%       sclkdp(nrec), where `nrec' is the number of pointing records,
%       the error SPICE(INVALIDDESCRTIME) is signaled by a routine in
%       the call tree of this routine.
%
%   7)  If the name of the reference frame is not one of those
%       supported by the Mice routine cspice_namfrm, the error
%       SPICE(INVALIDREFFRAME) is signaled by a routine in the call
%       tree of this routine.
%
%   8)  If `nrec', the number of pointing records, is less than or
%       equal to 0, the error SPICE(INVALIDNUMREC) is signaled by a
%       routine in the call tree of this routine.
%
%   9)  If `nints', the number of interpolation intervals, is less
%       than or equal to 0, the error SPICE(INVALIDNUMINT) is signaled
%       by a routine in the call tree of this routine.
%
%   10) If the encoded SCLK interval start times are not strictly
%       increasing, the error SPICE(TIMESOUTOFORDER) is signaled by a
%       routine in the call tree of this routine.
%
%   11) If an interval start time does not coincide with a time for
%       which there is an actual pointing instance in the segment, the
%       error SPICE(INVALIDSTARTTIME) is signaled by a routine in the
%       call tree of this routine.
%
%   12) This routine assumes that the rotation between adjacent
%       quaternions that are stored in the same interval has a
%       rotation angle of `theta' radians, where
%
%          0  <=  theta  <  pi.
%
%       The routines that evaluate the data in the segment produced
%       by this routine cannot distinguish between rotations of `theta'
%       radians, where `theta' is in the interval [0, pi), and
%       rotations of
%
%          theta   +   2 * k * pi
%
%       radians, where k is any integer. These `large' rotations
%       will yield invalid results when interpolated. You must
%       ensure that the data stored in the segment will not be
%       subject to this sort of ambiguity.
%
%   13) If any quaternion has magnitude zero, the error
%       SPICE(ZEROQUATERNION) is signaled by a routine in the call
%       tree of this routine.
%
%   14) If the start time of the first interval and the time of the
%       first pointing instance are not the same, the error
%       SPICE(TIMESDONTMATCH) is signaled by a routine in the call
%       tree of this routine.
%
%   15) If any of the input arguments, `handle', `begtim', `endtim',
%       `inst', `ref', `avflag', `segid', `sclkdp', `quats', `avvs' or
%       `starts', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   16) If any of the input arguments, `handle', `begtim', `endtim',
%       `inst', `ref', `avflag', `segid', `sclkdp', `quats', `avvs' or
%       `starts', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%   17) If the input vector arguments `sclkdp', `quats' and `avvs' do
%       not have the same dimension (N), an error is signaled by the
%       Mice interface.
%
%-Files
%
%   This routine adds a type 3 segment to a C-kernel. The C-kernel
%   may be either a new one or an existing one opened for writing.
%
%-Restrictions
%
%   1)  The creator of the segment is given the responsibility for
%       determining whether it is reasonable to interpolate between
%       two given pointing values.
%
%   2)  This routine assumes that the rotation between adjacent
%       quaternions that are stored in the same interval has a
%       rotation angle of `theta' radians, where
%
%           0  <=  theta  <  pi.
%
%       The routines that evaluate the data in the segment produced
%       by this routine cannot distinguish between rotations of `theta'
%       radians, where `theta' is in the interval [0, pi), and
%       rotations of
%
%           theta   +   2 * k * pi
%
%       radians, where k is any integer. These `large' rotations will
%       yield invalid results when interpolated. You must ensure that
%       the data stored in the segment will not be subject to this
%       sort of ambiguity.
%
%   3)  All pointing instances in the segment must belong to one and
%       only one of the intervals.
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
%   -Mice Version 1.0.2, 29-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 11-JUL-2012 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 19-MAY-2006 (EDW)
%
%-Index_Entries
%
%   write CK type_3 pointing data segment
%
%-&

function cspice_ckw03( handle, ...
                       begtim, ...
                       endtim, ...
                       inst,   ...
                       ref,    ...
                       avflag, ...
                       segid,  ...
                       sclkdp, ...
                       quats,  ...
                       avvs,   ...
                       starts )

   switch nargin
      case 11

         handle = zzmice_int(handle);
         begtim = zzmice_dp(begtim);
         endtim = zzmice_dp(endtim);
         avflag = zzmice_int(avflag);
         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         quats  = zzmice_dp(quats);
         avvs   = zzmice_dp(avvs);
         starts = zzmice_dp(starts);

      otherwise

         error ( [ 'Usage: '                    ...
                   'cspice_ckw03( handle, '     ...
                                 'begtim, '     ...
                                 'endtim, '     ...
                                 'inst, '       ...
                                 'ref, '        ...
                                 'avflag, '     ...
                                 'segid, '      ...
                                 'sclkdp(N), '  ...
                                 'quats(4,N), ' ...
                                 'avvs(3,N), '  ...
                                 'starts(M))' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'ckw03_c', handle, begtim, endtim, inst, ref,  avflag,  ...
                       segid,  sclkdp,  quats, avvs, starts )
   catch spiceerr
      rethrow(spiceerr)
   end



