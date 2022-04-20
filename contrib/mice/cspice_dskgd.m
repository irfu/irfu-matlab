%-Abstract
%
%   CSPICE_DSKGD returns the DSK descriptor from a DSK segment identified
%   by a DAS handle and DLA descriptor.
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
%      handle   the handle of a DSK file that is open for read access.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      dladsc   the DLA segment descriptor corresponding to a DSK segment.
%
%               [SPICE_DLA_DSCSIZ,1]  = size(dladsc);
%                               int32 = class(dladsc)
%
%   the call:
%
%      [dskdsc] = cspice_dskgd( handle, dladsc )
%
%   returns:
%
%      dskdsc   the DSK segment descriptor of the segment designated by the
%               input handle and DLA descriptor.
%
%               [SPICE_DSK_DSCSIZ,1]  = size(dskdsc);
%                              double = class(dskdsc)
%
%-Parameters
%
%   See the parameter definitions file
%
%      MiceDLA.m
%
%   for declarations of DLA descriptor sizes and documentation of the
%   contents of DLA descriptors.
%
%   See the parameter definitions file
%
%      MiceDSK.m
%
%   for declarations of DSK descriptor sizes and documentation of the
%   contents of DSK descriptors.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Dump the DSK descriptor of the first segment of a DSK file.
%
%      Example code begins here.
%
%
%      function dskgd_ex1()
%
%         %
%         % Declare DSK Mice parameters for use in API calls.
%         %
%         MiceUser
%
%         %
%         % Prompt for the name of the DSK file.
%         %
%         dsk = input( 'Name of DSK file > ', 's' );
%
%         %
%         % Open the DSK file for read access.
%         % We use the DAS-level interface for
%         % this function.
%         %
%         handle  = cspice_dasopr( dsk );
%
%         %
%         % Begin a forward search through the
%         % kernel, treating the file as a DLA.
%         % In this example, it's a very short
%         % search.
%         %
%         [dladsc, found] = cspice_dlabfs( handle );
%
%         if ~found
%
%            %
%            % We arrive here only if the kernel
%            % contains no segments. This is
%            % unexpected, but we're prepared for it.
%            %
%            fprintf( 'No segments found in DSK file %s\n', dsk )
%            return
%
%         end
%
%         %
%         % If we made it this far, `dladsc' is the
%         % DLA descriptor of the first segment.
%         %
%         % Fetch the DSK descriptor of this segment.
%         %
%         dskdsc = cspice_dskgd( handle, dladsc );
%
%         fprintf( '\n' )
%         fprintf( 'DSK descriptor contents:\n\n' )
%         fprintf( '   Surface ID:              %d\n',                     ...
%                        dskdsc(SPICE_DSK_SRFIDX) )
%         fprintf( '   Center ID:               %d\n',                     ...
%                        dskdsc(SPICE_DSK_CTRIDX) )
%         fprintf( '   Data class:               %d\n',                    ...
%                        dskdsc(SPICE_DSK_CLSIDX) )
%         fprintf( '   Data type:                %d\n',                    ...
%                        dskdsc(SPICE_DSK_TYPIDX) )
%         fprintf( '   Frame ID:                 %d\n',                    ...
%                        dskdsc(SPICE_DSK_FRMIDX) )
%         fprintf( '   Coordinate system:        %d\n',                    ...
%                        dskdsc(SPICE_DSK_SYSIDX) )
%         fprintf( '   Parameters:       %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_PARIDX) )
%
%         for   i = 1:(SPICE_DSK_NSYPAR-1)
%
%            fprintf( '                     %15.6f\n', ...
%                           dskdsc(SPICE_DSK_PARIDX + i) )
%         end
%
%         fprintf( '   Coordinate 1 min: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MN1IDX) )
%         fprintf( '   Coordinate 1 max: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MX1IDX) )
%         fprintf( '   Coordinate 2 min: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MN2IDX) )
%         fprintf( '   Coordinate 2 max: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MX2IDX) )
%         fprintf( '   Coordinate 3 min: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MN3IDX) )
%         fprintf( '   Coordinate 3 max: %15.6f\n',                        ...
%                        dskdsc(SPICE_DSK_MX3IDX) )
%         fprintf( '   Start time:    %18.6f\n',                           ...
%                        dskdsc(SPICE_DSK_BTMIDX) )
%         fprintf( '   Stop time:     %18.6f\n',                           ...
%                        dskdsc(SPICE_DSK_ETMIDX) )
%
%         %
%         % Close file.
%         %
%         cspice_dascls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the DSK file named phobos512.bds, the output
%      was:
%
%
%      Name of DSK file > phobos512.bds
%
%      DSK descriptor contents:
%
%         Surface ID:              401
%         Center ID:               401
%         Data class:               1
%         Data type:                2
%         Frame ID:                 10021
%         Coordinate system:        1
%         Parameters:              0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%                                  0.000000
%         Coordinate 1 min:       -3.141593
%         Coordinate 1 max:        3.141593
%         Coordinate 2 min:       -1.570796
%         Coordinate 2 max:        1.570796
%         Coordinate 3 min:        8.049632
%         Coordinate 3 max:       13.940940
%         Start time:    -1577879958.816059
%         Stop time:      1577880066.183913
%
%
%   2) Again, dump the DSK descriptors of a DSK file, this time
%      interpreting the descriptor information and displaying
%      it in a user-friendly form. This display is a simplified
%      version of that created by the utility DSKBRIEF.
%
%      This program requests the name of an optional meta-kernel.
%      The meta-kernel can be used to define surface name-ID
%      associations. If no meta-kernel is needed, the user can
%      enter a carriage return at the prompt for this file.
%
%
%      Example code begins here.
%
%
%      function dskgd_ex2()
%
%         %
%         % Declare DSK Mice parameters for use in API calls.
%         %
%         MiceUser
%
%         %
%         % Initial values
%         %
%         clsnms = {'Single-valued surface', 'General surface'};
%
%         sysnms = {'Latitudinal', 'Cylindrical', 'Rectangular',           ...
%                  'Planetodetic'};
%
%         dsk    = input( 'Enter DSK name         > ', 's' );
%         meta   = input( 'Enter meta-kernel name > ', 's' );
%
%         if ( ~strcmp(meta, '') && ~strcmp(meta, ' ') )
%            cspice_furnsh( meta );
%         end
%
%         %
%         % Open the DLA file and begin a forward search
%         % for segments.
%         %
%         [handle] = cspice_dasopr( dsk );
%
%         segno    = 0;
%
%         [nxtdsc, found] = cspice_dlabfs( handle );
%
%         while found
%
%            segno = segno + 1;
%
%            %
%            % Make the DLA descriptor we just fetched
%            % the current one.
%            %
%            dladsc   = nxtdsc;
%
%            [dskdsc] = cspice_dskgd( handle, dladsc );
%
%            bodyid   = int64( dskdsc(SPICE_DSK_CTRIDX) );
%            surfid   = int64( dskdsc(SPICE_DSK_SRFIDX) );
%            framid   = int64( dskdsc(SPICE_DSK_FRMIDX) );
%            dtype    = int64( dskdsc(SPICE_DSK_TYPIDX) );
%            dclass   = int64( dskdsc(SPICE_DSK_CLSIDX) );
%
%            [bodnam]         = cspice_bodc2s( bodyid );
%            [srfnam, isname] = cspice_srfc2s( surfid, bodyid );
%            [frame]          = cspice_frmnam( framid );
%
%            if ( frame == ' ' )
%               sprintf( frame, '%d', framid );
%            end
%
%            [btime]  = cspice_etcal( dskdsc(SPICE_DSK_BTMIDX) );
%            [etime]  = cspice_etcal( dskdsc(SPICE_DSK_ETMIDX) );
%
%            corsys   = int64( dskdsc(SPICE_DSK_SYSIDX) );
%
%            sysnam   = sysnms( corsys );
%
%            fprintf( '====================================\n' )
%            fprintf( 'DSK descriptor for segment  %2d\n', segno )
%            fprintf( '  Body:               %s\n', bodnam )
%            fprintf( '  Surface:            %s\n', srfnam )
%            fprintf( '  Frame:              %s\n', frame )
%            fprintf( '  Start time (TDB):   %s\n', btime )
%            fprintf( '  Stop time  (TDB):   %s\n', etime )
%            fprintf( '  Data type:          %d\n', dtype )
%            fprintf( '  Data class:         %d   %s\n',                   ...
%                           dclass, char(clsnms(dclass)) )
%            fprintf( '  Coordinate system:  %s\n', char(sysnam) )
%
%            if ( corsys == SPICE_DSK_PDTSYS )
%
%               re = dskdsc(SPICE_DSK_PARIDX  );
%               f  = dskdsc(SPICE_DSK_PARIDX+1);
%               rp = re * ( 1.0 - f );
%
%               fprintf( '     Equatorial radius (km):  %f\n', re )
%               fprintf( '     Polar radius      (km):  %f\n', rp )
%
%            end
%
%            fprintf( '  Segment boundaries:\n' )
%
%            if ( corsys == SPICE_DSK_LATSYS )
%
%               fprintf( '   Longitude (deg):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN1IDX)* cspice_dpr,   ...
%                                  dskdsc(SPICE_DSK_MX1IDX)* cspice_dpr )
%               fprintf( '   Latitude  (deg):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN2IDX)* cspice_dpr,   ...
%                                  dskdsc(SPICE_DSK_MX2IDX)* cspice_dpr )
%               fprintf( '   Radius     (km):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN3IDX),               ...
%                                  dskdsc(SPICE_DSK_MX3IDX)         )
%
%            elseif ( corsys == SPICE_DSK_CYLSYS )
%
%               fprintf ( 'Coordinate system was Cylindrical\n' )
%               exit
%
%            elseif ( corsys == SPICE_DSK_RECSYS )
%
%               fprintf( '   X-coordinate (km):  %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN1IDX),               ...
%                                  dskdsc(SPICE_DSK_MX1IDX)         )
%               fprintf( '   Y-coordinate (km):  %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN2IDX),               ...
%                                  dskdsc(SPICE_DSK_MX2IDX)         )
%               fprintf( '   Z-coordinate (km):  %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN3IDX),               ...
%                                  dskdsc(SPICE_DSK_MX3IDX)         )
%
%            elseif ( corsys == SPICE_DSK_PDTSYS )
%
%               fprintf( '   Longitude (deg):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN1IDX)* cspice_dpr,   ...
%                                  dskdsc(SPICE_DSK_MX1IDX)* cspice_dpr )
%               fprintf( '   Latitude  (deg):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN2IDX)* cspice_dpr,   ...
%                                  dskdsc(SPICE_DSK_MX2IDX)* cspice_dpr )
%               fprintf( '   Altitude   (km):    %18.12f %18.12f\n',       ...
%                                  dskdsc(SPICE_DSK_MN3IDX),               ...
%                                  dskdsc(SPICE_DSK_MX3IDX)         )
%            end
%
%            %
%            % Find the next segment, if it exists.
%            %
%            [nxtdsc, found] = cspice_dlafns( handle, dladsc );
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
%      platform, using the DSK file named phobos512.bds and an empty
%      string instead of the meta-kernel name, the output was:
%
%
%      Enter DSK name         > phobos512.bds
%      Enter meta-kernel name >
%      ====================================
%      DSK descriptor for segment   1
%        Body:               PHOBOS
%        Surface:            401
%        Frame:              IAU_PHOBOS
%        Start time (TDB):   1950 JAN 01 00:00:41.183
%        Stop time  (TDB):   2050 JAN 01 00:01:06.183
%        Data type:          2
%        Data class:         1   Single-valued surface
%        Coordinate system:  Latitudinal
%        Segment boundaries:
%         Longitude (deg):     -180.000000000000   180.000000000000
%         Latitude  (deg):      -90.000000000000    90.000000000000
%         Radius     (km):        8.049632248722    13.940939832124
%
%
%   3) Again, dump the DSK descriptors of a DSK file, using the
%      program from example 2, but this time reading the DSK file
%
%         phobos_3_3_3seg.bds
%
%      which can be created by running an example program from
%      cspice_dskw02. Use the meta-kernel shown below to demonstrate surface
%      name-ID mapping.
%
%
%         KPL/MK
%
%         File: dskgd_ex3.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The file contents shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%
%         \begindata
%
%         NAIF_SURFACE_NAME += ( 'Phobos example surface 1',
%                                'Phobos example surface 2',
%                                'Phobos example surface 3' )
%         NAIF_SURFACE_CODE += (   1,   2,   3 )
%         NAIF_SURFACE_BODY += ( 401, 401, 401 )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      When Example #2 was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the DSK file named phobos_3_3_3seg.bds and the
%      meta-kernel dskgd_ex3.tm, the output was:
%
%
%      Enter DSK name         > phobos_3_3_3seg.bds
%      Enter meta-kernel name > dskgd_ex3.tm
%      ====================================
%      DSK descriptor for segment   1
%        Body:               PHOBOS
%        Surface:            Phobos example surface 1
%        Frame:              IAU_PHOBOS
%        Start time (TDB):   1950 JAN 01 00:00:00.000
%        Stop time  (TDB):   2050 JAN 01 00:00:00.000
%        Data type:          2
%        Data class:         2   General surface
%        Coordinate system:  Latitudinal
%        Segment boundaries:
%         Longitude (deg):     -180.000000000000   180.000000000000
%         Latitude  (deg):      -90.000000000000    90.000000000000
%         Radius     (km):        8.225298075974    14.011768145626
%      ====================================
%      DSK descriptor for segment   2
%        Body:               PHOBOS
%        Surface:            Phobos example surface 2
%        Frame:              IAU_PHOBOS
%        Start time (TDB):   1950 JAN 01 00:00:00.000
%        Stop time  (TDB):   2050 JAN 01 00:00:00.000
%        Data type:          2
%        Data class:         2   General surface
%        Coordinate system:  Rectangular
%        Segment boundaries:
%         X-coordinate (km):     -1.300000000000     1.310000000000
%         Y-coordinate (km):     -1.210000000000     1.200000000000
%         Z-coordinate (km):     -9.452932357788     9.638179779053
%      ====================================
%      DSK descriptor for segment   3
%        Body:               PHOBOS
%        Surface:            Phobos example surface 3
%        Frame:              IAU_PHOBOS
%        Start time (TDB):   1950 JAN 01 00:00:00.000
%        Stop time  (TDB):   2050 JAN 01 00:00:00.000
%        Data type:          2
%        Data class:         2   General surface
%        Coordinate system:  Planetodetic
%           Equatorial radius (km):  13.000000
%           Polar radius      (km):  9.100000
%        Segment boundaries:
%         Longitude (deg):     -180.000000000000   180.000000000000
%         Latitude  (deg):      -90.000000000000    90.000000000000
%         Altitude   (km):       -3.728668683604     1.372015791081
%
%
%-Particulars
%
%   This is a convenience routine intended for use by low-level
%   routines that read DSK segments. This routine may also be called
%   by user applications that must access DSK files at the segment
%   level.
%
%-Exceptions
%
%   1)  If the size of the double precision component of the segment
%       is smaller than that of a DSK descriptor, the error
%       SPICE(INVALIDFORMAT) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the input handle is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the input DLA descriptor is invalid, the effect of this
%       routine is undefined. The error *may* be diagnosed by
%       routines in the call tree of this routine, but there are no
%       guarantees.
%
%   4)  If any DAS read error is detected, the error is signaled by a
%       routine in the call tree of this routine.
%
%   5)  If any of the input arguments, `handle' or `dladsc', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `handle' or `dladsc', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   See input argument `handle'.
%
%-Restrictions
%
%   1)  See Exception #3.
%
%-Required_Reading
%
%   DAS.REQ
%   DSK.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Updated
%       code example to prompt for the input DSK file and added start and
%       stop times to the output. Added second and third examples.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 05-MAY-2014 (NJB) (EDW)
%
%-Index_Entries
%
%   return DSK segment descriptor
%
%-&

function [dskdsc] = cspice_dskgd( handle, dladsc )

   switch nargin
      case 2

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );

      otherwise

         error ( 'Usage: [dskdsc(24)] = cspice_dskgd( handle, dladsc )' )

   end

   %
   % Call the MEX library.
   %
   try
      [dskdsc] = mice('dskgd_c', handle, dladsc );
   catch spiceerr
      rethrow(spiceerr)
   end


