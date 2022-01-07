%-Abstract
%
%   CSPICE_GETFAT determines the file architecture and file type of most
%   SPICE kernel files.
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
%      file     the name of a SPICE kernel file whose architecture and type
%               are desired.
%
%               [1,c1] = size(file); char = class(file)
%
%                  or
%
%               [1,1] = size(file); cell = class(file)
%
%               This file must be closed when this routine is called.
%
%   the call:
%
%      [arch, kertyp] = cspice_getfat( file )
%
%   returns:
%
%      arch     the file architecture of the SPICE kernel file specified by
%               `file'.
%
%               [1,c2] = size(arch); char = class(arch)
%
%                  or
%
%               [1,1] = size(arch); cell = class(arch)
%
%               If the architecture cannot be determined or is not
%               recognized the value '?' is returned.
%
%               Architectures currently recognized are:
%
%                  DAF -- The file is based on the DAF architecture.
%                  DAS -- The file is based on the DAS architecture.
%                  XFR -- The file is in a SPICE transfer file format.
%                  DEC -- The file is an old SPICE decimal text file.
%                  ASC -- An ASCII text file.
%                  KPL -- Kernel Pool File (i.e., a text kernel)
%                  TXT -- An ASCII text file.
%                  TE1 -- Text E-Kernel type 1.
%                   ?  -- The architecture could not be determined.
%
%               This variable must be at least 3 characters long.
%
%      kertyp   the type of the SPICE kernel file.
%
%               [1,c3] = size(kertyp); char = class(kertyp)
%
%                  or
%
%               [1,1] = size(kertyp); cell = class(kertyp)
%
%               If the type can not be determined the value '?' is
%               returned.
%
%               Kernel file types may be any sequence of at most four
%               printing characters. NAIF has reserved for its use
%               types which contain all upper case letters.
%
%               A file type of 'PRE' means that the file is a
%               pre-release file.
%
%               This variable may be at most 4 characters long.
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
%   1) Determine the file architecture and file type of all the
%      kernels loaded through a meta-kernel and of a kernel in
%      transfer format.
%
%      Use the SPK kernel below to provide an example of a kernel in
%      transfer format.
%
%         earthstns_itrf93_050714.xsp
%
%
%      Use the meta-kernel shown below to load the other types of
%      SPICE kernels.
%
%
%         KPL/MK
%
%         File: getfat_ex1.tm
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
%            File name                        Contents
%            ---------                        --------
%            de430.bsp                        Planetary ephemeris
%            mar097.bsp                       Mars satellite ephemeris
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0011.tls                     Leapseconds
%            mgs_moc_v20.ti                   MGS MOC instrument
%                                             parameters
%            mgs_sclkscet_00061.tsc           MGS SCLK coefficients
%            mgs_sc_ext12.bc                  MGS s/c bus attitude
%            mgs_ext12_ipng_mgs95j.bsp        MGS ephemeris
%            megr90n000cb_plate.bds           Plate model based on
%                                             MEGDR DEM, resolution
%                                             4 pixels/degree.
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de430.bsp',
%                                'mar097.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls',
%                                'mgs_moc_v20.ti',
%                                'mgs_sclkscet_00061.tsc',
%                                'mgs_sc_ext12.bc',
%                                'mgs_ext12_ipng_mgs95j.bsp',
%                                'megr90n000cb_plate.bds'      )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function getfat_ex1()
%
%         %
%         % Check the file architecture and type of an SPK
%         % in transfer format.
%         %
%         fname1 = 'earthstns_itrf93_050714.xsp';
%         [arch, ktype] = cspice_getfat( fname1 );
%         fprintf( 'File name     :  %s\n', fname1 )
%         fprintf( '  Architecture:  %s\n', arch )
%         fprintf( '  Kernel type :  %s\n', ktype )
%         fprintf( '\n' )
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'getfat_ex1.tm' );
%
%         %
%         % Get the file architecture and kernel type for each of
%         % the kernels in the kernel pool.
%         %
%         [count] = cspice_ktotal( 'ALL' );
%
%         for i=1:count
%            [fname,  ktype, source,                          ...
%             handle, found]         = cspice_kdata( i, 'ALL' );
%
%            [arch, ktype] = cspice_getfat( fname );
%
%            fprintf( 'File name     :  %s\n', fname )
%            fprintf( '  Source      :  %s\n', source )
%            fprintf( '  Architecture:  %s\n', arch )
%            fprintf( '  Kernel type :  %s\n', ktype )
%            fprintf( '\n' )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      File name     :  earthstns_itrf93_050714.xsp
%        Architecture:  XFR
%        Kernel type :  DAF
%
%      File name     :  getfat_ex1.tm
%        Source      :
%        Architecture:  KPL
%        Kernel type :  MK
%
%      File name     :  de430.bsp
%        Source      :  getfat_ex1.tm
%        Architecture:  DAF
%        Kernel type :  SPK
%
%      File name     :  mar097.bsp
%        Source      :  getfat_ex1.tm
%        Architecture:  DAF
%        Kernel type :  SPK
%
%      File name     :  pck00010.tpc
%        Source      :  getfat_ex1.tm
%        Architecture:  KPL
%        Kernel type :  PCK
%
%      File name     :  naif0011.tls
%        Source      :  getfat_ex1.tm
%        Architecture:  KPL
%        Kernel type :  LSK
%
%      File name     :  mgs_moc_v20.ti
%        Source      :  getfat_ex1.tm
%        Architecture:  KPL
%        Kernel type :  IK
%
%      File name     :  mgs_sclkscet_00061.tsc
%        Source      :  getfat_ex1.tm
%        Architecture:  KPL
%        Kernel type :  SCLK
%
%      File name     :  mgs_sc_ext12.bc
%        Source      :  getfat_ex1.tm
%        Architecture:  DAF
%        Kernel type :  CK
%
%      File name     :  mgs_ext12_ipng_mgs95j.bsp
%        Source      :  getfat_ex1.tm
%        Architecture:  DAF
%        Kernel type :  SPK
%
%      File name     :  megr90n000cb_plate.bds
%        Source      :  getfat_ex1.tm
%        Architecture:  DAS
%        Kernel type :  DSK
%
%
%-Particulars
%
%   This function is a support utility routine that determines the
%   architecture and type of a SPICE kernel file.
%
%-Exceptions
%
%   1)  If the filename specified is blank, the error
%       SPICE(BLANKFILENAME) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any inquire on the filename specified by `file', required to
%       obtain information about the physical file, fails for some
%       reason, the error SPICE(INQUIREERROR) is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If the file specified by `file' does not exist, the error
%       SPICE(FILENOTFOUND) is signaled by a routine in the call tree
%       of this routine.
%
%   4)  If the file specified by `file' is already open but not through
%       SPICE interfaces, the error SPICE(EXTERNALOPEN) is signaled by
%       a routine in the call tree of this routine.
%
%   5)  If an attempt to open the file specified by `file' fails when
%       this routine requires that it succeed, the error
%       SPICE(FILEOPENFAILED) is signaled by a routine in the call
%       tree of this routine.
%
%   6)  If an attempt to read the file specified by `file' fails when
%       this routine requires that it succeed, the error
%       SPICE(FILEREADFAILED) is signaled by a routine in the call
%       tree of this routine.
%
%   7)  If an issue is detected during the opening the input file or
%       the process to determine its architecture and type, an error
%       is signaled by a routine in the call tree of this routine.
%
%   8)  If the ID word in a DAF based kernel is 'NAIF/DAF', then the
%       algorithm cspice_getfat uses to distinguish between CK and SPK
%       kernels may result in an indeterminate `kertyp' if the SPK or
%       CK files have invalid first segments.
%
%   9)  If the input argument `file' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   10) If the input argument `file' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   The SPICE kernel file specified by `file' is opened and then
%   closed by this routine to determine its file architecture and
%   type. Filenames of open files should not be passed to this
%   routine.
%
%-Restrictions
%
%   1)  In order to properly determine the type of DAF based binary
%       kernels, the routine requires that their first segments and
%       the meta data necessary to address them are valid.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   determine the architecture and type of a kernel file
%
%-&
function [arch, kertyp] = cspice_getfat( file )

   switch nargin
      case 1

         file = zzmice_str(file);

      otherwise

         error ( 'Usage: [`arch`, `kertyp`] = cspice_getfat( `file` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [arch, kertyp] = mice('getfat_c', file);
   catch spiceerr
      rethrow(spiceerr)
   end
