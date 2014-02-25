function mms_sdc_cdfwrite( varargin )
% MMS_SDC_CDFWRITE writes the input data to the corresponding CDF file in
% the output dir defined in ENVIR.DROPBOX_ROOT.
%   Three different types of data files can be created SITL and QUICKLOOK
%   DCE files and USC potential files.
%	Example:
%      MMS_SDC_CDFWRITE( filename.cdf, scID, 'sitl', xyz_pgse_data, ...
%      xyz_dsl_data, bitmask);
%   will produce a SITL data file called filename.cdf.
%      MMS_SDC_CDFWRITE( filename.cdf, scID, 'ql', xyz_pgse_data, ...
%      xyz_dsl_data, bitmask, quality);
%   will create a QuickLook data file called filename.cdf and subsequently
%      MMS_SDC_CDFWRITE( filename.cdf, scID, 'usc', ESCP, PSP, Delta, ...
%      PSP_P, bitmask);
%   will create a Spacecraft potential data file called filename.cdf.
%
%   If the corresponding mms_sdc_cdfwritec.mexa64 file has not yet been
%   created this script will also build it using mex, the C source file and
%   linking it to CDF code found in ENVIR.CDF_BASE.
%
%	Note 1: It is assumed that other SDC processing scripts will move the 
%   created output file to its final destination (from /ENIVR.DROPBOX_ROOT/ 
%   to /path/as/defined/by/	SDC Developer Guide, v1.7).
%	Note 2: 20140204, FIELDS consider SITL to be a dataLevel even if the
%   SDC Developer Guide is unclear about that issue. Therefor when running
%   MMS SDP SITL processing this script will check corresponding folder
%   structure and write filenames accordingly.
%
% 	See also MMS_CDF_WRITING, MMS_BITMASKING, MMS_INIT.


% Check number of inputs
narginchk(7,9);

global ENVIR;

if(ischar(varargin{1}))
		filename_output = lower(varargin{1});
else
	irf.log('critical','Matlab:MMS_SDC_CDFWRITE: String input required as filename.');
	error('Matlab:MMS_SDC_CDFWRITE:String input required.');
end

if( isnumeric(varargin{2}) && varargin{2}<=4 && varargin{2}>=1 )
    scID = int8(varargin{2}); % Make sure it is of right size.
else
    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: SC id must be number between 1 and 4.');
    error('Matlab:MMS_SDC_CDFWRITE: SC id must be number between 1 and 4.');
end

% Validate number of inputs for each type.
if(ischar(varargin{3}))
    calledBy = lower(varargin{3});
    switch(calledBy)
        case('sitl')
            if(nargin~=7)
                irf.log('critical','Matlab:MMS_SDC_CDFWRITE: requires 7 input arguments for STIL: Filename, scId, "sitl", Epoch, PGSE data, DSL data, bitmask.');
                error('Matlab:MMS_SDC_CDFWRITE: requires 7 input arguments for STIL.');
            else
                % 7 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==3))
                    dce_xyz_pgse = varargin{5};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: SITL dce_xyz_pgse must be a numerical matrix with three columns.');
                    error('Matlab:MMS_SDC_CDFWRITE: SITL dce_xyz_pgse must be a numerical matrix with three columns.');
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==3))
                    dce_xyz_dsl = varargin{6};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: SITL dce_xyz_dsl must be a numerical matrix with three columns.');
                    error('Matlab:MMS_SDC_CDFWRITE: SITL dce_xyz_dsl must be a numerical matrix with three columns.');
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    bitmask = uint16(varargin{7});
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: SITL bitmask must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: SITL bitmask must be a numerical vector.');
                end
            end
            
        case('ql')
            if(nargin~=8)
                irf.log('critical','Matlab:MMS_SDC_CDFWRITE: requires 8 input arguments for Quicklook: Filename, "ql", scId, Epoch, PGSE data, DSL data, bitmask, quality.');
                error('Matlab:MMS_SDC_CDFWRITE: requires 8 input arguments for Quicklook.');
            else
                % 8 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==3))
                    dce_xyz_pgse = varargin{5};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: QuickLook dce_xyz_pgse must be a numerical matrix with three columns.');
                    error('Matlab:MMS_SDC_CDFWRITE: QuickLook dce_xyz_pgse must be a numerical matrix with three columns.');
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==3))
                    dce_xyz_dsl = varargin{6};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: QuickLook dce_xyz_dsl must be a numerical matrix with three columns.');
                    error('Matlab:MMS_SDC_CDFWRITE: QuickLook dce_xyz_dsl must be a numerical matrix with three columns.');
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    bitmask = uint16(varargin{7}); % make sure it is stored with correct size
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: QuickLook bitmask must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: QuickLook bitmask must be a numerical vector.');
                end
                if((isnumeric(varargin{8}) && size(varargin{8},2)==1))
                    qualityMark = uint16(varargin{8}); % make sure it is store with correct size
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: QuickLook quality must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: QuickLook quality must be a numerical vector.');
                end
            end
        case('usc')
            if(nargin~=9)
                irf.log('critical','Matlab:MMS_SDC_CDFWRITE: requires 9 input arguments for Usc: Filename, scId, "usc", Epoch, PGSE data, DSL data, bitmask, quality.');
                error('Matlab:MMS_SDC_CDFWRITE: requires 9 input arguments for Usc.');
            else
                % 9 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==1))
                    escp = varargin{5};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Usc ESCP must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: Usc ESCP must be a numerical vector.');
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==1))
                    psp = varargin{6};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Usc PSP must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: Usc PSP must be a numerical vector.');
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    delta = varargin{7};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Usc DELTA must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: Usc DELTA must be a numerical vector.');
                end
                if((isnumeric(varargin{8}) && size(varargin{8},2)==6))
                    psp_p = varargin{8};
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Usc PSP_P must be a numerical matrix with six columns.');
                    error('Matlab:MMS_SDC_CDFWRITE: Usc PSP_P must be a numerical matrix with six columns.');
                end
                if((isnumeric(varargin{9}) && size(varargin{9},2)==1))
                    bitmask = uint16(varargin{9}); % make sure it is stored with correct size
                else
                    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Usc bitmask must be a numerical vector.');
                    error('Matlab:MMS_SDC_CDFWRITE: Usc bitmask must be a numerical vector.');
                end
            end
    otherwise
        irf.log('critical','MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid for types.');
        error('MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid for types.');        
    end
else
    irf.log('critical','MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid types and must be sent as string.');
    error('MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid for types.');
end

if( isnumeric(varargin{4}) )
    EpochTimes = int64(varargin{4}); % Make sure it is of right size.
else
    irf.log('critical','Matlab:MMS_SDC_CDFWRITE: Epoch times must be numerical.');
    error('Matlab:MMS_SDC_CDFWRITE: Epoch times must be numberical.');
end


% Check to see if compiled file exists if not build it.
if(~exist('mms_sdc_cdfwritec.mexa64','file'))
    
    % Provide warning to users logfile to ensure they are know.
    irf.log('warning','MMS_SDC_CDFWRITEC.MEXA64 file does not exist. Building it.');
    
    % Locate the path to source and output directory.
    pathToSrc = which('mms_sdc_cdfwritec.c');
    [s, filename, ext] = fileparts(pathToSrc);
    
    % Call mex with the required arguments, if user has set up some other
    % mexopts.sh the following values will be used instead. If user runs on
    % a very old Matlab (7.11) there might be issues with linking to
    % Matlab's own lib.
    mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -D_GNU_SOURCE -pthread -fexceptions'], ...
        ['-I',ENVIR.CDF_BASE,'/include/'] ,['-L',ENVIR.CDF_BASE,'/lib/'], ['-Wl,-rpath,',ENVIR.CDF_BASE,'/lib/'],...
        '-outdir',[s,'/'],'-lcdf', [s,'/mms_sdc_cdfwritec.c']);
end

% Actually call the mexa64 file.
switch(calledBy)
    
    case('sitl')
        
        try
            mms_sdc_cdfwritec(filename_output, scID, 'sitl', EpochTimes, dce_xyz_pgse', dce_xyz_dsl', bitmask);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if (strcmp(err.identifier,'MATLAB:mms_sdc_cdfwrite:filename_output:exists'))
                % If our cdfwrite code resulted in error write proper log message.
                irf.log('critical',err.message);
                % Then end with MATLAB:SDCcode and numberical error code to
                % be fetched by bash script.
                error('MATLAB:SDCcode', '183');
            else
                % Display any other errors as usual.
                rethrow(err);
            end
        end % End of try
        
    case('ql')
        
        try
            mms_sdc_cdfwritec(filename_output, scID, 'ql', EpochTimes, dce_xyz_pgse', dce_xyz_dsl', bitmask, qualityMark);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if (strcmp(err.identifier,'MATLAB:mms_sdc_cdfwrite:filename_output:exists'))
                % If our cdfwrite code resulted in error write proper log message.
                irf.log('critical',err.message);
                % Then end with MATLAB:SDCcode and numberical error code to
                % be fetched by bash script.
                error('MATLAB:SDCcode', '183');
            else
                % Display any other errors as usual.
                rethrow(err);
            end
        end % End of try
        
    case('usc')
        
        try
            mms_sdc_cdfwritec(filename_output, scID, 'usc', EpochTimes, escp, psp, delta, psp_p', bitmask);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if (strcmp(err.identifier,'MATLAB:mms_sdc_cdfwrite:filename_output:exists'))
                % If our cdfwrite code resulted in error write proper log message.
                irf.log('critical',err.message);
                % Then end with MATLAB:SDCcode and numberical error code to
                % be fetched by bash script.
                error('MATLAB:SDCcode', '183');
            else
                % Display any other errors as usual.
                rethrow(err);
            end
        end % End of try
    otherwise
        irf.log('critical','MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid types.');
        error('MATLAB:MMS_SDC_CDFWRITE: Only "ql", "sitl" or "usc" are valid for types.');
end


end

