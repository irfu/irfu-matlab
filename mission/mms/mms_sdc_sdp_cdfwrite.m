function mms_sdc_sdp_cdfwrite( varargin )
% MMS_SDC_SDP_CDFWRITE writes the input data to the corresponding CDF file in
% the output dir defined in ENVIR.DROPBOX_ROOT.
%   Three different types of data files can be created SITL and QUICKLOOK
%   DCE files and USC potential files.
%	Example:
%      MMS_SDC_SDP_CDFWRITE( filename.cdf, scID, 'sitl', xyz_pgse_data, ...
%      xyz_dsl_data, bitmask);
%   will produce a SITL data file called filename.cdf.
%      MMS_SDC_SDP_CDFWRITE( filename.cdf, scID, 'ql', xyz_pgse_data, ...
%      xyz_dsl_data, bitmask, quality);
%   will create a QuickLook data file called filename.cdf and subsequently
%      MMS_SDC_SDP_CDFWRITE( filename.cdf, scID, 'usc', ESCP, PSP, Delta, ...
%      PSP_P, bitmask);
%   will create a Spacecraft potential data file called filename.cdf.
%
%   If the corresponding mms_sdc_sdp_cdfwritec.mexa64 file has not yet been
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
% 	See also MMS_SDC_SDP_CDF_WRITING, MMS_SDC_SDP_INIT.


% Check number of inputs
narginchk(7,9);

global ENVIR;

if(ischar(varargin{1}))
    filename_output = lower(varargin{1});
else
    err_str = 'MMS_SDC_SDP_CDFWRITE filename required as first argument.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
end

if( isnumeric(varargin{2}) && varargin{2}<=4 && varargin{2}>=1 )
    scID = int8(varargin{2}); % Make sure it is of right size.
else
    err_str = 'MMS_SDC_SDP_CDFWRITE SC id must be number between 1 and 4.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
end

% Validate number of inputs for each type.
if(ischar(varargin{3}))
    calledBy = lower(varargin{3});
    switch(calledBy)
        case('sitl')
            if(nargin~=7)
                err_str = 'MMS_SDC_SDP_CDFWRITE requires 7 input arguments for STIL: Filename, scId, "sitl", Epoch, PGSE data, DSL data, bitmask.';
                irf.log('critical', err_str);
                error('MATLAB:MMS_SDC_SDP_CDFWRITE', err_str);
            else
                % 7 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==3))
                    dce_xyz_pgse = varargin{5};
                else
                    err_str = 'SITL dce_xyz_pgse must be a numerical matrix with three columns.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==3))
                    dce_xyz_dsl = varargin{6};
                else
                    err_str = 'SITL dce_xyz_dsl must be a numerical matrix with three columns.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    bitmask = uint16(varargin{7});
                else
                    err_str = 'SITL bitmask must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
            end
            
        case('ql')
            if(nargin~=8)
                err_str = 'MMS_SDC_SDP_CDFWRITE: requires 8 input arguments for Quicklook: Filename, "ql", scId, Epoch, PGSE data, DSL data, bitmask, quality.';
                irf.log('critical', err_str);
                error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
            else
                % 8 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==3))
                    dce_xyz_pgse = varargin{5};
                else
                    err_str = 'QuickLook dce_xyz_pgse must be a numerical matrix with three columns.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==3))
                    dce_xyz_dsl = varargin{6};
                else
                    err_str = 'QuickLook dce_xyz_dsl must be a numerical matrix with three columns.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    bitmask = uint16(varargin{7}); % make sure it is stored with correct size
                else
                    err_str = 'QuickLook bitmask must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{8}) && size(varargin{8},2)==1))
                    qualityMark = uint16(varargin{8}); % make sure it is store with correct size
                else
                    err_str = 'QuickLook quality must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
            end
            
        case('usc')
            if(nargin~=9)
                err_str='MMS_SDC_SDP_CDFWRITE: requires 9 input arguments for Usc: Filename, scId, "usc", Epoch, PGSE data, DSL data, bitmask, quality.';
                irf.log('critical', err_str);
                error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
            else
                % 9 Arguments Ok.             
                if((isnumeric(varargin{5}) && size(varargin{5},2)==1))
                    escp = varargin{5};
                else
                    err_str='Usc ESCP must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{6}) && size(varargin{6},2)==1))
                    psp = varargin{6};
                else
                    err_str='Usc PSP must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{7}) && size(varargin{7},2)==1))
                    delta = varargin{7};
                else
                    err_str='Usc DELTA must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{8}) && size(varargin{8},2)==6))
                    psp_p = varargin{8};
                else
                    err_str='Usc PSP_P must be a numerical matrix with six columns.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
                if((isnumeric(varargin{9}) && size(varargin{9},2)==1))
                    bitmask = uint16(varargin{9}); % make sure it is stored with correct size
                else
                    err_str='Usc bitmask must be a numerical vector.';
                    irf.log('critical', err_str);
                    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
                end
            end
            
    otherwise
        err_str='MMS_SDC_SDP_CDFWRITE: Only "ql", "sitl" or "usc" are valid for types.';
        irf.log('critical', err_str);
        error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
        
    end
    
else
    err_str='MATLAB:MMS_SDC_SDP_CDFWRITE: Only "ql", "sitl" or "usc" are valid types and must be sent as string.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
    
end

if( isnumeric(varargin{4}) )
    EpochTimes = int64(varargin{4}); % Make sure it is of right size.
    
else
    err_str='MMS_SDC_SDP_CDFWRITE: Epoch times must be numerical.';
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);
    
end


% Check to see if compiled file exists if not build it.
if( ~exist(['mms_sdc_sdp_cdfwritec.', mexext],'file') )
    
    % Provide warning to users logfile to ensure they are know.
    irf.log('warning',['MMS_SDC_SDP_CDFWRITEC.', mexext,...
        ' file does not exist. Building it.']);
    
    % Locate the path to source and output directory.
    pathToSrc = which('mms_sdc_sdp_cdfwritec.c');
    [s, ~, ~] = fileparts(pathToSrc);
    
    % Call mex with the required arguments, if user has set up some other
    % mexopts.sh the following values will be used instead. If user runs on
    % a very old Matlab (7.11) there might be issues with linking to
    % Matlab's own lib.
    try
        mex(['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -D_GNU_SOURCE -pthread -fexceptions'],...
        ['-I',ENVIR.CDF_BASE,filesep,'include',filesep],...
        ['-L',ENVIR.CDF_BASE,filesep,'lib',filesep],...
        '-outdir',[s,filesep],'-lcdf', [s,filesep,'mms_sdc_sdp_cdfwritec.c']);
    catch err
        % If running at SDC, there will be NO GCC installed and this will
        % fail. MEXA64 to be built at IRFU or on dedicated machine at SDC
        % after filing a ticket for it at SDC.
        irf.log('critical', 'Mex failed. Are we running at SDC?');
        irf.log('critical', err.message);
        rethrow(err);
    end
    
end

% Actually call the mexa64 file.
switch(calledBy)
    case('sitl')
        try
            mms_sdc_sdp_cdfwritec(filename_output, scID, 'sitl', ...
                EpochTimes, dce_xyz_pgse', dce_xyz_dsl', bitmask);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if( strcmp(err.identifier,...
                    'MATLAB:mms_sdc_sdp_cdfwrite:filename_output:exists') )
                % If our cdfwrite code resulted in error write proper log.
                irf.log('critical', err.message);
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
            mms_sdc_sdp_cdfwritec(filename_output, scID, 'ql', ...
                EpochTimes, dce_xyz_pgse', dce_xyz_dsl', bitmask, ...
                qualityMark);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if( strcmp(err.identifier, ...
                    'MATLAB:mms_sdc_sdp_cdfwrite:filename_output:exists') )
                % If our cdfwrite code resulted in error write proper log.
                irf.log('critical', err.message);
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
            mms_sdc_sdp_cdfwritec(filename_output, scID, 'usc', ...
                EpochTimes, escp, psp, delta, psp_p', bitmask);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if( strcmp(err.identifier, ...
                    'MATLAB:mms_sdc_sdp_cdfwrite:filename_output:exists') )
                % If our cdfwrite code resulted in error write proper log.
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
        err_str='MATLAB:MMS_SDC_SDP_CDFWRITE: Only "ql", "sitl" or "usc" are valid types and must be sent as string.';
        irf.log('critical', err_str);
        error('Matlab:MMS_SDC_SDP_CDFWRITE:INPUT', err_str);

end

end
