function [ filename_output ] = mms_cdf_writing( dataOut, bitmask, HeaderInfo, qualityMark )
% MMS_CDF_WRITING writes the dataobj dataOut to the corresponding CDF file.
%	[filename_output] = MMS_CDF_WRITING( dataOut, bitmask, HeaderInfo, qualityMark) will
%       write an MMS CDF file containing the dataobj dataOut's data to a temporary output folder
%	defined by ENVIR.DROPBOX_ROOT. The struct HeaderInfo will help determine which output
%	file is to be created.
%
%	Example:
%		filename_output = mms_cdf_writing( dataOut, bitmask, HeaderInfo, qualityMark);
%
%	Note 1: It is assumed that other SDC processing scripts will move the created output
%	file to its final destination (from /ENIVR.DROPBOX_ROOT/ to /path/as/defined/by/
%	SDC Developer Guide, v1.7).
%	Note 2: 20140204, FIELDS consider SITL to be a dataLevel even if the
%	SDC Developer Guide is unclear about that issue. Therefor when running
%	MMS SDP SITL processing this script will check corresponding folder structure
%	and write filenames accordingly.
%
%       NOTE 3: The number of input arguments must be altered when the processing actually takes
%       place in previous functions. For now it will only multiply the input dataObj data to fill
%       CDF files of various types.
%
% 	See also MMS_CDF_IN_PROCESS, MMS_BITMASKING, MMS_INIT.

% Verify that we have all information requried.
narginchk(3,4);

global ENVIR;
global MMS_CONST;

% HeaderInfo contains information about what is to be written, as if it is
% SDP_L1B to SITL_DCE or SDP_L1B to Quicklook or DCV_L1B to Usc. It also
% contains information about source files used ("Parents") and what data
% type (fast / slow / burst) etc.


switch(HeaderInfo.calledBy)
    case('sitl_dce')
        %disp('sitl_dce');
        % List existing files. SITL is has it bottommost level as monthly directories.
        preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, '/science/', HeaderInfo.scId, ...
            '/',HeaderInfo.instrumentId, '/', 'sitl/', HeaderInfo.dataLevel, '/'...
            HeaderInfo.startTime(1:4), '/', HeaderInfo.startTime(5:6), '/', ...
	    HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_sitl_',HeaderInfo.dataLevel,'_dce2d_', ...
            HeaderInfo.startTime, '_v*.cdf']);

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(size(preExistingFiles,1), '%u'));
        else
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(MMS_CONST.Version.Z));
        end

        % Create the new output filename. (excl extension).
        filename_output = [HeaderInfo.scId, '_', HeaderInfo.instrumentId, ...
            '_sitl_',HeaderInfo.dataLevel,'_dce2d_', HeaderInfo.startTime, ...
            '_', versionNum];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYImms2_sdp_slow_l2_uscdcv_20150410000000_v0.0.0NG TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.

        cd(ENVIR.DROPBOX_ROOT);

        % FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = getv(dataOut, dataOut.vars{1,1}); % The epoch times
        data1 = getv(dataOut, dataOut.vars{5,1}); % The data
        
        irf.log('debug',['MATLAB:mms_cdf_writing:sitl Ready to write data to temporary file in DROPBOX_ROOT/', filename_output,'.cdf']);
        try
            mms_sdc_cdfwrite(filename_output, int8(str2num(HeaderInfo.scId(end))), 'sitl', epochTT.data, data1.data', data1.data', uint16(bitmask));
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
                % It was an unkown error..
                irf.log('critical',['Error when writing CDF file', err.identifier, err.message]);
                rethrow(err);
            end
        end % End of try
        
        
    case('ql')
        % List existing files.
        % Quicklook may have its bottommost level as monthly or daily directories depending on which dataproduct it is. For SRVY and/or SITL bottommost are monthly, for brst (and fast and slow?) they are daily. Note SITL is run with sitl_dce as HeaderInfo.calledby and should not end up here.
	if(strcmp(HeaderInfo.dataMode,'srvy'))
           preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, '/science/', HeaderInfo.scId, ...
            '/',HeaderInfo.instrumentId, '/', HeaderInfo.dataMode, '/', ...
	    HeaderInfo.calledBy ,'/', ...
            HeaderInfo.startTime(1:4), '/', HeaderInfo.startTime(5:6), '/', ...
	    HeaderInfo.scId, '_', HeaderInfo.instrumentId, '_', ...
	    HeaderInfo.dataMode, '_ql_dce2d_', ...
            HeaderInfo.startTime, '_v*.cdf']);
	else
           preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, '/science/', HeaderInfo.scId, ...
            '/',HeaderInfo.instrumentId, '/', HeaderInfo.dataMode, '/', ...
	    HeaderInfo.calledBy ,'/', ...
            HeaderInfo.startTime(1:4), '/', HeaderInfo.startTime(5:6), '/', ...
            HeaderInfo.startTime(7:8), '/', HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, '_ql_dce2d_', ...
            HeaderInfo.startTime, '_v*.cdf']);
	end

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(size(preExistingFiles,1), '%u'));
        else
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(MMS_CONST.Version.Z));
        end

        % Create the new output filename. (excl extension).
        filename_output = [HeaderInfo.scId, '_', HeaderInfo.instrumentId, ...
            '_', HeaderInfo.dataMode, '_ql_dce2d_', HeaderInfo.startTime, ...
            '_', versionNum];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.

        cd(ENVIR.DROPBOX_ROOT);

        % FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = getv(dataOut, dataOut.vars{1,1}); % The epoch times
        data1 = getv(dataOut, dataOut.vars{5,1}); % The data
        
        irf.log('debug',['MATLAB:mms_cdf_writing:ql Ready to write data to temporary file in DROPBOX_ROOT/', filename_output,'.cdf']);
        try
            mms_sdc_cdfwrite(filename_output, int8(str2num(HeaderInfo.scId(end))), 'ql', epochTT.data, data1.data', data1.data', uint16(bitmask), uint16(qualityMark));
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
        %disp('usc');
        % List existing files.
	% Survey have the bottommost level as month while Burst (fast, slow?) have their as day directories. Sitl should not end up here.
        if(strcmp(HeaderInfo.dataMode,'srvy'))
          preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, '/science/', HeaderInfo.scId, ...
            '/',HeaderInfo.instrumentId, '/', HeaderInfo.dataMode, '/l2/', ...
            HeaderInfo.startTime(1:4), '/', HeaderInfo.startTime(5:6), '/', ...
	    HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, '_l2_uscdcv_', ...
            HeaderInfo.startTime, '_v*.cdf']);
	else
          preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, '/science/', HeaderInfo.scId, ...
            '/',HeaderInfo.instrumentId, '/', HeaderInfo.dataMode, '/l2/', ...
            HeaderInfo.startTime(1:4), '/', HeaderInfo.startTime(5:6), '/', ...
            HeaderInfo.startTime(7:8), '/', HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, '_l2_uscdcv_', ...
            HeaderInfo.startTime, '_v*.cdf']);
	end

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(size(preExistingFiles,1), '%u'));
        else
            versionNum = strcat('v',num2str(MMS_CONST.Version.X),'.',num2str(MMS_CONST.Version.Y),'.',num2str(MMS_CONST.Version.Z));
        end

        % Create the new output filename. (excl extension).
        filename_output = [HeaderInfo.scId, '_', HeaderInfo.instrumentId, ...
            '_', HeaderInfo.dataMode, '_l2_uscdcv_', HeaderInfo.startTime, ...
            '_', versionNum];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.

        cd(ENVIR.DROPBOX_ROOT);

        % FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = getv(dataOut, dataOut.vars{1,1}); % The epoch times
        data1 = getv(dataOut, dataOut.vars{5,1}); % The data
     %%%%%   
     
        psp_p = [data1.data, data1.data];
        irf.log('debug',['MATLAB:mms_cdf_writing:usc Ready to write data to temporary file in DROPBOX_ROOT/', filename_output,'.cdf']);
        try
            mms_sdc_cdfwrite(filename_output, int8(str2num(HeaderInfo.scId(end))), 'usc', epochTT.data, data1.data(:,1), data1.data(:,2), data1.data(:,3), psp_p', uint16(bitmask));
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
        irf.log('warning','MATLAB:mms_cdf_writing:HeaderInfo.calledBy unknown or not implemented yet. As of now, "ql", "sitl" and "usc" exists.');
end


% % % % % %
% DISABLED: AS MathWork's MatLab 2013b does not handle TT2000 CDF files 
%  properly using this has been disabled.
% % % % % %
% Create a new output file from the skeleton source file. This way all
% global parameters etc. are already created as expected. Catch status and
% cmdoutput in order for it not to write to screen, which would be at SDC.
%
%[status, cmdout] = unix([ENVIR.CDF_BASE,'/bin/skeletoncdf -cdf',' ', filename_output,'.cdf ',ENVIR.DATA_PATH_ROOT,'/MMS_sdc_skt-20131029/mms1_sdp_dce.skt']);
% Verify it was created correctly: %INFO_skeleton = cdfinfo(filename_output);
% Correct format for number of records, but not for time variable.
%cdfwrite(filename_output, {'Epoch', num2cell(dataoutCDFobj.data.mms1_sdp_epoch_dce.data), dataoutCDFobj.vars{2}, num2cell(dataoutCDFobj.data.mms1_sdp_timetag_dce.data), dataoutCDFobj.vars{3}, num2cell(dataoutCDFobj.data.mms1_sdp_samplerate_dce.data), dataoutCDFobj.vars{4}, dataoutCDFobj.data.DCE_LABL_1.data, dataoutCDFobj.vars{5}, num2cell(dataoutCDFobj.data.mms1_sdp_dce_sensor.data,2)}, 'TT2000', true);





% Update some of the global parameters that are not static.

% Generation date is today (when script runs).
GATTRIB.Generation_date = {0, 'CDF_CHAR', datestr(now,'yyyymmdd')};

% Data version is the version number. Version number should be "X.Y.Z"
GATTRIB.Data_version = {0, 'CDF_CHAR', versionNum};

% FIXME: ADD Generated by as gitversion of software.
irfVersion = irf('version');
GATTRIB.Generated_by = {0, 'CDF_CHAR',irfVersion};

% Parents is the source file logical id, if multiple sources add subsequent
% entries for each source file. 
% Multiple.Parents={0, 'CDF_CHAR', [['CDF>',source1CDF.GloablAttributes.Logigal_file_id{1,1}]; ['CDF>',source2CDF.GloablAttributes.Logigal_file_id{1,1}]]}
%GATTRIB.Parents = {0, 'CDF_CHAR', ['CDF>',dataOut.GlobalAttributes.Logical_file_id{1,1}]};


if(HeaderInfo.numberOfSources == 1)
    GATTRIB.Parents = {0, 'CDF_CHAR', ['CDF>',HeaderInfo.parents_1]};
elseif(HeaderInfo.numberOfSources == 2)
    diffLength = length(HeaderInfo.parents_1)-length(HeaderInfo.parents_2);
    if(diffLength>0)
        GATTRIB.Parents = {0, 'CDF_CHAR', [['CDF>',HeaderInfo.parents_1]; ['CDF>',HeaderInfo.parents_2, blanks(diffLength)]]};
    elseif(diffLength<0)
        GATTRIB.Parents = {0, 'CDF_CHAR', [['CDF>',HeaderInfo.parents_1, blanks(abs(diffLength))]; ['CDF>',HeaderInfo.parents_2]]};
    else
        GATTRIB.Parents = {0, 'CDF_CHAR', [['CDF>',HeaderInfo.parents_1]; ['CDF>',HeaderInfo.parents_2]]};
    end
end

% Update all the new values to GlobalAttributes
irf.log('debug','MATLAB:mms_cdf_writing:UpdatingGlobalAttributes');
cdfupdate(filename_output,'GlobalAttributes',GATTRIB);


end
