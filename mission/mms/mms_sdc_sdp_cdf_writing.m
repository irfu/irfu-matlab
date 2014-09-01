function [ filename_output ] = mms_sdc_sdp_cdf_writing( HeaderInfo )
% MMS_SDC_SDP_CDF_WRITING writes the data to the corresponding CDF file.
%
%	filename_output = MMS_SDC_SDP_CDF_WRITING( HeaderInfo)
%   will write an MMS CDF file containing the data stored to a temporary 
%   output folder defined by ENVIR.DROPBOX_ROOT. The struct HeaderInfo will
%   help determine which output file is to be created.
%
%   Example:
%   filename_output = mms_sdc_sdp_cdf_writing( HeaderInfo);
%
%	Note 1: It is assumed that other SDC processing scripts will move the
%   created output file to its final destination (from /ENIVR.DROPBOX_ROOT/
%   to /path/as/defined/by/	SDC Developer Guide, v1.7).
%
% 	See also MMS_SDC_SDP_CDF_IN_PROCESS, MMS_SDC_SDP_INIT.

% Verify that we have all information requried.
narginchk(1,1);

global ENVIR;
global MMS_CONST;

% Simply recall all data from memory.
global DataInMemory; % TO BE CHANGED LATER ON! ! ! !

% HeaderInfo contains information about what is to be written, as if it is
% SDP_L1B to SITL_DCE or SDP_L1B to Quicklook or DCV_L1B to Usc. It also
% contains information about source files used ("Parents") and what data
% type (fast / slow / burst) etc.

oldDir = pwd;

switch(HeaderInfo.calledBy)
    case('sitl')
        % List existing files.
        % SITL is has it bottom most level as monthly directories.
        preExistingFiles = dir([ENVIR.DATA_PATH_ROOT, filesep,'science',...
            filesep, HeaderInfo.scId, filesep, HeaderInfo.instrumentId, ...
            filesep, HeaderInfo.dataMode, filesep, 'sitl', filesep, ...
            HeaderInfo.startTime(1:4), filesep, ...
            HeaderInfo.startTime(5:6), filesep, HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
            '_sitl_dce2d_', HeaderInfo.startTime, '_v*.cdf']);

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(size(preExistingFiles,1), '%u') ];
        else
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(MMS_CONST.Version.Z) ];
        end

        % Create the new output filename. (excl extension).
        filename_output = [ HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
            '_sitl_dce2d_', HeaderInfo.startTime, '_', versionNum ];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.
        
        cd(ENVIR.DROPBOX_ROOT);

        %% FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = DataInMemory.dce.time;
        data1(:,1) = DataInMemory.dce.e12.data;
        data1(:,2) = DataInMemory.dce.e34.data;
        data1(:,3) = DataInMemory.dce.e56.data;
        bitmask = DataInMemory.dce.e12.bitmask;
        
        %%%
        irf.log('debug',['MATLAB:mms_sdc_sdp_cdf_writing:sitl Ready to',...
            ' write data to temporary file in DROPBOX_ROOT/', ...
            filename_output,'.cdf']);
        mms_sdc_sdp_cdfwrite( filename_output, ...
            int8(str2double(HeaderInfo.scId(end))), 'sitl', epochTT, ...
            data1, data1, uint16(bitmask) );
        
        
    case('ql')
        % List existing files.
        % Quicklook may have its bottommost level as monthly or daily 
        % directories depending on which dataproduct it is. For SRVY and/or
        % SITL bottommost are monthly, for brst (and fast and slow?) they 
        % are daily. Note SITL is run with sitl_dce as HeaderInfo.calledby
        % and should not end up here.
        if( strcmp(HeaderInfo.dataMode,'srvy') )
           preExistingFiles = dir( [ENVIR.DATA_PATH_ROOT, filesep, ...
               'science', filesep, HeaderInfo.scId, filesep, ...
               HeaderInfo.instrumentId, filesep, HeaderInfo.dataMode, ...
               filesep, HeaderInfo.calledBy, filesep, ...
               HeaderInfo.startTime(1:4), filesep, ...
               HeaderInfo.startTime(5:6), filesep, HeaderInfo.scId, ...
               '_', HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
               '_ql_dce2d_', HeaderInfo.startTime, '_v*.cdf'] );
        else
           preExistingFiles = dir( [ENVIR.DATA_PATH_ROOT, filesep, ...
               'science', filesep, HeaderInfo.scId, filesep, ...
               HeaderInfo.instrumentId, filesep, HeaderInfo.dataMode, ...
               filesep, HeaderInfo.calledBy, filesep, ...
               HeaderInfo.startTime(1:4), filesep, ...
               HeaderInfo.startTime(5:6), filesep, ...
               HeaderInfo.startTime(7:8), filesep, HeaderInfo.scId, ...
               '_', HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
               '_ql_dce2d_', HeaderInfo.startTime, '_v*.cdf'] );
        end

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(size(preExistingFiles,1), '%u') ];
        else
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(MMS_CONST.Version.Z) ];
        end

        % Create the new output filename. (excl extension).
        filename_output = [ HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
            '_ql_dce2d_', HeaderInfo.startTime, '_', versionNum ];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.

        cd(ENVIR.DROPBOX_ROOT);

        %% FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = DataInMemory.dce.time;
        data1(:,1) = DataInMemory.dce.e12.data;
        data1(:,2) = DataInMemory.dce.e34.data;
        data1(:,3) = DataInMemory.dce.e56.data;
        bitmask = DataInMemory.dce.e12.bitmask;
        %%%
        irf.log('debug',['MATLAB:mms_sdc_sdp_cdf_writing:ql Ready to',...
            ' write data to temporary file in DROPBOX_ROOT/', ...
            filename_output,'.cdf']);
        
        mms_sdc_sdp_cdfwrite( filename_output, ...
            int8(str2double(HeaderInfo.scId(end))), 'ql', epochTT, ...
            data1, data1, uint16(bitmask), ...
            uint16(mms_sdc_sdp_bitmask2quality('e',bitmask)) );
        
        
    case('usc')
        %disp('usc');
        % List existing files.
        % Survey have the bottommost level as month while Burst (fast, 
        % slow?) have their as day directories.
        if(strcmp(HeaderInfo.dataMode,'srvy'))
            preExistingFiles = dir( [ENVIR.DATA_PATH_ROOT, filesep, ...
                'science', filesep, HeaderInfo.scId, filesep, ...
                HeaderInfo.instrumentId, filesep, HeaderInfo.dataMode, ...
                filesep, 'l2', filesep, HeaderInfo.startTime(1:4), ...
                filesep, HeaderInfo.startTime(5:6), filesep, ...
                HeaderInfo.scId, '_', HeaderInfo.instrumentId, '_', ...
                HeaderInfo.dataMode, '_l2_uscdcv_', ...
                HeaderInfo.startTime, '_v*.cdf'] );
        else
            preExistingFiles = dir( [ENVIR.DATA_PATH_ROOT, filesep, ...
                'science', filesep, HeaderInfo.scId, filesep, ...
                HeaderInfo.instrumentId, filesep, HeaderInfo.dataMode, ...
                filesep, 'l2', filesep, HeaderInfo.startTime(1:4), ...
                filesep, HeaderInfo.startTime(5:6), filesep, ...
                HeaderInfo.startTime(7:8), filesep, HeaderInfo.scId, ...
                '_', HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
                '_l2_uscdcv_', HeaderInfo.startTime, '_v*.cdf'] );
        end

        if(size(preExistingFiles,1)>0)
            % Next version number is number of preExistingFiles + 1 (one).
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(size(preExistingFiles,1), '%u') ];
        else
            versionNum = [ 'v', num2str(MMS_CONST.Version.X), '.', ...
                num2str(MMS_CONST.Version.Y), '.', ...
                num2str(MMS_CONST.Version.Z) ];
        end

        % Create the new output filename. (excl extension).
        filename_output = [HeaderInfo.scId, '_', ...
            HeaderInfo.instrumentId, '_', HeaderInfo.dataMode, ...
            '_l2_uscdcv_', HeaderInfo.startTime, '_', versionNum];
        
        % NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS 
        % CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING 
        % DIRECTORY WHEN EXECUTING.
        
        cd(ENVIR.DROPBOX_ROOT);

        %% FIXME: DUMMY DATA FOR NOW.
        % For now store data temporarly
        epochTT = DataInMemory.dcv.time;
        data1(:,1) = DataInMemory.dcv.v1.data;
        data1(:,2) = DataInMemory.dcv.v3.data;
        data1(:,3) = DataInMemory.dcv.v5.data;     
        psp_p = [data1, data1];
        bitmask = DataInMemory.dcv.v1.bitmask;
        %%%%
        
        irf.log('debug', ['MATLAB:mms_sdc_sdp_cdf_writing:usc Ready to',...
            ' write data to temporary file in DROPBOX_ROOT/', ...
            filename_output,'.cdf']);
        
        mms_sdc_sdp_cdfwrite( filename_output, ...
            int8(str2double(HeaderInfo.scId(end))), 'usc', epochTT, ...
            data1(:,1), data1(:,2), data1(:,3), psp_p, uint16(bitmask) );    
    
  otherwise
    errStr = 'unrecognized HeaderInfo.calledBy';
    irf.log('critical', errStr);
    error(errStr)
end


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
irf.log('debug','MATLAB:mms_sdc_sdp_cdf_writing:UpdatingGlobalAttributes');
cdfupdate(filename_output,'GlobalAttributes',GATTRIB);

% Return to previous working directory.
cd(oldDir);

end
