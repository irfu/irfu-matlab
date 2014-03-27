classdef mms_sdc_sdp_db
  %MMS_DSC_SDP_DB Class containing raw data 
  %   Detailed explanation goes here
  
  properties
    % All these are cell arrays, as generally there can be multiple files
    % spanning a burst. For survey and HK data these must be daily files.
    bdcv  % Burst files
    bdce
    dcv   % Survey files fast/slow
    dce
    hk101 % sun pulse, required for SITL/QL
    hk105 % FLDHK4, e.g. sweep
    hk10d % FLDHK8
    hk10e % FLDHK64, e.g. bias, guard settings
  end
  properties (SetAccess = private)
    scId
    fastSurvey % fast=1, slow=0
  end
  properties (Dependent, SetAccess = private)
    tint
  end
  
  methods
    function dbObj = mms_sdc_sdp_db(dcvFile,dceFile)
    % Constructor
      narginchk(2,2)
      dbObj.scId = []; dbObj.fastSurvey = [];
      modeTmp = '';
      props = {'dcv','dce'}; files = {dcvFile,dceFile};
      for i=1:length(props)
        if ~isempty(files{i})
          try
            [~, fileName, ~] = fileparts(files{i});
            [ok,mmsId,mode] = validate_file_name(props{i});
            if ~ok
              errMsg = sprintf('invalid %s filename: %s',props{i},fileName);
              irf.log('critical',errMsg)
              error(errMsg) %#ok<SPERR>
            end
            if isempty(dbObj.scId), dbObj.scId = mmsId;
            else
              if mmsId~=dbObj.scId
                errMsg = sprintf('invalid %s filename: %s, expecting : mms%d',...
                  props{i},fileName,dbObj.scId);
                irf.log('critical',errMsg)
                error(errMsg) %#ok<SPERR>
              end
            end
            % XXX: implement modes for fast/slow + burst
            if isempty(modeTmp), modeTmp = mode;
            else
              if ~strcmp(modeTmp,mode)
                errMsg = sprintf('invalid %s filename: %s, expecting : %s',...
                  props{i},fileName,modeTmp);
                irf.log('critical',errMsg)
                error(errMsg) %#ok<SPERR>
              end
            end
            dbObj.(props{i}).data = dataobj(files{i});
          catch err
            errMsg = sprintf('failed to load %s file: %s',props{i},files{i});
            irf.log('critical',errMsg)
            rethrow(err)
          end
        end
      end
      function [ok,mmsId,mode] = validate_file_name(prod)
        ok = false; mmsId = 0; mode = '';
        toks = tokenize(fileName,'_');
        if ~strcmp(toks{1}(1:3),'mms'), return, end
        mmsId =  str2double(toks{1}(4));
        if mmsId>4 || mmsId<1, return, end
        if ~strcmp(toks{2},'sdp'), return, end
        mode = toks{3};
        if isempty(intersect(mode,{'slow','fast','burst'}))
          return
        end
        if ~strcmp(toks{4},'l1b'), return, end
        if ~strcmp(toks{5},prod), return, end
      end
    end
    function res = get.tint(obj)
    % Return time interval spanning the files
      res = [];
      if ~isempty(obj.dcv)
        res = obj.dcv.data.(sprintf('mms%s_sdp_epoch_dce',obj.scId)).data([1 end]);
      end
      if ~isempty(obj.dcv)
        tmpInt = obj.dcv.data.(sprintf('mms%s_sdp_epoch_dce',obj.scId)).data([1 end]);
        if isempty(res), res = tmpInt;
        else
          if tmpInt(1)<res(1), res(1) =  tmpInt(1); end
          if tmpInt(2)>res(2), res(2) =  tmpInt(2); end
        end
      end
    end
  end
  
end

