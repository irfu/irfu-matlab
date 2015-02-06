classdef mms_sdc_sdp_db < handle
  %MMS_DSC_SDP_DB Class containing raw data 
  %   Detailed explanation goes here
  
  properties (Constant)
    % Files of valid/supported files
    files = {'bdcv','bdce','dcv','dce','hk101','hk105','hk10d','hk10e'};
  end
  properties
    % All these are cell arrays, as generally there can be multiple files
    % spanning a burst. For survey and HK data these must be daily files.
    % Empty variables mean no attempt to load, NaN means failed load (no
    % file).
    bdcv  % Burst files
    bdce
    dcv   % Survey files fast/slow
    dce
  end
  properties (SetAccess = private)
    hk101 % sun pulse, required for SITL/QL
    hk105 % FLDHK4, e.g. sweep
    hk10d % FLDHK8
    hk10e % FLDHK64, e.g. bias, guard settings
    scId
    fastSurvey % fast=1, slow=0
    tint % time interval defined by input file(s)
  end
  
  methods
    function obj = mms_sdc_sdp_db()
    % Constructor
      narginchk(0,0)
      obj.scId = []; 
      obj.fastSurvey = [];
      obj.tint = [];
      for f=mms_sdc_sdp_db.files
        obj.(f{:}) = {};
      end
    end
    
    function [varargout] = subsref(obj,idx)
      %SUBSREF handle indexing
        switch idx(1).type
          % Use the built-in subsref for dot notation
          case '.'
            % For exmpty fileds try to load
            if ~isempty(intersect(idx(1).subs,mms_sdc_sdp_db.files)) && ...
                isempty(obj.(idx(1).subs))
              load_files(obj,idx(1).subs);
            end
            [varargout{1:nargout}] = builtin('subsref',obj,idx);
          case '()'
            [varargout{1:nargout}] = builtin('subsref',obj,idx);
          case '{}'
            [varargout{1:nargout}] = builtin('subsref',obj,idx);
        end
    end
    
    function load_files(obj,fileId)
      % Load file Id
      if ~isempty(obj.tint)
        fprintf('loading %s',fileId)
        obj.(fileId) = NaN;
      end
    end
    
    function set.bdcv(obj,fileName)
      obj.bdcv = set_file_private(obj,'bdcv',fileName);
    end
    function set.bdce(obj,fileName)
      obj.bdce = set_file_private(obj,'bdce',fileName);
    end
    function set.dcv(obj,fileName)
      obj.dcv = set_file_private(obj,'dcv',fileName);
    end
    function set.dce(obj,fileName)
      obj.dce = set_file_private(obj,'dce',fileName);
    end
  
    function ok = validate_file_name(obj,prod,fileName)
      % validate L1b CDF file names
      % HK: mms1_fields_hk_l1b_101_20130904_v0.0.0.cdf	
      % B:  mms1_sdp_brst_l1b_dcv_20130904212500_v0.0.2.cdf
      % or possibly mms1_edp_fast_l1b_dce_20160101_v2.0.1.cdf (sdp->edp)
      if isempty(intersect(prod,mms_sdc_sdp_db.files))
        irf.log('critical','Invalid fileId')
        error('irf:mms_sdc_sdp_db:validate_file_name','Invalid fileId')
      end
      ok = false;
      toks = tokenize(fileName,'_');
      if ~strcmp(toks{1}(1:3),'mms'), return, end
      mmsId =  str2double(toks{1}(4));
      if mmsId>4 || mmsId<1, return, end
      if ~isempty(obj.scId) && obj.scId~=mmsId
        irf.log('warning',['sc id mismatch, expecting ' num2str(obj.scId)])
        return
      else obj.scId=mmsId;
      end
      if strcmp(prod(1:2),'hk')
        if ~strcmp(toks{2},'fields'), return, end
        if ~strcmp(toks{3},'hk'), return, end
        if ~strcmp(toks{4},'l1b'), return, end
        if ~strcmp(toks{5},prod(3:5)), return, end
        ok = true; return
      end
      % now DCE & DCV
      if ~strcmp(toks{2},'sdp') && ~strcmp(toks{2},'edp'), return, end
      if ~strcmp(toks{4},'l1b'), return, end
      
      if strcmp(prod(1:3),'bdc')
        if ~strcmp(toks{3},'brst'), return, end
        if ~strcmp(toks{5},prod(2:4)), return, end
      elseif strcmp(prod(1:2),'dc')
        if isempty(intersect(toks{3},{'slow','fast'}))
          return
        end
        if isempty(obj.fastSurvey) && strcmp(toks{3},'fast')
          obj.fastSurvey = true;
        else obj.fastSurvey = false;
        end
        if ~strcmp(toks{5},prod), return, end
      else % not DCE | DCV
        error('irf:mms_sdc_sdp_db:validate_file_name: should not be here')
      end
      ok = true;
    end
  end
  
  methods (Access=private)
    function res = set_file_private(obj,fileId,fullFileName)
      res = {};
      if isempty(fullFileName), return, end
      if ~isempty(obj.(fileId))
        irf.log('critical',[fileId ' already set'])
        error('irf:mms_sdc_sdp_db:set','already set')
      end
      % create dataobj from full file name 
      if ~ischar(fullFileName)
        errStr = 'Expecting string input (path+filename)';
        irf.log('critical',errStr);
        error('irf:mms_sdc_sdp_db:set_file_private',errStr)
      end
      [~, fileName, ~] = fileparts(fullFileName);
      if ~validate_file_name(obj,fileId,fileName)
        errStr = 'Invalid filename)';
        irf.log('critical',errStr);
        error('irf:mms_sdc_sdp_db:set_file_private',errStr)
      end
        
      try 
        % Create dataobj
        res = {struct('fileName',fileName,'data',dataobj(fullFileName, 'KeepTT2000'))};
        % Update time
        if fileId(1)=='b'
          epochVar = sprintf('mms%d_sdp_epoch_dc%s',obj.scId,fileId(end));
        else epochVar = 'Epoch';
        end
        tintTmp = res{1}.data.data.(epochVar).data([1 end]);
        if isempty(obj.tint), obj.tint = tintTmp;
        else
          if obj.tint(1) > tintTmp(1), obj.tint(1) = tintTmp(1); end
          if obj.tint(2) < tintTmp(2), obj.tint(2) = tintTmp(2); end
        end
      catch err
        errMsg = sprintf('failed to load %s file: %s',fileId,fullFileName);
        irf.log('critical',errMsg)
        rethrow(err)
      end
    end
  end
end

