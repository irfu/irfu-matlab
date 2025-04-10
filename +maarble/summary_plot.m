function summary_plot(fname,dataPath)
%SUMMARY_PLOT  Make a summary plot from a CAA-CDF file
%
% summary_plot(fileName,dataPath)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
%
% This software was developed as part of the MAARBLE (Monitoring,
% Analyzing and Assessing Radiation Belt Energization and Loss)
% collaborative research project which has received funding from the
% European Community's Seventh Framework Programme (FP7-SPACE-2011-1)
% under grant agreement n. 284520.

persistent Rth

%% Load data
if 0
  cd /Users/yuri/Dropbox/Projects/MAARBLE/WP3/HeadersAndExamples/CAA-Test-Files-Cluster-ULF/C1_CP_AUX_MAARBLE_ULF_PC12  %#ok<UNRCH>
  fname='C1_CP_AUX_MAARBLE_ULF_PC12__20101013_120000_20101013_150000_V130628.cdf';
end

dataDir='/data/themis';

if nargin<2, dataPath = '.'; end

iSep=strfind(fname,'__');
if isempty(iSep)
  irf.log('critical','invalid file name')
  error('invalid file name')
end


productName = fname(1:iSep-1);
d = dataobj([dataPath filesep fname]);

%% Construct ebsp

ebsp = struct('t',[],'f',[],'flagFac',1,...
  'bb_xxyyzzss',[],'ee_xxyyzzss',[],'ee_ss',[],...
  'pf_xyz',[],'pf_rtp',[],...
  'dop',[],'dop2d',[],'planarity',[],'ellipticity',[],'k_tp',[],...
  'fullB',[],'B0',[],'r',[]);

fields = {{'t','Time'},...
  {'f','Frequency'},...
  {'df','Frequency_BHW'},...
  {'bb_xxyyzzss','BB_xxyyzz_fac'},...
  {'k_tp','KSVD_fac'},...
  {'ee_ss','ESUM'},...
  {'dop','DOP'},...
  {'planarity','PLANSVD'},...
  {'ellipticity','ELLSVD'},...
  {'pf_rtp','PV_fac'},...
  {'B0','BMAG'}
  };
notValid = 0;
for iField = 1:length(fields)
  fieldName = fields{iField}{1}; varName = fields{iField}{2};
  tmpVar = get_variable(d,[varName '__' productName]);
  if isempty(tmpVar)
    % Specially treat GB data which are not in FAC system
    if strcmp(fieldName,'bb_xxyyzzss')
      tmpVar = get_variable(d,['BB_xxyyzz__' productName]);
      if isempty(tmpVar), continue
      else, ebsp.flagFac = 0;
      end
    else, continue
    end
  end
  if isnumeric(tmpVar.FILLVAL)
    tmpVar.data(tmpVar.data == tmpVar.FILLVAL) = NaN;
  end
  ebsp.(fieldName) = struct('data',tmpVar.data,'units',tmpVar.UNITS,'scale','');
  if isfield(tmpVar,'SCALETYP') && ischar(tmpVar.SCALETYP)
    ebsp.(fieldName).scale = tmpVar.SCALETYP;
  end
  if strcmp(fieldName,'bb_xxyyzzss')
    if numel(size(ebsp.(fieldName).data))==2 && all(all(isnan(ebsp.(fieldName).data)))
      % For no data Qtran produces one record with fill values
      irf.log('warning',sprintf('No B data for %s',fname))
      return
    else
      ebsp.(fieldName).data(:,:,4) = sum(ebsp.(fieldName).data(:,:,1:3),3);
    end
  end
  notValid = notValid + validate_data();
end
if notValid
  irf.log('critical',sprintf('Data failed validation for %s',fname))
  return
end

% Handling of position
if regexp(productName,'^C[1-4]_CP')==1
  ebsp.r = local.c_read(['R' productName(2)],ebsp.t.data([1 end]));
elseif regexp(productName,'^CC_CP_AUX_MAARBLE_TH[A-E]_[U,V]LF')==1
  if isempty(Rth)
    Rth = load(sprintf('%s%smRth.mat',dataDir,filesep), '-mat');
  end
  ebsp.r =Rth.(['Rth' lower(productName(21))]);
elseif regexp(productName,'^CC_CP_AUX_MAARBLE_G1[1-2]_ULF')==1
  try
    dd=dataobj([dataPath '/../../FACMATR/CDF/CC_CP_AUX_MAARBLE_G1' fname(21) '_ULF_FACMATR__' fname(33:47) '*.cdf']);
    ebsp.r = getmat(dd,['sc_pos_xyz_GSE__CC_CP_AUX_MAARBLE_G1' fname(21) '_ULF_FACMATR']);
  catch
    irf.log('warning','No FACMART file found')
  end
end

flagNoE = 0;
flagNoPolarization = 0;
if isempty(ebsp.ee_ss), flagNoE = 1; end
if isempty(ebsp.dop), flagNoPolarization = 1; end

%% plot
if flagNoE && flagNoPolarization
  h = irf_pl_ebsp(ebsp,{{'bb_xxyyzzss',1:4}});
elseif flagNoE
  limByDopStruct = struct('type','low','val',0.6,'param','dop','comp',1);
  limByPlanarityStruct = struct('type','low','val',0.6,'param','planarity','comp',1);
  params = {{'bb_xxyyzzss',4},...
    {'dop'},{'planarity'},...
    {'ellipticity',[],{limByDopStruct,limByPlanarityStruct}},...
    {'k_tp',[],{limByDopStruct,limByPlanarityStruct}},...
    };
  h = irf_pl_ebsp(ebsp,params);
else
  h = irf_pl_ebsp(ebsp);
end

ht=title(h(1),fname(1:end-4));
set(ht,'FontSize',9)

%% Save
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
set(gcf,'InvertHardCopy','off');
print('-r600','-dpng',fname(1:end-4))

%% Help functions
  function res = validate_data
    res = 0;
    switch fieldName
      case {'bb_xxyyzzss','ee_ss'}
        if any(any(any( ebsp.(fieldName).data < 0 )))
          res = 1;
          irf.log('critical',['Ivalid ' varName ' : negative values' ])
          return
        end
      case {'k_tp','pf_rtp'}
        c = tokenize(fieldName,'_'); c = c{2};
        idx = find(c=='r');
        if ~isempty(idx) && any(any( ebsp.(fieldName).data(:,:,idx) < 0 ))
          res = 1;
          irf.log('critical',['Ivalid ' varName ' : negative values of R' ])
        end
        idx = find(c=='t');
        if any( any( ebsp.(fieldName).data(:,:,idx) > 90 ) | ...
            any(ebsp.(fieldName).data(:,:,idx) < 0 ) )
          res = 1;
          irf.log('critical',...
            ['Ivalid ' varName ' : THETA outside range, expected 0<theta<90'])
        end
        idx = find(c=='p');
        if any( any( ebsp.(fieldName).data(:,:,idx) > 180 ) | ...
            any(ebsp.(fieldName).data(:,:,idx) < -180 ) )
          res = 1;
          irf.log('critical',...
            ['Ivalid ' varName ' : PHI outside range -180<phi<180'])
        end
        if res, return, end
      case {'dop','planarity'}
        if any(any( ebsp.(fieldName).data < 0 | ebsp.(fieldName).data > 1 ))
          res = 1;
          irf.log('critical',...
            ['Ivalid ' varName ' : values outside range, expected 0<val<1' ])
          return
        end
      case 'ellipticity'
        if any(any( ebsp.(fieldName).data < -1 | ebsp.(fieldName).data > 1 ))
          res = 1;
          irf.log('critical',...
            ['Ivalid ' varName ' : values outside range, expected -1<val<1' ])
          return
        end
      case 'B0'
        if any( ebsp.(fieldName).data < 0)
          res = 1;
          irf.log('critical',...
            ['Ivalid ' varName ' : negative values' ])
          return
        end
      otherwise
        return
    end
  end
end