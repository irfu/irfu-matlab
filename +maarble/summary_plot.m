function summary_plot(fname,dataPath)
%SUMMARY_PLOT  Make a summary plot from a CAA-CDF file
%
% summary_plot(fileName,dataPath)

% ----------------------------------------------------------------------------
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

%% Load data
if 0
  cd /Users/yuri/Dropbox/Projects/MAARBLE/WP3/HeadersAndExamples/CAA-Test-Files-Cluster-ULF/C1_CP_AUX_MAARBLE_ULF_PC12
  fname='C1_CP_AUX_MAARBLE_ULF_PC12__20101013_120000_20101013_150000_V130628.cdf';
end

if nargin<2, dataPath = '.'; end

iSep=strfind(fname,'__');
if isempty(iSep)
  irf_log('error','invalid file name')
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

for iField = 1:length(fields)
  fieldName = fields{iField}{1};
  varName = fields{iField}{2};
  tmpVar = get_variable(d,[varName '__' productName]);
  if isempty(tmpVar)
    % Specially treat GB data which are not in FAC system
    if strcmp(fieldName,'bb_xxyyzzss')
      tmpVar = get_variable(d,['BB_xxyyzz__' productName]);
      if isempty(tmpVar), continue
      else ebsp.flagFac = 0;
      end
    else continue
    end
  end
  if isnumeric(tmpVar.FILLVAL)
    tmpVar.data(tmpVar.data == tmpVar.FILLVAL) = NaN;
  end
  ebsp.(fieldName) = struct('data',tmpVar.data,'units',tmpVar.UNITS);
  if strcmp(fieldName,'bb_xxyyzzss')
    ebsp.(fieldName).data(:,:,4) = sum(ebsp.(fieldName).data(:,:,1:3),3);
  end
end

%XXX TODO: add handling of position

flagNoE = 0;
flagNoPolarization = 0;
if isempty(ebsp.ee_ss), flagNoE = 1; end
if isempty(ebsp.dop), flagNoPolarization = 1; end

%% plot
if flagNoE && flagNoPolarization
  h = irf_pl_ebsp(ebsp,{{'bb_xxyyzzss',1:4}});
elseif flagNoE
  params = {{'bb_xxyyzzss',4},...
    {'dop'},{'planarity'},...
    {'ellipticity',[],{limByDopStruct,limByPlanarityStruct}},...
    {'k_tp',[],{limByDopStruct,limByPlanarityStruct}},...
    };
  h = irf_pl_ebsp(ebsp,params);
else
  h = irf_pl_ebsp(ebsp);
end

ht=title(h(1),fname);
set(ht,'FontSize',8)

%% Save
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print('-r600','-dpng',fname(1:end-4))