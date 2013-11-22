%function res = summary_plot(fname)

%% Load data
cd /Users/yuri/Dropbox/Projects/MAARBLE/WP3/HeadersAndExamples/CAA-Test-Files-Cluster-ULF/C1_CP_AUX_MAARBLE_ULF_PC12
fname='C1_CP_AUX_MAARBLE_ULF_PC12__20101013_120000_20101013_150000_V130628.cdf';

iSep=strfind(fname,'__');
if isempty(iSep)
  irf_log('error','invalid file name')
  error('invalid file name')
end


productName = fname(1:iSep-1);
d = dataobj(fname);

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
  if isempty(tmpVar), continue, end
  if isnumeric(tmpVar.FILLVAL)
    tmpVar.data(tmpVar.data == tmpVar.FILLVAL) = NaN;
  end
  ebsp.(fieldName) = struct('data',tmpVar.data,'units',tmpVar.UNITS);
  if strcmp(fieldName,'bb_xxyyzzss')
    ebsp.(fieldName).data(:,:,4) = sum(ebsp.(fieldName).data(:,:,1:3),3);
  end
end

%% plot

h = irf_pl_ebsp(ebsp);