function Out = wave_psd_ell_thetak(TT)
%MAARBLE.WAVE_PSD_ELL_THETAK  Extract PSD, ELL, thetaK for a TT
%
% Out = wave_psd_ell_thetak(TT)
%
% TT is a MAARBLE TimeTable containing frequency limits for a banded
% wave element in entry description
%
% Out is a structure having the following fileds:
%    time            - time vector
%    fPeak           - peak frequency
%    psdBpeak        - PSD B at wave peak
%    ellipticity     - ellipticity at wave peak
%    thetaK          - wave normal angle (deg) at wave peak
%    gseR, mlt, mLat - position

ULF_PATH = '/data/caa/MAARBLE/WaveDatabase/ULF/';
cl_s=''; th_s = '';
chk_input()
if ~isempty(th_s) % load THEMIS positions 1min resolution
  Rth = load('/data/caalocal/THEMIS/mRth',['Rth' lower(th_s)]);
end

Out = struct('time',[],'fPeak',[],'psdBpeak',[],'ellipticity',[],'thetaK',[],...
  'gseR',[],'mlt',[],'mLat',[],...
  'header',TT.Header{1},'created',irf_time(now,'date>utc'));

oldPwd = pwd;
for ievent=1:numel(TT)
  fprintf('Event #%d (out of %d) : %s\n',ievent,numel(TT),...
    irf_disp_iso_range(TT.TimeInterval(ievent,:),1));
  lowerFreqBound = str2double(TT.Description{ievent}{1});
  upperFreqBound = str2double(TT.Description{ievent}{2});
  flagSkipPC12 = false; flagSkipPC35 = false;
  if upperFreqBound<=0.1, flagSkipPC12 = true; end
  if lowerFreqBound>=0.1, flagSkipPC35 = true; end

  tint = TT.TimeInterval(ievent,:);
  tintDl = tint + 300*[-1 1]; %5 minutes are included in the timetable on either side of the EMIC wave
  if isempty(th_s) % Cluster
    gseR = local.c_read(['R' cl_s],tintDl);
    if isempty(gseR) || size(gseR,1)<=3
      gseR = c_caa_var_get(['sc_r_xyz_gse__C' cl_s '_CP_AUX_POSGSE_1M'],...
        'mat','tint',tintDl);
    end
    if size(gseR,1)<=3, continue, end
  else % THEMIS
    gseR = irf_tlim(Rth.(['Rth' lower(th_s)]),tintDl);
  end

  if flagSkipPC35, ebspPC35 = []; else, ebspPC35 = get_ebsp('PC35'); end
  if flagSkipPC12, ebspPC12 = []; else, ebspPC12 = get_ebsp('PC12'); end

  if isempty(ebspPC35) && isempty(ebspPC12)
    irf.log('warning','No data')
    continue
  end

  %% Combine PC12  &  PC35
  % Frequency limits
  if flagSkipPC35, freq = tocolumn(ebspPC12.f.data);
  elseif flagSkipPC12, freq = tocolumn(ebspPC35.f.data);
  else
    freq = [tocolumn(ebspPC35.f.data); tocolumn(ebspPC12.f.data)];
  end
  loweridx = find(abs(freq-lowerFreqBound)<.001);
  upperidx = find(abs(freq-upperFreqBound)<.001);
  idxFreq = loweridx:upperidx;

  if flagSkipPC35
    tint1min = [floor(tint(1)/60) ceil(tint(2)/60)]*60;
    outTime = ((tint1min(1)):60:(tint1min(end)))';
  else
    [outTime, idxPC35] = irf_tlim(ebspPC35.t.data,tint);
    BBssPC35 = double(ebspPC35.bb_xxyyzzss.data(idxPC35,:,4));
    ellPC35  = double(ebspPC35.ellipticity.data(idxPC35,:));
    ktPC35   = double(ebspPC35.k_tp.data(idxPC35,:,1));
  end
  if flagSkipPC12, BBss = BBssPC35; ell = ellPC35; thK = ktPC35;
  else
    [inTime, idxPC12] = irf_tlim(ebspPC12.t.data, tint + 30*[-1 1]); % 30 sec on each side for PC12 for propoper averaging
    if isempty(idxPC12), error('No PC12 data'), end
    BBssPC12 = double(ebspPC12.bb_xxyyzzss.data(idxPC12,:,4));
    BBssPC12av = AverageData(BBssPC12,inTime,outTime);
    ellPC12 = double(ebspPC12.ellipticity.data(idxPC12,:));
    ellPC12av = AverageData(ellPC12,inTime,outTime);
    ktPC12 = double(ebspPC12.k_tp.data(idxPC12,:,1));
    ktPC12av = AverageData(ktPC12,inTime,outTime);
    if flagSkipPC35, BBss = BBssPC12av; ell = ellPC12av; thK = ktPC12av;
    else % combined PC12 & PC35
      BBss = [BBssPC35, BBssPC12av];
      ell = [ellPC35, ellPC12av];
      thK = [ktPC35, ktPC12av];
    end
  end
  BBss=BBss(:,idxFreq);
  [psdBpeak,idxPeakTmp]=max(BBss,[],2);
  fPeak = double(freq(idxFreq(1)-1+idxPeakTmp));
  idxPeak = ((1:size(BBss,1))'-1)*size(BBss,2)+idxPeakTmp;
  ell = ell(:,idxFreq)'; % transpose to get the idx right
  ellPeak = ell(idxPeak);
  thK=thK(:,idxFreq)'; % transpose to get the idx right
  thKpeak=thK(idxPeak);

  % Mag coordinates
  gseR = irf_resamp(gseR,outTime);
  smR = irf.geocentric_coordinate_transformation(gseR,'gse>sm');
  [azimuth,elevation,~] = cart2sph(smR(:,2),smR(:,3),smR(:,4));
  mLat = elevation*180/pi;
  mlt = mod(azimuth*180/pi+180,360)/360*24;

  % Save output
  Out.time = [Out.time; outTime];
  Out.fPeak = [Out.fPeak; fPeak];
  Out.psdBpeak = [Out.psdBpeak; psdBpeak];
  Out.ellipticity = [Out.ellipticity; ellPeak];
  Out.thetaK = [Out.thetaK; thKpeak];
  Out.gseR = [Out.gseR; gseR];
  Out.mlt = [Out.mlt; mlt];
  Out.mLat = [Out.mLat; mLat];
end
cd (oldPwd)

  function chk_input()
    header = TT.Header{1};
    fprintf('Table : %s\n',header);
    if header(1) ~= 'C', error('TT Header must start with C'), end
    cl_s = header(2); th_s = '';
    if cl_s=='C' % THEMIS
      if header(4)~= 'T' || header(5)~= 'H'
        error('TT Header must start with CC_TH')
      end
      th_s = header(6);
      if ~any(th_s=='ABCDE')
        error('TT THEMIS header must start with CC_THA..E')
      end
    elseif ~any(cl_s=='1234'), error('TT Cluster header must start with C1..4')
    end
  end
  function res = get_ebsp(range)
    res = [];
    fname = find_file();
    if isempty(fname), return, end
    dobj=dataobj(fname);

    res = struct('t',[],'f',[],'flagFac',1,...
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

    iSep=strfind(fname,'__');
    productName = fname(1:iSep-1);

    for iField = 1:length(fields)
      fieldName = fields{iField}{1}; varName = fields{iField}{2};
      tmpVar = get_variable(dobj,[varName '__' productName]);
      if isempty(tmpVar)
        % Specially treat GB data which are not in FAC system
        if strcmp(fieldName,'bb_xxyyzzss')
          tmpVar = get_variable(dobj,['BB_xxyyzz__' productName]);
          if isempty(tmpVar), continue
          else, res.flagFac = 0;
          end
        else, continue
        end
      end
      if isnumeric(tmpVar.FILLVAL)
        tmpVar.data(tmpVar.data == tmpVar.FILLVAL) = NaN;
      end
      res.(fieldName) = struct('data',tmpVar.data,'units',tmpVar.UNITS);
      if strcmp(fieldName,'bb_xxyyzzss')
        if numel(size(res.(fieldName).data))==2 && all(all(isnan(res.(fieldName).data)))
          % For no data Qtran produces one record with fill values
          irf.log('warning',sprintf('No B data for %s',fname))
          return
        else
          res.(fieldName).data(:,:,4) = sum(res.(fieldName).data(:,:,1:3),3);
        end
      end
    end
    function fn = find_file()
      fn = '';
      if isempty(th_s) % Cluster
        indexDir = [ULF_PATH 'Cluster/C' cl_s '/' range '/CDF'];
      else % THEMIS
        indexDir = [ULF_PATH 'THEMIS/TH' th_s '/' range '/CDF'];
      end
      tStartStr = irf_fname(tintDl(1),3);
      tStopStr = irf_fname(tintDl(end),3);
      cd(indexDir);
      d = dir(['C' cl_s '*_' tStartStr '_*.cdf']);
      if ~strcmpi(tStopStr,tStartStr)
        d1 = dir(['C' cl_s '*_' tStopStr '_*.cdf']);
        d = [d; d1];
      end
      m=0; i=1;
      while i <= length(d) && m == 0
        if isempty(th_s), tend = d(i).name(45:59);
        else, tend = d(i).name(49:63);
        end
        tenddate = datenum(tend,'yyyymmdd_HHMMSS');
        if tenddate > epoch2date(tintDl(1)), m=i; end
        i=i+1;
      end
      if m, fn = d(m).name;
      else, irf.log('warning', ['no files found for ' range]);
      end
    end
  end % GET_EBSP
end

function out = AverageData(data,x,y,avWindow,flagSerial)
% average data with time x to time y using window
dtx = median(diff(x)); dty = median(diff(y));
if nargin<4, avWindow = dty; end
if nargin<5, flagSerial = 0; end
dt2 = avWindow/2;
ndataOut = length(y);

% Pad data with NaNs from each side
nPointToAdd = ceil(dt2/dtx);
padNan = zeros(nPointToAdd,size(data,2))*NaN;
data = [padNan; data; padNan];
padTime = dtx*(1:nPointToAdd);
x = [x(1)-fliplr(padTime)'; x; x(end)+padTime'];

out = zeros(ndataOut,size(data,2));
if flagSerial % Serial execution
  for i=1:length(y)
    out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
  end
else % Parallel execution
  parfor i=1:length(y)
    out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
  end
end
end
function m = FastNanMean(x,idx)
% Faster version of nanmean()
xx = x(idx,:);
% Find NaNs and set them to zero
nans = isnan(xx); xx(nans) = 0;
% Count up non-NaNs.
n = sum(~nans,1);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(xx,1) ./ n;
m(n<size(xx,1)*0.75) = NaN; % minDataFrac = .075
end
