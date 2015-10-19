function [ampBtotal,elltotal,thetaKtotal]=wave_psd_ell_thetak(TT)
%MAARBLE.WAVE_PSD_ELL_THETAK  Extract PSD, ELL, thetaK for a TT
%
% [ampBtotal,elltotal,thetaKtotal]=wave_psd_ell_thetak(TT)

ampBtotal = []; elltotal = []; thetaKtotal = []; 

header12 = TT.Header{1}(1:2);
if header12(1) ~= 'C', error('TT HEader must start with C'), end
cl_s = header12(2);
if ~any(cl_s=='1234'), error('TT HEader must start with C1..4'), end

oldPwd = pwd;
for ievent=1:numel(TT),
    fprintf('Event #%d (%d) : %s',ievent,numel(TT),...
      irf_disp_iso_range(TT.TimeInterval(ievent,:),1));
    tint = TT.TimeInterval(ievent,:);
    tintDl = tint + 300*[-1 1]; %5 minutes are included in the timetable on either side of the EMIC wave
    gseR = c_caa_var_get(['sc_r_xyz_gse__C' cl_s '_CP_AUX_POSGSE_1M'],...
        'mat','tint',tintDl);
    
    if size(gseR,1)<=3, continue, end
    indexDir = ['/data/caa/MAARBLE/WaveDatabase/ULF/Cluster/C' cl_s '/PC12/CDF'];
    cd(indexDir);
    d = dir(['C' cl_s '*.cdf']);
    m=0;
    i=1;
    while i <= length(d) && m == 0,
      tend = d(i).name(45:59);
      tenddate = datenum(tend,'yyyymmdd_HHMMSS');
      if tenddate > epoch2date(tintDl(1)),
        m=i;
      end
      i=i+1;
    end
    fname=d(m).name;
    dobj=dataobj(fname);
    indexDir2 = ['/data/caa/MAARBLE/WaveDatabase/ULF/Cluster/C' cl_s '/PC35/CDF'];
    cd(indexDir2);
    iSep=strfind(fname,'__');
    productName = fname(1:iSep-1);
    PCsep=strfind(fname,'PC12');
    fnamePC35=[fname(1:PCsep-1) 'PC35' fname(PCsep+4:end)];
    dobjPC35=dataobj(fnamePC35);
    productNamePC35 = fnamePC35(1:iSep-1);
    
    ebspPC12 = struct('t',[],'f',[],'flagFac',1,...
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
      fieldName = fields{iField}{1}; varName = fields{iField}{2};
      tmpVar = get_variable(dobj,[varName '__' productName]);
      tmpVarPC35 = get_variable(dobjPC35,[varName '__' productNamePC35]);
      if isempty(tmpVar)
        % Specially treat GB data which are not in FAC system
        if strcmp(fieldName,'bb_xxyyzzss')
          tmpVar = get_variable(dobj,['BB_xxyyzz__' productName]);
          if isempty(tmpVar), continue
          else ebspPC12.flagFac = 0;
          end
        else continue
        end
      end
      if isempty(tmpVarPC35)
        % Specially treat GB data which are not in FAC system
        if strcmp(fieldName,'bb_xxyyzzss')
          tmpVarPC35 = get_variable(dobjPC35,['BB_xxyyzz__' productName]);
          if isempty(tmpVarPC35), continue
          else ebspPC35.flagFac = 0;
          end
        else continue
        end
      end
      if isnumeric(tmpVar.FILLVAL)
        tmpVar.data(tmpVar.data == tmpVar.FILLVAL) = NaN;
      end
      if isnumeric(tmpVarPC35.FILLVAL)
        tmpVarPC35.data(tmpVarPC35.data == tmpVarPC35.FILLVAL) = NaN;
      end
      ebspPC12.(fieldName) = struct('data',tmpVar.data,'units',tmpVar.UNITS);
      ebspPC35.(fieldName) = struct('data',tmpVarPC35.data,'units',tmpVarPC35.UNITS);
      if strcmp(fieldName,'bb_xxyyzzss')
        if numel(size(ebspPC12.(fieldName).data))==2 && all(all(isnan(ebspPC12.(fieldName).data)))
          % For no data Qtran produces one record with fill values
          irf.log('warning',sprintf('No B data for %s',fname))
          return
        else
          ebspPC12.(fieldName).data(:,:,4) = sum(ebspPC12.(fieldName).data(:,:,1:3),3);
        end
      end
      if strcmp(fieldName,'bb_xxyyzzss')
        if numel(size(ebspPC35.(fieldName).data))==2 && all(all(isnan(ebspPC35.(fieldName).data)))
          % For no data Qtran produces one record with fill values
          irf.log('warning',sprintf('No B data for %s',fname))
          return
        else
          ebspPC35.(fieldName).data(:,:,4) = sum(ebspPC35.(fieldName).data(:,:,1:3),3);
        end
      end
    end
    
    %% Construct ebsp
    
    [outTime, idxPC35] = irf_tlim(ebspPC35.t.data,tint);
    BBssPC35 = double(ebspPC35.bb_xxyyzzss.data(idxPC35,:,4));
    
    [inTime, idxPC12] = irf_tlim(ebspPC12.t.data, tint + 30*[-1 1]); % 30 sec on each side for PC12 for propoper averaging
    BBssPC12 = double(ebspPC12.bb_xxyyzzss.data(idxPC12,:,4));
    BBssPC12av = AverageData(BBssPC12,inTime,outTime);
    BBss = [BBssPC35, BBssPC12av]; % combined BBss for PC12 & 35
    
    f = [tocolumn(ebspPC35.f.data); tocolumn(ebspPC12.f.data)];
    lowerFreqBound = str2double(TT.Description{ievent}{1});
    upperFreqBound = str2double(TT.Description{ievent}{2});
    loweridx = find(abs(f-lowerFreqBound)<.001);
    upperidx = find(abs(f-upperFreqBound)<.001);
    idxFreq = loweridx:upperidx;
    
    BBss=BBss(:,idxFreq);
    [ampBpeak,idxPeakTmp]=max(BBss,[],2);
    idxPeak = ((1:size(BBss,1))'-1)*size(BBss,2)+idxPeakTmp;
    
    ellPC35 = double(ebspPC35.ellipticity.data(idxPC35,:));
    ellPC12 = double(ebspPC12.ellipticity.data(idxPC12,:));
    ellPC12 = AverageData(ellPC12,inTime,outTime);
    ell = [ellPC35, ellPC12]; % combined ellipticity for PC12 & 35
    ell = ell(:,idxFreq);
    ellpeak = ell(idxPeak);
    
    ktPC12=double(ebspPC12.k_tp.data(idxPC12,:,1));
    ktPC35=double(ebspPC35.k_tp.data(idxPC35,:,1));
    ktPC12 = AverageData(ktPC12,inTime,outTime);
    thetaK = [ktPC35, ktPC12];
    thetaK=thetaK(:,idxFreq);
    thetaKpeak=thetaK(idxPeak);
      
    ampBtotal = [ampBtotal;ampBpeak]; %#ok<AGROW>
    elltotal = [elltotal;ellpeak]; %#ok<AGROW>
    thetaKtotal = [thetaKtotal;thetaKpeak]; %#ok<AGROW>
end

cd (oldPwd)

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
