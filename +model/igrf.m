function [out1,out2]=igrf(t,flag)
% MODEL.IGRF return IGRF model
%	[lambda,phi]=MODEL.IGRF(t,'dipole') return magnetic dipole latitude and longitude

persistent hIGRF gIGRF yearsIGRF

if isempty(hIGRF)
  fileIGRF = [fileparts(which('irf.m')), filesep, '+model', filesep, 'igrf14coeffs.txt'];
  irf.log('warning',['Reading IGRF coefficients from file:' fileIGRF]);
  %file reading
  fid = fopen(fileIGRF);
  out = textscan(fid, '%s', 'delimiter',sprintf('\n')); %#ok<SPRINTFN> % cell array with lines
  out = out{1};
  fclose(fid);
  % construct IGRF coefficient matrices
  % g(i,[n m years]) 0
  % h(i,[n,m,years]) 1
  nCoefIGRF = numel(out) - 4;
  a = textscan(out{4},'%s');
  a = a{1};
  nYears      = numel(a)-3;
  formatYears = repmat('%f',1,nYears);
  yearsIGRF   = textscan(out{4},['%*s%*s%*s' formatYears]);
  yearsIGRF   = horzcat(yearsIGRF{:});
  yearsIGRF(end)=yearsIGRF(end-1)+5;
  % read in all IGRF coefficients from file
  iIGRF=zeros(nCoefIGRF,nYears+3);
  ii=1;
  for iLine=5:numel(out)
    d=textscan(out{iLine},['%s%f%f' formatYears]);
    iIGRF(ii,2:end)=horzcat(d{2:end});
    if d{1}{1}-'h'==0, iIGRF(ii,1)=1;end
    ii=ii+1;
  end
  % the last column is derivative, make to value for +5 year
  iIGRF(:,end) = iIGRF(:,end-1)+5*iIGRF(:,end);
  hIGRF=iIGRF(iIGRF(:,1)==1,2:end);
  gIGRF=iIGRF(iIGRF(:,1)==0,2:end);
end

timeVec = irf_time(t,'vector');
yearRef = timeVec(:,1);
if min(yearRef) < min(yearsIGRF)
  irf.log('warning',...
    ['requested time is earlier than the first available IGRF model from ' ...
    num2str(min(yearRef)) ', extrapolating in past.. ']);
end

year = yearRef + ...
  (t - irf_time([yearRef repmat([1 1 0 0 0],numel(yearRef),1)],'vector>epoch'))...
  /(365.25*86400);

switch flag
  case 'dipole'
    g01 = interp1(yearsIGRF,gIGRF(1,3:end),year,'linear','extrap');
    g11 = interp1(yearsIGRF,gIGRF(2,3:end),year,'linear','extrap');
    h11 = interp1(yearsIGRF,hIGRF(1,3:end),year,'linear','extrap');
    lambda = atand(h11./g11);
    % There is a correction (Hapgood 1997) to the Hapgood 1992 model, which
    % replaces the previous arcsin with an arctan.
    phi = 90-atand((g11.*cosd(lambda)+h11.*sind(lambda))./g01); % latitude
    if nargout == 2
      out1 = lambda;
      out2 = phi;
    else
      disp(['lamda = ' num2str(lambda) ' deg']);
      disp(['  phi = ' num2str(phi) ' deg']);
    end
  otherwise
    irf.log('critical','input flag is not recognized');
    error('model.igrf() input flag not recognized');
end
