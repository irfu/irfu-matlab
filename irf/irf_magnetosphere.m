function [x,y,omni] = irf_magnetosphere(model,Dp,Bz)
%IRF_MAGNETOSPHERE Return the location of magnetopause
%
%  [X,Y]=IRF_MAGNETOSPHERE(model,Dp,Bz)
%     X,Y - vectors with X and Y coordinates of magnetopause location
%
% INPUT:
%       model - model to use. Implemented - 'mp_shue1998','bs'
%       Dp    - dynamic pressure
%       Bz    - IMF Bz GSM
%
%  [X,Y]=IRF_MAGNETOSPHERE(model,time)
%       get solar wind parameters online form OMNI  database
%       if no OMNI data available, return empty dataset
%
%  [X,Y]=IRF_MAGNETOSPHERE(model)
%       use default parameters, Dp=2nPa, Bz=0nT
%
%  [X,Y,OMNI]=IRF_MAGNETOSPHERE(model,time)
%       return OMNI data used in structure OMNI (omni.Dp and omni.Bx, By, Bz in GSM)
%
% Examples:
%  [x,y] = irf_magnetosphere('mp_shue1998',10,-2)
%  [x,y] = irf_magnetosphere('mp_shue1998',irf_time([2001 10 01 18 0 0]))
%

persistent dpbz

if nargout>0 % default return empty variables
  x=[];y=[];omni=[];
end

if nargin == 1 % use default solar wind values
  Dp=2;
  Bx=0;By=0;Bz=0;
  M=4;
elseif nargin == 2 % IRF_MAGNETOPAUSE(model, time)
  t=Dp;
  tint=t + [-2 2]*3600;
  if isempty(dpbz) || t<dpbz(1,1) || t>dpbz(end,1)
    dpbz=irf_get_data(tint,'P,bzgsm,bx,bygsm,Ma','omni_min');
  end
  if isempty(dpbz) % no OMNI data, return empty
    return
  else
    dpbz_t=irf_resamp(dpbz,t);
    Dp = dpbz_t(2);
    if isnan(Dp)
      dp_tmp = irf_get_data(tint, 'P','omni');
      dp_tmp = irf_resamp(dp_tmp, t);
      Dp = dp_tmp(2);
      irf_log('fcal', 'Dp uses 1-hour omni data. ');
    end
    Bz = dpbz_t(3);
    if isnan(Bz)
      bz_tmp = irf_get_data(tint, 'bzgsm', 'omni');
      bz_tmp = irf_resamp(bz_tmp, t);
      Bz = bz_tmp(2);
      irf_log('fcal', 'Bz uses 1-hour omni data. ');
    end
    Bx=dpbz_t(4);
    if isnan(Bx)
      bx_tmp = irf_get_data(tint, 'bx', 'omni');
      bx_tmp = irf_resamp(bx_tmp, t);
      Bx = bx_tmp(2);
      irf_log('fcal', 'Bx uses 1-hour omni data. ');
    end
    By=dpbz_t(5);
    if isnan(By)
      by_tmp = irf_get_data(tint, 'bygsm', 'omni');
      by_tmp = irf_resamp(by_tmp, t);
      By = by_tmp(2);
      irf_log('fcal', 'By uses 1-hour omni data. ');
    end
    M=dpbz_t(6);
    if isnan(M)
      m_tmp = irf_get_data(tint, 'Ma', 'omni');
      m_tmp = irf_resamp(m_tmp, t);
      M = m_tmp(2);
      irf_log('fcal', 'Ma uses 1-hour omni data. ');
    end
    if isnan(Dp) || isnan(Bz) || isempty(Dp) || isempty(Bz)
      return
    end
  end
elseif nargin==3 % specified Dp and Bz
  Bx=0;By=0;M=4;
end
omni.Dp=Dp;omni.Bz=Bz;omni.Bx=Bx;omni.By=By;

switch lower(model)
  case 'mp_shue1998'
    % Reference: Shue et al 1998
    %  Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
    %  Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
    % Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
    % Default values: Dp=2nPa, Bz=0nT
    %
    theta=0:0.1:pi;
    rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6);
    alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
    r=rzero*(2./(1+cos(theta))).^alpha;
    x=r.*cos(theta);
    y=r.*sin(theta);
    ii=find(abs(x)>100);
    x(ii)=[];
    y(ii)=[];
  case 'bs'
    % 'bs'
    %  standoff distance (Farris and Russell 1994)
    %  rstandoff=rmp*(1+1.1*((gamma-1)*M^2+2)/((gamma+1)*(M^2-1)))
    [xmp,~] = irf_magnetosphere('mp_shue1998',Dp,Bz);
    gamma=5/3;
    rmp=xmp(1);
    rstandoff=rmp*(1+1.1*((gamma-1)*M^2+2)/((gamma+1)*(M^2-1)));
    %x=rstandoff:-0.5:-100;
    x=rstandoff-logspace(log10(0.1),log10(100+rstandoff),300)+0.1; % Smaller increments at the subsolar point
    rho=sqrt(0.04*(x-rstandoff).^2-45.3*(x-rstandoff)); % original F/G model adds rstandoff^2=645
    y=rho;
  otherwise
    irf_log('fcal','Unknown model.');
    return;
end

