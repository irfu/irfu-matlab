function [x,y] = irf_magnetosphere(model,Dp,Bz)
%IRF_MAGNETOPAUSE Return the location of magnetopause
%  
%  IRF_MAGNETOPAUSE(model,Dp,Bz)
%
% INPUT: 
%       model - model to use. Implemented - "mp_shue1998"
%       Dp    - dynamic pressure
%       Bz    - IMF Bz GSM
%
%  [X,Y] = IRF_MAGNETOPAUSE('mp_shue1998',Dp,Bz)
%     X,Y - vectors with X and Y coordinates of magnetopause location
%

% Reference: Shue et al 1998
%  Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
%  Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
% Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
% Default values: Dp=2nPa, Bz=0nT
% 

% 'bs'
%  standoff distance (Farris and Russell 1994)
%  rstandoff=rmp*(1+1.1*((gamma-1)*M^2+2)/((gamma+1)*(M^2-1)))

if nargin == 1, % use default solar wind values
    Dp=2;
    Bz=0;
end

switch lower(model)
    case 'mp_shue1998'
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
        [xmp,~] = irf_magnetosphere('mp_shue1998',Dp,Bz);
        gamma=5/3;
        M=4;
        rmp=xmp(1);
        rstandoff=rmp*(1+1.1*((gamma-1)*M^2+2)/((gamma+1)*(M^2-1)));
        x=rstandoff:-0.5:-100;
        rho=sqrt(0.04*(x-rstandoff).^2-45.3*(x-rstandoff)); % original F/G model adds rstandoff^2=645
        y=rho;
    otherwise
        irf_log('fcal','Unknown model.');
        return;
end

