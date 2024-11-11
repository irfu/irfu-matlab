function [varargout] = solar_wind(model,varargin)
%SOLAR_WIND return solar wind parameters
%
%  [..]=MODEL.SOLAR_WIND(model,{input_parameters})
%
% INPUT:
%       model - model to use. Implemented
%             'stverak2015_electrons_slow_wind' (doi:10.1002/2015JA021368)
%             'stverak2015_electrons_fast_wind'
%
%  [n_m3,Tpar_K,Tperp_K,T_K,qpar_Wm2,qperp_Wm2,q_Wm2]=...
%       MODEL.SOLAR_WIND('stverak2015_electrons_slow_wind',distanceToSunInAU)
%%
% Examples:
%  [n,Tpar,Tperp,T,qpar,qperp,q]=model.solar_wind('stverak2015_electrons_slow_wind',1);

toPrintOutResult = false;
if nargout==0 % default return empty variables, printout values
  toPrintOutResult = true;
end

if nargin == 0
  help model.solar_wind;
  return;
elseif nargin == 1 % use default solar wind values
  distanceToSunInAU = 1;
elseif nargin == 2 % IRF_MAGNETOPAUSE(model, time)
  distanceToSunInAU = varargin{1};
end

switch lower(model)
  case 'stverak2015_electrons_slow_wind'
    % Reference: doi:10.1002/2015JA021368 Table 1
    % Based on Helios 1 and 2

    n_m3    = 9.32e6 * distanceToSunInAU.^(-2.03);
    Tpar_K  = 1.68e5 * distanceToSunInAU.^(-0.67);
    Tperp_K = 1.37e5 * distanceToSunInAU.^(-0.53);
    T_K     = 1.48e5 * distanceToSunInAU.^(-0.59);
    Qpar_Wm2  = 3.23e-5 * distanceToSunInAU.^(-2.95);
    Qperp_Wm2 = 7.05e-6 * distanceToSunInAU.^(-2.40);
    Q_Wm2     = 1.52e-5 * distanceToSunInAU.^(-2.84);
    varoutNames = {'n_m3','Tpar_K','Tperp_K','Qpar_Wm2','Qperp_Wm2','Q_Wm2'};
    varargout={n_m3,Tpar_K,Tperp_K,T_K,Qpar_Wm2,Qperp_Wm2,Q_Wm2};
  case 'stverak2015_electrons_fast_wind'
    % Reference: doi:10.1002/2015JA021368 Table 1

    n_m3    = 5.61e6 * distanceToSunInAU.^(-1.83);
    Tpar_K  = 1.74e5 * distanceToSunInAU.^(-0.35);
    Tperp_K = 9.92e4 * distanceToSunInAU.^(-0.27);
    T_K     = 1.24e5 * distanceToSunInAU.^(-0.31);
    Qpar_Wm2  = 2.86e-5 * distanceToSunInAU.^(-2.47);
    Qperp_Wm2 = 2.97e-6 * distanceToSunInAU.^(-2.28);
    Q_Wm2     = 1.15e-5 * distanceToSunInAU.^(-2.44);
    varoutNames = {'n_m3','Tpar_K','Tperp_K','Qpar_Wm2','Qperp_Wm2','Q_Wm2'};
    varargout={n_m3,Tpar_K,Tperp_K,T_K,Qpar_Wm2,Qperp_Wm2,Q_Wm2};
  otherwise
    irf_log('error','Unknown model.');
    return;
end

if toPrintOutResult
  for iVar = 1:numel(varoutNames)
    disp([varoutNames{iVar} ': ' num2str(varargout{iVar}(1),'%10.2e\n') ]);
  end
end

if nargout == 0
  clear varargout;
end
