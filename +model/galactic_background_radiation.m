function [S,f]=galactic_background_radiation(varargin)
% MODEL.GALACTIC_BACKGROUND_RADIATION 
%	[S,f]=MODEL.GALACTIC_BACKGROUND_RADIATION() returns the galactic
%	radiation in (V^2/m^2/Hz) and frequency vector in Hz
%
%	[S]=MODEL.GALACTIC_BACKGROUND_RADIATION(f) returns the galactic
%	radiation at the given frequency f in (V^2/m^2/Hz)
%
%	[S]=MODEL.GALACTIC_BACKGROUND_RADIATION(f,'W') returns the galactic
%	radiation at the given frequency f in (W/m^2/Hz/sr)
%
% 	MODEL.GALACTIC_BACKGROUND_RADIATION('plot') plots the spectra in new
% 	window with units  (V^2/m^2/Hz)
%
% 	MODEL.GALACTIC_BACKGROUND_RADIATION('plot','W') plots the spectra in new
% 	window with units  (W/m^2/Hz/sr)
%
%	Example:
%		S=MODEL.GALACTIC_BACKGROUND_RADIATION(1e6); % radiation at 1MHz
%


%% based on Dulk et al., AAP, 2001. https://dx.doi.org/10.1051/0004-6361:20000006
%
% For the direction of the galactic poles, in particular the south galactic pole
% (SGP), the resulting spectrum of specific intensity in units of W m?2 Hz?1
% sr?1 is well described by the equation:
% 
% $$I_? = I_g ?^{?0.52} \frac{1 ? exp[??(?)]}{?(?)} + I_{eg} ?^{?0.80} exp[??(?)]$$ 
% 
% where ? is in MHz, the first and second terms are the galactic and
% extragalactic contributions respectively, and the opacity in the polar
% direction is ?(?). The values of the parameters are:
% 
% $$I_g = 2.48 10^{?20}, I_{eg} = 1.06 10^{?20}, ?(?) = 5.0 ?^{?2.1} $$
% 
% Cane (1977, 1978) also produced 10 and 30 MHz maps of apparent brightness
% temperature, Tb, for essentially the whole sky. We note that I? and Tb are
% related by the Rayleigh-Jeans Law:
% 
% $$kT_b	= I_? c^2/2?^2$$
% 
% Over large solid angles away from the galactic plane the brightness
% temperature is fairly constant and Eq. (1) is applicable. However, the
% brightness temperature rises by a factor of five to ten near the galactic
% plane and displays much more structure in angle. It would be rather difficult
% to accurately predict the response of most antennas to galactic emission near
% the plane, but fortunately the high brightness structures subtend a much
% smaller solid angle than the mid-to-high latitude regions; thus the total
% response of a low-gain antenna can be estimated with reasonable accuracy, i.e.
% to a factor of 2 or better.
% 
% 
% ### Converting to (V/m)^2/Hz
% 
% * Assume isotropic emission 
% * The galactic radiation energy density Eg=I_?* 4 pi / c 
% * The energy is half in E, half in B (e.m. emission) 
% 		Eg=2*eps0E^2/2=epso E^2
% * We have three electric field components thus
%		Ex^2=1/3 Eg/eps0 = (4 pi /3 ) I_? / c / eps0

%% Check input
isWUnits=1;		% default return in units (W/m^2/Hz/sr)
isPlot=0;		% default not to plot
fDefault=10.^(-1:.03:log10(30))*1e6; % default frequency vector

	% [S,f]=IRF_MODEL_GALACTIC_BACKGROUND_RADIATION() 
if nargin == 0 && nargout == 2 
	isWUnits=0; 
	f=fDefault;		% 0.1 - 30 MHz
	%
	% [S]=IRF_MODEL_GALACTIC_BACKGROUND_RADIATION(f) 
elseif nargin == 1 && isnumeric(varargin{:}) && nargout == 1 
	isWUnits=0;
	f=varargin{1};
	%
	% [S]=IRF_MODEL_GALACTIC_BACKGROUND_RADIATION(f,'W')
elseif nargin == 1 && isnumeric(varargin{:}) && nargout == 2 && strcmp(varargin{2},'W')
	f=varargin{1};
	%
	% IRF_MODEL_GALACTIC_BACKGROUND_RADIATION('plot')
elseif nargin == 1 && strcmp(varargin{1},'plot')
	f=fDefault;
	isWUnits=0;
	isPlot=1;
	%
	% IRF_MODEL_GALACTIC_BACKGROUND_RADIATION('plot','W')
elseif nargin == 2 && strcmp(varargin{1},'plot') && ...
		strcmp(varargin{2},'W')
	f=fDefault;
	isPlot=1;
else
	disp('WARNING!!!')
	disp('Check syntax of using function irf_model_galactic_background_radiation!');
	disp('Nothing returned');
	return;
end

%% Calculating emission

Ig	=2.48e-20;
Ieg =1.06e-20;
tau = @(f) 5*f.^(-2.1);
s_galactic = @(f) Ig*f.^(-0.52).*(1-exp(-tau(f)))./(tau(f));
s_extragalactic = @(f) Ieg*f.^(-0.80).*exp(-tau(f));
s_total = @(f) s_galactic(f)+s_extragalactic(f);
Units=irf_units;
s_in_volts_per_meter = @(f) s_total(f)*4/3*pi/Units.c/Units.eps0;

frequencyMHz=f./1e6;
if isWUnits
	S=s_total(frequencyMHz);
else
	S=s_in_volts_per_meter(frequencyMHz);
end

%% Plotting if required

if isPlot
	h=irf_plot(1,'newfigure');
    plot(h,f/1e6,S,'-');
	grid(h,'on');
	set(h,'xscale','log');
	set(h,'yscale','log');
	title(h,'Galactic background radiation');
	xlabel(h,'frequency [MHz]');
	if isWUnits
		ylabel(h,'S [W/m^2/Hz/s]');
	else
		ylabel(h,'S [V^2/m^2/Hz]');
	end
	irf_legend(h,'References: Dulk et al., A&A, 2001',[0.98,0.02]);
	clear S f; % do not reutrn anything
end
