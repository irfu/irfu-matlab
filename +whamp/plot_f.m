function [varargout]=plot_f(varargin)
%function [h,varargout]=plot_f(n,m,t,vd,d,a1,a2,pitchangles,plotoption,title_option)
% WHAMP.PLOT_F plotting distribution functions used in WHAMP calculations
%
% WHAMP.PLOT_F(n,m,t,vd,d,a1,a2) Plot distribution functions
%     in the interval +-3 thermal velocities of the first component. 
%     n,m,.. can be vectors (if more than one plasma component) or input
%     can also be one matrix size Nx7.
%
% Input:
%     n  - density [m^-3]
%     m  - particle mass, 0-electrons, 1-protons, 16-oxygen 
%     t  - temperature [keV]
%     vd - V_drift/V_term
%     d  - loss cone parameter (1 when no loss cone)
%     a1 - T_perp/T_par
%     a2 - loss cone parameter (0 when no loss cone)
%
% WHAMP.PLOT_F(Plasmamodel) Plasma model can be specified as Plasmamodel
%     structure or Species structure, see for defintion WHAMP.RUN.
%
%           WHAMP.PLOT_F(H,..) Plot in axis with handle H   
%         [h]=WHAMP.PLOT_F(..) Returns handle to the plot axis   
% [h,f,vp,vz]=WHAMP.PLOT_F(..) Returns perp and par velocities 
%  [h,f,vtot]=WHAMP.PLOT_F(..) Returns total velocity or energy depending
%                              on 'PSDvs..' flag
%
% WHAMP.PLOT_F(..,'pitchangles',pitchangles) If pitchangles are given PSD
%     vs velocity are plotted for the given angles.
%
% WHAMP.PLOT_F(..,'PSDvsV') plot PSD vs velocity in [m/s] (default)
% WHAMP.PLOT_F(..,'PSDvsE') plot PSD vs energy in [eV]
% WHAMP.PLOT_F(..,'FredvsVz') plot reduce distribution function in sum_j m_1/m_j*F_j(V_z)
% WHAMP.PLOT_F(..,'notitle') do not print title (default is print plasma
%      parameters)
% WHAMP.PLOT_F(..,'km/s') use km/s instead of default m/s
%
% Output: 	h - handle to plot
%		f - phase space density [s^3/m^6] corresponding to:
%		vp, vz - matrices for plotting countour plots
%		vtot - nxm matrix where n-# of pitchangles, m-length of vtot [m/s].
%		Etot - nxm matrix where n-# of pitchangles, m-length of vtot [eV].
%
% Examples:
%	whamp.plot_f(4e6,1,0.3,0.9,1,1,0);
%	whamp.plot_f(4e6,1,0.3,0.9,1,1,0,'pitchangles',[0 45 90]);
%  Oxygen = struct('m',16,'n',1,'t',10,'a',5,'vd',1,'d',1,'b',0);
% [h,f,vp,vz] = whamp.plot_f(Oxygen,'km/s');
%
% WHAMP quick reference: https://github.com/irfu/whamp/blob/master/whamp_manual.pdf?raw=true
% original WHAMP code: http://www.tp.umu.se/forskning/space/WHAMP/

Units = irf_units;
Me = Units.me; % electron mass
Mp = Units.mp; % proton mass
e  = Units.e; % elementary charge

%% Flags that can be overwritten by input parameters
plotPSDvsV      = false;
plotPSDvsE      = false;
plotFredvsVz    = false;
printTitle      = true;
plotPitchangles = false;
unitsVelocity   = 1; % default use SI unit [m/s]
labelVelocity   = 'm/s';

%% Check input 
[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax)
	figure;
	h=gca;
else
	h=ax;
end

x=args{1};
if isempty(x) % nothing to plot, first input parameter empty
  return;
end

isPlasmaModelDefined = false; % default
while ~isempty(args)
	if 	isstruct(args{1}) && isfield(args{1},'Species') % whamp.run_f(PlasmaModel)
		plasmaModel = args{1};
		for iSpecies = numel(plasmaModel.Species):-1:1
			species = plasmaModel.Species{iSpecies};
			n(iSpecies)  = species.n;
			m(iSpecies)  = species.m;
			t(iSpecies)  = species.t;
			vd(iSpecies) = species.vd;
			d(iSpecies)  = species.d;
			a(iSpecies) = species.a;
			b(iSpecies) = species.b;
		end
		isPlasmaModelDefined = true;
		args(1)=[];
	elseif isstruct(args{1}) && isfield(args{1},'m')  % whamp.run_f(Species)
		species = args{1};
		n  = species.n;
		m  = species.m;
		t  = species.t;
		vd = species.vd;
		d  = species.d;
		a = species.a;
		b = species.b;
		isPlasmaModelDefined = true;
		args(1)=[];
	elseif nargs < 1 && isnumeric(args{1}) &&	size(args{1},2)==7   % whamp.run_f([i x 7]) matrix with all values
		species = args{1};
		n =species(:,1);
		m =species(:,2);
		t =species(:,3);
		vd=species(:,4);
		d =species(:,5);
		a=species(:,6);
		b=species(:,7);
		isPlasmaModelDefined = true;
		args(1)=[];
	elseif nargs >= 7 && isnumeric(args{1})   % whamp.run_f(n,m,t,vd,d,a1,a2,...) matrix with all values
		n =args{1};
		m =args{2};
		t =args{3};
		vd=args{4};
		d =args{5};
		a=args{6};
		b=args{7};
		isPlasmaModelDefined = true;
		args(1:7)=[];
	elseif ischar(args{1}) && strcmpi(args{1},'notitle')
		printTitle = false;
		args(1) = [];
	elseif ischar(args{1}) && strcmpi(args{1},'PSDvsV')
		plotPSDvsV = true;
		args(1) = [];
	elseif ischar(args{1}) && strcmpi(args{1},'PSDvsE')
		plotPSDvsE = true;
		args(1) = [];
	elseif ischar(args{1}) && strcmpi(args{1},'FredvsVz')
		plotFredvsVz = true;
		args(1) = [];
	elseif ischar(args{1}) && strcmpi(args{1},'km/s')
		unitsVelocity = 1e3;
		labelVelocity = 'km/s';
		args(1) = [];
	elseif ischar(args{1}) && strcmpi(args{1},'pitchangles')
		plotPitchangles = true;
		args(1) = [];
		if numel(args) > 0 && isnumeric(args{1})
			pitchangles = args{1};
			args(1) = [];
		else
			irf.log('critical','ERROR: pitchangles not specified, see help!');
			return;
		end
	else
		irf.log('critical','ERROR: WHAMP input unrecognized, see help!');
		help whamp.plot_f;
		return;
	end
	nargs = numel(args);
end

if ~isPlasmaModelDefined
	irf.log('critical','ERROR: WHAMP plasma model not specified!');
  help whamp.plot_f;
	return
end

%% Check which plot to make
if plotPSDvsV
	plotPSDvsE      = false;
	plotFredvsVz    = false;
elseif plotPSDvsE
	plotFredvsVz    = false;
elseif ~plotFredvsVz
	plotPSDvsV = true;
end
%% Calculate distribution function

for j=length(n):-1:1
  % estimate thermal velocity in m/s
  if m(j)==0, mm(j)=Me;else, mm(j)=Mp*m(j);end
  vt(j)=sqrt(2*e * 1000*t(j)/mm(j)); % [m/s] (t in keV therefore *1000)
end


if plotPitchangles
	pitchangles=pitchangles*pi/180;
	for I=length(pitchangles):-1:1
		vvec(I,:)=0:vt(1)/20:10*vt(1);
		vpvec(I,:)=sin(pitchangles(I)).*vvec(I,:);
		vzvec(I,:)=cos(pitchangles(I)).*vvec(I,:);
	end
	vp=vpvec;
	vz=vzvec;
	for I=length(pitchangles):-1:1
		s = num2str(180/pi*pitchangles(I));
		M(I)={s};
	end
else
	vpvec=-3*vt(1):vt(1)/20:3*vt(1);
  vzvec=min([-3*vt(1) -2*vt(1)+vd(1)]):vt(1)/20:max([vt(1)*3 vt(1)*3+vd(1)]);
  [vp,vz]=meshgrid(vpvec,vzvec);
end


f=vp.*0;
vz_reduced=min([-3*vt(1) -3*vt(1)+vd(1)]):vt(1)/20:max([vt(1)*3 vt(1)*3+vd(1)]);
F_reduced=vz_reduced.*0; % reduced distribution function

for j=1:length(n)
  if a(j)==0 
    ea1=vp.*0;
  else
    ea1=exp(-1*(vp.^2./a(j)./vt(j)./vt(j)));
  end
  if b(j)==0 
    ea2=vp.*0;
  else
    ea2=exp(-1*(vp.^2./b(j)./vt(j)./vt(j)));
  end
  if a(j)==b(j) || d(j)==1
    K=0;
  else
    K=(1-d(j))/(a(j)-b(j));
  end
  ff=exp(-1*(vz./vt(j)-vd(j)).^2);
  %f=f+ff.*n(j).*(d(j)/a1(j).*ea1+K.*(ea1-ea2));	%old
  f=f+1/(pi^(3/2)*(vt(j))^3)*ff.*n(j).*(d(j)/a(j).*ea1+K.*(ea1-ea2));	%added normalization /DS
  										%f units [s^3/m^6]
										%f/ntot=f in whamp.
  F_reduced=F_reduced+(mm(1)/mm(j))/(pi^(1/2)*vt(j))*exp(-1*(vz_reduced./vt(j)-vd(j)).^2).*n(j);
end

f = f*unitsVelocity^6;

if plotFredvsVz
        % E_reduced=(1/e)*mm(1)/2*vz_reduced.^2;	%normalized to first species, Etot[eV]
        semilogy(h,vz_reduced/unitsVelocity,F_reduced);
        grid(h,'on');
        xlabel(h,['vz_{reduced} [' labelVelocity ']'])
        ylabel(h,['F_{reduced} [[' labelVelocity ']^{-1} m^{-3}]'])
elseif plotPitchangles
    vtot=sqrt(vp.^2+vz.^2);
    Etot=(1/e)*mm(1)/2*vtot.^2;	%normalized to first species, Etot[eV]
    if plotPSDvsE
        loglog(h,Etot',f');
        grid(h,'on');
        xlabel(h,'Etot [eV]')
        ylabel(h,['PSD [s^3/' labelVelocity	'^6]'])
        legend(h,M)
		elseif plotPSDvsV
        plot(h,vtot'/unitsVelocity,f');
        grid(h,'on');
        xlabel(h,['Vtot [' labelVelocity	']'])
				ylabel(h,['PSD [s^3/' labelVelocity '^6]'])
        legend(h,M)
    end
else
	contour3(h,vp/unitsVelocity,vz/unitsVelocity,f,30);view(0,90);
	axis(h,'equal');
	grid(h,'on');
	xlabel(h,['V_{\perp} [' labelVelocity ']']);
	ylabel(h,['V_{||} [' labelVelocity ']']);
end
if printTitle
	title(h,['1st comp: T=' num2str(t(1),3) 'keV'...
		', vd/vt=' num2str(vd(1))...
		', Tperp/Tpar=' num2str(a(1)) ', d=' num2str(d(1)) ', b=' num2str(b(1)) ]);
end

%% Output definition

if nargout >= 1 
	varargout{1} = h;
end
if nargout == 2
	if plotPSDvsE
		varargout{2}=f;
		varargout{3}=vtot/unitsVelocity;
	elseif plotPSDvsE
		varargout{2}=f;
		varargout{3}=Etot;
	elseif plotFredvsVz
		varargout{2}=F_reduced;
		varargout{3}=vz_reduced/unitsVelocity;
	end
elseif nargout ==3
	varargout{2}=f;
	varargout{3}=vp/unitsVelocity;
	varargout{4}=vz/unitsVelocity;
end
	
