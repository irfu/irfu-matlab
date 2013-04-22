function [h,varargout]=plot_f_f(n,m,t,vd,d,a1,a2,pitchangles,plotoption,title_option)
% Usage:
% [h,varargout]=whamp.plot_f_f(n,m,t,vd,d,a1,a2,[pitchangles],[plotoption]);
% [h,f,vp,vz]=whamp.plot_f_f(n,m,t,vd,d,a1,a2);
% [h,f,vtot]=whamp.plot_f_f(n,m,t,vd,d,a1,a2,[pitchangles],[plotoption],[title_option]);
%
% plot the distribution function, parameters as defined in whamp 
%
% Input:
%     n  - density [m^-3]
%     m  - particle mass, 0-electrons, 1-protons, 16-oxygen 
%     t  - temperature [keV]
%     vd - V_drift/V_term
%     d  - loss cone parameter, 1.0 = no loss cone)
%     a1 - T_perp/T_par
%     a2 - ??? (use 0)
%
% n,m,.. can be vectors (if more than one plasma component) 
% input can also be one matrix size Nx7 (if plotting countour plot)
% the plot is +-thermal velocities of first component around the mass center of first component
% 
% If pitchangles are given PSD vs velocity are plotted for the given angles.
% plotoption (optional) 0 - PSD vs V [m/s]
%			1 - PSD vs E [eV]
%			2 - F_reduced vs V_z [m/s] - sum_j m_1/m_j*F_j(V_z)
% title_option: 1 - default title (info on plasma parameters)
%               0 - no title
%               string - string as title 
%
% Output: 	h - handle to plot
%		f - phase space density [s^3/m^6] corresponding to:
%		vp, vz - matrices for plotting countour plots
%		vtot - nxm matrix where n-# of pitchangles, m-length of vtot [m/s].
%		Etot - nxm matrix where n-# of pitchangles, m-length of vtot [eV].
%
% Examples:
%	whamp.plot_f_f(4e6,1,0.3,0.9,1,1,0);
%	[h,f,vp,vz]=whamp.plot_f_f(4e6,1,0.3,0.9,1,1,0);
%	[h,f,vp,vz]=whamp.plot_f_f([4e6 1 0.3 0.9 1 1 0]);
%	[h,f,vtot]=whamp.plot_f_f(20e6,0,0.025,0,1,1,0,[0 45 90 135 180]);
%	[h,f,Etot]=whamp.plot_f_f(20e6,0,0.025,0,1,1,0,[0 45 90 135 180],1);
%
% short WHAMP manual: http://www.space.irfu.se/~andris/whamp/whamp_manual.pdf
% original WHAMP code: http://www.tp.umu.se/forskning/space/WHAMP/

Me=9.1094e-31; % electron mass
Mp=1.6726e-27; % proton mass
e=1.6022e-19; % elementary charge


if nargin < 1 , 
  help whamp.plot_f_f;
  return
end
if nargin == 1,
  if size(n,2)==7, 
	m=n(:,2);
	t=n(:,3);
	vd=n(:,4);
	d=n(:,5);
	a1=n(:,6);
	a2=n(:,7);
	n(:,2:end)=[];
  else
	help whamp.plot_f_f;
	return;
  end
end

if nargin<7 && nargin>1,return;end
for j=1:length(n),
  % estimate thermal velocity in m/s
  if m(j)==0, mm(j)=Me;else mm(j)=Mp*m(j);end
  vt(j)=sqrt(2*e * 1000*t(j)/mm(j)); % [m/s] (t in keV therefore *1000)
end

if nargin<8
  vpvec=-3*vt(1):vt(1)/20:3*vt(1);
  vzvec=min([-3*vt(1) -2*vt(1)+vd(1)]):vt(1)/20:max([vt(1)*3 vt(1)*3+vd(1)]);
  [vp,vz]=meshgrid(vpvec,vzvec);
%added pitchangles option /DS
elseif nargin>=8
	if nargin==8
		plotoption=0;
	end
	pitchangles=pitchangles*pi/180;
	for I=1:length(pitchangles)
		vvec(I,:)=0:vt(1)/20:10*vt(1);
		vpvec(I,:)=sin(pitchangles(I)).*vvec(I,:);
		vzvec(I,:)=cos(pitchangles(I)).*vvec(I,:);
    end
    vp=vpvec;
    vz=vzvec;
    for I=1:length(pitchangles)
        s = num2str(180/pi*pitchangles(I));
        M(I)={s};
    end
end


f=vp.*0;
vz_reduced=min([-3*vt(1) -2*vt(1)+vd(1)]):vt(1)/20:max([vt(1)*3 vt(1)*3+vd(1)]);
F_reduced=vz_reduced.*0; % reduced distribution function

ntot=sum(n)/length(n);	%added DS

for j=1:length(n),
  if a1(j)==0, 
    ea1=vp.*0;
  else
    ea1=exp(-1*(vp.^2./a1(j)./vt(j)./vt(j)));
  end
  if a2(j)==0, 
    ea2=vp.*0;
  else
    ea2=exp(-1*(vp.^2./a2(j)./vt(j)./vt(j)));
  end
  if a1(j)==a2(j) || d(j)==1,
    K=0;
  else
    K=(1-d(j))/(a1(j)-a2(j));
  end
  ff=exp(-1*(vz./vt(j)-vd(j)).^2);
  %f=f+ff.*n(j).*(d(j)/a1(j).*ea1+K.*(ea1-ea2));	%old
  f=f+1/(pi^(3/2)*(vt(j))^3)*ff.*n(j).*(d(j)/a1(j).*ea1+K.*(ea1-ea2));	%added normalization /DS
  										%f units [s^3/m^6]
										%f/ntot=f in whamp.
  F_reduced=F_reduced+(mm(1)/mm(j))/(pi^(1/2)*vt(j))*exp(-1*(vz_reduced./vt(j)-vd(j)).^2).*n(j);
end

QJAS_UNITS=1;
if QJAS_UNITS,  f = f*1e18; end

if nargin<8
	h=contour3(vp,vz,log10(f),30);view(0,90);
	axis equal;
	grid on;
	xlabel('Vperp');ylabel('Vpar');
% added PSD vs vtot for pitchangles /DS
elseif nargin>=8
    vtot=sqrt(vp.^2+vz.^2);
%	Etot=1e3*1/3.517388e14*1*vtot.^2;	%normalized to first species, Etot[eV]
    Etot=(1/e)*mm(1)/2*vtot.^2;	%normalized to first species, Etot[eV]
    if plotoption==1
        h=loglog(Etot',f');
        grid on
        xlabel('Etot [eV]')
        if QJAS_UNITS, ylabel('PSD [s^3/km^6]')
        else ylabel('PSD [s^3/m^6]')
        end
        legend(M)
    elseif plotoption==2
        % E_reduced=(1/e)*mm(1)/2*vz_reduced.^2;	%normalized to first species, Etot[eV]
        h=semilogy(vz_reduced,F_reduced);
        grid on
        xlabel('vz_{reduced} [m/s]')
        ylabel('F_{reduced} [[m/s]^{-1} m^{-3}]')
    else % default f(v)
        if QJAS_UNITS, vtot = vtot*1e-3; end
        h=plot(vtot',f');
        grid on
        if QJAS_UNITS
            xlabel('Vtot [km/s]')
            ylabel('PSD [s^3/km^6]')
        else
            xlabel('Vtot [m/s]')
            ylabel('PSD [s^3/m^6]')
        end
        legend(M)
    end
end
if exist('title_option','var'),
  if ischar(title_option),
    title(title_option);
  elseif title_option==1,
    title(['vt(1)=' num2str(vt(1)) ', vd(1)=' num2str(vd(1)) ', d(1)=' num2str(d(1)) ', a1(1)=' num2str(a1(1)) ', a2(1)=' num2str(a2(1)) ]);
  end
else
  title(['vt(1)=' num2str(vt(1)) ', vd(1)=' num2str(vd(1)) ', d(1)=' num2str(d(1)) ', a1(1)=' num2str(a1(1)) ', a2(1)=' num2str(a2(1)) ]);
end

%added varoutput /DS

if nargin<8
	varargout(1)={f};
	varargout(2)={vp};
	varargout(3)={vz};
elseif nargin==8
	varargout(1)={f};
	varargout(2)={vtot};
elseif nargin==9
	varargout(1)={f};
	varargout(2)={Etot};
end
	
