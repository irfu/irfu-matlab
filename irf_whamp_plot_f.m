function [h,varargout]=irf_whamp_plot_f(n,m,t,vd,d,a1,a2,pitchangles,plotoption);
% Usage:
% [h,varargout]=irf_whamp_plot_f(n,m,t,vd,d,a1,a2,[pitchangles],[plotoption]);
% [h,f,vp,vz]=irf_whamp_plot_f(n,m,t,vd,d,a1,a2);
% [h,f,vtot]=irf_whamp_plot_f(n,m,t,vd,d,a1,a2,[pitchangles],[plotoption]);
%
% plot the distribution function, parameters as defined in whamp 
% n,m,.. can be vectors (if more than one plasma component) 
% input can also be one matrix size Nx7 (if plotting countour plot)
% the plot is +-thermal velocities of first component around the mass center of first component
% 
% if pitchangles are given PSD vs velocity are plotted for the given angles.
% plotoption (optional) 0 - PSD vs V [m/s]
%			1 - PSD vs E [eV]
%
% Output: 	h - handle to plot
%		f - phase space density [s^3/m^6] corresponding to:
%		vp, vz - matrices for plotting countour plots
%		vtot - nxm matrix where n-# of pitchangles, m-length of vtot [m/s].
%		Etot - nxm matrix where n-# of pitchangles, m-length of vtot [eV].
% Examples:
%	irf_whamp_plot_f(4e6,1,0.3,0.9,1,1,0);
%	[h,f,vp,vz]=irf_whamp_plot_f(4e6,1,0.3,0.9,1,1,0);
%	[h,f,vp,vz]=irf_whamp_plot_f([4e6 1 0.3 0.9 1 1 0]);
%	[h,f,vtot]=irf_whamp_plot_f(20e6,0,0.025,0,1,1,0,[0 45 90 135 180]);
%	[h,f,Etot]=irf_whamp_plot_f(20e6,0,0.025,0,1,1,0,[0 45 90 135 180],1);
%
% Last change 20050415
%
%
% short WHAMP manual: http://www.space.irfu.se/~andris/whamp/whamp_manual_v1.1.pdf
% original WHAMP code: http://www.tp.umu.se/forskning/space/WHAMP/


if nargin < 1 , 
  help av_whamp_plot_f;
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
  else,
	help av_whamp_plot_f;
	return;
  end
end

if nargin<7 & nargin>1,return;end
for j=1:length(n),
  % estimate thermal velocity
  if m(j)==0, mm=1;else mm=1836.1*m(j);end
  vt(j)=sqrt(3.517388e14 * t(j)/mm);
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
		vpvec(I,:)=0:vt(1)/20:3*vt(1);
		vzvec(I,:)=cos(pitchangles(I)).*vpvec(I,:);
	end
	vp=vpvec;
	vz=vzvec;
end


f=vp.*0;

ntot=sum(n)/length(n);	%added DS

for j=1:length(n),
  if a1(j)==0, 
    ea1=vp.*0;;
  else
    ea1=exp(-1*(vp.^2./a1(j)./vt(j)./vt(j)));
  end
  if a2(j)==0, 
    ea2=vp.*0;;
  else
    ea2=exp(-1*(vp.^2./a2(j)./vt(j)./vt(j)));
  end
  if a1(j)==a2(j),
    K=0;
  else
    K=(1-d(j))/(a1(j)-a2(j));
  end
  ff=exp(-1*(vz./vt(j)-vd(j)).^2);
  %f=f+ff.*n(j).*(d(j)/a1(j).*ea1+K.*(ea1-ea2));	%old
  f=f+1/(pi^(3/2)*(vt(j))^3)*ff.*n(j).*(d(j)/a1(j).*ea1+K.*(ea1-ea2));	%added normalization /DS
  										%f units [s^3/m^6]
										%f/ntot=f in whamp.
end

if nargin<8
	h=contour3(vp,vz,log10(f),30);view(0,90);
	axis equal;
	grid on;
	title(['vt(1)=' num2str(vt(1)) ', vd(1)=' num2str(vd(1)) ', d(1)=' num2str(d(1)) ', a1(1)=' num2str(a1(1)) ', a2(1)=' num2str(a2(1)) ]);
	xlabel('Vperp');ylabel('Vpar');
% added PSD vs vtot for pitchangles /DS
elseif nargin>=8
	vtot=sqrt(vp.^2+vz.^2);
	Etot=1e3*1/3.517388e14*1*vtot.^2;	%normalized to first species, Etot[eV]
	if plotoption==1
		h=loglog(Etot',f');
		grid on
		xlabel('Etot [eV]')
	else
		h=plot(vtot',f');
		xlabel('Vtot [m/s]')
	end
	title(['vt(1)=' num2str(vt(1)) ', vd(1)=' num2str(vd(1)) ', d(1)=' num2str(d(1)) ', a1(1)=' num2str(a1(1)) ', a2(1)=' num2str(a2(1)) ]);
	ylabel('PSD [s^3/m^6]')
	for I=1:length(pitchangles)
		s=[num2str(180/pi*pitchangles(I))];
		M(I)={s};
	end
	legend(M)
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
	
