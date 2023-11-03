


function specFlux=c_rapid_pad_plot_bf2004(varargin)
% C_RAPID_PAD_PLOT_BF2004 plot the electron pitch angle distribution of RAPID as function of time
%
% C_RAPID_PAD_PLOT_BF2004(dataObjectEPADEX,dataObjectEPITCH,energyChannel);
%	Input dataobjects are cdf files downloaded from CAA and loaded with CAA_LOAD
%		dataObjectEPADEX - caa data object C?_CP_RAP_EPADEX?
%		dataObjectEPITCH - caa data object C?_CP_RAP_EPITCH
%		energyChannel    - channel 1 is 40.7 keV; channel 2 is 68.1keV
%
% C_RAPID_PAD_PLOT_BF2004(AX,...)
%		plot in axis with handle AX
%		if AX not given plot in current axis GCA
%
%---------------------------------------
%Example:
%c_rapid_pad_plot_bf2004(C1_CP_RAP_EPADEX2, C1_CP_RAP_EPITCH, 1);
%---------------------------------------


[ax,args,nargs] = axescheck(varargin{:});
if nargs==0 % nothing to plot
  return;
end
if isempty(ax) % if no axis input, then plot in current axis
  ax=gca;
end
RAPflux_obj=args{1};
RAPpitang_obj=args{2};
E_channel=args{3};


varsFlux=RAPflux_obj.Variables;
Flux_str=varsFlux{3,1};
%Azm_str=varsFlux{6,1};  % NEVER USED
Time_str=varsFlux{2,1};
varsPitang=RAPpitang_obj.Variables;
Pitang_str=varsPitang{3,1};

RapFlux=getmat(RAPflux_obj, Flux_str);
RapPitang=getmat(RAPpitang_obj,Pitang_str);
%RapAzmang=getmat(RAPflux_obj, Azm_str); % NEVER USED
Raptime=getmat(RAPflux_obj, Time_str);
Fluxunits=getunits(RAPflux_obj,Flux_str);

aaa=RapFlux.data;
bbb=RapPitang.data;
ttt=Raptime(:,1);


Pitang=[10 30 50 70 90 110 130 150 170];
angtemp=[0 20 40 60 80 100 120 140 160 180];
Flux=zeros(length(ttt),length(Pitang));
count=zeros(length(ttt),length(Pitang));

for kk=1:length(ttt)
  Flux_each=squeeze(aaa(kk,E_channel,:,:));
  Pitang_each=squeeze(bbb(kk,:,:));
  for ii=1:16
    for jj=1:9
      if isnan(Flux_each(ii,jj))==0 && isnan(Pitang_each(ii,jj))==0
        ind=find(angtemp>=Pitang_each(ii,jj)); index=ind(1)-1;
        if index==0; index=1; end

        Flux(kk,index)=Flux(kk,index)+Flux_each(ii,jj);
        count(kk,index)=count(kk,index)+1;
      end
    end
  end
  for jj=1:length(Pitang)
    if count(kk,jj)==0
      Flux(kk,jj)=NaN;
    else
      Flux(kk,jj)=Flux(kk,jj)/count(kk,jj);
    end
  end
end

specFlux=struct('t',ttt,'f',Pitang,'p',Flux,'f_unit','Pitch angle [deg]');
if nargout == 0 % no output, just plot the spectrogram
  irf_spectrogram(ax, specFlux);
  if isa(ax,'handle'), hcb = colorbar(ax); % HG2
  else, hcb = colorbar('peer',ax);
  end
  colormap(hcb,jet);
  ylabel(hcb,Fluxunits);
  clear specFlux; % do not return anything
end
end





