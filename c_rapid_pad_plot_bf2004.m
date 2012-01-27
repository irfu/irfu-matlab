


function c=c_rapid_pad_plot_bf2004(varargin)
% c_rapid_pad_plot_bf2004 plot the electron pitch angle distribution of RAPID as function of time
%
%Input parameter (3 in total, downloaded from CAA website and then constructed into the object format)
%1. C?_CP_RAP_EPADEX?
%2. C?_CP_RAP_EPITCH
%3. E_channel=1;             channel 1 is 40.7 keV; channel 2 is 68.1keV
%
%---------------------------------------
%Example: 
%c_rapid_pad_plot_bf2004(C1_CP_RAP_EPADEX2, C1_CP_RAP_EPITCH, 1);
%---------------------------------------


[ax,args,nargs] = axescheck(varargin{:});
RAPflux_obj=args{1};
RAPpitang_obj=args{2};
E_channel=args{3};


    varsFlux=RAPflux_obj.Variables;
    Flux_str=varsFlux{3,1};  Azm_str=varsFlux{6,1};  Time_str=varsFlux{2,1}; 
    varsPitang=RAPpitang_obj.Variables;
    Pitang_str=varsPitang{3,1};

    RapFlux=getmat(RAPflux_obj, Flux_str);
    RapPitang=getmat(RAPpitang_obj,Pitang_str);
    RapAzmang=getmat(RAPflux_obj, Azm_str);
    Raptime=getmat(RAPflux_obj, Time_str);
    
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
                ind=find(angtemp>=Pitang_each(ii,jj)); index=ind(1)-1;
                   if index==0; index=1; end
                   if isnan(Flux_each(ii,jj))==0; 
                       Flux(kk,index)=Flux(kk,index)+Flux_each(ii,jj);
                       count(kk,index)=count(kk,index)+1;
                   end
                
            end
        end
        for jj=1:length(Pitang)
            if count(kk,jj)==0,
               Flux(kk,jj)=NaN;
            else
               Flux(kk,jj)=Flux(kk,jj)/count(kk,jj);
            end
        end
    end
    
    
    specFlux=struct('t',ttt,'f',Pitang,'p',Flux,'f_unit','Pitch angle [deg]');
    
    irf_spectrogram(gca, specFlux);
    colorbar;
    colormap(jet);


end




    
   