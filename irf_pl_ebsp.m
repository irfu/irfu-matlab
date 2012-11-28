function h=irf_pl_ebsp(varargin)
%IRF_PL_EBSP   Plot E wavelet spectra, B wavelet spectra, Poynting flux,
% E/B, ellipticity, polarization
%
% h=irf_pl_ebsp(e,b,pos,arguments)
% modified from irf_pl_ebs
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
% e = wave electric field, columns (t ex ey ez)
% b = wave magnetic field, columns (t bx by bz)
% pos = position vector of spacecraft in GSE, columns (t x y z)
% choose which parameters to plot in arguments, options include:
%  'e', 'b', 's', 'eb', 'pol', 'ellip' %can be included in any order
% to set the frequency range, add the optional arguments
%  'freq', [fmin fmax]   %NOTE these must be the last two arguments
%  (Or use default [.01 5])
% 
% Returns in h handles to subplots
%
% Examples:
%    h=irf_pl_ebsp(e,b,xyz,'b','ellip','pol','e','s','freq',[.01 1.5]);
%    h=irf_pl_ebsp(e,b,xyz,'b','pol','ellip');


% %% Check input 
%   [ax,args,nargs] = axescheck(varargin{:});
%   %x=args{1};
%   if isempty(args), % nothing to plot, first input parameter empty
%     return;
%   end
%   
%   original_args=args;
%   var_desc{1} = '';
%   flag_subplot = 0;
%   have_options = 1;
% 
%   if nargs > 1, have_options = 1; end
%   
%   % TODO: this needs to be determined from data ?
%   sampl=22.5;
%   
%     %% Default values that can be overridden by options
% % dt = 0;
% % flag_yy = 0;
% % scaleyy = 1;
% % plot_type = '';
% % marker = '-';
% % flag_plot_all_data=1;
% flag_colorbar=1;
% colorbar_scale=1;
% 
%   %% Set some important parameters
%   freq_int=[.01 5];
%   freq_number=25;
%   Morlet_width=5.36;
% 
%   if isnumeric(args{end})
%       freq_int=args{end};
%       args=args(1:end-2);
%   elseif ismember({'freq'},args)
%       disp('frequency interval values missing. using default')
%       args=args(1:end-1);
%   end
%   
%   amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
%   a=logspace(amin,amax,anumber);
%   w0=sampl/2; % The maximum frequency
%   newfreq=w0./a;

%   %% Calculations
%   included_params=ismember({'e' 'b' 's' 'eb' 'pol' 'ellip'},lower(args));
%   if any(included_params(1:4) == 1) && all(included_params(5:6) == 0),
%     [e1,b1,powerEx_plot,powerEy_plot,powerEz_plot,power2E_plot,powerBx_plot,powerBy_plot,powerBz_plot,power2B_plot,Spar_plot,EtoB_plot,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_ebs(e1,b1,xyz,freq_int);
% %    [e1,power2E_plot,Spar_plot,power2B_SM_plot]=irf_ebs(e1,b1,xyz);
%     disp('irf_ebs ... calculations');
%   elseif all(included_params(1:4) == 0) && any(included_params(5:6) == 1),
%     [e1,b1,Lp,Ls5]=irf_pol(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5]=irf_pol(e1,b1,xyz);
%     disp('irf_pol ... calculations');
%   elseif all(included_params([1 3:4]) == 0) && any(included_params([2 5:6]) == 1),
%     [e1,b1,Lp,Ls5,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_bpol(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5,power2B_SM_plot]=irf_bpol(e1,b1,xyz);
%     disp('irf_bpol ... calculations');
%   else    
%     %for everything
%     [e1,b1,Lp,Ls5,powerEx_plot,powerEy_plot,powerEz_plot,power2E_plot,powerBx_plot,powerBy_plot,powerBz_plot,power2B_plot,Spar_plot,EtoB_plot,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_ebsp(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5,power2E_plot,Spar_plot,power2B_SM_plot]=irf_ebsp(e1,b1,xyz);
%     disp('irf_ebsp ... calculations');
%   end
  
  
  

% %% Check input options
% while have_options
%     l = 1;
%     switch(lower(args{1})) 
%         case 'linestyle'
%             marker = args{2};
%             l = 2;
%         case 'nocolorbar'
%             flag_colorbar=0;
%             l = 1;
%         case 
%         case 'freq'
%             if nargs>1
%                 if isnumeric(args{2})
%                     dt = args{2};
%                     l = 2;
%                 else irf_log('fcal,','wrongArgType : dt must be numeric')
%                 end
%             else irf_log('fcal,','wrongArgType : dt value is missing')
%             end
%         
%         otherwise
%             %irf_log('fcal',['Assuming ''' args{1} ''' is a LineStyle'])
%             marker = args{1};
%             args = args(2:end);
%             break
%     end
%     args = args(l+1:end);
%     if isempty(args), break, end
% end

%save('2009May22.mat', 'e1','b1','power2E_plot','power2B_SM_plot','Spar_plot','EtoB_plot');

  [ax,args,nargs] = axescheck(varargin{:});

%% plot everything
clear h;
npl=length(args)-3;clf;
ipl=1;
colorbar_scale=1;
if colorbar_scale==1,
   it2 = 0:.018:.9; it2=it2';
   it3 = ones(size(it2)); 
   it3(end)=.9;
%   it=0:.02:1;it=it'; xcm=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
   it=0:.02:1;it=it'; 
   xcm=[ [0*it flipud(it) it];[it2 it2 it3];...
       [flipud(it3) flipud(it2) flipud(it2)]; [flipud(it) 0*it 0*it]];
   xcm2=[ [it2 it2 it3];[flipud(it3) flipud(it2) flipud(it2)]];
   it4 = 0:.01:.5; it4=it4';
   it5 = 0:.008:.4; it5=it5';
   it6 = .4:.004:.6; it6=it6';
   it7 = .6:.008:1; it7=it7';
   it9 = 1:1:10;
   it10 = 1:1:9;
   xcm3=[ [it5 it./3 it4]; [it6 it./3+.33 0*it+.5]; [it7 it./3+.66 it4+.5]];
   %it8 = [it9*0 0:.04:.6 fliplr(0:.04:.6) it10*0]; it8=it8';
   it8 = [0:.024:.6 fliplr(.3+.015:.015:.6) fliplr(0+.06:.06:.3)]; it8=it8';
   it11 = [fliplr(.1:.016:.5) .1:.02:.5 .5:.02:.6-.04]; it11=it11';
   it12 = [0:.00375:.15 .15:.063:.78-.063]; it12=it12';
   %xcm4=[it8 it./1.3 it11];
   %xcm4=[it8 it12 it11];
   it13 = [0:.009:.405 .405:.02:.505-.02]; it13=it13';
   %xcm4=[it8+.1 it13+.1 it11+.1];
   it14 = .1:.01125:1; it14=it14';
   it15 = .1:.0025:.3; it15=it15';
   it16 = .2:.005:.6; it16=it16';
   it17 = .2:.035:.9; it17=it17';
   xcm4=[[it14 it15*0+.2 flipud(it16)]; [ones(size(it17)) it17 it17]];
   clear it; clear it2; clear it3;
else
      colormap('default');xcm=colormap;
end

load caa/cmap.mat % default map
cmapStandard=cmap;
it=0:.02:1;it=it(:);
cmapPoynting=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
%cmapCombo=[cmapStandard;cmapPoynting];
cmapCombo=[cmapStandard;xcm2];

   %[ax,args,nargs] = axescheck(varargin{:});
    %t=e(:,1);
    timeVector=args{1};
    t=timeVector;
    t_start_epoch=get_t_start_epoch(t(1,1));

    frequencyVector=args{2};
    newfreq=frequencyVector;
    BVector=args{3};
    args=args(4:end);

%     power2B_SM_plot = BB_xxyyzz_fac;
%     power2E_plot = EESum_xxyyzz_ISR2;
%     %EE_xxyyzz_FAC = ;
%     Spar_plot = Poynting_xyz_FAC;
%     %Poynting_rThetaPhi_FAC = ;
%     polarizationSign = polarization;
%     %polarizationEllipseRatio.*polarizationSign = ellipticity;
    power2B_SM_plot = args{1};
    power2E_plot = args{2};
    EESum_xxyyzz_ISR2 = args{2};
    EE_xxyyzz_FAC = args{3};
    Spar_plot = args{4};
    Poynting_rThetaPhi_FAC = args{5};
    k_thphSVD_fac = args{6};
    polSVD_fac = args{7};
    %polarizationEllipseRatio.*polarizationSign = ellipticity;
    ellipticity = args{8};


    %%%%%%%%% E spectra %%%%%%%%%%%%
    
    if nargs==8 || nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        %pcolor(t-t_start_epoch,newfreq,log10(abs(power2E_plot.'))) % With edge effects removed
%         if max(newfreq)>.01,
%             freqIndex=find(newfreq>.1);
%             freqIndex2=find(newfreq<=.1);
%             if min(newfreq)>.01, 
%                 pcolor(t-t_start_epoch,newfreq,log10(abs(EE_xxyyzz_FAC.')))
%             else    
%                 EEcombined=EE_xxyyzz_FAC;
%                 EEcombined(:,freqIndex2)=EESum_xxyyzz_ISR2(:,freqIndex2);
%                 %pcolor(t-t_start_epoch,newfreq(freqIndex),log10(abs(EE_xxyyzz_FAC(:,freqIndex).'))) % With edge effects removed
%                 %pcolor(t-t_start_epoch,newfreq(freqIndex2),log10(abs(EESum_xxyyzz_ISR2(:,freqIndex2).'))) % With edge effects removed
%                 pcolor(t-t_start_epoch,newfreq,log10(abs(EEcombined.'))) % With edge effects removed
%             end
%         else 
%             pcolor(t-t_start_epoch,newfreq,log10(abs(EESum_xxyyzz_ISR2.'))) % With edge effects removed
%         end    
         %pcolor(t-t_start_epoch,newfreq,log10(abs(EE_xxyyzz_FAC.'))) % With edge effects removed
        pcolor(t-t_start_epoch,newfreq,log10(abs(EESum_xxyyzz_ISR2.'))) % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log','tickdir','out');
        if any(isnan(power2E_plot)),
            caxis([-3.5 3.5]);
        else
            %cmean=nanmean(nanmean(log10(abs(power2E_plot))));
            cmean=nanmean(nanmean(log10(abs(EE_xxyyzz_FAC))));
            caxis(floor(cmean)+[-3.5 3.5]);
        end
        %caxis([-3.5 3.5]);
        colormap(xcm);
        %freezeColors;
        hca = colorbar;
        %set(h, 'ylim', [0 257]);
        %hca=colorbar('peer',h(1));
        ylabel(hca,{'log(E)'; '[(mV/m)^2/Hz]'});
        %irf_colormap(hca,'standard');
        %cbfreeze(hca)
        
%         h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
%         %pcolor(t-t_start_epoch,newfreq,log10(abs(power2E_plot.'))) % With edge effects removed  
%         pcolor(t-t_start_epoch,newfreq,log10(abs(EE_xxyyzz_FAC.'))) % With edge effects removed
%         shading flat
%         ylabel('f [Hz]')
%         set(gca,'yscale','log','tickdir','out');
%         if any(isnan(power2E_plot)),
%             caxis([-3.5 3.5]);
%         else
%             %cmean=nanmean(nanmean(log10(abs(power2E_plot))));
%             cmean=nanmean(nanmean(log10(abs(EE_xxyyzz_FAC))));
%             caxis(floor(cmean)+[-3.5 3.5]);
%         end
%         %caxis([-3.5 3.5]);
%         %colormap(xcm4);
%         %freezeColors;
%         hca = colorbar;
%         ylabel(hca,{'log(E)'; '[(mV/m)^2/Hz]'});
%         %cbfreeze(hca)
        
    end

    %%%%% B spectra from spectral matrix
   
    h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
    hold on
    pcolor(t-t_start_epoch,newfreq,log10(abs(power2B_SM_plot(:,:,4).'))) 
%    pcolor(t-t_start_epoch,newfreq,log10(abs(power2B_plot.'))) 
    shading flat
    ylabel('f [Hz]')
    set(gca,'yscale','log','tickdir','out');
  %  set(gca,'tickdir','out');
    cmean=nanmean(nanmean(log10(abs(power2B_SM_plot(:,:,4)))));
    %hold on
    %plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1.6e-19./1.67e-27);
    plot(BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi,'color','k');
    plot(BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4,'--','color','k');
    plot(BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16,'-.','color','k');
    axis([min(t-t_start_epoch) max(t-t_start_epoch) min(newfreq) max(newfreq)]);
    hold off
    caxis(floor(cmean)+[-3.5 3.5]);
  %  caxis(floor(cmean)+[-2.5 2.5]);
%    caxis([-6.5 0.5]);
    %colormap(cmapCombo);
    %freezeColors;
    %hca2 = colorbar;
    hca2=colorbar('peer',h(2));
    %set(h, 'ylim', [0 257]);
    %irf_colormap(hca2,'standard');
    ylabel(hca2,{'log(B)'; '[nT^2/Hz]'});
    %cbfreeze(hca)
    
    %%%%%% propogation direction from SVD
    
    if nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        pcolor(t-t_start_epoch,newfreq,180./pi.*k_thphSVD_fac(:,:,1).') % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log');set(gca,'tickdir','out');
    %    set(gca,'tickdir','out');
        caxis([0 90]);
        %colormap(cmapCombo);
        %freezeColors;
        hca3 = colorbar;
        %hca3=colorbar('peer',h(3));
        %set(h, 'ylim', [0 257]);
        ylabel(hca3,'k (theta)');
        %irf_colormap(hca3,'standard');
        %cbfreeze(hca)

        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        pcolor(t-t_start_epoch,newfreq,180./pi.*k_thphSVD_fac(:,:,2).') % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log');set(gca,'tickdir','out');
    %    set(gca,'tickdir','out');
        caxis([-180 180]);
        %colormap(cmapCombo);
        %freezeColors;
        %hca4 = colorbar;
        hca4=colorbar('peer',h(4));
        %set(h, 'ylim', [257 461]);
        ylabel(hca4,'k (phi)');
        %irf_colormap(hca4,'poynting');
        %cbfreeze(hca)
    end


      %%%%% Degree of Polarization
    
    if nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        pcolor(t-t_start_epoch,newfreq,polSVD_fac.') % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log');set(gca,'tickdir','out');
    %    set(gca,'tickdir','out');
        caxis([0 1]);
        %colormap(cmapCombo);
        %freezeColors;
        %hca5 = colorbar;
        hca5=colorbar('peer',h(5));
        %set(h, 'ylim', [0 257]);
%         colorbarHandles=findobj(gca,'Tag','Colorbar') 
%         for i=1:length(colorbarHandles) 
%             if colorbarHandles(i)~=hca 
%                 cBarAxis1=colorbarHandles(i); 
%             end 
%         end
        
        ylabel(hca5,{'Degree of'; 'Polarization'});
        %irf_colormap(hca5,'standard');
        %cbfreeze(hca)
    end


    %%%%%% Ellipticity
    
    if nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
    %    pcolor(t-t_start_epoch,newfreq,log10(Pol.')) 
        %pcolor(t-t_start_epoch,newfreq,(polarizationEllipseRatio.*polarizationSign).') 
        pcolor(t-t_start_epoch,newfreq,(ellipticity).') 
    %    pcolor(t-t_start_epoch,newfreq,(Lp.*Phase_dif./abs(Phase_dif)).') 
    %    pcolor(t-t_start_epoch,newfreq,Lp.') 
        shading flat
        ylabel('f [Hz]')
        %ht=text(0,0,'log10(E/B) [(1000 km/s)]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
        set(gca,'yscale','log');set(gca,'tickdir','out');
    %    set(gca,'tickdir','out');
        %hca6 = colorbar;
        %hca6=colorbar('peer',h(6));
        caxis([-1 1]);
        %cLim=[257,461];
    %    caxis([-2 2]);
        %colormap(cmapCombo);
        %freezeColors;
        hca6=colorbar('peer',h(6));
        %set(hca6, 'ylim', [257 462]);
        ylabel(hca6,'Ellipticity');
        %irf_colormap(hca6,'poynting');
        %freezeColors;
        %hca = colorbar;
        
        %cbfreeze(hca)
    %    ylabel(hca,'log10(E/B) [(1e3 km/s)]');
        
        %ylabel(hca,'Ellipticity');

        %from irf_colormap.m
%         hcb = hca;
%         hy=get(hcb,'ylabel');
%         ylabel_string=get(hy,'string');
%         ylabel_fontsize=get(hy,'fontsize');
%         new_hcb = cbfreeze(hcb);
%         new_hy=get(new_hcb,'ylabel');
%         set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
%         
        
        %hca=cbfreeze(hca)
    end


      %%%%%%%%% S spectra %%%%%%%%%%%%
    if nargs==8 || nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        %pcolor(t-t_start_epoch,newfreq,(sign(Spar_plot(:,:,1)+Spar_plot(:,:,2)).*sqrt(abs(Spar_plot(:,:,1)+Spar_plot(:,:,2)))).') % With edge effects removed
        pcolor(t-t_start_epoch,newfreq,(sign(Poynting_rThetaPhi_FAC(:,:,1)).*sqrt(abs(Poynting_rThetaPhi_FAC(:,:,1)))).') % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log');set(gca,'tickdir','out');

        if any(isnan(Spar_plot)),
%            caxis([-10 10]);
            %caxis([-50 50]);
            caxis([-30 30]);
        else
            cc = [-max(max(sqrt(abs(Spar_plot)))) max(max(sqrt(abs(Spar_plot))))];
            caxis(cc);
        end
    %    caxis([-10 10]);
        %colormap(cmapCombo);
        %freezeColors;
        %hca7 = colorbar;
        hca7=colorbar('peer',h(7));
        %set(h, 'ylim', [257 461]);
        %ylabel(hca,{'S_{perp}'; '[\mu W/m^2Hz]^{1/2}'});
        ylabel(hca7,{'S_{r}'; '[\mu W/m^2Hz]^{1/2}'});
        %cbfreeze(hca)
        %unfreezeColors;
        %irf_colormap(hca7,'poynting');
    end
    if nargs==8 || nargs==11,
        h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
        pcolor(t-t_start_epoch,newfreq,(sign(Spar_plot(:,:,3)).*sqrt(abs(Spar_plot(:,:,3)))).') % With edge effects removed
        shading flat
        ylabel('f [Hz]')
        set(gca,'yscale','log');set(gca,'tickdir','out');

        if any(isnan(Spar_plot)),
%            caxis([-10 10]);
            %caxis([-50 50]);
            caxis([-30 30]);
        else
            cc = [-max(max(sqrt(abs(Spar_plot)))) max(max(sqrt(abs(Spar_plot))))];
            caxis(cc);
        end
    %    caxis([-10 10]);
        %colormap(cmapCombo);
        %freezeColors;
        %hca8 = colorbar;
        hca8=colorbar('peer',h(8));
        %set(h, 'ylim', [257 461]);
        ylabel(hca8,{'S_{II}'; '[\mu W/m^2Hz]^{1/2}'});
        %cbfreeze(hca)
        %unfreezeColors;
        %irf_colormap(hca8,'standard');
        %colorbarHandles=findobj(gca,'Tag','Colorbar');
        %length(colorbarHandles)
        
%         irf_colormap(h(1),'standard');
%         irf_colormap(h(2),'standard');
%         irf_colormap(h(3),'standard');
%         irf_colormap(h(4),'poynting');
%         irf_colormap(h(5),'standard');
%         irf_colormap(h(6),'poynting');
%         irf_colormap(h(7),'poynting');
%         %ylabel(h(8),{'S_{II}'; '[\mu W/m^2Hz]^{1/2}'});
%         irf_colormap(h(8),'standard');
        
    end
    


%     %%%%%%%% E/B spectra %%%%%%%%%%%%
%     elseif strcmp(plot_param,'eb'),
%         pcolor(t-t_start_epoch,newfreq,log10(abs(EtoB_plot.'))) % With edge effects removed
%         shading flat
%         ylabel('f [Hz]')
%         %ht=text(0,0,'log10(E/B) [(1000 km/s)]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
%         set(gca,'yscale','log');set(gca,'tickdir','out');
%         caxis([-1 4]);
%         colormap(xcm);
%         hca = colorbar;
%         ylabel(hca,'log10(E/B) [(1e3 km/s)]');
% 
%     else
%         
%     end


    %% Add figure menu
    irf_figmenu;


      axes(h(1));
      %title(['Width Morlet wavelet = ' num2str(Morlet_width)]);
      %ht=irf_pl_info([mfilename '  ' datestr(now)]); set(ht,'interpreter','none'); % add information to the plot
      irf_zoom(h,'x',[min(t) max(t)]);
      irf_legend(h(1),'Figure reference',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);
      irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
      irf_timeaxis(h)
    %irf_legend('C4', 'color','cluster');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_start_epoch=get_t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud=get(gcf,'userdata');
ii = find(~isnan(t));
if ii,
  valid_time_stamp=t(ii(1));
else
  valid_time_stamp=[];
end

if isfield(ud,'t_start_epoch'),
  t_start_epoch=ud.t_start_epoch;
elseif valid_time_stamp,
  if valid_time_stamp > 1e8, % set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
    t_start_epoch=valid_time_stamp;
    ud.t_start_epoch=t_start_epoch;
    set(gcf,'userdata',ud);
    irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
  else
    t_start_epoch=0;
  end
else
  t_start_epoch=0;
end
end

function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   Revision: 1.1.8.1   Date: 2010/03/16 00:15:50 

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end
end

