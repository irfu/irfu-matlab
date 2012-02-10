function irf_polar_distr(varargin)
% IRF_POLAR_DISTR  Plot polar distribution
%
% Plots a polar distribution of particle data.
% First argument is time interval in epoch.
% Can plot PEACE PITCH and 3DXPH products.
% 
% Example: irf_polar_distr(tint,'C4_CP_PEA_PITCH_3DXH_DEFlux')
%          irf_polar_distr(tint,'C4_CP_PEA_3DXPH_PSD')

% For PEACE only PITCH and 3DXPH because 3DXPH is the only one included 
% in c_caa_construct_subspin_res_data. The angular resolution also depends
% on this function.

% If two products are given, plot together in 
% opposite side of plot. Not yet implemented.

% Matlab function surf is used. The grid is defined by the energy and pitch
% angle intervals, i.e. one value extra per dimension. The color of each
% square/interval is taken from the data.

set(0,'defaultAxesFontSize',14);
set(0,'defaultTextFontSize',14);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');

time=varargin{1};
switch size(varargin,2) % One product given
    case 2
        product=varargin{2};
        
        % Used for colorbar label
        if any(strfind(product,'PSD')); distr='PSD';
        elseif any(strfind(product,'DEFlux')); distr='DEFlux';
        elseif any(strfind(product,'DPFlux')); distr='DPFlux';
        end
        
        if any(strfind(product,'PITCH'))
            %%
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',product]);
            [caaSEDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', product]);
            [caaSEDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product]);
            
            theta=Data.dep_x{2}.data(1,:);
            t=Data.t;
            en=[Data.dep_x{3}.data(1,:)-SEDL(1,2:end) Data.dep_x{3}.data(1,end)+SEDU(1,end)];
            nan_en=isnan(en);en(nan_en)=[];
            phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
            
            dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            
            % Define subspin time form angle phi.
            tt=subspintime(dataobject,phi);
            [~,ind]=irf_tlim(tt,time);
            
            specrec.f=theta;
            specrec.en=flipdim(en,2);
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' Data_units ']'];
            specrec.f_label='Pitch angle';
            specrec.data=double(data(ind,:,:));
                 
            r=double(specrec.en)';
            rlog = log10( r ); % energy levels in eV
            r_man = rlog-rlog(1);
            r0=rlog(1)-r_man(1);
            theta = double(specrec.f)-90-(specrec.f(2)-specrec.f(1))/2; % pitch angles
            theta=[theta theta(end)+(specrec.f(2)-specrec.f(1))]; % Kolla vinklarna idag
            X = r_man*cosd(theta); % xy-coordinates
            Y = r_man*sind(theta); % xy-coordinates

            resultat=log10(specrec.data(:,:,:));
            C=nanmean(resultat,1);
            CC=flipdim(squeeze(C),2);
            CC=CC';

            Xplot=[-flipdim(X(2:end,:),1); X];
            Yplot=[flipdim(Y(2:end,:),1); Y];
            Cplot=[flipdim(CC(1:end,:),1); CC(1:end,:)];

            figure;
            h=surf(Xplot,Yplot,Xplot*0,Cplot(1:end,:));
            view(2)
            axis equal tight
            shading flat 
            grid off
            cb=colorbar;
            ylabel(cb,specrec.p_label)
            xlabel('Energy  [keV]')
            ylabel('Energy  [keV]')
            
            t1str=datestr(epoch2date(time(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(time(2)),'HH:MM:SS.FFF');
            title([product,'    ',t1str,'-',t2str,'UT'])   
            
            % Ticks
            xticks=log10([1e-1 1e0 1e1 1e10 1e100])+r0;
            xticks=[xticks(find(xticks>0))];
            xticks=[xticks(find(xticks<rlog(end)))];
            xticklabels=cell(size(xticks));
            for k=1:length(xticklabels)
            xticklabels{k}=num2str(10.^(xticks(k)-r0),'%.1f');
            end         
            xticks=[-flipdim(xticks,2) xticks];
            xticklabels=[flipdim(xticklabels,2) xticklabels];
            yticks=xticks;
            yticklabels=xticklabels; 
            set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
                'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)
        elseif any(strfind(product,'3DXPH'))
            %%
            res=c_caa_construct_subspin_res_data(['Data__', product]);
            [caaSEDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', product]);
            [caaSEDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product]);
            SEDL=flipdim(SEDL(1,2:end),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
            SEDU=flipdim(SEDU(1,2:end),2); en_nan=isnan(SEDU);SEDU(en_nan)=[];                    
            [~,ind]=irf_tlim(res.tt,time);
            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            specrec.en=res.en(:);
            specrec.en=[res.en-SEDL res.en(end)+SEDU(end)];
            specrec.data=res.data(ind,:,:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];
           
            thi=1:length(res.theta);                       
            r=double(specrec.en)';
            rlog = log10( r ); % energy levels in eV
            r_man = rlog-rlog(1);
            r0=rlog(1)-r_man(1);
            theta = double(specrec.f(thi))-90-(specrec.f(2)-specrec.f(1))/2; % pitch angles
            theta=[theta theta(end)+(specrec.f(2)-specrec.f(1))]; % Kolla vinklarna idag
            X = r_man*cosd(theta); % xy-coordinates
            Y = r_man*sind(theta); % xy-coordinates

            resultat=log10(specrec.data(:,thi,:));
            C=nanmean(resultat,1);
            CC=squeeze(C);
            CC=CC';

            Xplot=[-flipdim(X(2:end,:),1); X];
            Yplot=[flipdim(Y(2:end,:),1); Y];
            Cplot=[flipdim(CC(1:end,:),1); CC(1:end,:)];

            figure;
            h=surf(Xplot,Yplot,Xplot*0,Cplot(1:end,:));
            view(2); axis equal tight; shading flat; grid off
            cb=colorbar;
            ylabel(cb,specrec.p_label)
            xlabel('Energy  [keV]')
            ylabel('Energy  [keV]')
            
            % Ticks
            xticks=log10([1e-1 1e0 1e1 1e10 1e100])+r0;
            xticks=[xticks(find(xticks>0))];
            xticks=[xticks(find(xticks<rlog(end)))];
            xticklabels=cell(size(xticks));
            for k=1:length(xticklabels)
            xticklabels{k}=num2str(10.^(xticks(k)-r0),'%.1f');
            end
            xticks=[-flipdim(xticks,2) xticks];
            xticklabels=[flipdim(xticklabels,2) xticklabels];
            yticks=xticks;
            yticklabels=xticklabels;            
            set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
                'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)
            
            t1str=datestr(epoch2date(time(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(time(2)),'HH:MM:SS.FFF');
            title([product,'    ',t1str,'-',t2str,'UT'])
            
        end
    case 3 % Plot two distributions in opposite sides
        product1=varargin{2};
        product2=varargin{3};
        if any(strfind(product1,'PSD')); distr1='PSD';
        elseif any(strfind(product1,'DEFlux')); distr1='DEFlux';
        elseif any(strfind(product1,'DPFlux')); distr1='DPFlux';
        end
        if any(strfind(product2,'PSD')); distr2='PSD';
        elseif any(strfind(product2,'DEFlux')); distr2='DEFlux';
        elseif any(strfind(product2,'DPFlux')); distr2='DPFlux';
        end
        
        if any(strfind(product,'3DXPH'))
        elseif any(strfind(product,'PITCH_3DXH')) || any(strfind(product,'PITCH_3DRH')) || any(strfind(product,'PITCH_3DRL'))
            %%
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',product]);
            theta=Data.dep_x{2}.data(1,:);
            t=Data.t;
            en=Data.dep_x{3}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
            phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
            dataraw=Data.data;
            dataraw(:,:,:,nan_en)=[];
            dataraw(:,nan_phi,:,:)=[];

            %dataraw(dataraw==caaData.FILLVAL)=NaN;
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            
            tt=subspintime(dataobject,phi);
            [~,ind]=irf_tlim(tt,time);
            
            repmat(t(:),1,length(phi));
            specrec.f=theta;
            specrec.data=data;
            specrec.en=en;
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' ']'];
            specrec.f_label='Pitch angle';
            specrec.data=double(specrec.data(ind,:,:));
            
            thi=1:length(theta);
            [caaSweepEnergyDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', product]);
            [caaSweepEnergyDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product]);
            SEDU
            SEDL
            r=double(specrec.en)';
            r=[r(1)-(r(2)-r(1))/2; r(1:end-1)+(r(2:end)-r(1:end-1))/2 ;r(end)+(r(end)-r(end-1))/2];
            rlog = log10( r ); % energy levels in eV
            r_man = rlog-rlog(1);
            r0=rlog(1)-r_man(1);
            %[r_man; 2*r_man(end)-r_man(end-1)];
            theta = double(specrec.f(thi))-90-(specrec.f(2)-specrec.f(1))/2; % pitch angles
            theta=[theta theta(end)+(specrec.f(2)-specrec.f(1))]; % Kolla vinklarna idag
            X = r_man*cosd(theta); % xy-coordinates
            Y = r_man*sind(theta); % xy-coordinates

            resultat=log10(specrec.data(:,thi,:));
            C=nanmean(resultat,1);
            CC=squeeze(C);
            CC=CC';

            Xplot=[-flipdim(X(2:end,:),1); X];
            Yplot=[flipdim(Y(2:end,:),1); Y];
            Cplot=[flipdim(CC(1:end,:),1); CC(1:end,:)];

            figure;
            h=surf(Xplot,Yplot,Xplot*0,Cplot(1:end,:));
            view(2)
            axis equal tight
            shading flat 
            grid off
            cb=colorbar;
            ylabel(cb,specrec.p_label)
            xlabel('Energy  [keV]')
            ylabel('Energy  [keV]')
            t1str=datestr(epoch2date(tin(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(tin(2)),'MM:SS.FFF');
            title([product,'    ',t1str,'-',t2str,'UT'])   
            
            % Ticks
            xticks=log10([1e-1 1e0 1e1 1e10 1e100])+r0;
            xticks=[xticks(find(xticks>0))];
            xticks=[xticks(find(xticks<rlog(end)))];
            xticklabels=cell(size(xticks));
            for k=1:length(xticklabels)
            xticklabels{k}=num2str(10.^(xticks(k)-r0),'%.1f');
            end
            
            xticks=[-flipdim(xticks,2) xticks];
            xticklabels=[flipdim(xticklabels,2) xticklabels];
            yticks=xticks;
            yticklabels=xticklabels;
        elseif any(strfind(product,'3DXPH'))


        end
end

% Add pitch angle labels
text(0,rlog(end)-r0-0.5,0,'0^o')
text(0,-rlog(end)+r0+0.5,0,'180^o')
text(-rlog(end)+r0+0.5,0,0,'90^o')
text(rlog(end)-r0-0.5,0,0,'90^o')

end

function [tt,dtsampling]=subspintime(dataobject,phi)
% construct subspin time vector
% phi are azimuthal angles (spin period is divided in the number of azimuth
% angles)
timevar=getv(dataobject,dataobject.VariableAttributes.DEPEND_0{1,2});
tt=timevar.data(:);
tt=repmat(tt,1,length(phi));

if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
    if ischar(timevar.DELTA_PLUS)
        deltaplus= getv(dataobject,timevar.DELTA_PLUS);
        dtplus=deltaplus.data(1,:);
        if dtplus>5, % temporary solution for CIS problems
            if (dtplus/2>3.5) && (dtplus/2 < 4.5), dtplus=dtplus/2;
            elseif (dtplus/3>3.5) && (dtplus/3 < 4.5), dtplus=dtplus/3;
            elseif (dtplus/4>3.5) && (dtplus/4 < 4.5), dtplus=dtplus/4;
            end
        end
    elseif isnumeric(timevar.DELTA_PLUS)
        dtplus=timevar.DELTA_PLUS;
    end
    if ischar(timevar.DELTA_MINUS)
        deltaminus= getv(dataobject,timevar.DELTA_MINUS);
        dtminus=deltaplus.data(1,:);
    elseif isnumeric(timevar.DELTA_MINUS)
        dtminus=timevar.DELTA_MINUS;
    end
else
    dtplus=2;
    dtminus=2;
end
spin_period=dtplus+dtminus;
dtsampling=spin_period/length(phi);
for j=length(phi):-1:1,
    tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
end
tt=reshape(tt',numel(tt),1);
end
