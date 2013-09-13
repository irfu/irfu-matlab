function c_caa_distribution_function(varargin)
% C_CAA_DISTRIBUTION_FUNCTION  Plot particle distribution.
%
% Plots a polar, cross-section distribution of particle data.
% Time interval is specified as 'tint',[tint_in_epoch].
% Axis handle can be specified.
% It is possible to define the origion by option 'emin', can be
% 0<emin<log10(minE)
% The type of plot is specified as an option: 
%   'polar' - Default, can plot one or two products at a time.
%   'cross-section' - Plots cross-sections of the distribution 
%                     function at 0, 90 and 180 degrees as well
%                     as zero-count level. Can plot one product.
%                     at a time.
% 
% Supported products: PEA: PITCH, 3DXPH 
%                     CIS: CODIF_LH, CODIF_HS, HIA_HS_MAG
%                     RAP: PAD
% 
% Example: c_caa_distributon_function(h,'C4_CP_RAP_PAD_E3DD','tint',tint)
%          c_caa_distributon_function('C4_CP_PEA_3DXPH_PSD','tint',tint)
%          c_caa_distributon_function(h,'C4_CP_PEA_PITCH_3DXH_DPFlux','tint',tint,'cross-section')
%          c_caa_distribution_function('tint',tint,'C3_CP_PEA_PITCH_3DXH_DEFlux','C3_CP_CIS_HIA_HS_MAG_IONS_PEF','polar')

% Current issues: o 3DRL: problem creating subspin time since sampling
% period seems very irregular. dt varies between 20 and 12 whereas it is
% constantly ~4 for 3DXH. 
% o Colorbar: While it is possible to plot two different products in the
% same plot, the color coding can not be separated.
% o There's something strange with background level for PITCH products

% For PEACE only PITCH and 3DXPH because 3DXPH is the only one included 
% in c_caa_construct_subspin_res_data. The angular resolution also depends
% on this function.

% Matlab function surf is used. The grid is defined by the energy and pitch
% angle intervals, i.e. one value extra per dimension. The color of each
% square/interval is taken from the data.

% General figure options
set(0,'defaultAxesFontSize',14);
set(0,'defaultTextFontSize',14);
set(0,'defaultAxesFontUnits','pixels');
set(0,'defaultTextFontUnits','pixels');

% Read input
[ax,args,nargs] = axescheck(varargin{:});
original_args=args;

% Default values that one can override with options.
plot_type='polar';
n_products=0;
emin=[]; % Defines the origin
products={};

% Check for products
l=1;
while l<=nargs
    if any(strfind(args{l},'RAP'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='RAP';
    elseif any(strfind(args{l},'CIS'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='CIS'; 
    elseif any(strfind(args{l},'PEA'))
       n_products=n_products+1; 
       products{n_products}=args{l};
       product_type{n_products}='PEA';       
    else
        switch lower(args{l})
            case 'tint'
                l=l+1;
                tint=args{l};
            case 'polar'
                plot_type='polar';
            case 'cross-section'
                plot_type='cross-section';              
            case 'emin'
                l=l+1;
                emin=args{l};
        end
    end
l=l+1;
end

% Return if not enough input is given.
% Reduce products if too much input is given.
if isempty(tint)
    disp('No time interval was given!')
    return;
end
switch n_products
    case 0
        disp('No particle product was given!')
        return;
    case 1
    case 2
        if strcmp(plot_type,'cross-section')
            disp('Cross-section distributions are only plotted for single products.')
            disp('Plotting first product given.')
            products=products(1);
        end
    otherwise
        switch plot_type
            case 'polar'
                disp('Polar distributions are plotted for at most two products.')
                disp('Plotting two first products given.')
                products=products(1:2);
                n_products=2;
            case 'cross-section'
                disp('Cross-section distributions are only plotted for single products.')
                disp('Plotting first product given.')
                products=products(1);
                n_products=1;
        end
end

% If no axes is given, initialize figure.
if isempty(ax) 
    ax=irf_plot(1);
end

% Plot distributions
switch plot_type
    case 'polar' % Plot polar distribution             
        for k=1:n_products % Prepare matrices for surf
            % For colorbar label
            distr{k}=products{k}(max(strfind(products{k},'_'))+1:end);

            % Preparing different sets of data: theta,
            % specrec.data,specrec.en
            if any(strfind(products{k},'PITCH'))
                [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',products{k}]);

                % Set up grid
                theta_plus=Data.dep_x{2}.df.plus;
                theta_minus=Data.dep_x{2}.df.minus;
                theta=Data.dep_x{2}.data(1,:);
                nan_theta=isnan(theta);theta(nan_theta)=[]; % for data
                theta=[theta-theta_minus theta(end)+theta_plus(end)]; % grid

                en_plus=Data.dep_x{3}.df.plus;
                en_minus=Data.dep_x{3}.df.minus;
                en=Data.dep_x{3}.data(1,:);
                nan_en=isnan(en);en(nan_en)=[]; % for data
                en_plus(nan_en)=[];en_minus(nan_en)=[];                    
                en=[en+en_plus en(end)-en_minus(end)];

                phi=Data.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];

                % Construct subspin reoslution data
                dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
                dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
                data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
                
                [tt dtsampling]=subspintime(dataobject,phi); % Define subspin time from angle phi.
                [~,ind]=irf_tlim(tt,tint);

                specrec.f=theta;
                specrec.en=flipdim(en,2); % stigande
                specrec.t=tt;
                specrec.p_label{k}=['Log ' distr{k} ' [' Data_units ']'];
                specrec.f_label='Pitch angle';
                specrec.data=double(data(ind,:,:));
                specrec.data=flipdim(specrec.data,3);
                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN;                  
            end
            if any(strfind(products{k},'PAD'))
                [caaData,dataobject,Data,Data_units]=c_caa_var_get(['PAD_Electron_Dif_flux__',products{k}]);

                % Set up grid
                theta_plus=Data.dep_x{2}.df.plus;
                theta_minus=Data.dep_x{2}.df.minus;
                theta=Data.dep_x{2}.data(1,:);
                nan_theta=isnan(theta);theta(nan_theta)=[]; % for data
                theta=[theta-theta_minus theta(end)+theta_plus(end)]; % grid

                en_plus=Data.dep_x{1}.df.plus;
                en_minus=Data.dep_x{1}.df.minus;
                en=Data.dep_x{1}.data(1,:);
                nan_en=isnan(en);en(nan_en)=[]; % for data
                en_plus(nan_en)=[];en_minus(nan_en)=[];                    
                en=[en+en_plus en(end)-en_plus(1)];

                dataraw=Data.data; dataraw(:,nan_en,:)=[];
                
                % Pick out defined time interval
                dtsampling=(Data.t(2)-Data.t(1))/2;
                tmin=Data.t-dtsampling;
                tmax=Data.t+dtsampling;
                tind=find(tmin<tint(2));
                tind=find(tmax(tind)>tint(1));

                specrec.f=theta;
                specrec.en=en; % stigande
                specrec.t=Data.t(tind);
                specrec.p_label{k}=['Log ' distr{k} ' [' Data_units ']'];
                specrec.f_label='Pitch angle';
                specrec.data=permute(double(dataraw(tind,:,:)),[1 3 2]); % order: time, pitch angle, energy
                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN;                  
            end
            if any([strfind(products{k},'3DXPH'),...
                    strfind(products{k},'CODIF_HS'),...
                    strfind(products{k},'CODIF_LS'),...
                    strfind(products{k},'HIA_HS'),...
                    strfind(products{k},'HIA_LS')])
                % Loads data differently depending on product
                if strfind(products{k},'3DXPH') % PEACE 
                    res=c_caa_construct_subspin_res_data(['Data__', products{k}]);
                    [caaSEDL,~,SEDL]=c_caa_var_get(['Sweep_Energy_DeltaLower__', products{k}]);
                    [caaSEDU,~,SEDU]=c_caa_var_get(['Sweep_Energy_DeltaUpper__', products{k}]);
                    SEDL=flipdim(SEDL(1,2:end),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
                    SEDU=flipdim(SEDU(1,2:end),2); en_nan=isnan(SEDU);SEDU(en_nan)=[];                    
                else % CIS products
                    products{k}(strfind(products{k},'-'))='_';
                    res=c_caa_construct_subspin_res_data(['x3d_ions__', products{k}]);
                    caaSEDL=c_caa_var_get(['delta_plus_energy_table__', products{k}]);
                    caaSEDU=c_caa_var_get(['delta_minus_energy_table__', products{k}]);                            
                    SEDL=flipdim(caaSEDL.data(1,:),2); en_nan=isnan(SEDL);SEDL(en_nan)=[];
                    SEDU=flipdim(caaSEDU.data(1,:),2); en_nan=isnan(SEDU);SEDU(en_nan)=[]; 
                end
            
                [~,ind]=irf_tlim(res.tt,tint);
                specrec.t=res.tt;
                specrec.dt=res.dtsampling/2;
                dtheta=(res.theta(2)-res.theta(1))/2;
                specrec.f=[res.theta-dtheta res.theta(end)+dtheta];
                specrec.f_label='Pitch angle';
                specrec.p=res.pitch_angle(ind,:);
                specrec.en=res.en(:);
                specrec.en=[res.en-SEDL res.en(end)+SEDU(end)];
                specrec.data=res.data(ind,:,:);
                specrec.p_label{k}={['Log ' distr{k} ' [' res.dataunits ']']};

                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN; 
                
            end
            % Save data in cell array 
            data_save{k}=specrec.data;
            % Set up r
            r=double(specrec.en)'; % In eV
            rlog{k}=log10(r); % log eV
            % Pitch angles, turn so that pitch angle 0 is on top
            theta_save{k} = double(specrec.f)+90;             
        end
        if n_products==1
            rlog{2}=rlog{1};
        end
        % Take out r0
        if isempty(emin) || log10(emin)>min([min(rlog{1}) min(rlog{2})])
            r0 = min([min(rlog{1}) min(rlog{2})]);
        else
            r0 = log10(emin);
        end        
        % Create surf grids    
        for k=1:n_products
            r_man{k}=rlog{k}-r0;
            X{k} = r_man{k}*cosd(theta_save{k});
            Y{k} = r_man{k}*sind(theta_save{k});
            C{k} = squeeze(nanmean(log10(data_save{k}),1))';
        end
        if n_products==1
            r_man{2}=r_man{1};
        end
        if n_products==1 % Duplicate matrices
            X{2}=X{1};Y{2}=Y{1};C{2}=C{1};
        end            
        if 1 % Plot data
            surf(ax,X{1},Y{1},X{1}*0,C{1}); hold(ax,'on');
            surf(ax,-flipdim(X{2},2),Y{2},X{2}*0,C{2});         
            view(ax,2); axis(ax,'equal','tight'); shading(ax,'flat'); grid(ax,'off');
            cb=colorbar('peer',ax); 
        end
        if 1 % Title and labels
            t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');     
            if n_products==1 
                str_product=products{1};
                str_product(strfind(str_product,'_'))=' ';
                title_str={[t1str,'-',t2str,' UT'],str_product};
                ylabel(cb,specrec.p_label{1})
            elseif n_products==2
                str_product1=products{1};
                str_product1(strfind(str_product1,'_'))=' ';
                str_product2=products{2};
                str_product2(strfind(str_product2,'_'))=' ';
                title_str={[t1str,'-',t2str,' UT'],...
                    ['Left: ',str_product1],['Right: ', str_product2]};
                ylabel(cb,[char(specrec.p_label{1}), '  /  ', char(specrec.p_label{2})])
            end
            title(ax,title_str)      
        end
        if 1 % Ticks, something not quite right with these
            xticks=log10([1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]/1000)-r0;
            %x_ind=find(xticks>0);
            xticks=[xticks(find(xticks>0))];
            xticks=[xticks(find(xticks<max([max(r_man{1}) max(r_man{2})])))];
            xticklabels=cell(size(xticks));
            for k=1:length(xticklabels)
                xticklabels{k}=num2str(10.^(xticks(k)+r0)/1000);
            end         
            xticks=[-flipdim(xticks,2) xticks];
            xticklabels=[flipdim(xticklabels,2) xticklabels];
            yticks=xticks;
            yticklabels=xticklabels; 
            set(gca,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
                'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)  
            xlabel(ax,'Energy  [keV]'); ylabel(ax,'Energy  [keV]')
        end
        if 1 % Pitch angle labels
            rmax=max([max(r_man{1}) max(r_man{2})]);
            text(0-0.2,rmax-0.5,0,'0^o')
            text(0-0.2,-rmax+0.5,0,'180^o')
            text(-0.2-rmax+0.5,0,0,'90^o')
            text(-0.2+rmax-0.5,0,0,'90^o')
        end        
    case 'cross-section' % Plot cross-section distribution
         % Only one product that will be plotted on both sides          
        product=products{1};
        % For colorbar label
        distr=product(max(strfind(product,'_'))+1:end);
        % Prepare data
        if any(strfind(product,'PITCH')) % Something fishy with background level
            disp('Something fishy with background level...')
            % Load pitch angle data
            [caaData,dataobject,Data,Data_units]=c_caa_var_get(['Data__',product]);
            [caabg,~,bg,bg_units]=c_caa_var_get(['BackgroundLevel__',product]);

            % Pitch angles
            theta=Data.dep_x{2}.data(1,:);
            nan_theta=isnan(theta);theta(nan_theta)=[];
            % Energy levels
            en=Data.dep_x{3}.data(1,:);
            nan_en=isnan(en);en(nan_en)=[];
            % Sub spin angles
            phi=Data.dep_x{1}.data(1,:);
            nan_phi=isnan(phi);phi(nan_phi)=[];

            % Construct subspin reoslution data
            dataraw=Data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
            dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
            tt=subspintime(dataobject,phi); % Define subspin time from angle phi.
            [~,ind]=irf_tlim(tt,tint);
            % Same for background level
            bgraw=bg.data; bgraw(:,:,:,nan_en)=[];bgraw(:,nan_phi,:,:)=[];
            bgraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
            bg=reshape(bgraw,size(bgraw,1)*size(bgraw,2),size(bgraw,3),size(bgraw,4));

            specrec.f=theta;
            %specrec.en=flipdim(en,2);
            [specrec.en en_order]=sort(en);
            specrec.t=tt;
            specrec.p_label=['Log ' distr ' [' Data_units ']'];
            specrec.f_label='Pitch angle';
            specrec.bg=squeeze(nanmean(double(bg(ind,:,en_order)),1)); % average over time
            specrec.data=squeeze(nanmean(double(data(ind,:,en_order)),1)); % average over time
            %specrec.bg=specrec.bg(:,en_order);
            %specrec.data=specrec.data(:,en_order);
        end
        if any([strfind(product,'3DXPH')])
            % Load pitch angle data
            res=c_caa_construct_subspin_res_data(['Data__', product]);
            
            bg=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);

            [~,ind]=irf_tlim(res.tt,tint);
            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            specrec.en=res.en(:);
            specrec.data=res.data(ind,:,:);
            specrec.bg=bg.data(ind,:,:);
            specrec.bg_en=bg.en(:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   

            specrec.data=squeeze(nanmean(specrec.data,1));
            specrec.bg=squeeze(nanmean(specrec.bg,1));
            energy=specrec.en;
        end
        if any([strfind(product,'CODIF_HS'),...
                strfind(product,'CODIF_LS'),...
                strfind(product,'HIA_HS_MAG')])
            % Load pitch angle data
            res=c_caa_construct_subspin_res_data(['Data__', product]);
            
            % No background level data?
            [~,ind]=irf_tlim(res.tt,tint);
            specrec.t=res.tt;
            specrec.dt=res.dtsampling/2;
            specrec.f=res.theta;
            specrec.f_label='Pitch angle';
            specrec.p=res.pitch_angle(ind,:);
            specrec.en=res.en(:);
            specrec.data=res.data(ind,:,:);
            specrec.bg_en=bg.en(:);
            specrec.p_label=['Log ' distr ' [' res.dataunits ']'];                   

            specrec.data=squeeze(nanmean(specrec.data,1));
            specrec.bg=squeeze(nanmean(specrec.bg,1));
            specrec.bg=specrec.data;
            specrec.bg(:)=NaN;
            energy=specrec.en;
        end
        if any(strfind(product,'PAD'))
                [caaData,dataobject,Data,Data_units]=c_caa_var_get(['PAD_Electron_Dif_flux__',product]);

                theta=Data.dep_x{2}.data(1,:);
                nan_theta=isnan(theta);theta(nan_theta)=[]; % for data
                
                en=Data.dep_x{1}.data(1,:);
                nan_en=isnan(en);en(nan_en)=[]; % for data
                
                % Construct subspin reoslution data
                dataraw=Data.data; dataraw(:,nan_en,:)=[];
                
                % Pick out defined time interval
                dtsampling=(Data.t(2)-Data.t(1))/2;
                tmin=Data.t-dtsampling;
                tmax=Data.t+dtsampling;
                tind=find(tmin<tint(2));
                tind=find(tmax(tind)>tint(1));

                specrec.f=theta;
                specrec.en=en'; % stigande
                specrec.t=Data.t(tind);
                specrec.p_label=['Log ' distr ' [' Data_units ']'];
                specrec.f_label='Pitch angle';
                specrec.data=squeeze(nanmean(permute(double(dataraw(tind,:,:)),[1 3 2]))); % order: time, pitch angle, energy
                data_0=find(specrec.data==0);
                specrec.data(data_0)=NaN;  
                specrec.bg=specrec.data;
                specrec.bg(:)=NaN;
            end
        if 1 % Select data to plot
            PA0=nanmean(specrec.data(1,:),1)';
            if mod(length(specrec.f),2)
                PA90=nanmean(specrec.data(find(length(specrec.f)),:),1)';
            else
                PA90=nanmean(specrec.data([length(specrec.f)/2 length(specrec.f)/2+1],:),1)';
            end
                PA180=specrec.data(end,:)';
            PAbg=nanmean(specrec.bg(:,:),1)'; % One common zero-count level for all levels
        end
        if 1 % Plotting data
            loglog(ax,specrec.en,PA0,specrec.en,PA90,specrec.en,PA180,specrec.en,PAbg,'--');
            %set(ax,'xlim',[specrec.en(1) specrec.en(end)])
            t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
            t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
            %prod_str{1}=[t1str,'-',t2str,'UT'];
            %prod_str{2}=product;
            prodstr=[t1str,'-',t2str,' UT    ',product];
            prodstr(strfind(prodstr,'_'))=' ';
            irf_legend(ax,{'0^o','90^o','180^o','-- Zero count'},[0.94 0.94]) 
            ylabel(ax,specrec.p_label)
            xlabel(ax,'Energy  [eV]')
            grid(ax,'off');
            title(ax,prodstr)
         end
end
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