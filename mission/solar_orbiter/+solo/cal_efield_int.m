function [DCE_SRF, params, E_vxb] = cal_efield_int(tt,win,plot_cali,output_res,vxb_dat)
% cali_efield_int calibrate SolO E-field for an specific interval, used to
% check e-field calibration perfomrance and also test new calibration
% methods.
%
% [DCE_SRF, params, E_vxb] = cali_efield_int(tt,win,plot_cali,output_res,vxb_dat)
%
% INPUTS:
%   tt is the center of window to calibrate in epochtt
%   win is the total length of the window in hours for calibration, if empty then a 3 hour window will be used with tt as the center
%   plot_cali = 1 will plot the calibration result
%   output_res = 1 (output full resolution), = 0 (use 1 second)
%   vxb_dat = 1 (L2), 0 (LL), 3 (L2 if available else use LL), 4 (don't compute vxb)
%
% OUTPUTS
%   DCE_SRF is the calibrated E-field in SRF
%   params is a matrix of the cailbataion parameters and comparison with E_vxb
%     row 1 = [k23 d23];
%     row 2 = [k123 d123];
%     row 3 = [cc_y,cc_z]; correlation coefficients for y and z efield with vxb
%   E_vxb is the convecitve electric field if PAS and MAG L2 or LL are available
%
% Example
% [DCE_SRF, params, E_vxb] = solo.cal_efield_int(irf_time('2022-07-25T06:00:00.000Z','utc>epochtt'),6,1,1,1)

if length(tt)==1
    if isempty(win)
        win = 3;
        irf.log('warning','Calibration window set to default 6 hours')
    end
    tint = tt+[-0.5*win*60*60 0.5*win*60*60];
else
    tint = tt;
    if ~isa(tint,'GenericTimeArray')
        error('TINT must be of GenericTimeArray type');
    elseif tint.stop-tint.start<=0
        error('TINT duration is zero or negative');
    end
end

vdc = solo.get_data('vdc_1sec',tint);
QUAL = solo.get_data('vdc_1sec_qual',tint);

if plot_cali
    psp = solo.get_data('scpot',tint);
else
    psp=[];
end

if ~isempty(vdc)

    if output_res
        vdc2 = solo.get_data('vdc',tint);
        QUAL2 = solo.get_data('vdc_qual',tint);
    else
        vdc2 = vdc;
        QUAL2 = QUAL;
    end

    iiBad = QUAL.data <2; vdc.data(iiBad,:) = NaN;
    iiBad2 = QUAL2.data <2; vdc2.data(iiBad2,:) = NaN;

    v3 = double(vdc.z.data); % V3
    idxNan = ~isnan(v3);
    v1 = double(vdc.x.data(idxNan)); % V1
    v2 = double(vdc.y.data(idxNan)); % V2
    v3 = v3(idxNan); % V3 with no nans

    if length(v2)>100

        %% calibration parameters - performed on the 1 second data
        [p23,~]=polyfit(v2,v3,1);
        k23  = p23(1);
        d23 = p23(2);

        v2_scaled = double(v2).*k23 +double(d23); %Remove potential offset between 2,3
        v23 = (v2_scaled+v3)/2; % Corresponding to a measurement point between the two antennas.
        %v23 = (v2+v3)/2; % Corresponding to a measurement point between the two antennas.

        [k123,d123] = lsqfitgm(v23,v1);

        %% Do the calibration on full resolution data
        V2_scaled = vdc2.y.data.*k23 + double(d23); %Remove potential offset between 2,3
        V_d23 = V2_scaled-vdc2.z.data; %Fixed V2-V3.
        Ey_SRF = -V_d23*1e3/6.99;

        V23 = (V2_scaled+vdc2.z.data)/2; % Corresponding to a measurement point between the two antennas.

        Ez_SRF = (V23.*k123+ d123) - vdc2.x.data; % Scale V23.
        Ez_SRF = Ez_SRF*1e3/11.2;

        corr_common_mode=corrcoef(v2_scaled-v3,v2_scaled+v3);
        corr_common_mode=abs(corr_common_mode(1,2));

        if corr_common_mode>0.7
            warning('Common mode detected')
        end

        DCE_SRF = irf.ts_vec_xyz(vdc2.time,[Ey_SRF*0 Ey_SRF Ez_SRF]);
        DCE_SRF.units = 'mV/m';
        DCE_SRF.coordinateSystem = 'SRF';

        params = [k23 d23 ; k123 d123];


if vxb_dat~=4


        if vxb_dat == 1 || vxb_dat == 3
            B = solo.get_data('b_srf_norm',tint);
        else
            B = [];
        end
        if isempty(B) % use LL data if L2 not available
            B = solo.get_data('LL_B_SRF',tint);
            if vxb_dat ~=0
                irf.log('warning','MAG data not found, usinf LL')
            end
            if isempty(B)
                irf.log('warning','MAG LL data not found')
            end
        end

        if vxb_dat == 1 || vxb_dat == 3
            V = solo.get_data('Vi_srf',tint);
        else
            V = [];
        end
        if isempty(V) % use LL data if L2 not available
            V = solo.get_data('LL_V_SRF',tint);
            if vxb_dat ~=0
                irf.log('warning','PAS data not found, usinf LL')
            end
            if isempty(V)
                irf.log('warning','PAS LL data not found')
            end
        end

        if ~isempty(B) && ~isempty(V)
            bvint = intersect(B.resample(V).time.epoch,V.time.epoch);
            if isempty(bvint)
                check_ol = 0;
            else
                check_ol = 1;
            end
        else
            check_ol = 0;
        end

        if check_ol
            B = B.resample(V);
            E_vxb = irf_e_vxb(V,B); % Compute convective E-field
            if isfield(B.userData,'GlobalAttributes')
                B_str = cell2mat(B.userData.GlobalAttributes.Data_type);
                B_str = B_str(1:2);
            else
                B_str = 'LL';
            end

            if isfield(V.userData,'GlobalAttributes')
                V_str = cell2mat(V.userData.GlobalAttributes.Data_type);
                V_str = V_str(1:2);
            else
                V_str = 'LL';
            end

            if ~isempty(intersect(DCE_SRF.resample(E_vxb).time.epoch,E_vxb.time.epoch))
                corr_y=corrcoef(fillmissing(DCE_SRF.resample(E_vxb).y.data,'nearest'), fillmissing(E_vxb.y.data,'nearest'));
                cc_y=corr_y(1,2);


                corr_z=corrcoef(fillmissing(DCE_SRF.resample(E_vxb).z.data,'nearest'), fillmissing(E_vxb.z.data,'nearest'));
                cc_z=corr_z(1,2);


                params = [params;cc_y,cc_z];
            else
                E_vxb = irf.ts_vec_xyz(DCE_SRF.time,DCE_SRF.data*0);
                B_str = 'n/a';
                V_str = 'n/a';
                cc_y = nan;
                cc_z = nan;

                params = [params;cc_y,cc_z];
            end


        else
            E_vxb = irf.ts_vec_xyz(DCE_SRF.time,DCE_SRF.data*0);
            B_str = 'n/a';
            V_str = 'n/a';
            cc_y = nan;
            cc_z = nan;

            params = [params;cc_y,cc_z];
        end
else
    params = [params;nan nan];
    plot_cali=0;
    E_vxb=[];
end

        if plot_cali

            % this is used to find the y-axis range for the plot since auto is not
            % ideal due to outliers/spikes.
            if range(E_vxb.y.data)~=0
                dat_y = rmoutliers([DCE_SRF.y.data;E_vxb.y.data]);
                dat_z = rmoutliers([DCE_SRF.z.data;E_vxb.z.data]);
            else
                dat_y = rmoutliers(DCE_SRF.y.data);
                dat_z = rmoutliers(DCE_SRF.z.data);
            end
            scal_val = 2;

            if ~plot_cali
                figure('Visible','off')
            end

            subplot(6,2,1)
            plot(v2,v3,'kx')
            hold
            xlabel('v2 [V]')
            ylabel('v3 [V]')
            axis equal
            set(gca,'linewidth',1.8,'fontsize',20)
            plot(v2,v2*k23 + d23,'r','linewidth',2)
            title(['k23 = ' num2str(k23) ', d23 = ' num2str(d23)]);
            leg_11 = legend('VDC','cal');
            leg_11.Box = 0;
            irf_legend(gca,'(a)',[0.01 0.98],'color','k','fontsize',25)
            grid

            subplot(6,2,2)
            plot(v23,v1,'kx')
            hold on
            xlabel('v2_{scaled}3 [V]')
            ylabel('v1 [V]')
            axis equal
            set(gca,'linewidth',1.8,'fontsize',20)
            plot(v23,v23*k123 + d123,'r','linewidth',2)
            title(['k123 = ' num2str(k123) ', d23 = ' num2str(d123)]);
            irf_legend(gca,'(b)',[0.01 0.98],'color','k','fontsize',25)
            grid

            hca = subplot(6,2,[3,4]);
            irf_plot(hca,DCE_SRF.y,'color','k','linewidth',1.5)
            hold(hca)
            irf_plot(hca,E_vxb.y,'color','r','linewidth',1.5)
            set(gca,'linewidth',1.8,'fontsize',20)
            ylabel('EDC-Y [mV/m]')
            ylim(hca,'auto')
            leg1 = legend('E_{RPW}',['E_{-VxB} (V-' V_str ' B-' B_str ')']);
            irf_zoom(hca,'x',tint)
            leg1.Box = 0;
            leg1.Orientation = "horizontal";
            ylim(hca,[min(dat_y)-abs(min(dat_y)*scal_val) max(dat_y)+abs(min(dat_y)*scal_val)])
            leg1.Position = [0.3768 0.7827 0.2652 0.0257];
            xlabel(hca,'')
            irf_legend(gca,'(c)',[0.01 0.98],'color','k','fontsize',25)
            irf_legend(gca,['corr coef = ' num2str(round(cc_y,2))],[0.5 0.98],'color','k','fontsize',25)

            hca = subplot(6,2,[5,6]);
            irf_plot(hca,DCE_SRF.z,'color','k','linewidth',1.5)
            hold(hca)
            irf_plot(hca,E_vxb.z,'color','r','linewidth',1.5)
            set(gca,'linewidth',1.8,'fontsize',20)
            ylabel('EDC-Z [mV/m]')
            ylim(hca,'auto')
            irf_zoom(hca,'x',tint);
            ylim(hca,[min(dat_z)-abs(min(dat_z)*scal_val) max(dat_z)+abs(min(dat_z)*scal_val)])
            irf_legend(gca,'(d)',[0.01 0.98],'color','k','fontsize',25)
            irf_legend(gca,['corr coef = ' num2str(round(cc_z,2))],[0.5 0.98],'color','k','fontsize',25)
            xlabel('')

            if ~isempty(B)
                hca = subplot(6,2,[7,8]);
                irf_plot(hca,B,'linewidth',1.5)
                set(hca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]; [0 1 0]]);
                set(gca,'linewidth',1.8,'fontsize',20)
                ylabel('B [nT]')
                irf_zoom(hca,'x',tint);
                leg3 = legend('Bx','By','Bz');
                leg3.Box=0;
                leg3.Position = [0.9138 0.4354 0.0628 0.0525];
                irf_legend(gca,'(e)',[0.01 0.98],'color','k','fontsize',25)
                xlabel('')
            else
                hca = subplot(6,2,[7,8]);
                irf_legend(hca,'(e)',[0.01 0.98],'color','k','fontsize',25)
                irf_legend(hca,'(No data)',[0.5 0.98],'color','k','fontsize',25)
            end


            if ~isempty(V)
                hca = subplot(6,2,[9,10]);
                yyaxis left
                irf_plot(hca,V.x,'color','k','linewidth',1.5)
                ylabel('Vx [km/s]')
                hca.YColor = [0 0 0];
                yyaxis right
                irf_plot(hca,V.y,'color','b','linewidth',1.5)
                hold(hca)
                irf_plot(hca,V.z,'color','r','linewidth',1.5)
                set(gca,'linewidth',1.8,'fontsize',20)
                irf_zoom(hca,'x',tint);
                ylim(hca,'auto')
                hca.YColor = [0 0 0];
                ylabel('Vyz [km/s]')
                leg4 = legend('Vx','Vy','Vz');
                leg4.Box=0;
                leg4.Orientation="horizontal";
                leg4.Position=[0.1416 0.2128 0.1830 0.0192];
                irf_legend(gca,'(f)',[0.01 0.98],'color','k','fontsize',25)
                xlabel('')
            else
                hca = subplot(6,2,[9,10]);
                irf_legend(hca,'(f)',[0.01 0.98],'color','k','fontsize',25)
                irf_legend(hca,'(No data)',[0.5 0.98],'color','k','fontsize',25)
                set(gca,'linewidth',1.8,'fontsize',20)
            end


            hca = subplot(6,2,[11,12]);
            irf_plot(hca,psp,'color','k','linewidth',1.5)
            ylabel(hca,'Av. PSP [V]')
            irf_zoom(hca,'x',tint);
            irf_legend(gca,'(g)',[0.01 0.98],'color','k','fontsize',25)
            set(gca,'linewidth',1.8,'fontsize',20)


            set(gcf,'color','w','position',[399 1 1237 1325])

        end

    else
        DCE_SRF = [];
        params = [];
        E_vxb = [];
        irf.log('critical','Not enough points for calibration')
    end

else
    DCE_SRF = [];
    params = [];
    E_vxb = [];
    irf.log('warning','VDC data was not found')
end
