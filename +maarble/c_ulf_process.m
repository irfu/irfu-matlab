freqRange = 'pc35';cl_id = 1; tint=iso2epoch('2010-10-13T12:00:00Z') + [0 3*3600]; % PC3-5 example
%freqRange = 'pc12';cl_id = 1;tint=iso2epoch('2007-01-03T16:00:00Z') + [0 0.5*3600]; % PC1-2 example

outDir = '.';
plotFlag = 1;

wantPC12 = 0;
wantPC35 = 0;
switch lower(freqRange)
    case 'all'
        wantPC12 = 1;
        wantPC35 = 1;
    case 'pc12'
        wantPC12 = 1;
    case 'pc35'
        wantPC35 = 1;
    otherwise
        error('Invalid value for freqRange')
end


% Round time interval to minutes
tint = round(tint/60)*60;
t_1min = (tint(1):60:tint(end))';

cl_s = int2str(cl_id);

%% Download data
if 0
    caa_download(tint+[-30 30],['C' cl_s '_CP_FGM_5VPS_ISR2'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_EFW_L3_E'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_AUX_POSGSE_1M'],'nowildcard');
    caa_download(tint+[-30 30],'CL_SP_AUX','nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_FGM_FULL_ISR2'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_EFW_L2_E'],'nowildcard');
end

%% Load
% Construct 1 min-B0 (lowpassed at 1/600 Hz)
gseB_5VPS = c_caa_var_get(['B_vec_xyz_gse__C' cl_s '_CP_FGM_5VPS'],...
    'mat','tint',tint+[-30 30]);
gseR = c_caa_var_get(['sc_pos_xyz_gse__C' cl_s '_CP_FGM_5VPS'],...
    'mat','tint',tint+[-30 30]);
gseV = c_caa_var_get(['sc_v_xyz_gse__C' cl_s '_CP_AUX_POSGSE_1M'],...
    'mat','tint',tint+[-30 30]);
SAXlat = c_caa_var_get(['sc_at' cl_s '_lat__CL_SP_AUX'],...
    'mat','tint',tint+[-30 30]);
SAXlong = c_caa_var_get(['sc_at' cl_s '_long__CL_SP_AUX'],...
    'mat','tint',tint+[-30 30]);
[xspin,yspin,zspin] = sph2cart(mean(SAXlong(:,2))*pi/180,...
    mean(SAXlat(:,2))*pi/180,1); SAX = [xspin yspin zspin];
R = c_coord_trans('gse','isr2',gseR,'SAX',SAX);
V = c_coord_trans('gse','isr2',gseV,'SAX',SAX);
B_5VPS = c_coord_trans('gse','isr2',gseB_5VPS,'SAX',SAX);
                
if wantPC35
    E_4SEC = c_caa_var_get(['E_Vec_xy_ISR2__C' cl_s '_CP_EFW_L3_E'],...
        'mat','tint',tint+[-30 30]);
end
if wantPC12
    B_FULL = c_caa_var_get(['B_vec_xyz_isr2__C' cl_s '_CP_FGM_FULL_ISR2'],...
        'mat','tint',tint+[-1 1]);
    % XXX
    % TODO: We need to check the bitmask here and not use any data with low
    % quality.
    E_L2 = c_caa_var_get(['E_Vec_xy_ISR2__C' cl_s '_CP_EFW_L2_E'],...
        'mat','tint',tint+[-1 1]);
end

%% Calculate and plot
bf = irf_filt(B_5VPS,0,1/600,1/5,5);
B0_1MIN = irf_resamp(bf,t_1min); clear bf

if wantPC35
    t_4SEC = ((tint(1)+2):4:tint(end))';
    B_4SEC = irf_resamp(B_5VPS,t_4SEC);
    E3D_4SEC = irf_edb(irf_resamp(E_4SEC,t_4SEC),B_4SEC,15,'Eperp+NaN'); % Ez
    
    % Construct the inertial frame
    evxb = irf_tappl(irf_cross(B_4SEC,irf_resamp(V,t_4SEC)),'*1e-3*(-1)');
    iE3D_4SEC = E3D_4SEC;
    iE3D_4SEC(:,2:4) = iE3D_4SEC(:,2:4) - evxb(:,2:4);
    
    tic; [timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,polSVD_fac,planarity,ellipticity,k_thphSVD_fac] = ...
        irf_ebsp(iE3D_4SEC,B_4SEC,[],B0_1MIN,R,'pc35','noresamp','fullB=dB');
    toc
    BMAG = irf_abs(B0_1MIN); BMAG(:,2:4) = []; BMAG = irf_resamp(BMAG,timeVector);
    h=irf_pl_ebsp(cl_id,R,timeVector,'pc35',BMAG,BB_xxyyzz_fac,...
        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
        k_thphSVD_fac,polSVD_fac,ellipticity);
end
if wantPC12
    fFGM = 22.5;
    t_BASE = (fix(min(B_FULL(1,1),E_L2(1,1))):2.0/fFGM:ceil(max(B_FULL(end,1),E_L2(end,1))))';
    B_BASE = irf_resamp(B_FULL,t_BASE);
    E3D_BASE = irf_edb(irf_resamp(E_L2,t_BASE),B_BASE,15,'Eperp+NaN'); % Ez
    
    % Construct the inertial frame
    evxb = irf_tappl(irf_cross(B_BASE,irf_resamp(V,t_BASE)),'*1e-3*(-1)');
    iE3D_BASE = E3D_BASE;
    iE3D_BASE(:,2:4) = iE3D_BASE(:,2:4) - evxb(:,2:4);
    
    [timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,polSVD_fac,planarity,ellipticity,k_thphSVD_fac] = ...
        irf_ebsp(iE3D_BASE,B_BASE,[],B0_1MIN,R,'pc12','noresamp','fullB=dB','dedotb=0');
    BMAG = irf_abs(B0_1MIN); BMAG(:,2:4) = []; BMAG = irf_resamp(BMAG,timeVector);
    h=irf_pl_ebsp(cl_id,R,timeVector,'pc12',BMAG,BB_xxyyzz_fac,...
        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
        k_thphSVD_fac,polSVD_fac,ellipticity);
end