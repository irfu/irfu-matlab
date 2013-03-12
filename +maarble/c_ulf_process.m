cl_id = 2;
tint=iso2epoch('2010-10-13T12:00:00Z') + [0 3*3600];
freqRange = 'all';
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

%% Download data && Load
if 0
    caa_download(tint+[-30 30],['C' cl_s '_CP_FGM_5VPS_ISR2'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_EFW_L3_E'],'nowildcard');
end
% Construct 1 min-B0 (lowpassed at 1/600 Hz)
B_5VPS = c_caa_var_get(['B_vec_xyz_isr2__C' cl_s '_CP_FGM_5VPS_ISR2'],...
    'mat','tint',tint+[-30 30]);
R = c_caa_var_get(['sc_pos_xyz_isr2__C' cl_s '_CP_FGM_5VPS_ISR2'],...
    'mat','tint',tint+[-30 30]);
if wantPC35
    E_4SEC = c_caa_var_get(['E_Vec_xy_ISR2__C' cl_s '_CP_EFW_L3_E'],...
        'mat','tint',tint+[-30 30]);
end
if wantPC12
    %xxx
end

%% Construct FAC
bf = irf_filt(B_5VPS,0,1/600,1/5,5);
B0_1MIN = irf_resamp(bf,t_1min); clear bf

%% PC3-5
if wantPC35
    t_4SEC = ((tint(1)+2):4:tint(end))';
    E3D_4SEC = irf_edb(irf_resamp(E_4SEC,t_4SEC),B_5VPS,15); % Ez
    B_4SEC = irf_resamp(B_5VPS,t_4SEC);
    
    %XXX: to construct inertial frame
    
    [timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity] = ...
        irf_ebsp(E3D_4SEC,B_4SEC,B0_1MIN,R,'pc35','noresamp');
    BMAG = irf_abs(B0_1MIN); BMAG(:,2:4) = []; BMAG = irf_resamp(BMAG,timeVector);
    h=irf_pl_ebsp(cl_id,R,timeVector,'pc35',BMAG,BB_xxyyzz_fac,...
        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
        k_thphSVD_fac,polSVD_fac,ellipticity);
end