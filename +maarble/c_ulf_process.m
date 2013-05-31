freqRange = 'pc35';cl_id = 1; tint=iso2epoch('2010-10-13T12:00:00Z') + [0 3*3600]; % PC3-5 example
%freqRange = 'pc12';cl_id = 1;tint=iso2epoch('2007-01-03T16:00:00Z') + [0 0.5*3600]; % PC1-2 example
%freqRange = 'pc12';cl_id = 1;tint=iso2epoch('2011-11-01T20:13:00Z') + [0 25*60]; % PC1-2 example
%freqRange = [10 180]; cl_id = 4;tint=iso2epoch('2001-02-26T05:18:00Z') + [0 60]; % VLF example

outDir = '.';
plotFlag = 1;

wantPC12 = 0;
wantPC35 = 0; wantSCM = 0;
if ischar(freqRange)
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
else
    if freqRange(1) > 1, wantSCM = 1; end
end


% Round time interval to minutes
tint = round(tint/60)*60;
t_1min = (tint(1):60:tint(end))';

cl_s = int2str(cl_id);

%% Download data
if 0
    caa_download(tint+[-30 30],['C' cl_s '_CP_FGM_5VPS'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_EFW_L3_E'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_AUX_POSGSE_1M'],'nowildcard');
    caa_download(tint+[-30 30],'CL_SP_AUX','nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_FGM_FULL_ISR2'],'nowildcard');
    caa_download(tint+[-30 30],['C' cl_s '_CP_EFW_L2_E'],'nowildcard');
    if wantSCM
      caa_download(tint+[-30 30],['C' cl_s '_CP_STA_CWF_HBR_ISR2'],'nowildcard');
    end
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
elseif wantPC12
    B_FULL = c_caa_var_get(['B_vec_xyz_isr2__C' cl_s '_CP_FGM_FULL_ISR2'],...
        'mat','tint',tint+[-1 1]);
    % XXX
    % TODO: We need to check the bitmask here and not use any data with low
    % quality.
    E_L2 = c_caa_var_get(['E_Vec_xy_ISR2__C' cl_s '_CP_EFW_L2_E'],...
        'mat','tint',tint+[-1 1]);
else
    E_L2 = c_caa_var_get(['E_Vec_xy_ISR2__C' cl_s '_CP_EFW_L2_E'],...
        'mat','tint',tint+[-1 1]);
    B_FULL = c_caa_var_get(['B_vec_xyz_isr2__C' cl_s '_CP_FGM_FULL_ISR2'],...
        'mat','tint',tint+[-1 1]);
    if wantSCM
        BSC = c_caa_var_get('B_vec_xyz_Instrument__C4_CP_STA_CWF_HBR_ISR2','mat');
    end
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
    
    tic; ebsp = ...
        irf_ebsp(iE3D_4SEC,B_4SEC,[],B0_1MIN,R,'pc35',...
        'fac','polarization','noresamp','fullB=dB');
    toc
elseif wantPC12
    fFGM = 22.5;
    t_BASE = (fix(min(B_FULL(1,1),E_L2(1,1))):2.0/fFGM:ceil(max(B_FULL(end,1),E_L2(end,1))))';
    B_BASE = irf_resamp(B_FULL,t_BASE);
    E3D_BASE = irf_edb(irf_resamp(E_L2,t_BASE),B_BASE,15,'Eperp+NaN'); % Ez
    
    % Construct the inertial frame
    evxb = irf_tappl(irf_cross(B_BASE,irf_resamp(V,t_BASE)),'*1e-3*(-1)');
    iE3D_BASE = E3D_BASE;
    iE3D_BASE(:,2:4) = iE3D_BASE(:,2:4) - evxb(:,2:4);
    
    tic
    ebsp = irf_ebsp(iE3D_BASE,B_BASE,[],B0_1MIN,R,'pc12',...
        'fac','polarization','noresamp','fullB=dB','dedotb=0','nav',12);
    toc
else
  tic
  if wantSCM
    ebsp = irf_ebsp(E_L2,BSC,B_FULL,B_5VPS,R,freqRange,...
      'fac','polarization','dedotb=0');
  else
    ebsp = irf_ebsp(E_L2,B_FULL,[],B_5VPS,R,freqRange,...
      'fac','polarization','dedotb=0','fullB=dB');
  end
  toc
end
h = irf_pl_ebsp(ebsp);
irf_zoom(h,'x',tint)
title(h(1),['Cluster ' cl_s ', ' irf_disp_iso_range(tint,1)])