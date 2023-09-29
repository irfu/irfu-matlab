function res = get_data(varStr,Tint)
% res = solo.get_data(varStr, Tint)
%  
% Read Solar Orbiter data from the IRF server
%
% varStr is one of the following (note aliases can also be used):
%
% MAG:
%   'L2_mag-srf-normal' (alias: b_srf_norm), 'L2_mag_rtn_normal' (alias: b_rtn_norm)
%   'L2_mag-srf-normal-1-minute' (alias: b_srf_norm_1min), 'L2_mag_rtn_normal-1-minute' (alias: b_rtn_norm_1min))
%   'L2_mag-srf-burst' (alias:  b_srf_brst), 'L2_mag-rtn-burst' (alias: b_rtn_brst)
%
% RPW:
%   'L3_rpw-bia-scpot' (alias: scpot)
%   'L3_rpw-bia-efield_srf' (alias: e_srf)
%   'L3_rpw-bia-efield_rtn' (alias: e_rtn)
%   'L3_rpw-bia-efield-10-seconds_srf' (alias: e_srf_10sec)
%   'L3_rpw-bia-efield-10-seconds_rtn' (alias: e_rtn_10sec)
%   'L3_rpw-bia-density' (alias: nescpot)
%   'L3_rpw-bia-density-10-seconds' (alias: nescpot_10sec)
%   'L2_rpw-lfr-surv-cwf-b-cdag_srf' (Search coil) (alias: b_scm_srf)
%   'L2_rpw-lfr-surv-cwf-b-cdag_rtn' (alias: b_scm_rtn)
%   'L2_rpw-lfr-surv-cwf-e-1-second' (alias: vdc_1sec)
%   'L2_rpw-lfr-surv-cwf-e' (alias: vdc)
%   'L2_rpw-lfr-surv-cwf-e-1-second_qual' (alias: vdc_1sec_qual)
%   'L2_rpw-lfr-surv-cwf-e_qual' (alias: vdc_qual)
%    Snapshots and other products need to be added!
%
% SWA-PAS:
%   'L2_swa-pas-eflux (alias: pas_eflux)',
%   'L2_swa-pas-grnd-mom_V_RTN' (alias: vi_rtn)
%   'L2_swa-pas-grnd-mom_V_SRF' (alias: vi_srf), 'L2_swa-pas-grnd-mom_N' (alias: ni)
%   'L2_swa-pas-grnd-mom_T' (total T) (alias: ti),
%   'L2_swa-pas-grnd-mom_TxTyTz_SRF' (alias: ti_xyz_srf)
%   'L2_swa-pas-grnd-mom_TxTyTz_RTN' (alias: ti_xyz_rtn)
%   'L2_swa-pas-grnd-mom_Tani' ([Tpar, Tperp1 Tperp2]) (alias: ti_fac)
%   'L2_swa-pas-grnd-mom_Pi_SRF' (pressure tensor) (alias: pi_srf),
%   'L2_swa-pas-grnd-mom_Pi_RTN' (alias: pi_rtn)
%   'L2_swa-pas-vdf' (alias: pas_vdf) 'L2_swa-pas-quality_factor' (alias: pas_qf)
%
% LOW LATENCY (NOT FOR SCIENCE!)
%   'LL_B_RTN', 'LL_B_SRF', 'LL_V_RTN', 'LL_V_SRF', 'LL_N'
%
% EPHEMERIS
%   'pos_rtn' - from solo.get_position
%
% Example
%   Tint = irf.tint('2020-07-29T10:20:00.000Z/2020-07-29T11:20:00.000Z');
%   V = solo.get_data('L2_swa-pas-grnd-mom_V_RTN',Tint) % Solar wind speed in RTN


Units = irf_units;

if ~isa(Tint,'GenericTimeArray')
    error('TINT must be of GenericTimeArray type');
elseif Tint.stop-Tint.start<=0
    error('TINT duration is zero or negative');
end

% List of full variable names. Alias names not inlcluded but changed to
% proper variable name below.
vars = {'L2_mag-srf-normal','L2_mag-srf-normal-1-minute','L2_mag-rtn-normal','L2_mag-rtn-normal-1-minute', 'L2_mag-srf-burst', 'L2_mag-rtn-burst', ...
    'L3_rpw-bia-scpot', 'L3_rpw-bia-efield_srf', 'L3_rpw-bia-efield_rtn','L3_rpw-bia-efield-10-seconds_srf','L3_rpw-bia-efield_rtn'...
    'L3_rpw-bia-scpot', 'L3_rpw-bia-density', 'L3_rpw-bia-density-10-seconds', 'L2_rpw-lfr-surv-cwf-b-cdag_srf', 'L2_rpw-lfr-surv-cwf-b-cdag_rtn', ...
    'L2_rpw-lfr-surv-cwf-e-1-second', 'L2_swa-pas-eflux', 'L2_swa-pas-grnd-mom_V_RTN', 'L2_swa-pas-grnd-mom_V_SRF', 'L2_swa-pas-grnd-mom_N', ...
    'L2_swa-pas-grnd-mom_T', 'L2_swa-pas-grnd-mom_TxTyTz_SRF', 'L2_swa-pas-grnd-mom_TxTyTz_RTN', 'L2_rpw-lfr-surv-cwf-e','L2_rpw-lfr-surv-cwf-e-1-second_qual',...
    'L2_swa-pas-grnd-mom_Tani','L2_swa-pas-grnd-mom_P_SRF', 'L2_swa-pas-grnd-mom_P_RTN', 'L2_swa-pas-vdf', 'L2_rpw-lfr-surv-cwf-e_qual',...
    'pos_rtn','L2_swa-pas-quality_factor', 'LL_B_RTN', 'LL_B_SRF', 'LL_V_RTN', 'LL_V_SRF', 'LL_N','L2_rpw-tds-surv-stat'};

%% check if alias is used and change to full variable name
if ~ismember(varStr, vars)
    switch lower(varStr) % effectivly ignore letter case
        % short alias and full variable names
        case 'b_rtn_brst',      varStrNew = 'L2_mag-rtn-burst';
        case 'b_rtn_norm',      varStrNew = 'L2_mag-rtn-normal';
        case 'b_rtn_norm_1min', varStrNew = 'L2_mag-rtn-normal-1-minute';
        case 'b_srf_brst',      varStrNew = 'L2_mag-srf-burst';
        case 'b_srf_norm',      varStrNew = 'L2_mag-srf-normal';
        case 'b_srf_norm_1min', varStrNew = 'L2_mag-srf-normal-1-minute';
        case 'vi_rtn',          varStrNew = 'L2_swa-pas-grnd-mom_V_RTN';
        case 'vi_srf',          varStrNew = 'L2_swa-pas-grnd-mom_V_SRF';
        case 'ni',              varStrNew = 'L2_swa-pas-grnd-mom_N';
        case 'ti',              varStrNew = 'L2_swa-pas-grnd-mom_T';
        case 'pas_vdf',         varStrNew = 'L2_swa-pas-vdf';
        case 'pas_qf',          varStrNew = 'L2_swa-pas-quality_factor';
        case 'pas_eflux',       varStrNew = 'L2_swa-pas-eflux';
        case 'ti_xyz_srf',      varStrNew = 'L2_swa-pas-grnd-mom_TxTyTz_SRF';
        case 'ti_xyz_rtn',      varStrNew = 'L2_swa-pas-grnd-mom_TxTyTz_RTN';
        case 'ti_fac',          varStrNew = 'L2_swa-pas-grnd-mom_Tani';
        case 'pi_srf',          varStrNew = 'L2_swa-pas-grnd-mom_P_SRF';
        case 'pi_rtn',          varStrNew = 'L2_swa-pas-grnd-mom_P_RTN';
        case 'scpot',           varStrNew = 'L3_rpw-bia-scpot';
        case 'vdc_1sec',        varStrNew = 'L2_rpw-lfr-surv-cwf-e-1-second';
        case 'vdc',             varStrNew = 'L2_rpw-lfr-surv-cwf-e';
        case 'vdc_1sec_qual',   varStrNew = 'L2_rpw-lfr-surv-cwf-e-1-second_qual';
        case 'vdc_qual',        varStrNew = 'L2_rpw-lfr-surv-cwf-e_qual';
        case 'e_srf',           varStrNew = 'L3_rpw-bia-efield_srf';
        case 'e_rtn',           varStrNew = 'L3_rpw-bia-efield_rtn';
        case 'e_srf_10sec',     varStrNew = 'L3_rpw-bia-efield-10-seconds_srf';
        case 'e_rtn_10sec',     varStrNew = 'L3_rpw-bia-efield-10-seconds_rtn';
        case 'nescpot',         varStrNew = 'L3_rpw-bia-density';
        case 'nescpot_10sec',   varStrNew = 'L3_rpw-bia-density-10-seconds';
        case 'b_scm_srf',       varStrNew = 'L2_rpw-lfr-surv-cwf-b-cdag_srf';
        case 'b_scm_rtn',       varStrNew = 'L2_rpw-lfr-surv-cwf-b-cdag_rtn';
        case 'tds_stat',        varStrNew = 'L2_rpw-tds-surv-stat';
        otherwise
            % fallback, it was not a full variable name nor short alias
            errStr = ['"varStr":', varStr, ' incorrect alias used.'];
            irf.log('critical', errStr);
            error(errStr);
    end
    % Print what alias has been changed to
    irf.log('debug', ['Alias used: ', varStr, ' changed to ', varStrNew]);
    % replace alias with the full variable name
    varStr = varStrNew;
end

% Check if the varStr matches the list of acceptable variables
if ~ismember(varStr, vars)
    errStr = ['"varStr":', varStr, ' is not a recognized paramter.'];
    irf.log('critical', errStr);
    error(errStr);
end

C = strsplit(varStr,'_');

if strcmp(varStr(1),'L') && ~strcmp(varStr(2),'L') % check if request L2/3 data
    switch C{2}(1:3)
        case 'mag'
            % MAG variables
            C2 = strsplit(C{2},'-');
            res = solo.db_get_ts(['solo_', varStr], ['B_', upper(C2{2})], Tint);
        case 'rpw'
            % RPW variables
            C2 = strsplit(C{2},'-');
            switch C2{3}
                case 'scpot'
                    % Probe-spacecraft potential
                    res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'PSP', Tint);
                case 'density'
                    % Electron density from probe-spacecraft potential
                    res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'DENSITY', Tint);
                case 'efield'
                    % E-field
                    EDC_SRF = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'EDC_SRF', Tint);
                    if strcmp(C{3},'rtn') && ~isempty(EDC_SRF)
                        EDC_RTN = EDC_SRF; EDC_RTN.data(:,1:2) = -EDC_RTN.data(:,1:2);
                        EDC_RTN.name = 'EDC_RTN';
                        res=EDC_RTN;
                    else
                        res=EDC_SRF;
                    end
                case 'surv'
                    if strcmp(C2{4},'stat')
                        res = solo.db_get_ts(['solo_', C{1}, '_', C{2},'-cdag'], 'WA_MED_FREQ', Tint);
                    else
                        switch C2{5}
                            case 'b' % search-coil
                                BSCM = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'B_RTN', Tint);
                                if strcmp(C{3},'srf') && ~isempty(BSCM)
                                	res = solo.srf2rtn(BSCM, -1);
                                    res.name = 'B_SRF';
                                else
                                    res = BSCM;
                                end
                            case 'e'  % VDC
                                if length(C)>2
                                    if strcmp(C{3},'qual')
                                        res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'QUALITY_FLAG', Tint);
                                    end
                                else
                                res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'VDC', Tint);
                                end
                        end
                    end
                otherwise
                    errStr = 'Not yet defined';
                    irf.log('critical', errStr);
                    error(errStr);
            end
        case 'swa'
            % SWA-PAS variables
            C2 = strsplit(C{2},'-');
            switch C2{2}
                case 'pas'
                    switch C2{3}
                        case 'eflux'
                            % omni energy flux
                            ieflux = solo.db_get_ts(['solo_',varStr],'eflux',Tint);
                            if ~isempty(ieflux)
                                efulx_file = solo.db_list_files(['solo_',varStr],Tint);
                                iEnergy = spdfcdfread([efulx_file(1).path, filesep, efulx_file(1).name],'variables','Energy');
                                res = struct('t', ieflux.time.epochUnix);
                                res.p = ieflux.data;
                                res.p_label='dEF';
                                res.f = repmat(iEnergy,1,numel(res.t))';
                            else
                                res = [];
                            end
                        case 'grnd'
                            % ground-moments
                            switch C{3}
                                case {'P', 'V'}
                                    res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], [C{3}, '_', C{4}], Tint);
                                case {'N', 'T'}
                                    res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], C{3}, Tint);
                                case 'Tani'
                                    % requires MAG data and also density
                                    PSRF = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','P_SRF',Tint);
                                    if isempty(PSRF)
                                        errStr = 'Pressure tensor not loaded';
                                        irf.log('critical', errStr);
                                        error(errStr);
                                    end
                                    NPAS = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N',Tint);
                                    if isempty(PSRF)
                                        errStr='Density not loaded';
                                        irf.log('critical', errStr);
                                        error(errStr);
                                    end
                                    BSRF = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF',Tint);
                                    if isempty(BSRF)
                                        errStr='MAG data not loaded';
                                        irf.log('critical', errStr);
                                        error(errStr);
                                    end
                                    PXX = PSRF.data(:,1);
                                    PYY = PSRF.data(:,2);
                                    PZZ = PSRF.data(:,3);
                                    PXY = PSRF.data(:,4);
                                    PXZ = PSRF.data(:,5);
                                    PYZ = PSRF.data(:,6);
                                    Pten = TSeries(PSRF.time,[PXX,PXY,PXZ,PYY,PYZ,PZZ]); % to be consistent with mms.rotate_tensor
                                    B0 = irf_filt(BSRF,0,0.1,[],3);
                                    PfacT = mms.rotate_tensor(Pten,'fac',B0.resample(Pten));
                                    Pfac = TSeries(PfacT.time,[PfacT.data(:,1,1), PfacT.data(:,2,2),PfacT.data(:,3,3)]);
                                    res = TSeries(PfacT.time,(Pfac.data./(NPAS.resample(Pfac).data.*Units.kB))./(Units.eV/Units.kB));
                                    res.userData = {'par','perp1','perp2'};
                                case 'TxTyTz'
                                    res = solo.db_get_ts(['solo_', C{1}, '_', C{2}], ['TxTyTz_', C{4}], Tint);
                                otherwise
                                    errStr = 'Not yet defined';
                                    irf.log('critical', errStr);
                                    error(errStr);
                            end
                        case 'quality'
                            res = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','quality_factor',Tint);
                        case 'vdf'
                            vdf_files  = solo.db_list_files(['solo_',varStr],Tint);
                            if ~isempty(vdf_files)
                                for k=1:length(vdf_files)
                                    tmpDataObj = dataobj([vdf_files(k).path, filesep, vdf_files(k).name]);
                                    PDout = solo.make_pdist(tmpDataObj);
                                    res.(irf_time(Tint(k),'epochtt>utc_Tyyyymmdd')) = PDout.tlim(Tint);
                                    clear PDout
                                end
                            else
                                res = [];
                            end
                        otherwise
                            errStr = 'Not yet defined';
                            irf.log('critical', errStr);
                            error(errStr);
                    end
                case 'eas'
                    % Electron data - can be added at a later date
                    warning('Electron data to be added in the future')
                otherwise
                    errStr = 'Not yet defined';
                    irf.log('critical', errStr);
                    error(errStr);
            end
    end
end
if strcmp(C{1},'pos')
    % return the solo position - predicted, should we use flown?
    res = solo.get_position(Tint, 'predicted', 'frame', 'SOLO_SUN_RTN');
end
if strcmp(C{1},'LL')
    res = solo.read_LL(varStr(4:end),Tint);
end
