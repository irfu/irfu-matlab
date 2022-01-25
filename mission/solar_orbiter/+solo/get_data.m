function res = get_data(varStr,Tint)
% res = solo.get_data(varStr, Tint)
%
% varStr is one of:
%
% MAG:
%   'L2_mag-srf-normal' (alias: B_srf_norm), 'L2_mag_rtn_normal' (alias: B_rtn_norm)
%   'L2_mag-srf-burst' (alias:  B_srf_brst), 'L2_mag-rtn-burst' (alias: B_rtn_brst)
%
% RPW:
%   'L3_rpw-bia-scpot' (alias: scpot)
%   'L3_rpw-bia-efield_srf' (alias: E_srf)
%   'L3_rpw-bia-efield_rtn' (alias: E_rtn)
%   'L3_rpw-bia-density' (alias: Nescpot)
%   'L2_rpw-lfr-surv-cwf-b-cdag_srf' (Search coil) (alias: B_scm_srf)
%   'L2_rpw-lfr-surv-cwf-b-cdag_rtn' (alias: B_scm_rtn)
%    Snapshots and other products need to be added!
%
% SWA-PAS:
%   'L2_swa-pas-eflux (alias: pas_eflux)',
%   'L2_swa-pas-grnd-mom_V_RTN' (alias: Vi_rtn)
%   'L2_swa-pas-grnd-mom_V_SRF' (alias: Vi_srf), 'L2_swa-pas-grnd-mom_N' (alias: Ni)
%   'L2_swa-pas-grnd-mom_T' (total T) (alias: Ti),
%   'L2_swa-pas-grnd-mom_TxTyTz_SRF' (alias: Ti_xyz_srf)
%   'L2_swa-pas-grnd-mom_TxTyTz_RTN' (alias: Ti_xyz_rtn)
%   'L2_swa-pas-grnd-mom_Tani' ([Tpar, Tperp1 Tperp2]) (alias: Ti_fac)
%   'L2_swa-pas-grnd-mom_Pi_SRF' (pressure tensor) (alias: Pi_srf),
%   'L2_swa-pas-grnd-mom_Pi_RTN' (alias: Pi_rtn)
%   'L2_swa-pas-vdf' (alias: pas_vdf) 'L2_swa-pas-quality_factor' (alias: pas_qf)
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
vars = {'L2_mag-srf-normal','L2_mag-rtn-normal', 'L2_mag-srf-burst', 'L2_mag-rtn-burst', ...
  'L3_rpw-bia-scpot', 'L3_rpw-bia-efield_srf', 'L3_rpw-bia-efield_rtn', 'L3_rpw-bia-scpot', 'L3_rpw-bia-density', ...
  'L2_rpw-lfr-surv-cwf-b-cdag_srf', 'L2_rpw-lfr-surv-cwf-b-cdag_rtn', ...
  'L2_swa-pas-eflux', 'L2_swa-pas-grnd-mom_V_RTN', 'L2_swa-pas-grnd-mom_V_SRF', 'L2_swa-pas-grnd-mom_N', ...
  'L2_swa-pas-grnd-mom_T', 'L2_swa-pas-grnd-mom_TxTyTz_SRF', 'L2_swa-pas-grnd-mom_TxTyTz_RTN', ...
  'L2_swa-pas-grnd-mom_Tani','L2_swa-pas-grnd-mom_P_SRF', 'L2_swa-pas-grnd-mom_P_RTN', 'L2_swa-pas-vdf', ...
  'pos_rtn','L2_swa-pas-quality_factor'};

%%% check if alias is used and change to full variable name
if isempty(intersect(varStr,vars))
  var_old = varStr;

  % replace alias with the full variable name
  if strcmpi(varStr,'B_rtn_brst');   varStr = 'L2_mag-rtn-burst'; end
  if strcmpi(varStr,'B_rtn_norm');   varStr = 'L2_mag-rtn-normal'; end
  if strcmpi(varStr,'B_srf_brst');   varStr = 'L2_mag-srf-burst'; end
  if strcmpi(varStr,'B_srf_norm');   varStr = 'L2_mag-srf-normal'; end
  if strcmpi(varStr,'Vi_rtn');       varStr = 'L2_swa-pas-grnd-mom_V_RTN'; end
  if strcmpi(varStr,'Vi_srf');       varStr = 'L2_swa-pas-grnd-mom_V_SRF'; end
  if strcmpi(varStr,'Ni');           varStr = 'L2_swa-pas-grnd-mom_N'; end
  if strcmpi(varStr,'Ti');           varStr = 'L2_swa-pas-grnd-mom_T'; end
  if strcmpi(varStr,'pas_vdf');      varStr = 'L2_swa-pas-vdf'; end
  if strcmpi(varStr,'pas_qf');       varStr = 'L2_swa-pas-quality_factor'; end
  if strcmpi(varStr,'pas_eflux');    varStr = 'L2_swa-pas-eflux'; end
  if strcmpi(varStr,'Ti_xyz_srf');   varStr = 'L2_swa-pas-grnd-mom_TxTyTz_SRF'; end
  if strcmpi(varStr,'Ti_xyz_rtn');   varStr = 'L2_swa-pas-grnd-mom_TxTyTz_RTN'; end
  if strcmpi(varStr,'Ti_fac');       varStr = 'L2_swa-pas-grnd-mom_Tani'; end
  if strcmpi(varStr,'Pi_srf');       varStr = 'L2_swa-pas-grnd-mom_P_SRF'; end
  if strcmpi(varStr,'Pi_rtn');       varStr = 'L2_swa-pas-grnd-mom_P_RTN'; end
  if strcmpi(varStr,'scpot');        varStr = 'L3_rpw-bia-scpot'; end
  if strcmpi(varStr,'E_srf');        varStr = 'L3_rpw-bia-efield_srf'; end
  if strcmpi(varStr,'E_rtn');        varStr = 'L3_rpw-bia-efield_rtn'; end
  if strcmpi(varStr,'Nescpot');      varStr = 'L3_rpw-bia-density'; end
  if strcmpi(varStr,'B_scm_srf');    varStr = 'L2_rpw-lfr-surv-cwf-b-cdag_srf'; end
  if strcmpi(varStr,'B_scm_rtn');    varStr = 'L2_rpw-lfr-surv-cwf-b-cdag_rtn'; end

  % Sanity check that we have actually changed the variable
  if strcmp(var_old,varStr)
    errStr = ['"varStr":', varStr, ' incorrect alias used.'];
    irf.log('critical', errStr);
    error(errStr);
  else
    % Print what alias has been changed to
    fprintf(['Alias used: ' var_old ' changed to ' varStr '\n'])
  end
end

% Check if the varStr matches the list of acceptable variables
if ~ismember(varStr, vars)
  errStr = ['"varStr":', varStr, ' is not a recognized paramter.'];
  irf.log('critical', errStr);
  error(errStr);
end

C = strsplit(varStr,'_');

if strcmp(varStr(1),'L') % check if request L2/3 data
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
          res = solo.db_get_ts('solo_L3_rpw-bia-density', 'DENSITY', Tint);
        case 'efield'
          % E-field
          EDC_SRF = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'EDC_SRF', Tint);
          if strcmp(C{3},'rtn')
            EDC_RTN = EDC_SRF; EDC_RTN.data(:,1:2) = -EDC_RTN.data(:,1:2);
            EDC_RTN.name = 'EDC_RTN';
            res=EDC_RTN;
          else
            res=EDC_SRF;
          end
        case 'surv' % we might need to change this case when we add more RPW products!!!
          % search-coil
          BSCM = solo.db_get_ts(['solo_', C{1}, '_', C{2}], 'B_RTN', Tint);
          if strcmp(C{3},'srf')
            res = solo.srf2rtn(BSCM, -1);
            res.name = 'B_SRF';
          else
            res = BSCM;
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
              efulx_file = solo.db_list_files(['solo_',varStr],Tint);
              iEnergy = spdfcdfread([efulx_file(1).path, filesep, efulx_file(1).name],'variables','Energy');
              res = struct('t', ieflux.time.epochUnix);
              res.p = ieflux.data;
              res.p_label='dEF';
              res.f = repmat(iEnergy,1,numel(res.t))';
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
              % PAS ion VDFs: Daniel G. please check!
              t1 = irf_time(Tint(1),'epochtt>vector');
              t2 = irf_time(Tint(2),'epochtt>vector');
              dl = datenum(t1(1:3)):1:datenum(t2(1:3));
              for k=1:length(dl)
                vdf_fname = ['/Volumes/solo/soar/swa/L2/swa-pas-vdf/' datestr(dl(k),'yyyy') '/' datestr(dl(k),'mm') '/solo_L2_swa-pas-vdf_' datestr(dl(k),'yyyy') datestr(dl(k),'mm') datestr(dl(k),'dd') '_V02.cdf'];
                if exist(vdf_fname,'file')==2
                  tmpDataObj = dataobj(vdf_fname);
                  PDout = solo.make_pdist(tmpDataObj);
                  clear tmpDataObj
                  if k==1
                    tlim1 = irf.tint([irf_time(Tint(1),'epochtt>utc') '/' irf_time(PDout.time(end),'epochtt>utc')]);
                    PDout = PDout.tlim(tlim1);
                  end
                  if k==length(dl)
                    tlim2 = irf.tint([irf_time(PDout.time(1),'epochtt>utc') '/' irf_time(Tint(2),'epochtt>utc')]);
                    PDout = PDout.tlim(tlim2);
                  end
                else
                  PDout = [];
                end
                res.(datestr(dl(k),'Tyyyymmdd')) = PDout;
                clear PDout
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
elseif strcmp(C{1},'pos')
  % return the solo position - predicted, should be use flown?
  res = solo.get_position(Tint, 'predicted', 'frame', 'SOLO_SUN_RTN');
end