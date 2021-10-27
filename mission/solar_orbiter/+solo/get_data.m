function res = get_data(varStr,Tint)

% res = solo.get_data(varStr, Tint)
%
% varStr is one of:
%
% MAG:
%   'L2_mag-srf-normal', 'L2_mag-rtn-normal',
%   'L2_mag-srf-burst', 'L2_mag-rtn-burst'
%
% RPW:
%   'L3_rpw-bia-scpot' ,'L3_rpw-bia-efield_srf'
%   'L3_rpw-bia-efield_rtn', 'L3_rpw-bia-scpot'
%   'L3_rpw-bia-density', 'L2_rpw-lfr-surv-cwf-b-cdag_srf' (Search coil)
%   'L2_rpw-lfr-surv-cwf-b-cdag_rtn'
%    Snapshots and other products need to be added!
%
% SWA-PAS:
%   'L2_swa-pas-eflux', 'L2_swa-pas-grnd-mom_V_RTN'
%   'L2_swa-pas-grnd-mom_V_SRF', 'L2_swa-pas-grnd-mom_N'
%   'L2_swa-pas-grnd-mom_T' (total T), 'L2_swa-pas-grnd-mom_TxTyTz_SRF'
%   'L2_swa-pas-grnd-mom_TxTyTz_RTN', 'L2_swa-pas-grnd-mom_Tani' ([Tpar, Tperp])
%   'L2_swa-pas-grnd-mom_P_SRF' (pressure tensor), 'L2_swa-pas-grnd-mom_P_RTN'
%   'L2_swa-pas-vdf'
%
% EPHEMERIS
%   'pos_rtn' - from solo.get_position
%
% Example
%   Tint = irf.tint('2020-07-29T10:20:00.000Z/2020-07-29T11:20:00.000Z');
%   V = solo.get_data('L2_swa-pas-grnd-mom_V_RTN',Tint) % Solar wind speed in RTN
%


if ~isa(Tint,'GenericTimeArray')
  error('TINT must be of GenericTimeArray type');
  elseif Tint.stop-Tint.start<=0
  error('TINT duration is zero or negative');
end


vars = {'L2_mag-srf-normal','L2_mag-rtn-normal','L2_mag-srf-burst','L2_mag-rtn-burst','L3_rpw-bia-scpot',...
  'L3_rpw-bia-efield_srf','L3_rpw-bia-efield_rtn','L3_rpw-bia-scpot','L3_rpw-bia-density','L2_rpw-lfr-surv-cwf-b-cdag_srf','L2_rpw-lfr-surv-cwf-b-cdag_rtn','L2_swa-pas-eflux','L2_swa-pas-grnd-mom_V_RTN',...
  'L2_swa-pas-grnd-mom_V_SRF','L2_swa-pas-grnd-mom_N','L2_swa-pas-grnd-mom_T','L2_swa-pas-grnd-mom_TxTyTz_SRF','L2_swa-pas-grnd-mom_TxTyTz_RTN','L2_swa-pas-grnd-mom_Tani','L2_swa-pas-grnd-mom_P_SRF',...
  'L2_swa-pas-grnd-mom_P_RTN','L2_swa-pas-vdf','pos_rtn'};

% Check if the varStr matches the list of acceptable variables
var_match = 0;
for k=1:length(vars)
  if strcmp(vars{k},varStr)
    var_match=1;
  end
end
if var_match==0;error('variable not valid'); end % error if variable not recognised

C = strsplit(varStr,'_');

if strcmp(varStr(1),'L') % check if request L2/3 data
  
  switch C{2}(1:3)
    
    
    
        
    case 'mag' % MAG variables
      C2 = strsplit(C{2},'-');
      res = solo.db_get_ts(['solo_',varStr],['B_',upper(C2{2})],Tint);
      
      
      
      
    case 'rpw' % RPW variables
      C2 = strsplit(C{2},'-');
      
      % Probe-spacecraft potential
      if strcmp(C2{3},'scpot')
        res = solo.db_get_ts(['solo_',C{1},'_',C{2}], 'PSP', Tint);
      end
      
      % Electron density from probe-spacecraft potential
      if strcmp(C2{3},'density')
        res = solo.db_get_ts('solo_L3_rpw-bia-density','DENSITY',Tint);
      end
      
      % E-field
      if strcmp(C2{3},'efield')
        EDC_SRF = solo.db_get_ts(['solo_',C{1},'_',C{2}],'EDC_SRF', Tint);
        if strcmp(C{3},'rtn')
          EDC_RTN = EDC_SRF; EDC_RTN.data(:,1:2) = -EDC_RTN.data(:,1:2);
          EDC_RTN.name = 'EDC_RTN';
          res=EDC_RTN;
        else
          res=EDC_SRF;
        end
      end
      
      % search-coil
      if strcmp(C{2},'rpw-lfr-surv-cwf-b-cdag')
        BSCM = solo.db_get_ts(['solo_',C{1},'_',C{2}],'B_RTN', Tint);
        if strcmp(C{3},'srf')
          res = solo.srf2rtn(BSCM, -1);
          res.name = 'B_SRF';
        else
          res = BSCM;
        end
      end
      
      
      
      
    case 'swa' % SWA-PAS variables
      C2 = strsplit(C{2},'-');
      if strcmp(C2{2},'pas')
        
        if strcmp(C2{3},'eflux') % omni energy flux
          ieflux = solo.db_get_ts(['solo_',varStr],'eflux',Tint);
          efulx_file=solo.db_list_files(['solo_',varStr],Tint);
          iEnergy = spdfcdfread([efulx_file(1).path '/' efulx_file(1).name],'variables','Energy');
          res   = struct('t', ieflux.time.epochUnix);
          res.p = ieflux.data;
          res.p_label='dEF';
          res.f = repmat(iEnergy,1,numel(res.t))';
        end
        
        if strcmp(C2{3},'grnd') % ground-moments
          
          if strcmp(C{3},'V')
            res = solo.db_get_ts(['solo_',C{1},'_',C{2}],['V_',C{4}],Tint);
          end
          
          if strcmp(C{3},'N')
            res = solo.db_get_ts(['solo_',C{1},'_',C{2}],'N',Tint);
          end
          
          if strcmp(C{3},'T')
            res = solo.db_get_ts(['solo_',C{1},'_',C{2}],'T',Tint);
          end
          
          if strcmp(C{3},'P')
            res = solo.db_get_ts(['solo_',C{1},'_',C{2}],['P_',C{4}],Tint);
          end
          
          if strcmp(C{3},'Tani') % requires MAG data and also density
            PSRF = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','P_SRF',Tint);
            if isempty(PSRF); error('Pressure tensor not loaded'); end
            NPAS = solo.db_get_ts('solo_L2_swa-pas-grnd-mom','N',Tint);
            if isempty(PSRF); error('Density not loaded'); end
            BSRF = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF',Tint);
            if isempty(BSRF); error('MAG data not loaded'); end
            PXX = PSRF.data(:,1);
            PYY = PSRF.data(:,2);
            PZZ = PSRF.data(:,3);
            PXY = PSRF.data(:,4);
            PXZ = PSRF.data(:,5);
            PYZ = PSRF.data(:,6);
            Pten = TSeries(PSRF.time,[PXX,PXY,PXZ,PYY,PYZ,PZZ]); % to be consistent with mms.rotate_tensor
            B0 = irf_filt(BSRF,0,0.1,[],3);
            PfacT = mms.rotate_tensor(Pten,'fac',B0.resample(Pten),'pp');
            Pfac = TSeries(PfacT.time,[PfacT.data(:,1,1), PfacT.data(:,2,2)]);
            res = TSeries(PfacT.time,(Pfac.data./(NPAS.resample(Pfac).data.*(1.38064852*10^(-23))))./1.160451812e4);
            res.userData = {'par','perp'};
          end
          
          if strcmp(C{3},'TxTyTz')
            res = solo.db_get_ts(['solo_',C{1},'_',C{2}],['TxTyTz_',C{4}],Tint);
          end
          
        end
        
        % PAS ion VDFs: Daniel G. will add later
        if strcmp(C2{3},'vdf')
          error('VDFs not added yet')
        end
        
        
      end
      
      if strcmp(C2{2},'eas') % Electron data - can be added at a later date
        warning('Electron data to be added in the future')
      end
      
  end
  
  
end


if strcmp(C{1},'pos') % return the solo position - predicted, should be use flown?
  res = solo.get_position(Tint,'predicted','frame','SOLO_SUN_RTN');
end