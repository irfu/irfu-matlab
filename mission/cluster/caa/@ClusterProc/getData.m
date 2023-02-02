function data = getData(cp,cl_id,quantity,varargin)
%GETDATA(cp) produce Cluster level 2 and 3 data from the raw data
%
% data = getData(cp,cl_id,quantity,options)
%
% Input:
%   cp - ClusterProc object
%   cl_id - SC#
%   quantity - one of the following:
%
%   ec :   wcE{cl_id}p12/32, wcE{cl_id}p34 -> mERC // correct raw data
%          correct_sw_wake - correct wakes in the Solar Wind.
%   dies : diEs{cl_id}p12/32, diEs{cl_id}p34 -> mEDSI // spin fits [DSI]
%          also creates delta offsets D{cl_id}p12p34.
%          If the offset is real then it must be applied to p12/32,
%          if imaginary - to p34
%          Has the following OPTIONS:
%          sfit_ver - 0 (AIE c_efw_onesfit), 1 (BHN c_efw_c_efw_c_efw_spinfit_mx)
%          // default is to use the one specified in c_efw_sfit
%          probe_p - probe pair to use 12 or 34 [default 34]
%   die : diE{cl_id}p1234 -> mEDSIf // despun full res E [DSI]
%          also created ADC offsets Da{cl_id}p12 and Da{cl_id}p34
%   dief: diEF{cl_id}p1234, diESPEC{cl_id}p1234 -> mEDSIf
%          also created ADC offsets Da{cl_id}p12 and Da{cl_id}p34
%          // high-pass filtered full res E [DSI]
%   diespec: diESPEC{cl_id}p1234 -> mEDSI
%          // spectrum of full res E [DSI] +
%   idies, idie : idiEs{cl_id}p12, idiEs{cl_id}p34, idiE{cl_id}p1234 -> mEIDSI
%          Transform from SC to inertial frame
%   dieburst : dibE{cl_id}p1234 -> mEFWburst // despun i-burst E [DSI]
%   pburst : bP{cl_id}, bP{cl_id}_info, bNVps{cl_id}, P{freq}Hz{cl_id} -> mEFWburst
%          // i-burst P averaged from several probes
%
%   edbs,edb,iedb,iedbs : // Ez from E.B=0 [DSI+GSE]
%   E[s]{cl_id}, diE[s]{cl_id} -> mEdB    // SC frame
%   iE[s]{cl_id}, idiE[s]{cl_id} -> mEdBI  // Inertial frame
%          Has the following OPTIONS:
%          ang_limit - minimum angle(B,spin plane) [default 10 deg]
%          ang_blank - put Ez to NaN for points below ang_limit [default]
%          ang_fill - fill points below ang_limit with 1e27
%          ang_ez0 - use Ez=0 for points below ang_limit
%   p : P{cl_id}, P{cl_id}_info, NVps{cl_id},P10Hz{cl_id} -> mP	// P averaged
%   ps : Ps{cl_id} -> mP	// P spin resolution
%   whip: WHIP{cl_id} -> mEFW   // Whisper pulses present +1 precceding sec
%   bdump: DBUMP{cl_id} -> mEFW	// Burst dump present
%   efwa: EFWA{cl_id} -> mA	    // Phase from EFW FDM
%   sweep: SWEEP{cl_id} -> mEFW	// Sweep + dump present
%   badbias: BADBIASRESET{cl_id}, BADBIAS{cl_id}p[1..4] -> mEFW
%          // Bad bias settings
%   probesa: PROBESA{cl_id}p[1..4], SAASASE{cl_id}, SAASADI{cl_id} -> mEFW
%          // Probe saturation, including due to shadow for SAA=90 deg
%   hbiassa: HBIASSA{cl_id}p{12/32,34}, HBSATDSC{cl_id}p{12/32,34} -> mEFW
%          // Saturation due to high bias current
%   rawspec: RSPEC{cl_id}p{12/32,34} -> mEFW // Spectrum of raw signal (1,2..5 omega)
%   wake : PSWAKE{cl_id}p{12/34}, LOWAKE{cl_id}p{12/34} -> mEFW // wakes
%   edi : EDI{cl_id}, diEDI{cl_id} -> mEDI // EDI E in sc ref frame
%   br, brs : Br[s]{cl_id}, diBr[s]{cl_id} -> mBr // B resampled to E[s]
%   dibsc : diBSC{cl_id}, BSC{cl_id} -> mBSC // despun B STAFF SC [DSI+GSE]
%   dibscburst : wBSC4kHz{cl_id}, Atwo{cl_id} -> mEFWburst // despun B STAFF SC burst [DSI+GSE]
%   vedbs, vedb : VExB[s]{cl_id}, diVExB[s]{cl_id} -> mEdB // E.B=0 [DSI+GSE]
%   vce : VCE(p,h){cl_id},diVCE(p,h){cl_id} ->mCIS	// E CIS PP [GSE+DSI]
%   manproblems: whip|sweep|bdump|badbias|probesa|hbiassa|wake -> mEFW
%          // Reads manually-set problems from database.
%   hk: ibias puck guard HK{cl_id} -> mEFW // save housekeeping data
%
% Example:
%       getData(cp,4,'edbs','ang_fill','ang_limit',20,'probe_p',12)
%
% General options:
%       nosave : do no save on disk
%       rmwhip/withwhip : do/do not remove time intervals with Whisper pulses
%            (affects dies,die,dief,diespec,dieburst,hbiassa,rawspec,p,pburst)
%
% Specific options:
%       rmwake : remove wakes (affects edb,edbs,iedb,iedbs)
%       rmhbsa : remove saturation due to high bias current
%                 (affects dies,die,dief,diespec,dieburst)
%       no_saved_adc : do not use ADC offset computed from spinfits, e.g.
%           compute the offset (affects 'die')
%       no_adc_offset : do not correct ADC offset (affects 'die')
%       adc_offset_5points : compute ADC offset using 5-point averages
%       no_caa_delta : do not use the CAA Delta offset(c_efw_delta_off),
%           e.g. compute the offset from spinfits and use it afterwards
%           (affects 'dies','die')
%       check_caa_sh_interval - check for .caa_sh_interval/.caa_ms_interval.
%           'ec|vce' only processed if .caa_sh_interval found. correct_sw_wake flag set if needed.
%           'brs|edi|wake' only processed if .caa_ms_interval found.
%       ec_extraparams : extra paramters to pass to c_efw_swwake (affects ec)
%
% See also C_GET, C_CTL

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(3,15)
if nargin > 3, have_options = 1; args = varargin;
else, have_options = 0;
end

% default options
flag_save = 1;
flag_usesaved_adc_off = 1;
flag_usecaa_del_off = 1;
flag_use_adc_off = 1;
flag_edb = 1;
flag_rmwake = 0;
flag_rmhbsa = 0;
sfit_ver = -1;
correct_sw_wake = 0;
flag_wash_p32 = 1;
check_caa_sh_interval=0;
ec_extraparams = [];
nPointsADCOffset = 7;

flag_rmwhip = c_ctl(cl_id,'rm_whip');
if isempty(flag_rmwhip), flag_rmwhip = 1; end
flag_rmwhip_force = 0;
ang_limit = c_ctl(cl_id,'ang_lim');
if isempty(ang_limit), ang_limit = 10; end
probe_p = caa_sfit_probe(cl_id);
deltaof_sdev_max = c_ctl(cl_id, 'deltaof_sdev_max');
if isempty(deltaof_sdev_max), deltaof_sdev_max = 2; end
deltaof_max = c_ctl(cl_id, 'deltaof_max');
if isempty(deltaof_max), deltaof_max = 1.5; end
while have_options
  l = 1;
  switch(args{1})
    case 'nosave'
      flag_save = 0;
    case 'withwhip'
      flag_rmwhip = 0;
      flag_rmwhip_force = 0;
    case 'rmwhip'
      flag_rmwhip = 1;
      flag_rmwhip_force = 1;
    case 'rmhbsa'
      flag_rmhbsa = 1;
    case 'rmwake'
      flag_rmwake = 1;
    case 'no_caa_delta'
      flag_usecaa_del_off = 0;
    case 'no_saved_adc'
      flag_usesaved_adc_off = 0;
    case 'no_adc_offset'
      flag_use_adc_off = 0;
    case 'adc_offset_5points'
      nPointsADCOffset = 5;
    case 'ang_limit'
      if length(args)>1
        if isnumeric(args{2})
          ang_limit = args{2};
          l = 2;
        else, irf_log('fcal,','wrongArgType : ang_limit must be numeric')
        end
      else, irf_log('fcal,','wrongArgType : ang_limit value is missing')
      end
    case 'ang_blank'
      flag_edb = 1;	% [default]
    case 'ang_fill'
      flag_edb = 2;	% fill points below ang_limit with 1e27
      fill_val = 1e27;
    case 'ang_ez0'
      flag_edb = 0;	% use Ez=0 for points below ang_limit
    case 'probe_p'
      if length(args)>1
        if isnumeric(args{2})
          probe_p_tmp = args{2};
          l = 2;
        else, probe_p_tmp = str2double(args{2});
        end
        if (probe_p_tmp==12 || probe_p_tmp==34), probe_p = probe_p_tmp;
        else, irf_log('fcal,','wrongArgType : probe_p must be 12 or 34')
        end
      else, irf_log('fcal,','wrongArgType : probe_p value is missing')
      end
    case 'sfit_ver'
      if length(args)>1
        if isnumeric(args{2})
          l = 2;
          if	args{2}>=0 && args{2}<2, sfit_ver = args{2};
          else, irf_log('fcal,','wrongArgType : sfit_ver must be 0 or 1')
          end
        else, irf_log('fcal,','wrongArgType : sfit_ver must be numeric')
        end
      else, irf_log('fcal,','wrongArgType : sfit_ver value is missing')
      end
    case 'correct_sw_wake'
      correct_sw_wake = 1;
    case 'check_caa_sh_interval'
      check_caa_sh_interval = 1;
    case 'ec_extraparams'
      if length(args)>1
        if iscell(args{2})
          ec_extraparams = args{2};
          l = 2;
        else, irf_log('fcal,','wrongArgType : ec_extraparams must be a cell array')
        end
      else, irf_log('fcal,','wrongArgType : ec_extraparams value is missing')
      end
    otherwise
      irf_log('fcal,',['Option ''' args{1} '''not recognized'])
  end
  if length(args) > l, args = args(l+1:end);
  else, break
  end
end

save_list = '';

old_pwd = pwd;
cd(cp.sp) %enter the storage directory
if cp.sp~='.', irf_log('save',[quantity ' Storage directory is ' cp.sp]), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ec - correct raw Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'ec')
  save_file = './mERC.mat';
  
  if check_caa_sh_interval
    if exist('./.caa_ms_interval','file')
      irf_log('proc','Inside magnetosphere. No solar wind cleaning performed.')
      data = []; cd(old_pwd), return
    end
    if ~exist('./.caa_sh_interval','file')
      irf_log('proc','**** No .caa_XX_interval file! No solar wind cleaning performed.')
      data = []; cd(old_pwd), return
    end
    correct_sw_wake=1;
  end
  
  if ~any([correct_sw_wake]) %#ok<NBRAK> % List all possible methods here
    irf_log('proc','no cleaning method defined')
    data = []; cd(old_pwd), return
  end
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  [ok,pha] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',...
      irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  n_ok = 0;
  for probe = [12 34] % probe = [12 32 34]
    [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
    if ~ok || isempty(da)
      irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
      continue
    end
    
    fsamp = c_efw_fsample(da(:,1),'hx');
    
    problems = 'reset|bbias|probesa|probeld|sweep|bdump|nsops|whip'; %#ok<NASGU>
    nsops_errlist = [caa_str2errid('bad_bias') caa_str2errid('bad_hx')];%#ok<NASGU>
    signal = da; %#ok<NASGU>
    remove_problems
    da = res; %#ok<NODEF>
    clear res signal problems
    
    % Check if we have at least 1 spin of data left
    if length(find(~isnan(da(:,2)))) < 4*fsamp
      irf_log('proc',irf_ssub('No p? data after removals',probe))
      continue
    end
    
    % Load Whisper
    whip = [];
    if cl_id==2 || fsamp==450
      [ok,whip] = c_load('WHIP?',cl_id);
      if ~ok
        irf_log('load',	irf_ssub('Cannot load WHIP?',cl_id))
      end
      clear ok
    end
    
    modified = 0;
    if correct_sw_wake
      irf_log('proc',...
        sprintf('correcting solar wind wake on p%d',probe))
      [da, n_corrected,wake_dsc] = c_efw_swwake(da,probe,pha,whip,0,ec_extraparams); %#ok<ASGLU>
      
      if n_corrected>0
        eval(irf_ssub(...
          'WAKE?p!=wake_dsc;save_list=[save_list ''WAKE?p! ''];',...
          cl_id,probe));
        modified = 1;
      end
      clear n_corrected wake_dsc
    end
    
    if modified
      eval(irf_ssub('wcE?p!=da;save_list=[save_list ''wcE?p! ''];',...
        cl_id,probe));
      n_ok = n_ok + 1;
    end
    clear ok da
  end
  if ~n_ok, data = []; cd(old_pwd), return, end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dies, diehxs, dielxs - spin fiting of Electric field
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dies') || strcmp(quantity,'diehxs') || strcmp(quantity,'dielxs')
  save_file = './mEDSI.mat';
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  [ok,pha] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',...
      irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  if strcmp(quantity,'dies'), qlist = {'diehxs', 'dielxs'};
  else, qlist = {quantity};
  end
  
  for qq = qlist
    n_ok = 0; wE = {}; flag_lx = 0; p12 = ''; flag_have_p34 = 0;
    if strcmp(qq{:},'diehxs')
      flag_p12_loaded = 0;
      for probe = [12,32,34]
        flag_corr = 0;
        if probe == 32 && flag_p12_loaded, continue, end
        [ok,da] = c_load(irf_ssub('wcE?p!',cl_id,probe));
        if ~ok || isempty(da)
          [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
          if ~ok || isempty(da)
            irf_log('load', irf_ssub('No/empty wcE?p! and wE?p!',cl_id,probe));
            continue
          else, irf_log('load', irf_ssub('read wE?p!',cl_id,probe));
          end
          irf_log('load','using raw (not corrected) data')
        else
          irf_log('proc',irf_ssub('Using corrected data wcE?p!',cl_id,probe))
          flag_corr = 1;
        end
        if probe == 12, flag_p12_loaded = 1; p12 = 12;
        elseif probe == 32 && isempty(p12), p12 = 32;
        elseif probe == 34, flag_have_p34 = 1;
        end
        e = struct('probe',probe,'corr',flag_corr,'e',da);
        wE = [wE {e}]; n_ok = n_ok + 1; %#ok<AGROW>
        clear ok da e
      end
      clear flag_corr flag_p12_loaded
    elseif strcmp(qq{:},'dielxs')
      flag_lx = 1;
      p1 = []; p2 = []; p3 = []; p4 = [];
      for probe=1:4
        [ok,da] = c_load(irf_ssub('P10Hz?p!',cl_id,probe));
        if ~ok || isempty(da)
          irf_log('load', irf_ssub('No/empty P10Hz?p!',cl_id,probe));
          continue
        end
        c_eval('p?=da;',probe)
        clear ok da
      end
      
      p_sep = .088;
      if ~isempty(p1) && ~isempty(p2)
        if size(p1,1) ~= size(p2,1)
          [ii1, ii2]=irf_find_comm_idx(p1,p2);
          tt(:,1) = p2(ii2,1); tt(:,2) = ( p2(ii2,2) - p1(ii1,2) )/p_sep;
          clear ii1 ii2
        else
          tt(:,1) = p2(:,1); tt(:,2) = ( p2(:,2) - p1(:,2) )/p_sep;
        end
        e = struct('probe',12,'corr',0,'e',tt);
        wE = [wE {e}]; n_ok = n_ok + 1; %#ok<AGROW>
        clear tt e
        p12 = 12;
      end
      if ~isempty(p3) && ~isempty(p4)
        if size(p3,1) ~= size(p4,1)
          [ii3, ii4]=irf_find_comm_idx(p3,p4);
          tt(:,1) = p4(ii4,1); tt(:,2) = ( p4(ii4,2) - p3(ii3,2) )/p_sep;
          clear ii3 ii4
        else
          tt(:,1) = p4(:,1); tt(:,2) = ( p4(:,2) - p3(:,2) )/p_sep;
        end
        e = struct('probe',34,'corr',0,'e',tt);
        wE = [wE {e}]; n_ok = n_ok + 1; %#ok<AGROW>
        clear tt e
        flag_have_p34 = 1;
      end
      p_sep = .066;
      if ~isempty(p3) && ~isempty(p2)
        if size(p2,1) ~=  size(p3,1)
          [ii2,ii3]=irf_find_comm_idx(p2,p3);
          tt(:,1) = p2(ii2,1); tt(:,2) = ( p2(ii2,2) - p3(ii3,2) )/p_sep;
          clear ii2 ii3
        else
          tt(:,1) = p2(:,1); tt(:,2) = ( p2(:,2) - p3(:,2) )/p_sep;
        end
        e = struct('probe',32,'corr',0,'e',tt);
        wE = [wE {e}]; n_ok = n_ok + 1; %#ok<AGROW>
        clear tt e
        if isempty(p12), p12 = 32; end
      end
      if ~isempty(p2) && ~isempty(p4)
        if size(p2,1) ~=  size(p4,1)
          [ii2,ii4]=irf_find_comm_idx(p2,p4);
          tt(:,1) = p2(ii2,1); tt(:,2) = ( p2(ii2,2) - p4(ii4,2) )/p_sep;
          clear ii2 ii4
        else
          tt(:,1) = p2(:,1); tt(:,2) = ( p2(:,2) - p4(:,2) )/p_sep;
        end
        e = struct('probe',42,'corr',0,'e',tt);
        wE = [wE {e}]; n_ok = n_ok + 1; %#ok<AGROW>
        clear tt e
        if isempty(p12), p12 = 42; end
      end
      clear p1 p2 p3 p4 p_sep
    end
    if ~n_ok, data = []; continue, end
    
    if flag_lx, lx_str = 'LX'; else, lx_str=''; end
    for iPr=1:length(wE)
      if wE{iPr}.corr, ss = 'c'; else, ss = ''; end
      irf_log('proc',sprintf('Spin fit w%sE%dp%d %s -> diE%ss%dp%d',...
        ss,cl_id,wE{iPr}.probe,lx_str,lx_str,cl_id,wE{iPr}.probe))
      
      % Check if we have at least 3 data points to start with
      if length(wE{iPr}.e(:,1)) < 3
        irf_log('proc',irf_ssub('p? data too short',wE{iPr}.probe));
        if     wE{iPr}.probe==34, flag_have_p34 = 0;
        elseif wE{iPr}.probe==p12, p12=[];
        end
        continue
      end
      
      aa = c_phase(wE{iPr}.e(:,1),pha);
      if isempty(aa)
        irf_log('proc', sprintf('Empty phase, skipping w%sE%dp%d %s',...
          ss,cl_id,wE{iPr}.probe,lx_str))
        if wE{iPr}.probe==34, flag_have_p34 = 0;
        elseif wE{iPr}.probe==p12, p12=[];
        end
        continue
      end
      
      problems = 'reset|bbias|probesa|probeld|sweep|bdump|nsops|whip';
      if flag_lx
        ps = sprintf('%d',wE{iPr}.probe);
        nsops_errlist = [caa_str2errid('bad_bias') caa_str2errid('hxonly')...
          caa_str2errid('bad_lx')...
          caa_str2errid(['no_p' ps(1)]) caa_str2errid(['no_p' ps(2)])];%#ok<NASGU>
      else
        nsops_errlist = [caa_str2errid('bad_bias') caa_str2errid('bad_hx')...
          caa_str2errid(irf_ssub('no_p?',wE{iPr}.probe))];%#ok<NASGU>
      end
      % We remove Whisper only if explicitely asked for this by user
      if flag_rmwhip && flag_rmwhip_force, problems = [problems '|whip']; end %#ok<AGROW>
      if flag_rmhbsa, problems = [problems '|hbiassa']; end %#ok<AGROW,NASGU>
      signal = wE{iPr}.e; %#ok<NASGU>
      remove_problems
      wE{iPr}.e = res; %#ok<AGROW,NODEF>
      clear res signal problems
      
      noData = false;
      nDataPoints = length(find(~isnan(wE{iPr}.e(:,2))));
      if ~nDataPoints, noData = true;
      else
        if flag_lx, tmode = 'lx'; else, tmode = 'hx'; end
        fsamp = c_efw_fsample(wE{iPr}.e,tmode);
        if ~fsamp, error('no sampling frequency'),end
        
        % Check if we have at least 1 spin of data left
        if nDataPoints < 4*fsamp, noData = true; end
      end
      
      if noData
        irf_log('proc',...
          sprintf('No p%d%s data after removals',wE{iPr}.probe,lx_str));
        if     wE{iPr}.probe==34, flag_have_p34 = 0;
        elseif wE{iPr}.probe==p12, p12=[];
        end
        continue
      end
      
      % Run spin fitting
      sfit_ver_bak = sfit_ver; % save sfit_ver for next probe pair
      if sfit_ver>=0
        irf_log('proc',['using SFIT_VER=' num2str(sfit_ver)])
      else
        if wE{iPr}.probe==32 || wE{iPr}.probe==42, sfit_ver = 2;
        else, sfit_ver = 1;
        end
      end
      
      sp = c_efw_sfit(wE{iPr}.probe,3,10,20,wE{iPr}.e(:,1),wE{iPr}.e(:,2),...
        aa(:,1),aa(:,2),sfit_ver,tmode);
      sfit_ver = sfit_ver_bak;
      
      % Check if we have any data left
      if isempty(sp)
        irf_log('load',sprintf('No data left after spinfit for C%d p%d',...
          cl_id,wE{iPr}.probe))
        continue
      end
      
      % Remove erroneous points with zero time
      ind = find(sp(:,1)>0);
      if length(ind)<length(sp(:,1))
        irf_log('proc',...
          [num2str(length(sp(:,1))-length(ind)) ' spins removed (bad time)']);
        sp = sp(ind,:);
      end
      % Remove spins with bad spin fit (obtained E > 10000 mV/m)
      ind = find(abs(sp(:,3))>1e4); sp(ind,:) = [];
      if ind, disp([num2str(length(ind)) ' spins removed due to E>10000 mV/m']);end
      
      % ADC offsets, i.e. DC offset of the raw sinals
      % Warn about points with sdev>.8
      ii = find(sp(:,6)>.8);
      if length(ii)/size(sp,1)>.05
        irf_log('proc',[sprintf('%.1f',100*length(ii)/size(sp,1))...
          '% of spins have SDEV>.8 (ADC offsets)']);
      end
      
      % Remove whisper pulses for ADC offset and 2 Omega
      [ok,whip,msg] = c_load('WHIP?',cl_id);
      if ok
        if ~isempty(whip)
          irf_log('proc','blanking Whisper pulses in ADC offsets')
          if (wE{iPr}.probe == 32 || wE{iPr}.probe == 42) && size(sp,2) == 10
            SPCOLS = [4 9 10];
          else, SPCOLS = 4;
          end
          % Obtain time interval boundaries around each
          % (averaged) data point;
          data_time_lower = sp(:, 1) - 2;
          data_time_upper = sp(:, 1) + 2;
          for num = 1:size(whip, 1)
            problem_start = whip(num, 1);
            problem_stop = whip(num, 2);
            sp((problem_start >= data_time_lower & problem_start < data_time_upper) | ...
              (problem_stop  >  data_time_lower & problem_stop <= data_time_upper) | ...
              (problem_start <= data_time_lower & problem_stop >= data_time_upper),...
              SPCOLS) = NaN;
          end
        end
      else, irf_log('load',msg)
      end
      clear ok whip msg
      
      % Continue with ADC offsets
      adc_off = sp(:,[1 4]);
      % Don't include hbiassa intervals when calculating mean, variance
      max_off=adc_off(~isnan(adc_off(:,2)),:);
      [ok,hbsa] = c_load(irf_ssub('HBIASSA?p!',cl_id,12));
      if ok && ~isempty(hbsa),max_off = caa_rm_blankt(max_off,hbsa); end
      [ok,hbsa] = c_load(irf_ssub('HBIASSA?p!',cl_id,34));
      if ok && ~isempty(hbsa),max_off = caa_rm_blankt(max_off,hbsa); end
      [ok,hbsa] = c_load(irf_ssub('HBIASSA?p!',cl_id,32));
      if ok && ~isempty(hbsa),max_off = caa_rm_blankt(max_off,hbsa); end
      [ok,hbsa] = c_load(irf_ssub('HBIASSA?p!',cl_id,42));
      if ok && ~isempty(hbsa),max_off = caa_rm_blankt(max_off,hbsa); end
      max_off=max_off(~isnan(max_off(:,2)),:);
      adc_off_mean = mean(max_off(:,2));
      max_off = 3*std(max_off(:,2));
      % Fill NaNs with mean value
      adc_off(isnan(adc_off(:,2)),2) = adc_off_mean;
      % For C4p34, check for DAC problems if date is after 2009-01-01
      adc_despike=1;
      if cl_id == 4 && wE{iPr}.probe == 34 && adc_off(1,1) > 1230768000
        % Check for BADDAC set earlier (i.e. by manual forcing)
        [ok, badDAC] = c_load(irf_ssub('BADDAC?p!', cl_id,wE{iPr}.probe));
        if ok && ~isempty(badDAC)
          irf_log('proc', 'BADDAC set manually')
          adc_despike=0;
        end
        
        % Identify intervals
        idx = find( abs( adc_off(:,2)) > 1);
        if length(idx) > 150 % Longer interval of DAC problems
          adc_despike=0;
          badDAC=[adc_off(idx(1),1)-20 adc_off(idx(end),1)+20]; %#ok<NASGU>
          irf_log('proc','Long interval of DAC problems.');
        else
          if length(idx) > 5
            % Check for short interval at start
            idx2=find(idx < 150);
            if length(idx2) > 3 && idx2(end) < idx(idx2(end))*1.3
              adc_despike=0;
              if adc_despike,badDAC=[adc_off(1,1)-20 adc_off(idx(idx2(end)),1)+20]; %#ok<UNRCH>
              else, badDAC=[badDAC' [adc_off(1,1)-20 adc_off(idx(idx2(end)),1)+20]']';
              end
              irf_log('proc','Short interval of DAC problems at start of interval.');
            end
            % Check for short interval at end
            idx2=find(idx > length(adc_off(:,1))-150);
            if length(idx2) > 3 && idx2(end)-idx2(1) < (idx(idx2(end))-idx(idx2(1)))*1.3
              if adc_despike, badDAC=[adc_off(idx(idx2(1)),1)-20 adc_off(end,1)+20]; %#ok<NASGU>
              else, badDAC=[badDAC' [adc_off(idx(idx2(1)),1)-20 adc_off(end,1)+20]']'; %#ok<NASGU>
              end
              adc_despike=0;
              irf_log('proc','Short interval of DAC problems at end of interval.');
            end
            % Check for short interval in middle
            if adc_despike
              idx2=find(diff(idx)<3);
              if length(idx2) > 15
                badDAC=[adc_off(idx(1),1)-20 adc_off(idx(end),1)+20]; %#ok<NASGU>
                adc_despike=0;
                irf_log('proc','Short interval of DAC problems in middle of interval.');
              end
            end
          end
        end
        
        % Save problem to file mEFW.mat (not save_file=mEDSI.mat)
        if ~adc_despike
          badDACname=irf_ssub('BADDAC?p!',cl_id,wE{iPr}.probe);
          irf_log('save', [badDACname ' -> mEFW.mat']);
          eval([badDACname '=badDAC;']);
          if exist('mEFW.mat','file')
            eval(['save -append mEFW.mat ' badDACname]);
          else
            eval(['save mEFW.mat ' badDACname]);
          end
          % Tell export routines to use p12
          caa_sfit_probe(cl_id,p12);
        end
      end
      
      % Unless a suspect DAC interval was detected, take care of spikes
      if adc_despike
        % Take care of spikes: replace extreme values with mean value
        idx = find( abs( adc_off(:,2) - adc_off_mean ) > max_off );
        if ~isempty(idx)
          adc_off(idx,2) = 0;
          adc_off_mean = mean(adc_off(abs(adc_off(:,2))>0,2));
          adc_off(adc_off(:,2)==0,2) = adc_off_mean;
          irf_log('proc',sprintf('%d spikes (ADC offsets)',length(idx)));
        end
      end
      % Smooth the ADC offset signal
      adc_off = irf_waverage(adc_off,1/4,nPointsADCOffset); %#ok<NASGU>
      
      % Save 2 omega separately
      if wE{iPr}.probe == 32 && size(sp,2) == 10
        % Fill NaNs with zeros and smooth the signal
        for col = 9:10
          sp( isnan(sp(:,col)) ,col) = 0;
          sp(:,[1 col]) = irf_waverage(sp(:,[1 col]),1/4);
        end
        
        eval(irf_ssub(...
          ['w2W' lx_str '?p!=sp(:,[1 9 10]);save_list=[save_list ''w2W' lx_str '?p! ''];'],...
          cl_id,wE{iPr}.probe));
      end
      
      sp = sp(:,[1:4 6]);
      sp(:,4) = 0*sp(:,4); %#ok<NASGU> % Z component
      
      eval(irf_ssub(...
        ['diE' lx_str 's?p!=sp;Dadc' lx_str ...
        '?p!=adc_off;save_list=[save_list ''diE' lx_str 's?p! Dadc' lx_str '?p! ''];'],...
        cl_id,wE{iPr}.probe));
      clear tt sp adc_off
    end
    
    % Delta offsets (offsets between E DSI obtained from p12/32 and p34)
    if flag_have_p34 && ~isempty(p12)
      
      % Compute delta offsets, i.e. difference between DSI E obtained
      % from different probe pairs
      % Remove points which are > deltaof_sdev_max*sdev
      % as this must be a stable quantity
      eval(irf_ssub(['[ii1,ii2] = irf_find_comm_idx(diE' lx_str 's?p!,diE'...
        lx_str 's?p34);'],cl_id,p12))
      eval(irf_ssub(['df=diE' lx_str 's?p!(ii1,1:3);df(:,2:3)=diE' lx_str ...
        's?p!(ii1,2:3)-diE' lx_str 's?p34(ii2,2:3);'],cl_id,p12))
      clear ii1 ii2
      
      % Remove saturation due to too high bias current
      if ~flag_rmhbsa % Otherwise the saturation is already blanked
        for probe=[p12 34]
          [ok,hbsa,msg] = c_load(irf_ssub('HBIASSA?p!',cl_id,probe));
          if ok
            if ~isempty(hbsa)
              irf_log('proc','blanking HB saturation')
              df = caa_rm_blankt(df,hbsa);
            end
          else, irf_log('load',msg)
          end
          clear ok hbsa msg
        end
      end
      
      df = df (:,2:3);
      df(isnan(df(:,2)),:) = [];
      iia = [];
      if size(df,1)>2
        sdev = std(df);
        for comp = 1:2
          ii = find(abs(df(:,comp)-mean(df(:,comp))) > deltaof_sdev_max*sdev(comp));
          iia = [iia; ii]; %#ok<AGROW>
        end
        iia = sortrows(iia(:));
        iia(diff(iia)==0) = [];
        irf_log('calb',sprintf('%d points are removed for delta offsets',...
          length(iia)))
      end
      Del = [0 0];
      for comp = 1:2
        ddd = df(:,comp); ddd(iia) = [];
        Del(comp) = mean(ddd);
      end
      
      irf_log('calb',sprintf('%s delta offsets are: %.2f [x] %.2f [y]', ...
        lx_str, Del(1), Del(2)))
      
      % Check for unreallistically large Del.
      % If it is larger than deltaof_max, we set it to zero and
      % NOT doing any corrections.
      if ( abs(Del(1)) > deltaof_max ) || ( abs(Del(2)) > deltaof_max )
        irf_log('calb',...
          irf_ssub('DELTA OFFSET TOO BIG >!. Setting D?p12p34=[0 0]',...
          cl_id,deltaof_max))
        Del=[0 0]; %#ok<NASGU>
      else
        % Always correct p12/p32.
        % Deprecated behavior: real offset is applied to p12, imaginary to p34.
        irf_log('calb',irf_ssub('correcting p?',p12))
        eval(irf_ssub(['diE' lx_str 's?p!(:,2:3)=diE' lx_str ...
          's?p!(:,2:3)-ones(size(diE' lx_str 's?p!,1),1)*Del;'],cl_id,p12));
      end
      
      eval(irf_ssub(['D' lx_str '?p12p34=Del;'],cl_id))
      clear m12 m34 Del
      
      eval(irf_ssub(['save_list=[save_list ''D' lx_str '?p12p34 ''];'],cl_id));
    elseif ~flag_have_p34 && ~flag_lx % XXX: fix me
      % If we have only p12, we should use it for sfits&co
      if exist(irf_ssub('diEs?p!',cl_id,12), 'var')
        irf_log('calb',irf_ssub('Will use p? for spin res data',12))
        caa_sfit_probe(cl_id,12);
      elseif exist(irf_ssub('diEs?p!',cl_id,32), 'var')
        irf_log('calb',irf_ssub('Will use p? for spin res data',32))
        caa_sfit_probe(cl_id,32);
        %elseif exist(irf_ssub('diEs?p!',cl_id,42), 'var')
        %	irf_log('calb',irf_ssub('Will use p? for spin res data',42))
        %	caa_sfit_probe(cl_id,42);
      end
    end
  end
  if strcmp(quantity,'dies')
    % Decide on which probes to use
    if ~any(strfind(save_list,sprintf('diEs%dp',cl_id))) && ...
        any(strfind(save_list,sprintf('diELXs%dp',cl_id)))
      for pp = [34 12 32 42]
        if strfind(save_list,sprintf('diELXs%dp%d',cl_id,pp)) %#ok<STRIFCND>
          caa_sfit_probe(cl_id,pp*10);
          irf_log('calb',irf_ssub('Will use p?LX for spin res data',pp))
          break
        end
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % die - despin of full resolution data.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'die') || strcmp(quantity,'dief') || ...
    strcmp(quantity,'diespec') || strcmp(quantity,'dieburst') || ...
    strcmp(quantity,'dielx') || strcmp(quantity,'dielxspec')
  
  if strcmp(quantity,'dieburst'), do_burst = 1; else, do_burst = 0; end
  if strcmp(quantity,'dief'), do_filter = 1; else, do_filter = 0; end
  if strcmp(quantity,'dielx') || strcmp(quantity,'dielxspec'), flag_lx = 1;
  else, flag_lx = 0;
  end
  if do_burst
    save_file = './mEFWburst.mat';
    var1_name = 'dibE?p1234';
  else
    if strcmp(quantity,'diespec')
      save_file = './mEDSI.mat';
      var1_name = 'diESPEC?p1234';
    elseif strcmp(quantity,'dielxspec')
      save_file = './mEDSI.mat';
      var1_name = 'diELXSPEC?p1234';
    else
      save_file = './mEDSIf.mat';
      if do_filter, var1_name = 'diEF?p1234';
      elseif flag_lx, var1_name = 'diELX?p1234';
      else
        var1_name = 'diE?p1234';
      end
    end
  end
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  % Aux quantities: Dadc?
  [ok,pha] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',...
      irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  flag_bp = 0; % Flag for 8kHz BP filter internal burst data
  % Make electric field for the burst
  if do_burst || flag_lx
    p_ok = [];
    loaded = 0;
    e12 = []; p12 = 12;
    if ~flag_lx
      % First try to see if we have any V12H and V43H
      [ok,e12] = c_load('wE8kHz?p12',cl_id);
      if ~ok || isempty(e12)
        irf_log('load', irf_ssub('No/empty wE8kHz?p12',cl_id));
      else
        loaded = 1;
        p_ok = 12;
        flag_bp = 1;
      end
      [ok,e34] = c_load('wE8kHz?p34',cl_id);
      if ~ok || isempty(e34)
        irf_log('load', irf_ssub('No/empty wE8kHz?p34',cl_id));
      else
        loaded = 1;
        p_ok = [p_ok 34];
        flag_bp = 1;
      end
    end
    
    if flag_lx || ~loaded
      if flag_lx
        var_name = 'wELX?p!';
        param={'10Hz'};
      else
        var_name = 'wbE?p!';
        param={'180Hz','4kHz','32kHz'};
      end
      for k=1:length(param)
        for probe=[2 4]
          [ok,p2] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,probe));
          if ok
            loaded = 1;
            p_sep = .088;
            if probe==2
              [ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,1));
              p12 = 12;
              if ~ok
                [ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,3));
                p_sep = .066;
                p12 = 32;
              end
              if ~ok && flag_lx
                [ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,4));
                p_sep = .066;
                p12 = 42;
              end
              if ~ok
                if flag_lx, irf_log('load', irf_ssub(['No P' param{k} '?p1/3/4'],cl_id));
                else, irf_log('load', irf_ssub(['No P' param{k} '?p1/3'],cl_id));
                end
                continue
              end
              if size(p2,1) ~= size(p1,1)
                [ii1, ii2]=irf_find_comm_idx(p1,p2);
                e12(:,1) = p2(ii2,1);
                e12(:,2) = ( p2(ii2,2) - p1(ii1,2) )/p_sep;
              else
                e12(:,1) = p2(:,1);
                e12(:,2) = ( p2(:,2) - p1(:,2) )/p_sep;
              end
              BURST_MEAN_THRESH = 100; % mV/m
              if abs(mean(e12(~isnan(e12(:,2)),2)))>BURST_MEAN_THRESH % Sanity check for realistic electric field value
                e12 = [];
                irf_log('proc',sprintf('Disregarding p%d (mean>%d mV/m)',...
                  p12,BURST_MEAN_THRESH))
              else
                p_ok = [p_ok 12]; %#ok<AGROW>
                eval(irf_ssub('wbE?p!=e12;save_list=[save_list ''wbE?p! ''];',cl_id, p12));
              end
            else
              [ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,3));
              if ~ok
                irf_log('load', irf_ssub(['No P' param{k} '?p3'],cl_id));
                continue
              end
              
              if size(p2,1) ~= size(p1,1)
                [ii1, ii2]=irf_find_comm_idx(p1,p2);
                e34(:,1) = p2(ii2,1);
                e34(:,2) = ( p2(ii2,2) - p1(ii1,2) )/p_sep;
              else
                e34(:,1) = p2(:,1);
                e34(:,2) = ( p2(:,2) - p1(:,2) )/p_sep;
              end
              p_ok = [p_ok 34]; %#ok<AGROW>
              eval(irf_ssub('wbE?p!=e34;save_list=[save_list ''wbE?p! ''];',cl_id, 34));
            end
          else
            irf_log('load', irf_ssub(['No P' param{k} '?p!'],cl_id, probe));
          end
        end
        if loaded, break, end
      end
      if ~loaded && ~flag_lx
        p_ok = [];
        for probe = [12 32 34]
          [ok,da] = c_load(irf_ssub(var_name,cl_id,probe));
          if ~ok || isempty(da)
            irf_log('load', irf_ssub(['No/empty ' var_name],cl_id,probe));
            continue
          end
          
          if probe==32
            p12 = 32;
            e12 = da;
            p_ok = [p_ok 12]; %#ok<AGROW>
          else
            c_eval('e?=da;',probe)
            p_ok = [p_ok probe]; %#ok<AGROW>
          end
          clear ok da
        end
      end
    end
  else % Not burst
    p12 = 12; e12 = []; e34 =[];
    p_ok = []; p12_ok = 0;
    for probe = [12 32 34]
      if probe == 32 && p12_ok, continue, end
      [ok,da] = c_load(irf_ssub('wcE?p!',cl_id,probe));
      
      if ~ok || isempty(da)
        irf_log('load', irf_ssub('No/empty wcE?p!',cl_id,probe));
        if ~ok || isempty(da)
          [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
          if ~ok || isempty(da)
            irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
            continue
          end
          irf_log('load','using raw (not corrected) data')
        end
      end
      
      if ok && probe == 12, p12_ok = 1; end
      
      if probe==32
        p12 = 32;
        e12 = da;
        p_ok = [p_ok 12]; %#ok<AGROW>
      else
        c_eval('e?=da;',probe)
        p_ok = [p_ok probe]; %#ok<AGROW>
      end
      clear ok da
    end
  end
  if isempty(p_ok), data = []; cd(old_pwd), return, end
  clear p12_ok
  
  % Load ADC offsets
  for probe = p_ok
    if probe==12, ps=p12; else, ps=probe; end
    if flag_lx, wStr = 'DadcLX?p!';
    else, wStr = 'Dadc?p!';
    end
    
    [ok,dadc] = c_load(irf_ssub(wStr,cl_id,ps));
    if ~ok
      irf_log('load',irf_ssub(['Cannot load ' wStr],cl_id,ps))
      dadc = []; %#ok<NASGU>
    elseif isempty(dadc)
      irf_log('load',irf_ssub(['Empty ' wStr],cl_id,ps))
    end
    c_eval('dadc?=dadc;',probe)
    clear ok dadc
  end
  
  % calibration coefficients // see c_efw_despin
  coef=[[1 0 0];[1 0 0]];
  
  n_sig = 0;
  for p = p_ok
    if p==12
      ps = p12; tt = e12; dadc = dadc12;
    else
      ps = p; tt = e34; dadc = dadc34;
    end
    if do_burst
      % Correct ADC offset
      if  ~flag_bp
        if ~isempty(dadc)
          irf_log('calb','using saved ADC offsets')
          tmp_adc = irf_resamp(dadc,tt,'fsample',c_efw_fsample(tt,'ib'));
          tt(:,2) = tt(:,2) - tmp_adc(:,2);
          clear tmp_adc
        else, irf_log('calb','saved ADC offset empty')
        end
      end
      problems = 'reset|bbias|probesa|probeld|sweep|bdump|nsops|whip'; %#ok<NASGU>
      signal = tt; %#ok<NASGU>
      probe = ps; %#ok<NASGU>
      remove_problems
      tt = res; %#ok<NASGU,NODEF>
      clear res signal problems probe
      n_sig = n_sig + 1;
      c_eval('e?=tt;',p)
      
    else
      % Check if we have at least 2 data points to start with
      if size(tt,1)<2
        irf_log('proc',irf_ssub('p? data too short',ps))
        c_eval('e?=[];',p)
        continue
      end
      
      
      
      problems = 'reset|bbias|probesa|probeld|sweep|bdump|nsops';
      nsops_errlist = [caa_str2errid('bad_bias') caa_str2errid(irf_ssub('no_p?',ps))];
      if flag_lx
        nsops_errlist = [nsops_errlist caa_str2errid('bad_lx')]; %#ok<AGROW,NASGU>
      else
        nsops_errlist = [nsops_errlist caa_str2errid('bad_hx')]; %#ok<AGROW,NASGU>
        
        fsamp = c_efw_fsample(tt,'hx');
        if ~fsamp, error('no sampling frequency'),end
        % Always remove Whisper when we use 180Hz filter
        if (fsamp == 450) || ...
            ( cl_id == 2 && tt(1,1)>toepoch([2001 07 23 13 54 18]) ) || ...
            ( flag_rmwhip && flag_rmwhip_force )
          problems = [problems '|whip']; %#ok<AGROW>
        end
      end
      if flag_rmhbsa, problems = [problems '|hbiassa']; end %#ok<AGROW,NASGU>
      signal = tt; %#ok<NASGU>
      probe = ps; %#ok<NASGU>
      remove_problems
      tt = res; %#ok<NODEF>
      clear res signal problems probe
      
      noData = false;
      nDataPoints = length(find(~isnan(tt(:,2))));
      
      if ~nDataPoints, noData = true;
      else
        if flag_lx, tmode = 'lx'; else, tmode = 'hx'; end
        fsamp = c_efw_fsample(tt,tmode);
        if ~fsamp, error('no sampling frequency'),end
        
        % Check if we have at least 1 sec of data left
        if nDataPoints < fsamp, noData = true; end
      end
      
      if noData
        irf_log('proc',irf_ssub('No p? data after removals',ps))
        c_eval('e?=[];',p)
        continue
      end
      
      % Correct ADC offset
      if flag_use_adc_off && flag_usesaved_adc_off && ~flag_bp
        % Correct ADC offset
        if ~isempty(dadc)
          irf_log('calb','using saved ADC offset')
          if size(dadc,1)==1, tmp_adc = dadc;
          else, tmp_adc = irf_resamp(dadc,tt,'fsample',fsamp);
          end
          tt(:,2) = tt(:,2) - tmp_adc(:,2);
          clear tmp_adc
        else
          irf_log('calb','saved ADC offset empty')
          flag_usesaved_adc_off = 0;
        end
      end
      if flag_use_adc_off && ~flag_usesaved_adc_off  && ~flag_bp
        irf_log('calb','computing ADC offsets (simple averaging)')
        [tt,dadc] = caa_corof_adc(tt); %#ok<ASGLU>
        irf_log('calb',sprintf('Da%ddp%d : %.2f',cl_id,ps,dadc))
        eval(irf_ssub('Da?p!=dadc;',cl_id,ps));
      end
      n_sig = n_sig + 1;
      c_eval('e?=tt;',p)
    end
  end
  
  if n_sig==0
    irf_log('load','No raw data left for despin')
    data = []; cd(old_pwd); return
  end
  
  if n_sig==2
    if p12==32
      E_info.probe = '3234';
    else
      E_info.probe = '1234';
    end
    if size(e12,1)~=size(e34,1) || any(e12(:,1)-e34(:,1)~=0)
      % Different timelines. Need to correct
      irf_log('proc','making common timeline')
      [ii12,ii34] = irf_find_comm_idx(e12,e34);

      % If more than 50% of data missing on one probe
      % we base timeline on the probe pair which has more data
      if length(ii12) < .5*max(size(e12,1),size(e34,1))
        if size(e12,1)>size(e34,1)
          irf_log('proc',['Setting Ep34 to 0, except for '...
            num2str(length(ii12)) ' data points'])
          E_info.probe = num2str(p12);
          e34_tmp = e12;
          e34_tmp(~isnan(e12(:,2)),2) = 0;
          e34_tmp(ii12,2) = e34(ii34,2);
          e34 = e34_tmp;
          clear e34_tmp
        else
          irf_log('proc',['Setting Ep' num2str(p12)...
            ' to 0, except for ' num2str(length(ii12)) ' data points'])
          E_info.probe = '34';
          e12_tmp = e34;
          e12_tmp(~isnan(e34(:,2)),2) = 0;
          e12_tmp(ii34,2) = e12(ii12,2);
          e12 = e12_tmp;
          flag_wash_p32=0;
          clear e12_tmp
        end
        irf_log('proc','using one probe pair some part of the interval')
      else
        irf_log('proc',['Ep' num2str(p12) ' ' num2str(length(e12)) '->'...
          num2str(length(ii12)) ' data points'])
        e12 = e12(ii12,:);
        irf_log('proc',['Ep34 ' num2str(length(e34)) '->'...
          num2str(length(ii34)) ' data points'])
        e34 = e34(ii34,:);
      end
    end
    
    % Check for problem with one probe pair
    ii12 = find(~isnan(e12(:,2)));
    ii34 = find(~isnan(e34(:,2)));
    if do_burst, fsamp = c_efw_fsample(e12,'ib');
    else, fsamp = c_efw_fsample(e12,'hx');
    end
    % If more than 50% of data missing on one probe
    % we base timeline on the probe pair which has more data
    if abs(length(ii12)-length(ii34)) >...
        .5*( max(e12(end,1),e34(end,1)) - min(e12(1,1),e34(1,1)) )*fsamp
      if length(ii12)>length(ii34)
        ii = find(isnan(e34(:,2)) & ~isnan(e12(:,2)));
        e34(ii,2) = 0;
        E_info.probe = '12';
        irf_log('proc',['p34: NaN->0 for ' num2str(length(ii)) ' points'])
      else
        ii = find(isnan(e12(:,2)) & ~isnan(e34(:,2)));
        e12(ii,2) = 0;
        irf_log('proc',['p' num2str(p12) ': NaN->0 for ' ...
          num2str(length(ii)) ' points'])
        E_info.probe = '34';
      end
      irf_log('proc','using one probe pair some part of the interval')
    end
    
    % Use WEC coordinate system E=[t,0,p34,p12]
    full_e = zeros(length(e12),4);
    full_e(:,[1,4]) = e12;
    full_e(:,3) = e34(:,2);
    clear e12 e34
    
    % Load Delta offsets D?p12p34
    Del = [];
    if flag_usecaa_del_off  && ~flag_bp
      Del = c_efw_delta_off(full_e(1,1),cl_id);
    end
    if (~flag_usecaa_del_off || isempty(Del))  && ~flag_bp
      % Try saved offsets
      [ok,Del] = c_load('D?p12p34',cl_id);
      if ~ok || isempty(Del)
        irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id))
      end
      E_info.deltaoffauto = 1;
    else
      E_info.deltaoffauto = 0;
    end
    E_info.deltaoff = Del;
    
    if ~isempty(Del)
      if (E_info.deltaoffauto), dos = 'auto'; else, dos='caa'; end
      for comp=1:2
        if comp==1, x='x'; else, x='y'; end
        if isreal(Del(comp))
          % Real Del means we must correct p12/p32.
          irf_log('calb',['corr E' x ' delta offset(' dos ') on p' num2str(p12)])
          i_c = 1; % p12
        else
          % Correcting p34 is generally DEPRECATED, but can be used
          % when necessary (see c_efw_delta_off)
          irf_log('calb',['corr E' x ' delta offset(' dos ') on p34'])
          Del = imag(Del);
          i_c = 2; % p34
        end
        if comp==1, coef(i_c,3) = coef(i_c,3) + Del(1);
        else,       coef(i_c,3) = coef(i_c,3) - Del(2)*1j;
        end
      end
    end
    clear Del
  else
    
    % We have one probe pair
    if ~isempty(e12)
      pp = 12;
      E_info.probe = num2str(p12); %#ok<STRNU>
      EE = e12;
      clear e12
    else
      pp = 34;
      E_info.probe = '34'; %#ok<STRNU>
      EE = e34;
      clear e34
    end
    % Use WEC coordinate system E=[t,0,p34,p12]
    full_e = zeros(length(EE),4);
    full_e(:,1) = EE(:,1);
    if pp==12, full_e(:,4) = EE(:,2);
    else, full_e(:,3) = EE(:,2);
    end
    clear EE pp
  end
  
  c_eval([var1_name '_info=E_info;save_list=[save_list ''' var1_name '_info ''];'],cl_id);
  
  % Do actual despin
  aa = c_phase(full_e(:,1),pha);
  if isempty(aa)
    irf_log('proc','Can not do despin')
    data = []; cd(old_pwd); return
  end
  if p12==32 || p12==42
    
    % x(4)*cos(2*pha)+x(5)*sin(2*pha)
    if flag_wash_p32
      [ok,w2] = c_load('w2W?p32',cl_id);
      if ~ok || isempty(w2)
        irf_log('load', irf_ssub('No/empty w2W?p32',cl_id));
      else
        irf_log('proc','Washing asymmetric probe')
        w2(isnan(w2(:,2)),:) = [];
        w2r = irf_resamp(w2,full_e(:,1),'fsample',c_efw_fsample(full_e(:,1),'hx'));
        [i1,i2] = irf_find_comm_idx(aa(:,1),full_e(:,1));
        full_e(i2,4) = full_e(i2,4) + w2r(i2,2).*cos(2*aa(i1,2)*pi/180) + ...
          w2r(i2,3).*sin(2*aa(i1,2)*pi/180);
      end
    end
    
    full_e=c_efw_despin(full_e,aa,coef,'asym');
  else, full_e=c_efw_despin(full_e,aa,coef);
  end
  
  if strcmp(quantity,'diespec') || strcmp(quantity,'dielxspec')
    % Make a spectrum and save it.
    if flag_lx, sfreq = c_efw_fsample(full_e(:,1),'lx');
    else, sfreq = c_efw_fsample(full_e(:,1),'hx');
    end
    if ~sfreq, error('no sampling frequency'),end
    if sfreq < 6, nfft = 128; %#ok<NASGU>
    elseif sfreq < 26, nfft = 512; %#ok<NASGU>
    else, nfft = 8192; %#ok<NASGU>
    end
    c_eval([var1_name...
      '=irf_powerfft(full_e(:,1:3),nfft,sfreq);save_list=[save_list '''...
      var1_name '''];'],cl_id);
  else
    % HP-Filter
    if do_filter
      full_e = caa_filter_e(full_e(:,1:3));
      full_e(:,4) = 0;
    end
    
    % DS-> DSI
    full_e(:,3)=-full_e(:,3); %#ok<NASGU>
    
    c_eval([var1_name '=full_e; save_list=[save_list ''' var1_name '''];'],cl_id);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % idie, idies - DSI inertial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'idie') || strcmp(quantity,'idies')
  save_file = './mEIDSI.mat';
  
  if strcmp(quantity,'idie')
    var_s = {'diE?p1234'}; e_opt = 'die';
    var_b = 'diBr?'; b_opt ='br';
  else
    e_opt = 'dies';
    var_s = {'diEs?p12', 'diEs?p34'};
    var_b = 'diBrs?'; b_opt ='brs';
  end
  
  % Load resampled B
  [ok,diB] = c_load(var_b,cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub(['No ' var_b ' in mBr. Use getData(CP,cl_id,''' b_opt ''')'],cl_id))
    data = []; cd(old_pwd); return
  end
  
  [ok,diV] = c_load('diV?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No diV? in mR. Use getData(CDB,...,cl_id,''v'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  evxb = irf_tappl(irf_cross(diB,irf_resamp(diV,diB)),'*1e-3*(-1)');
  
  err_s = '';
  for k=1:length(var_s)
    [ok,diE] = c_load(var_s{k},cl_id);
    if ~ok
      if isempty(err_s), err_s = var_s{k};
      else, err_s = [err_s ', ' var_s{k}]; %#ok<AGROW>
      end
      continue
    end
    
    enew = diE;
    % We take only X and Y components. Z must remain zero.
    etmp = irf_resamp(evxb(:,1:3),enew(:,1),'fsample',c_efw_fsample(enew));
    enew(:,2:3) = diE(:,2:3) - etmp(:,2:3); %#ok<NASGU>
    clear etmp
    
    c_eval(['i' var_s{k} '= enew; clear enew'],cl_id)
    save_list=[save_list 'i' irf_ssub(var_s{k},cl_id) ' ']; %#ok<AGROW>
  end
  if ~isempty(err_s)
    irf_log('load',...
      irf_ssub(['No ' err_s ' in mEDSI(f). Use getData(CP,cl_id,''' e_opt ''')'],cl_id))
    data = []; cd(old_pwd); return
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % edb,edbs,iedb,iedbs - E.B=0 (sc,inertial)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edb') || strcmp(quantity,'edbs') || ...
    strcmp(quantity,'iedb') || strcmp(quantity,'iedbs')
  
  if strcmp(quantity,'iedb') || strcmp(quantity,'iedbs'), inert = 1;
  else, inert = 0;
  end
  
  if inert, save_file = './mEdBI.mat';
  else, save_file = './mEdB.mat';
  end
  
  if strcmp(quantity,'edb') || strcmp(quantity,'iedb')
    var_s = irf_ssub('diE?p1234',cl_id);
    varo_s = irf_ssub('E?',cl_id);
    var_b = 'diBr?';
  else
    switch probe_p
      case 12
        irf_log('proc','using p12')
        var_s = irf_ssub('diEs?p12',cl_id);
      case 32
        irf_log('proc','using p32')
        var_s = irf_ssub('diEs?p32',cl_id);
      case 34
        irf_log('proc','using p34')
        var_s = irf_ssub('diEs?p34',cl_id);
      case 420
        irf_log('proc','using p42LX')
        var_s = irf_ssub('diELXs?p42',cl_id);
      otherwise
        error(['Invalid probe pair ' num2str(probe_p)])
    end
    varo_s = irf_ssub('Es?',cl_id);
    var_b = 'diBrs?';
  end
  
  % Load resampled B
  [ok,diB,msg] = c_load(var_b,cl_id);
  if ~ok
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  % Load V if we need to do SC->Inertial transformation
  if inert
    [ok,diV,msg] = c_load('diV?',cl_id);
    if ~ok
      irf_log('load',msg)
      data = []; cd(old_pwd); return
    end
  end
  
  [ok,diE,msg] = c_load(var_s);
  if ~ok
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  if flag_usecaa_del_off && ( strcmp(quantity,'edbs') || strcmp(quantity,'iedbs') )
    % Delta offsets: remove automatic and apply CAA
    Del_caa = c_efw_delta_off(diE(1,1),cl_id);
    if ~isempty(Del_caa)
      [ok,Delauto] = c_load('D?p12p34',cl_id);
      if ~ok || isempty(Delauto)
        irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id))
      else
        diE = caa_corof_delta(diE,probe_p,Delauto,'undo');
        diE = caa_corof_delta(diE,probe_p,Del_caa,'apply');
      end
    end
  end
  
  % Save stand-deviation in the 6-th column
  if strcmp(quantity,'edbs') || strcmp(quantity,'iedbs'), diE(:,6) = diE(:,5); end
  
  % Remove wakes
  if flag_rmwake
    problems = 'wake'; %#ok<NASGU>
    signal = diE; %#ok<NASGU>
    probe = probe_p; %#ok<NASGU>
    remove_problems
    diE = res; %#ok<NODEF>
    clear res signal problems probe
  end
  
  dsiof = c_ctl(cl_id,'dsiof');
  if isempty(dsiof)
    st = diE(~isnan(diE(:,1)),1);
    if ~isempty(st), st = st(1); else, st = 0; end
    
    [ok,Ps,msg] = c_load('Ps?',cl_id);
    if ~ok, irf_log('load',msg), end
    [dsiof_def, dam_def] = c_efw_dsi_off(st,cl_id,Ps); clear st
    
    [ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
    [ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end
    
    if ok1 || ok2, irf_log('calb','Using saved DSI offsets')
    else, irf_log('calb','Using default DSI offsets')
    end
    clear dsiof_def dam_def
  else
    Ddsi = dsiof(1); Damp = dsiof(2);
    irf_log('calb','Using user specified DSI offsets')
  end
  clear dsiof
  
  diE = caa_corof_dsi(diE,Ddsi,Damp); clear Ddsi Damp
  
  irf_log('proc',['using angle limit of ' num2str(ang_limit) ' degrees'])
  [diE,ang]=irf_edb(diE,diB,ang_limit);
  diE(:,5) = ang; clear ang
  
  ii = find(abs(diE(:,5)) < ang_limit);
  if length(ii) > 1
    switch(flag_edb)
      case 0 % Ez=0, do nothing
        irf_log('proc','using Ez=0')
      case 1 % Remove points
        irf_log('proc','setting points < ang_limit to NaN')
        diE(ii,4) = diE(ii,4)*NaN;
      case 2 % Fill with fill_val
        irf_log('proc','setting points < ang_limit to 1e27')
        diE(ii,4) = ones(size(diE(ii,4)))*fill_val;
    end
  end
  
  % SC -> Inertial
  if inert
    evxb = irf_tappl(irf_cross(diB,irf_resamp(diV,diB)),'*1e-3*(-1)');
    evxb = irf_resamp(evxb,diE(:,1));
    diE(:,2:4) = diE(:,2:4) - evxb(:,2:4); clear evxb %#ok<NASGU>
    s = 'i';
  else, s = '';
  end
  
  % DSI->GSE
  if c_load('SAX?',cl_id)
    c_eval([s varo_s '=c_gse2dsi(diE(:,1:4),SAX?,-1);' s varo_s '(:,5)=diE(:,5);save_list=[save_list ''' s varo_s ' ''];'],cl_id);
  else
    irf_log('load',irf_ssub('No SAX? in mEPH. Use getData(CDB,...,cl_id,''sax'')',cl_id))
  end
  
  eval([ s 'di' varo_s '=diE;']); clear diE
  save_list=[save_list s 'di' varo_s];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % whip - Whisper present.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'whip')
  save_file = './mEFW.mat';
  
  [ok,fdm] = c_load('FDM?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No FDM?. Use getData(CDB,...,cl_id,''fdm'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  [t_s,t_e,fdm_r] = caa_efw_mode_tab(fdm, 'r');
  ii = find(fdm_r==1);
  
  if ~isempty(ii)
    % Add 1 sec before
    t_s = t_s(ii) - 1;
    t_e = t_e(ii);
    % Check for exceptionally long intervals with Whisper on.
    % 20 sec is the max reasonable duration for Whisper according to PAL.
    WHIMAX = 20;
    ii = find(t_e-t_s>WHIMAX);
    if ii
      for in=ii
        irf_log('dsrc', ['throwing away ' num2str(t_e(ii)-t_s(ii)) ...
          ' sec WHI pulse at ' epoch2iso(t_s(in),1)])
      end
      t_s(ii) = []; t_e(ii) = []; %#ok<NASGU>
      if isempty(t_s), data = []; cd(old_pwd); return, end
    end
    c_eval('WHIP?=[double(t_s)'' double(t_e)''];',cl_id);
    c_eval('save_list=[save_list '' WHIP? ''];',cl_id);
  else
    irf_log('dsrc',irf_ssub('No active whisper on C?',cl_id))
    c_eval('WHIP?=[];save_list=[save_list '' WHIP? ''];',cl_id);
  end
  clear t_s t_e fdm_r ii fdm ok
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % efwpha - sun angle from EFW FDM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'efwa')
  save_file = './mA.mat';
  
  [ok,fdm] = c_load('FDM?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No FDM?. Use getData(CDB,...,cl_id,''fdm'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  if isempty(fdm)
    irf_log('load',	irf_ssub('Empty FDM?)',cl_id))
    data = []; cd(old_pwd); return
  end
  
  % 0.002 sec is from EFW Command and telemetry Decription
  % 32 degrees is determined emprirically
  efwpha = [fdm(:,1)+0.002 fdm(:,4)*360.0/255.0-32.0];
  ii = find( efwpha(:,2) < 0 );
  efwpha(ii,2) = efwpha(ii,2) + 360.0; %#ok<NASGU>
  c_eval('EFWA?=efwpha;save_list=[save_list '' EFWA? ''];',cl_id);
  clear efwpha ii
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % bdump - Burst dump
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bdump')
  save_file = './mEFW.mat';
  
  [ok,fdm] = c_load('FDM?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No FDM?. Use getData(CDB,...,cl_id,''fdm'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  [t_s,t_e,fdm_px] = caa_efw_mode_tab(fdm, 'px');
  if isempty(fdm_px)
    irf_log('load','Insufficient FDM samples for bdump flagging.')
    data = []; cd(old_pwd); return
  end
  ii = find(fdm_px(:,1)==1 & fdm_px(:,2)==0);
  
  if ~isempty(ii)
    t_s = t_s(ii); t_e = t_e(ii);
    
    % If the gap between the two burst dumps is shorter the 2 seconds,
    % we consider it to be one interval
    for ii=1:length(t_s)-1
      if t_s(ii+1)-t_e(ii)<2.1, t_s(ii+1)=NaN; t_e(ii)=NaN; end
    end
    t_s(isnan(t_s)) = []; t_e(isnan(t_e)) = []; %#ok<NASGU>
    
    % We add one second to the end of the interval for safety
    t_e = t_e +1; %#ok<NASGU>
    c_eval('BDUMP?=[double(t_s)'' double(t_e)''];',cl_id);
    c_eval('save_list=[save_list '' BDUMP? ''];',cl_id);
  else
    irf_log('dsrc',irf_ssub('No burst dumps on C?',cl_id))
    c_eval('BDUMP?=[];save_list=[save_list '' BDUMP? ''];',cl_id);
  end
  clear t_s t_e fdm_px ii fdm ok
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sweep - Sweep present
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sweep')
  save_file = './mEFW.mat';
  
  [ok,fdm] = c_load('FDM?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No FDM?. Use getData(CDB,...,cl_id,''fdm'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  [t_s,t_e,fdm_w] = caa_efw_mode_tab(fdm, 'w');
  ii = find(fdm_w==1);
  
  [t_s_px,t_e_px,fdm_px] = caa_efw_mode_tab(fdm, 'px');
  if isempty(fdm_px)
    irf_log('load','Insufficient FDM samples for sweep flagging.')
    data = []; cd(old_pwd); return
  end
  ii_px = find(fdm_px(:,1)==1 & fdm_px(:,2)==1);
  if ~isempty(ii) || ~isempty(ii_px)
    if isempty(ii)
      % We hit sweep dump directly
      irf_log('dsrc','found loonely sweep dump')
      if length(ii_px)>1, irf_log('proc','WARNING: too many lonely sweep dumps'), end
      bdump = zeros(1,2);
      % We add one second at the start of the interval for safety
      bdump(1,1) = t_s_px(ii_px(1)) -1;
      % We add one second to the end of the interval for safety
      bdump(1,2) = t_e_px(ii_px(end)) +1; %#ok<NASGU>
    else
      bdump = zeros(length(ii),2);
      for k=1:length(ii)
        % We add one second at the start of the interval for safety
        bdump(k,1) = t_s(ii(k)) -1;
        % We look for dump of the sweep in the FDM which follows the sweep
        % or in the next one
        jj = find(t_s_px>=t_e(ii(k)) & t_s_px<t_e(ii(k))+1.1);
        % We add one second to the end of the interval for safety
        if isempty(jj)
          bdump(k,2) = t_e(ii(k)) +1;
          irf_log('dsrc','no dump after sweep')
        else, bdump(k,2) = t_e_px(jj(end)) +1;
        end
      end
    end
    c_eval('SWEEP?=bdump;save_list=[save_list '' SWEEP? ''];',cl_id);
  else
    irf_log('dsrc',irf_ssub('No sweeps on C?',cl_id))
    c_eval('SWEEP?=[];save_list=[save_list '' SWEEP? ''];',cl_id);
  end
  clear t_s t_e fdm_px ii fdm ok
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % badbias - Bad bias settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'badbias')
  save_file = './mEFW.mat';
  
  GOOD_BIAS = -130;
  
  DELTA_MINUS = 60; %#ok<NASGU>
  DELTA_PLUS = 3*60;
  
  % First we check for EFW resets which usually mean bas bias settings
  % in 2-3 minutes following the reset
  [efwtok,efwt] = c_load('EFWT?',cl_id);
  if efwtok
    c_eval(['BADBIASRESET?=[];'...
      'save_list=[save_list '' BADBIASRESET? ''];'],cl_id);
    
    ii = find(efwt(:,2)<DELTA_PLUS);
    if ~isempty(ii)
      t0 = efwt(ii(1),1) - efwt(ii(1),2);
      irf_log('proc', ['EFW reset at ' epoch2iso(t0,1)]);
      c_eval('BADBIASRESET?=[double(t0-DELTA_MINUS)'' double(t0+DELTA_PLUS)''];',cl_id);
    end
  else, irf_log('dsrc',irf_ssub('Cannot load EFWT?',cl_id))
  end
  clear t0 ii
  
  % The reason we remove 300 seconds (DELTA_MINUS) of data before a bad
  % bias is that we get rid of all the EFW resets in a clean way.
  % The 64 seconds (DELTA_PLUS) afterwards is because we have trouble deciding
  % the timing of the bias current, so we takek 2 x DSC interval.
  DELTA_MINUS = 300;
  DELTA_PLUS = 64;
  
  startTime = [];
  for pro=1:4
    [ok,ibias] = c_load(irf_ssub('IBIAS?p!',cl_id,pro));
    if ~ok
      irf_log('load',	irf_ssub('Cannot load IBIAS?p!',cl_id,pro))
      continue
    end
    if isempty(ibias)
      irf_log('load',	irf_ssub('Empty IBIAS?p!',cl_id,pro))
      continue
    end
    
    % Adjust GOOD_BIAS
    if isempty(startTime)
      if efwtok, startTime = efwt(1,1);
      else, startTime = ibias(1,1);
      end
      
      if (cl_id == 2) && (startTime>iso2epoch('2011-04-30T00:00:00Z'))
        GOOD_BIAS = -18;
      elseif startTime>iso2epoch('2006-06-16T00:00:00Z')
        GOOD_BIAS = -95;
      end
    end
    
    % Good & bad points
    ii_bad = find(ibias(:,2)>GOOD_BIAS);
    if isempty(ii_bad)
      irf_log('dsrc',irf_ssub('Bias current OK on C? p!',cl_id,pro))
      c_eval(['BADBIAS' num2str(cl_id) ...
        'p?=[];save_list=[save_list '' BADBIAS' num2str(cl_id) 'p? ''];'],pro);
    else
      ii_good = find(ibias(:,2)<=GOOD_BIAS);
      if isempty(ii_good)
        irf_log('dsrc','The whole interval has bad BIAS!!!')
        c_eval(['BADBIAS' num2str(cl_id) ...
          'p?=[double(ibias(1,1)-DELTA_MINUS)'' double(ibias(end,1)+DELTA_PLUS)''];'...
          'save_list=[save_list '' BADBIAS' num2str(cl_id) 'p? ''];'],pro);
      else
        ibias(ii_good,2) = 1;
        ibias(ii_bad,2) = 0;
        ii = irf_find_diff(ibias(:,2));
        if ibias(1,2)==0, ii = [1; ii]; end %#ok<AGROW>
        if ibias(end,2)==0, ii = [ii; length(ibias(:,2))]; end %#ok<AGROW>
        ii = reshape(ii,2,length(ii)/2);
        bb_st = ibias(ii(1,:));
        bb_et = ibias(ii(2,:));
        
        % Take care of bias sweeps. Bias is bad buring a sweep, but
        % the data will be removed as a sweep. Timing of bias current is
        % not so good, so we cannot do anything with bias current itself.
        % If we will not do this, we loose too much data due to
        % DELTA_MINUS+DELTA_PLUS
        [ok, sweep] = c_load('SWEEP?',cl_id);
        if ok
          if ~isempty(sweep)
            for in=1:size(sweep,1)
              % 40 sec is a number proposed by PAL to work around
              % bad timing
              ii = find(bb_st>sweep(in,1) & bb_st<sweep(in,2)+40);
              if ii
                irf_log('proc',['bad bias after sweep on p' ...
                  num2str(pro) ':'])
                for nn = ii
                  irf_log('proc',[epoch2iso(bb_st(nn),1) ' -- ' ...
                    epoch2iso(bb_et(nn),1)])
                end
                % Bad bias longer then 48 sec is not related to sweep.
                % 48 = 32 + 32/2
                if bb_et(ii) - bb_st(ii) >48
                  irf_log('proc','not related to sweep, leaving it')
                else
                  irf_log('proc','throwing it away')
                  bb_st(ii) = [];
                  bb_et(ii) = [];
                end
              end
            end
          end
        else
          irf_log('load',irf_ssub('cannot load SWEEP? needed for bad bias',cl_id))
        end
        if isempty(bb_st)
          irf_log('dsrc',irf_ssub('Bias current OK on C? p!',cl_id,pro))
        end
        
        res = [bb_st-DELTA_MINUS; bb_et+DELTA_PLUS]'; %#ok<NASGU>
        c_eval(['BADBIAS' num2str(cl_id) 'p?=res;'...
          'save_list=[save_list '' BADBIAS' num2str(cl_id) 'p? ''];'],pro);
        clear ii res
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % probesa - Probe saturation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'probesa')
  save_file = './mEFW.mat';
  
  if ~exist('mER.mat','file') && ~exist('mPR.mat','file')
    irf_log('load','No EFW data');data = []; cd(old_pwd); return
  end
  
  % Saturation level nA
  SA_LEVEL = 66; %#ok<NASGU>
  SA_LEV_POS = 66; % 66.0 V
  SA_LEV_POS_2 = 66; % 66.0 V
  
  % N_CONST sets the minimum number of points of constant potential
  % which we consider bad
  N_CONST = 4;
  % DT_PLUMIN is the interval by which we extend saturation
  % sintervals from each side
  DT_PLUMIN = 4;
  % Delta = .1 sec, half sampling interval for P.
  DELTA = .1;
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  [ok,pha,msg] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',msg);data = []; cd(old_pwd); return
  end
  [ok,SAX,msg] = c_load('SAX?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  [ok, sweep] = c_load('SWEEP?',cl_id);
  if ~ok
    sweep = [];
    irf_log('load',irf_ssub('cannot load SWEEP? needed for bad bias',cl_id))
  end
  ns_ops = c_ctl('get',cl_id,'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops',[c_ctl('get',5,'data_path') '/caa-control'])
    ns_ops = c_ctl('get',cl_id,'ns_ops');
  end
  [iso_t,dt] = caa_read_interval; start_time = iso2epoch(iso_t);
  
  % SAA saturation
  shadow_2 = atand(8/150); % angular width with 8 cm puck, 150 cm thin wire
  saa = atan2d(-SAX(3),SAX(1));
  saasa_se_all = []; saasa_di_all = [];
  if abs(saa-90) > shadow_2, irf_log('proc','no shadow (SAA)')
  else
    irf_log('proc',sprintf('Probe shadow (|90-%.1f(SAA)| > %.1f)',saa,shadow_2))
    tTmp = fix(start_time):0.5:(ceil(start_time)+ceil(dt));
    [phaTmp, tInts] = c_phase(tTmp', pha);
    if ~isempty(phaTmp)
      for idxTInt = 1:size(tInts,1)
        if size(tInts,1)>1
          phaTmpTmp = irf_tlim(phaTmp,tInts(idxTInt,:));
          irf_log('proc',...
            ['processing shadow : ' irf_disp_iso_range(tInts(idxTInt,:))])
        else, phaTmpTmp = phaTmp;
        end
        [spin_period,err_angle,err_angle_mean,phc_coef] = c_spin_period(phaTmpTmp,1);
        if isempty(spin_period)
          irf_log('proc','cannot find spin period'),data = [];cd(old_pwd); return
        end
        if err_angle>1 || err_angle_mean>1
          irf_log('proc','This is not yet implemented, need to do spin-by-spin');
          continue
        end
        fsamp = [];
        for probe = [12 32 34]
          [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
          if ~ok || isempty(da)
            irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
            continue
          end
          fsamp = c_efw_fsample(da,'hx');
          if ~isempty(fsamp), break; end
        end
        
        DT_MINUS_M = 2*spin_period*shadow_2/360;
        DT_PLUS_M = 1.6*DT_MINUS_M;
        DT_MINUS_L = 3*spin_period*shadow_2/360;
        DT_PLUS_L = 10/25;
        DT_MINUS_LX = 2/5; DT_PLUS_LX = 3/5;
        % Determine if we use 180Hz filter
        if cl_id == 2 && start_time>toepoch([2001 07 23 13 54 18])
          dt_lx = [DT_MINUS_M DT_PLUS_M];
          dt_hx = dt_lx;
        elseif fix(fsamp) == 450
          dt_lx = [DT_MINUS_LX DT_PLUS_LX];
          dt_hx = [DT_MINUS_M DT_PLUS_M];
        else % 10 Hz filter
          dt_lx = [DT_MINUS_LX DT_PLUS_LX];
          dt_hx = [DT_MINUS_L DT_PLUS_L];
        end
        
        nSpins=ceil(dt/spin_period)+2; spinN = (1:nSpins)-2;
        probes = [1 3 2 4]; saasa = zeros(nSpins,4);
        % Shadow times
        for iProbe=1:4
          saasa(:,probes(iProbe)) = ...
            (pi/4 + pi/2*(iProbe-1) + spinN*2*pi - phc_coef(2))/phc_coef(1)...
            + phaTmp(1,1);
        end
        % LX Shadows
        saasa_se = zeros(nSpins,8);
        for iProbe=1:4
          saasa_se(:,iProbe*2-1) = saasa(:,iProbe) - dt_lx(1);
          saasa_se(:,iProbe*2) = saasa(:,iProbe) + dt_lx(2);
        end
        % HX Shadows
        probes = [12 34 32 42]; saasa_di = zeros(nSpins*2,8);
        for iProbe=1:4
          pA = fix(probes(iProbe)/10); pB = probes(iProbe) - pA*10;
          tmpT = sort([saasa(:,pA); saasa(:,pB)]);
          saasa_di(:,iProbe*2-1) = tmpT - dt_hx(1);
          saasa_di(:,iProbe*2) = tmpT + dt_hx(2);
        end
        clear tmpT
        saasa_se_all = [saasa_se_all; saasa_se]; %#ok<AGROW>
        saasa_di_all = [saasa_di_all; saasa_di]; %#ok<AGROW>
      end
    end
  end % SAA saturation
  c_eval(['SAASASE?=saasa_se_all;SAASADI?=saasa_di_all;'...
    'save_list=[save_list '' SAASASE? '' '' SAASADI? '' ];'],cl_id);
  
  c_eval('sa_int_p?=[];')
  for pro=1:4
    [ok,p] = c_load(irf_ssub('P10Hz?p!',cl_id,pro));
    if ~isempty(start_time)
      if (pro==3 && start_time>toepoch([2001 07 23 00 00 00]) && cl_id==2) || ...
          (pro==1 && start_time>toepoch([2002 07 29 09 06 59]) && cl_id==3) || ...
          (pro==1 && start_time>toepoch([2002 12 28 03 02 57]) && cl_id==1) || ...
          (pro==1 && start_time>toepoch([2007 05 13 03 23 30]) && cl_id==2) || ...
          (pro==4 && start_time>toepoch([2009 10 14 03 23 30]) && cl_id==1) || ...
          (pro==4 && start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]) && cl_id==1)
        sa_int = []; %#ok<NASGU>
        if ~isempty(ns_ops)
          probestr=['p' num2str(pro)];
          sa_int = caa_get_ns_ops_int(start_time,dt,ns_ops,['no_' probestr]);
          if ~isempty(sa_int), irf_log('proc',['Found no_' probestr ' in NS_OPS']), end
        end
        
        irf_log('dsrc',...
          irf_ssub('Using fake PROBELD?p!',cl_id,pro));
        c_eval(['p?=[];sa_int_p?=sa_int;PROBELD' num2str(cl_id) ...
          'p?=[];save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
        continue
      end
    end
    if ~ok
      irf_log('load',	irf_ssub('Cannot load P10Hz?p!',cl_id,pro))
      c_eval('p?=[];',pro)
      continue
    end
    if isempty(p)
      irf_log('load',	irf_ssub('Empty P10Hz?p!',cl_id,pro))
    elseif isempty(start_time), start_time = p(1,1);
    end
    
    % Read in ns_ops
    if ~isempty(ns_ops)
      sa_int = caa_get_ns_ops_int(start_time,dt,ns_ops,['no_p' num2str(pro)]);
      if ~isempty(sa_int)
        irf_log('proc',['Found no_p' num2str(pro) ' in NS_OPS'])
        c_eval('sa_int_p?=sa_int;',pro)
        p = caa_rm_blankt(p,sa_int);
        p(isnan(p(:,2)),:) = [];
      end
    end
    
    % Clean up after evil sweeps
    if sweep
      irf_log('proc',['blanking sweeps on p' num2str(pro)])
      p = caa_rm_blankt(p,sweep);
      p(isnan(p(:,2)),:) = [];
    end
    
    c_eval('p?=p;',pro)
  end
  
  % Points below SA_LEVEL should be flagged as PROBELD.
  % This is excluded from E (NaN + flag), but not from P
  % as they still contain valuable physical information.
  % Points with positive potential should be flagged HBIASSA,
  % and flagged in E (but not excluded).
  
  % Bad points are points below SA_LEVEL
  c_eval(['if isempty(p?),ii_bad?=[];ii_god?=[];else,'...
    'ii_bad?=find(p?(:,2)<-SA_LEVEL); ii_god? = find(p?(:,2)>=-SA_LEVEL);end'])
  
  for pro=1:4
    c_eval('p=p?;ii_bad=ii_bad?;ii_god=ii_god?;',pro)
    if isempty(p)
      c_eval(['if ~isempty(sa_int_p?),PROBESA' num2str(cl_id) 'p?=sa_int_p?;',...
        'save_list=[save_list '' PROBESA' num2str(cl_id) 'p? '']; end'],pro);
      continue
    end
    
    if isempty(ii_bad) %#ok<NODEF>
      c_eval(['PROBELD' num2str(cl_id) ...
        'p?=[];save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
      continue
    end
    
    % We check if there is a saturation at the same time for other probes
    if pro==1,     oprop = [2 3 4];
    elseif pro==2, oprop = [1 3 4];
    elseif pro==3, oprop = [4 1 2];
    elseif pro==4, oprop = [3 1 2];
    end
    
    found = 0;
    for opro = oprop
      op = []; % silence M-lint
      c_eval('op=p?;oii_bad=ii_bad?;oii_god=ii_god?;',opro)
      if ~isempty(op)
        if isempty(oii_bad), found = 1; break
        else
          ii1 = irf_find_comm_idx(p(ii_bad,:),op(oii_bad,:));
          if isempty(ii1), found = 1; break, end
        end
      end
    end
    if found
      irf_log('proc',['found no corresponding saturation to p' ...
        num2str(pro) ' on p' num2str(opro)])
      c_eval(['PROBELD' num2str(cl_id) ...
        'p?=[];save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
      continue
    end
    
    if isempty(ii_god) %#ok<NODEF>
      c_eval(['PROBELD' num2str(cl_id) ...
        'p?=[double(p(1,1))'' double(p(end,1))''];'...
        'PROBESA' num2str(cl_id) 'p?=sa_int_p?;'...
        'save_list=[save_list '' PROBELD' num2str(cl_id) ...
        'p? PROBESA' num2str(cl_id) 'p? ''];'],pro);
    else
      p_tmp = p;
      p_tmp(ii_god,2) = 1;
      p_tmp(ii_bad,2) = 0;
      ii = irf_find_diff(p_tmp(:,2));
      if p_tmp(1,2)==0, ii = [1; ii]; end %#ok<AGROW>
      if p_tmp(end,2)==0, ii = [ii; length(p_tmp(:,2))]; end %#ok<AGROW>
      ii = reshape(ii,2,length(ii)/2);
      res = [p_tmp(ii(1,:))-DELTA; p_tmp(ii(2,:))+DELTA]'; %#ok<NASGU>
      c_eval(['PROBELD' num2str(cl_id) 'p?=res;'...
        'save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
      
      clear res ii
    end
    c_eval('p?=p;',pro)
    clear ii_god ii_bad p
  end
  
  % Points with constant potential (latched probe) or potential > 10 V should be flagged
  % PROBESA, and excluded from E and P.
  for pro=1:4
    c_eval(['p=p?;if ~isempty(p), ldsa = PROBELD' num2str(cl_id) 'p?; end'],pro)
    if isempty(p), continue, end
    
    % Leave only good points for further exploration
    if ~isempty(ldsa)
      irf_log('proc',['blanking LD saturation on p' num2str(pro)])
      p = caa_rm_blankt(p,ldsa);
      p(isnan(p(:,2)),:) = [];
      if isempty(p), continue, end
    end
    
    % Bad points are points with positive and/or constant potential
    ii_bad = find( p(:,2) >=SA_LEV_POS_2 );
    ii_god = find( p(:,2) < SA_LEV_POS_2 );
    if isempty(ii_god)
      c_eval(['PROBESA' num2str(cl_id) ...
        'p?=[sa_int_p?; double(p(1,1))'' double(p(end,1))''];'...
        'save_list=[save_list '' PROBESA' num2str(cl_id) 'p? ''];'],pro);
    else
      dd = diff(p(ii_god,2));
      p(ii_god,2) = 1;
      
      % Check for constant P, means probe is in a strange state
      ii = find(dd==0);
      if length(ii)>1
        % at least three consequetive points are the same
        kk = find(ii(1:end-1)-ii(2:end)==-1);
        if ~isempty(kk)
          jj = 1;
          while jj<=length(kk)
            bad_i = find(ii-ii(kk(jj))-(1:length(ii))'+kk(jj)==0);
            if length(bad_i)>N_CONST
              p([ii_god(ii(bad_i)); ii_god(ii(bad_i(end)))+1],2) = 0;
            end
            if isempty(bad_i), jj = jj + 1;
            else
              ll = find(kk>bad_i(end));
              if isempty(ll), break
              else, jj = ll(1);
              end
            end
          end
        end
      end
      p(ii_bad,2) = 0; %#ok<FNDSB>
      ii = irf_find_diff(p(:,2));
      if p(1,2)==0, ii = [1; ii]; end %#ok<AGROW>
      if p(end,2)==0, ii = [ii; length(p(:,2))]; end %#ok<AGROW>
      ii = reshape(ii,2,length(ii)/2);
      
      % We add DT_PLUMIN sec on each side and check for overlapping intervals
      res = [p(ii(1,:))-DT_PLUMIN; p(ii(2,:)-1)+DT_PLUMIN]';
      if length(ii(1,:))>1
        pos = 1;
        while 1
          if res(pos,2)+DT_PLUMIN>=res(pos+1,1)
            res(pos,2) = res(pos+1,2);
            res(pos+1,:) = [];
          end
          pos = pos + 1;
          if pos>=size(res,1), break, end
        end
      end
      c_eval(['PROBESA' num2str(cl_id) 'p?=[sa_int_p?; res];'...
        'save_list=[save_list '' PROBESA' num2str(cl_id) 'p? ''];'],pro);
      clear ii res
    end
  end
  
  % Points with positive potential should be flagged HBIASSA.
  % Note that this will be something like HBIASSA3p4, rather
  % than the HBIASSA3p34 that would be written by the hbiassa routines.
  for pro=1:4
    c_eval('p=p?;',pro)
    if isempty(p)
      c_eval(['HBIASSA' num2str(cl_id) 'p?=[];'...
        'save_list=[save_list '' HBIASSA' num2str(cl_id) 'p? ''];'],pro);
      continue
    end
    
    % Bad points are points with positive potential
    ii_bad = find( p(:,2) >=SA_LEV_POS );
    ii_god = find( p(:,2) < SA_LEV_POS );
    if isempty(ii_god)
      c_eval(['HBIASSA' num2str(cl_id) ...
        'p?=[double(p(1,1)) double(p(end,1))];'...
        'save_list=[save_list '' HBIASSA' num2str(cl_id) 'p? ''];'],pro);
    else
      p(ii_god,2) = 1;
      p(ii_bad,2) = 0; %#ok<FNDSB>
      ii = irf_find_diff(p(:,2));
      if p(1,2)==0, ii = [1; ii]; end %#ok<AGROW>
      if p(end,2)==0, ii = [ii; length(p(:,2))]; end %#ok<AGROW>
      ii = reshape(ii,2,length(ii)/2);
      
      % We add DT_PLUMIN sec on each side and check for overlapping intervals
      res = [p(ii(1,:))-DT_PLUMIN; p(ii(2,:)-1)+DT_PLUMIN]';
      if length(ii(1,:))>1
        pos = 1;
        while 1
          if res(pos,2)+DT_PLUMIN>=res(pos+1,1)
            res(pos,2) = res(pos+1,2);
            res(pos+1,:) = [];
          end
          pos = pos + 1;
          if pos>=size(res,1), break, end
        end
      end
      c_eval(['HBIASSA' num2str(cl_id) 'p?=res;'...
        'save_list=[save_list '' HBIASSA' num2str(cl_id) 'p? ''];'],pro);
      clear ii res
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % HBIASSA - saturation due to high bias current
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'hbiassa')
  save_file = './mEFW.mat';
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  [ok,pha,msg] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',msg);
    data = []; cd(old_pwd); return
  end
  
  evxb = [];
  [ok,r,msg] = c_load('R?',cl_id);
  if ~ok || isempty(r), irf_log('load',msg);
  else
    r = irf_abs(r); r = r(:,[1 5]);
    if ~isempty(r) && any( r(:,2)<4*6371 )
      % Load resampled B
      [ok,diB] = c_load('diBrs?',cl_id);
      if ~ok
        irf_log('load',...
          irf_ssub('No diBs? in mBr. Use getData(CP,cl_id,''brs'')',cl_id))
      else
        [ok,diV] = c_load('diV?',cl_id);
        if ~ok
          irf_log('load',...
            irf_ssub('No diV? in mR. Use getData(CDB,...,cl_id,''v'')',cl_id))
        end
        evxb = irf_tappl(irf_cross(diB,irf_resamp(diV,diB)),'*1e-3*(-1)');
      end
    end
  end
  
  p12_ok = 0; sf_probe = []; np = [];
  for probe = [12 32 34]
    if probe == 32 && p12_ok, continue, end
    [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
    if ~ok || isempty(da)
      irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
      continue
    end
    if probe == 12, p12_ok = 1; end
    
    fsamp = c_efw_fsample(da,'hx');
    
    problems = 'reset|bbias|probeld|sweep|bdump|saasa|nsops';
    nsops_errlist = [caa_str2errid('bad_bias') caa_str2errid('bad_hx')];%#ok<NASGU>
    
    % Always remove Whisper when we use 180Hz filter
    if (fsamp == 450) || ...
        ( cl_id == 2 && da(1,1)>toepoch([2001 07 23 13 54 18]) ) || ...
        ( flag_rmwhip && flag_rmwhip_force )
      problems = [problems '|whip'];  %#ok<AGROW,NASGU>
    end
    signal = da; %#ok<NASGU>
    remove_problems
    da = res; %#ok<NODEF>
    clear res signal problems
    
    % Check if we have at least 1 spin of data left
    if length(find(~isnan(da(:,2)))) < 4*fsamp
      irf_log('proc',irf_ssub('No p? data after removals',probe))
      continue
    end
    
    if ~isempty(evxb)
      % Remove the V_SC x B field before chechikng for saturation
      switch probe
        case 12, dPhi = 3*pi/4;
        case 32, dPhi = pi/2;
        case 42, dPhi = pi;
        case 34, dPhi = pi/4; % angles when phase =0
        otherwise, error('bad probe pair')
      end
      a_evxb = c_phase(da(:,1),pha);
      if size(a_evxb,1) ~= size(da,1)
        [iia,iid]=irf_find_comm_idx(a_evxb(:,1),da(:,1));
      else
        iia = 1:size(a_evxb,1); iid = iia;
      end
      evxb_r = irf_resamp(evxb,da(iid,1));
      evxb_spin = evxb_r(:,1:2);
      Phase = a_evxb(iia,2)/180*pi + dPhi;
      evxb_spin(:,2)=evxb_r(:,2).*cos(Phase) - evxb_r(:,3).*sin(Phase);
      da(iid,2) = da(iid,2) - evxb_spin(:,2);
      irf_log('proc','Removing Vsc x B at perigee')
    end
    
    [HBIASSA,wakedesc] = c_efw_hbias_satur(da,probe,pha);  %#ok<ASGLU>
    
    % Below 2RE we can have real large fields
    if ~isempty(HBIASSA) && ~isempty(r)
      rs = irf_resamp(r,da(:,1));
      ii_ok = ones(size(HBIASSA,1),1);
      for jInt = 1:size(HBIASSA,1)
        rs_i = irf_tlim(rs,HBIASSA(jInt,:));
        if any(rs_i(:,2)<2*6371)
          irf_log('proc',['Disregarding HBIASSA below 2 RE at ' irf_disp_iso_range(HBIASSA(jInt,:),1)])
          ii_ok(jInt) = 0;
        end
      end
      HBIASSA(ii_ok==0,:) = [];
    end
    
    % Append single-ended high bias saturation intervals from probesa
    for single_probe=[mod(probe,10) floor(probe/10)]
      sp_name=irf_ssub('HBIASSA?p!',cl_id,single_probe);
      [ok, sp_hbiassa] = c_load(sp_name);
      if ok && ~isempty(sp_hbiassa)
        irf_log('proc',['Appending ' sp_name ' to ' irf_ssub('HBIASSA?p!',cl_id,probe)])
        HBIASSA=[HBIASSA' sp_hbiassa']';
      end
    end
    
    % Check which of the probe pairs is more affected by saturation
    if ~isempty(HBIASSA)
      irf_log('proc','blanking HB saturation')
      da = caa_rm_blankt(da,HBIASSA);
    end
    np_tmp = length( da( ~isnan(da(:,2)) ,2 ) );
    if np_tmp > 0
      if ~isempty(np) && np_tmp == np
        % Same time lines: no preference for a particular pair
        sf_probe = [];
      elseif probe == 32
        if isempty(np) || np==0 || np_tmp > (np*1.2)
          sf_probe = probe;
        end
      else
        if isempty(np) || np==0 || np_tmp*1.2 > np
          sf_probe = probe;
        end
      end
    end
    np = np_tmp;
    clear np_tmp da
    
    if isempty(HBIASSA)
      eval(irf_ssub('HBIASSA?p!=[];save_list=[save_list ''HBIASSA?p! ''];',cl_id,probe));
    else
      eval(irf_ssub('HBIASSA?p!=HBIASSA;HBSATDSC?p!=wakedesc; save_list=[save_list ''HBIASSA?p! HBSATDSC?p! ''];',cl_id,probe));
    end
  end
  
  % Save probe to be used for spin resolution data
  if ~isempty(sf_probe)
    irf_log('calb',irf_ssub('Will use p? for spin res data',sf_probe))
    caa_sfit_probe(cl_id,sf_probe);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % rawspec - Spectrum of raw EFW signal
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'rawspec')
  save_file = './mEFW.mat';
  
  % Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
  [ok,pha] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',...
      irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id));
    data = []; cd(old_pwd); return
  else
    %Find zero phase
    ii = find(pha(:,2)==0);
    if isempty(ii)
      irf_log('proc',...
        irf_ssub('Cannot find zero phase in Atwo?',cl_id))
      data = []; cd(old_pwd); return
    end
    tpha0 = pha(ii(1),1);
  end
  
  p12 = 12; e12 = []; e34 =[];
  n_ok = 0; flag_lx = 0;
  tpharef = [];
  corrected_raw_data_p12 = 1;
  corrected_raw_data_p34 = 1;
  for probe = [12 32 34]
    [ok,da] = c_load(irf_ssub('wcE?p!',cl_id,probe));
    if ~ok || isempty(da)
      irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
      [ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
      if ~ok || isempty(da)
        irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
        continue
      end
      irf_log('load','using raw (not corrected) data')
      if probe==34, corrected_raw_data_p34 = 0;
      else, corrected_raw_data_p12 = 0;
      end
    end
    n_ok = n_ok + 1;
    
    if size(da,1) > length(tpharef), tpharef = da(:,1); end
    
    if probe==32
      p12 = 32;
      e12 = da;
    else, c_eval('e?=da;',probe)
    end
    clear ok da
  end
  if ~n_ok % Try with LX data p42
    [ok,p4] = c_load(irf_ssub('P10Hz?p!',cl_id,4));
    if ok && ~isempty(p4)
      [ok,p2] = c_load(irf_ssub('P10Hz?p!',cl_id,2));
      if ok && ~isempty(p2)
        p_sep = .066;
        if size(p2,1) ~=  size(p4,1)
          [ii2,ii4]=irf_find_comm_idx(p2,p4);
          e12(:,1) = p2(ii2,1); e12(:,2) = ( p2(ii2,2) - p4(ii4,2) )/p_sep;
        else
          e12(:,1) = p2(:,1); e12(:,2) = ( p2(:,2) - p4(:,2) )/p_sep;
        end
        p12 = 42; flag_lx = 1; n_ok = n_ok + 1;
        irf_log('proc','Using p42 LX')
        if size(e12,1) > length(tpharef), tpharef = e12(:,1); end
      else
        irf_log('load', irf_ssub('No/empty P10Hz?p!',cl_id,2));
      end
      clear p2 p4 p_sep
    else
      irf_log('load', irf_ssub('No/empty P10Hz?p!',cl_id,4));
    end
  end
  
  if ~n_ok, data = []; cd(old_pwd), return, end
  
  %Compute spin period
  ph = c_phase(tpharef,pha);
  if ~isempty(ph), ph(isnan(ph(:,2)),:) = []; end
  if isempty(ph)
    irf_log('proc','Phase is empty')
    data = []; cd(old_pwd), return
  end
  ph(:,1) = ph(:,1) - ph(1,1);
  phc = unwrap( ph(1:2,2)/180*pi );
  phc_coef = polyfit( ph(1:2,1),phc,1 );
  for k=1:floor( log10(length(ph(:,1))) )
    ii = 10^k;
    dp = ph(ii,2) - mod( polyval(phc_coef,ph(ii,1))*180/pi,360 );
    dpm = [dp dp-360 dp+360];
    dph = dpm(abs(dpm)<180);
    phc_coef(1) = phc_coef(1) + dph*pi/180/ph(ii,1);
  end
  spin_p = 2*pi/phc_coef(1);
  if spin_p > 4.4 || spin_p < 3.6
    irf_log('proc','using spin period == 4 sec')
    spin_p = 4.0;
  else, irf_log('proc',['spin period is ' num2str(spin_p) ' sec'])
  end
  
  for pr = [12 34]
    
    if pr==12
      tt = e12;
      probe = p12;
      corrected_raw_data = corrected_raw_data_p12;
    else
      tt = e34;
      probe = 34;
      corrected_raw_data = corrected_raw_data_p34;
    end
    if isempty(tt)
      irf_log('load',sprintf('No raw spectrum for C%d p%d',cl_id,probe))
      continue
    end
    
    if corrected_raw_data, ss = 'c';
    else, ss = '';
    end
    irf_log('proc',sprintf('Raw spectrum w%sE%dp%d -> RSPEC%dp%d',...
      ss,cl_id,probe,cl_id,probe))
    
    if flag_lx, fsamp = c_efw_fsample(tt,'lx');
    else, fsamp = c_efw_fsample(tt,'hx');
    end
    if ~fsamp, error('no sampling frequency'),end
    
    problems = 'reset|bbias|probesa|probeld|sweep|bdump|nsops';
    nsops_errlist = caa_str2errid('bad_bias');
    if ~flag_lx, nsops_errlist = [nsops_errlist caa_str2errid('bad_hx')]; end %#ok<AGROW,NASGU>
    if flag_rmwhip, problems = [problems '|whip']; end %#ok<NASGU,AGROW>
    signal = tt; %#ok<NASGU>
    remove_problems
    tt = res; %#ok<NODEF>
    clear res signal problems
    
    % Check if we have at least 1 spin of data left
    if length(find(~isnan(tt(:,2)))) < 4*fsamp
      irf_log('proc',irf_ssub('No p? data after removals',probe))
      continue
    end
    
    % Time when a particular probe pair was directed towards the Sun
    if probe == 12, tpha0probe = tpha0 - 3*spin_p/8;
    elseif probe == 32, tpha0probe = tpha0 - spin_p/4;
    elseif probe == 34, tpha0probe = tpha0 - spin_p/8;
    elseif probe == 42, tpha0probe = tpha0 - spin_p/2;
    end
    
    tstart = tpha0probe + fix( (tt(1,1) - tpha0probe)/spin_p )*spin_p;
    tend = tstart + ceil( (tt(end,1) - tstart)/spin_p )*spin_p;
    n = floor((tend - tstart)/spin_p) + 1;
    if n<0, error(['N is negative: ' num2str(n)]), end
    rspec = zeros(n,11)*NaN;
    
    % N_EMPTY .75 means that we use only spins with more then 75% points.
    N_EMPTY = .9;
    n_gap = 0;
    
    % Do it:
    for i=1:n
      t0 = tstart + (i-1)*spin_p;
      rspec(i,1) = t0;
      eind = find((tt(:,1) >= t0) & (tt(:,1) < t0+spin_p));
      
      % Check for data gaps inside one spin.
      if fsamp>0 && length(eind)<N_EMPTY*4*fsamp, eind = []; end
      
      % Need to check if we have any data to fit.
      if ~isempty(eind)
        %Handle EFW data
        IFFT_Diff = ifft( tt(eind,2) - mean(tt(eind,2)) );
        
        AmpSin = 2*imag(IFFT_Diff);
        AmpCos = 2*real(IFFT_Diff);
        
        rspec(i,2)  =  AmpCos(1+1); % 1 omega
        rspec(i,3)  = -AmpSin(1+1); % 2 omega
        rspec(i,4)  =  AmpCos(2+1); % 3 omega
        rspec(i,5)  = -AmpSin(2+1); % 4 omega
        rspec(i,6)  =  AmpCos(3+1); % 5 omega
        rspec(i,7)  = -AmpSin(3+1); % 6 omega
        rspec(i,8)  =  AmpCos(4+1); % 7 omega
        rspec(i,9)  = -AmpSin(4+1); % 8 omega
        rspec(i,10) =  AmpCos(5+1);
        rspec(i,11) = -AmpSin(5+1);
      else, n_gap = n_gap + 1;
      end
    end
    irf_log('proc',sprintf('%d spins processed, %d gaps found',n,n_gap))
    eval(irf_ssub('RSPEC?p!=rspec;save_list=[save_list ''RSPEC?p! ''];',cl_id,probe));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % wake (wakes in lobe and plasmasphere)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'wake')
  save_file = './mEFW.mat';
  
  if check_caa_sh_interval
    if ~exist('./.caa_ms_interval','file')
      irf_log('proc','Outside magnetosphere. No lobe wake removal performed.')
      data = []; cd(old_pwd), return
    end
  end
  
  [spinFits,msg] = caa_sfit_load(cl_id);
  if isempty(spinFits), irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,diBrs,msg] = c_load('diBrs?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,Ps,msg] = c_load('Ps?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,R,msg] = c_load('R?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,SAX,msg] = c_load('SAX?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,diV,msg] = c_load('diV?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  [ok,diEDI,msg] = c_load('diEDI?',cl_id);
  if ~ok, irf_log('load',msg), diEDI = []; end
  
  % Correct for DSI offsets
  [Ddsi,Damp] = c_efw_dsi_off(spinFits.diEs(1,1),cl_id,Ps);
  spinFits.diEs = caa_corof_dsi(spinFits.diEs,Ddsi,Damp);
  
  % Detect lobe and plasmaspheric wakes
  pswake = c_efw_corrot(cl_id,spinFits.diEs,diBrs,Ps,R,SAX,diV);
  spinFits.diEs = caa_rm_blankt(spinFits.diEs,pswake);
  [lowake,dEx] = c_efw_lobewake(cl_id,spinFits.diEs,diBrs,Ps,R,diEDI,1);
  
  % Detect nonsinusoidal wakes if x_GSE>0
  nonsinwake = []; %#ok<NASGU>
  if any(R(:,2) > 0)
    % Load required data
    [ok,pha] = c_load('Atwo?',cl_id);
    if ~ok, irf_log('load',msg), pha=[]; end
    sfpp=spinFits.probePair;
    if sfpp > 100
      sfpp = sfpp/10;
    end
    [ok,da] = c_load(irf_ssub('wE?p!',cl_id,sfpp));
    if ~ok, irf_log('load',msg), da=[]; end
    
    if ~isempty(da) && ~isempty(pha)
      % Blank intervals where GSE X position is > 0 or wakes already
      % detected
      if any(R(:,2)<0)
        x=irf_resamp(R(:,1:2),da);
        da(x<0,:)=NaN;
      end
      da=caa_rm_blankt(da,pswake);
      da=caa_rm_blankt(da,lowake);
      
      % Find nonsinusoidal wakes
      if any(isfinite(da(:,2)))
        nonsinwake = c_efw_nonsinwake(cl_id,spinFits.probePair,pha,da); %#ok<NASGU>
      end
    end
  end
  
  if ~isempty(dEx)
    cmd = 'if length(Ddsi)==1,DdsiX?=Ddsi+dEx;else DdsiX?=Ddsi;DdsiX?(:,2)=DdsiX?(:,2)+dEx; end; save mXTRA DdsiX?';
    if exist('./mXTRA.mat','file'), cmd = [cmd ' -append']; end
    c_eval(cmd,cl_id), clear cmd
  end
  eval(irf_ssub(...
    'PSWAKE?p!=pswake;LOWAKE?p!=lowake;NONSINWAKE?p!=nonsinwake;save_list=[save_list ''PSWAKE?p! LOWAKE?p! NONSINWAKE?p!''];',...
    cl_id,probe_p));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % edi (sc)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edi')
  
  save_file = './mEDI.mat';
  
  if check_caa_sh_interval
    if ~exist('./.caa_ms_interval','file')
      irf_log('proc','Outside magnetosphere. No edi processing performed.')
      data = []; cd(old_pwd), return
    end
  end
  
  var_s = 'iEDI?';
  varo_s = 'EDI?';
  
  % Load BPP. We use BPP for EDI as it must be a rather approximation
  [ok,B] = c_load('BPP?',cl_id);
  if ~ok
    [ok,B] = c_load('B?',cl_id);
    if ~ok
      irf_log('load',...
        irf_ssub('No B? and BPP?. Use getData(CDB,...,cl_id,''b'')',cl_id))
      data = []; cd(old_pwd); return
    end
  end
  
  % Load V if we need to do SC->Inertial transformation
  [ok,V,msg] = c_load('V?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  % Load E EDI (inertial)
  [ok,E,msg] = c_load(var_s,cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  % SC -> Inertial
  B = irf_resamp(B,E);
  evxb = irf_tappl(irf_cross(B,irf_resamp(V,B)),'*1e-3*(-1)');
  E(:,2:4) = E(:,2:4) + evxb(:,2:4); clear evxb %#ok<NASGU>
  
  % GSE->DSI
  if c_load('SAX?',cl_id)
    c_eval(['di' varo_s '=c_gse2dsi(E(:,1:4),SAX?);save_list=[save_list '' di' varo_s ' ''];'],cl_id);
  else
    irf_log('load',irf_ssub('No SAX? in mEPH. Use getData(CDB,...,cl_id,''sax'')',cl_id))
  end
  
  c_eval([varo_s '= E;'],cl_id); clear E
  save_list=[save_list irf_ssub(varo_s,cl_id) ' '];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Vedb,Vedbs = ExB with E.B=0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vedb') || strcmp(quantity,'vedbs')
  save_file = './mEdB.mat';
  
  if strcmp(quantity,'vedb')
    var_s = 'diE?'; e_opt = 'edb';
    varo_s = 'VExB?';
    var_b = 'diBr?';
  else
    var_s = 'diEs?'; e_opt = 'edbs';
    varo_s = 'VExBs?';
    var_b = 'diBrs?';
  end
  
  % Load resampled B
  [ok,diB,msg] = c_load(var_b,cl_id); %#ok<ASGLU>
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  
  % Load data and calculate ExB
  if c_load(var_s,cl_id)
    c_eval(['di' varo_s '=irf_e_vxb(' var_s '(:,1:4),diB,-1);di' varo_s '(:,5)=' var_s '(:,5);'],cl_id)
  else
    irf_log('load',...
      irf_ssub(['No ' var_s ' in mEdB. Use getData(CP,cl_id,''' e_opt ''')'],cl_id))
    data = []; cd(old_pwd); return
  end
  
  save_list=[save_list 'di' irf_ssub(varo_s,cl_id) ' '];
  
  % DSI->GSE
  if c_load('SAX?',cl_id)
    eval(irf_ssub([varo_s '=c_gse2dsi(di' varo_s '(:,1:4),SAX?,-1);' varo_s '(:,5)=di' varo_s '(:,5);save_list=[save_list ''' varo_s ' ''];'],cl_id));
  else
    irf_log('load',irf_ssub('No SAX? in mEPH. Use getData(CDB,...,cl_id,''sax'')',cl_id))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B resampled to E and Es
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'br') || strcmp(quantity,'brs')
  save_file = './mBr.mat';
  
  if check_caa_sh_interval && strcmp(quantity,'brs')
    if ~exist('./.caa_ms_interval','file')
      irf_log('proc','Outside magnetosphere. No brs resampling performed.')
      data = []; cd(old_pwd), return
    end
  end
  
  if strcmp(quantity,'br')
    var_b = 'Br?';
    [ok,E_tmp,msg] = c_load('diE?p1234',cl_id);
    if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
    e_sf = c_efw_fsample(E_tmp,'hx');
  else
    var_b = 'Brs?';
    spinFits = caa_sfit_load(cl_id);
    if isempty(spinFits)
      irf_log('load','Cannot load spinFits. ')
      data = []; cd(old_pwd); return
    end
    E_tmp = spinFits.diEs;
    e_sf = .25;
  end
  
  % Load B GSE, as it is level 0 FGM data for us.
  B_tmp = c_load('B?',cl_id,'var');
  
  % Check for data coverage
  % In the current approach we compute it from the sampling frequency of B.
  dt = E_tmp(end,1) - E_tmp(1,1);
  if ~isempty(B_tmp)
    Binfo = 'FR'; %#ok<NASGU>
    bad_coverage = 0;
    cover = 0;
    B_tmp = irf_tlim(B_tmp,E_tmp(1,1) + [-.5/e_sf (dt+.5/e_sf)]);
    if size(B_tmp,1) < 2, bad_coverage = 1;
    else
      fgm_sf = 1/(B_tmp(2,1)-B_tmp(1,1));
      del_f = 1.5;
      if (fgm_sf > 22.5 - del_f) && (fgm_sf < 22.5 + del_f), fgm_sf = 22.5;
      elseif (fgm_sf > 67.5 - del_f) && (fgm_sf < 67.5 + del_f), fgm_sf = 67.5;
      else, irf_log('proc','cannot guess sampling frequency for B')
      end
      cover = length(B_tmp(:,1))/((dt+1/e_sf)*fgm_sf);
      % We allow for 10% of data gaps. (should we??)
      if cover < .9, bad_coverage = 1; end
    end
  else, bad_coverage = 1; cover = 0;
  end
  
  % Try to use BPP as a backup
  if bad_coverage
    BPP_tmp = c_load('BPP?',cl_id,'var');
    if length(BPP_tmp) < 8
      % Use FR data if there is any (cover != 0)
      if cover==0
        irf_log('load','Cannot load B. Please load B FGM or B PP.')
        data = []; cd(old_pwd); return
      end
    else
      BPP_tmp = irf_tlim(BPP_tmp,E_tmp(1,1) + [-.5/e_sf (dt+.5/e_sf)]);
      if size(BPP_tmp,1) < 2
        irf_log('load','Cannot find any usefull B data. Please load B FGM or B PP.')
        data = []; cd(old_pwd); return
      end
      
      fgm_sf = 1/(BPP_tmp(2,1)-BPP_tmp(1,1));
      del_f = .1;
      if (fgm_sf > .25 - del_f) && (fgm_sf < .25 + del_f), fgm_sf = .25;
      else, irf_log('proc','cannot guess sampling frequency for B PP')
      end
      cover_pp = length(BPP_tmp(:,1))/((dt+1/e_sf)*fgm_sf);
      
      % If there is more PP data, then use it.
      % Take .99 to avoid marginal effects.
      if .99*cover_pp > cover
        B_tmp = BPP_tmp;
        Binfo = 'PP'; %#ok<NASGU>
        irf_log('proc','Using B PP to calculate Br')
      else, irf_log('proc',sprintf('Use B has %2.2f%% coverage',cover*100))
      end
    end
  end
  
  % Resample the data
  Br = irf_resamp(B_tmp,E_tmp,'fsample',e_sf); %#ok<NASGU>
  c_eval([ var_b '=Br;' var_b '_info=Binfo;save_list=[save_list ''' var_b ' '' '' ' var_b '_info '' ];'],cl_id)
  
  % DSI->GSE
  if c_load('SAX?',cl_id)
    eval(irf_ssub(['di' var_b '=c_gse2dsi(Br,SAX?); di' var_b '_info=Binfo;save_list=[save_list ''di' var_b ' '' '' di' var_b '_info '' ];'],cl_id));
  else
    irf_log('load',irf_ssub('No SAX? in mEPH. Use getData(CDB,...,cl_id,''sax'')',cl_id))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % P averaged from several probes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'p') || strcmp(quantity,'pburst')
  if strcmp(quantity,'pburst')
    do_burst = 1; save_file = './mEFWburst.mat';
    param={'180Hz','4kHz','32kHz'};
  else, do_burst = 0; save_file = './mP.mat'; param={'10Hz'};
  end
  P = struct('p1',[],'p2',[],'p3',[],'p4',[]);
  n_ok = 0; loaded = 0;
  for k=1:length(param)
    for iProbe=1:4
      [ok,p] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,iProbe));
      if ~ok || isempty(p)
        irf_log('load', irf_ssub(['No/empty P' param{k} '?p!'],cl_id,iProbe));
        continue
      end
      if ~loaded && ok, loaded=1; loaded_param = param{k}; end
      P.(probeS(iProbe)) = p;
      n_ok = n_ok + 1;
      clear ok p
    end
    if loaded, break, end
  end
  if ~n_ok, data = []; cd(old_pwd), return, end
  
  if do_burst, problems = 'reset|bbias|sweep|saasa|probesa|nsops|spike';
  else, problems = 'reset|bbias|sweep|saasa|probesa|nsops';
  end
  if flag_rmwhip, problems = [problems '|whip']; end %#ok<NASGU>
  
  n_ok = 0; res = [];
  for probe=1:4
    signal = P.(probeS(probe));
    if ~isempty(signal)
      if ~do_burst
        nsops_errlist = [caa_str2errid('hxonly') caa_str2errid('bad_lx')...
          caa_str2errid(irf_ssub('no_p?',probe))];%#ok<NASGU>
        %      else
        %        nsops_errlist = [caa_str2errid(irf_ssub('no_p?',probe))];
      end
      remove_problems
      P.(probeS(probe)) = res;
      n_ok = n_ok + 1;
    end
  end
  clear res signal problems probe
  if ~n_ok
    irf_log('proc','No P data remaining after blanking.')
    data = []; cd(old_pwd); return
  end
  
  %Check for problem with bad DAC
  [ok,badDAC] = c_load('BADDAC?p34',cl_id);
  if ok || ~isempty(badDAC)
    irf_log('proc',irf_ssub('Bad DAC C?p34',cl_id))
    P.p3=[];P.p4=[];
  end
  [ok,badDAC] = c_load('BADDAC?p12',cl_id);
  if ok || ~isempty(badDAC)
    irf_log('proc',irf_ssub('Bad DAC C?p12',cl_id))
    P.p1=[];P.p2=[];
  end
  
  if isempty([P.p1; P.p2; P.p3; P.p4])
    irf_log('dsrc','Cannot compute P'), data=[]; cd(old_pwd); return
  end
  tComb = []; pList = [];
  for iProbe = 1:4
    ps = probeS(iProbe);
    if isempty(P.(ps)), continue, end
    tComb = [tComb; P.(ps)(:,1)]; pList = [pList; iProbe]; %#ok<AGROW>
  end
  if length(pList)==1
    Pinfo.probe = pList; ps = probeS(pList); %#ok<STRNU>
    irf_log('proc',['computing from ' ps])
    p = P.(ps); %#ok<NASGU>
  else
    tComb = sort(unique(tComb));
    for iProbe = 1:4
      ps = probeS(iProbe); pTmp = zeros(1,length(tComb))*NaN;
      if ~isempty(P.(ps))
        [~,ia,ib] = intersect(tComb,P.(ps)(:,1)); pTmp(ia) = P.(ps)(ib,2);
      end
      P.(ps) = pTmp;
    end
    clear pTmp ia ib
    % Fix periods of high bias saturation
    % from the two affected probes, we use only the one with the max value
    saProbes = [12 34 32];
    for iP = saProbes
      pr1 = fix(iP/10); pr2 = iP-pr1*10;
      if isempty(intersect([pr1 pr2],pList)), continue, end
      [ok,sa] = c_load(sprintf('HBIASSA%dp%d',cl_id,iP));
      if ~ok
        irf_log('load',sprintf('Cannot load HBIASSA%dp%d',cl_id,iP));
        continue
      end
      if isempty(sa), continue, end
      for j=1:size(sa,1)
        indx = tComb>=(sa(j,1)-10) & tComb<=(sa(j,2)+10);
        V1 = P.(probeS(pr1))(indx); V2 = P.(probeS(pr2))(indx);
        maxP = max([V1; V2],[],1); V1(V1~=maxP) = NaN; V2(V2~=maxP) = NaN;
        P.(probeS(pr1))(indx) = V1; P.(probeS(pr2))(indx) = V2;
      end
    end % hbiassa
    % We prefer to use the opposing ones
    if ~all(isnan(P.p1) == (isnan(P.p2) & isnan(P.p3) & isnan(P.p4)))
      ii12 = ~isnan(P.p1) & ~isnan(P.p2); ii34 = ~isnan(P.p3) & ~isnan(P.p4);
      fix12 = ii34 & ~ii12; fix34 = ii12 & ~ii34;
      P.p1(fix12) = NaN; P.p2(fix12) = NaN;
      P.p3(fix34) = NaN; P.p4(fix34) = NaN;
    end
    % Delete empty probes from the list
    for iProbe = pList'
      if all(isnan(P.(probeS(iProbe)))), pList(pList==iProbe) = []; end %#ok<AGROW>
    end
    if isempty(pList)
      irf_log('dsrc','Cannot compute P'), data=[]; cd(old_pwd); return
    end
    ps = sprintf('%d',pList);
    Pinfo.probe = str2double(ps); irf_log('proc',['computing from ' ps]) %#ok<STRNU>
    p = [tComb, irf.nanmean([P.p1; P.p2; P.p3; P.p4])']; %#ok<NASGU>
  end
  
  c_eval(['P' loaded_param '?=p;save_list=[save_list ''P' loaded_param '? ''];'],cl_id);
  if do_burst
    c_eval('bP?=p;bP?_info=Pinfo;bNVps?=c_efw_scp2ne(p);bNVps?(:,end+1)=p(:,2); save_list=[save_list '' bP? bP?_info bNVps?''];',cl_id)
  else
    c_eval('P?=p;P?_info=Pinfo;NVps?=c_efw_scp2ne(p);NVps?(:,end+1)=p(:,2); save_list=[save_list '' P? P?_info NVps?''];',cl_id)
  end
  clear p p1 p2 p3 p4 ptmp
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % P spin resolution
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'ps')
  save_file = './mP.mat';
  
  [ok,P_tmp,msg] = c_load('P?',cl_id);
  if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
  [ok,P_info] = c_load('P?_info',cl_id);
  if ~ok, P_info =[]; end
  
  P_tmp = P_tmp(~isnan(P_tmp(:,2)),:);
  if isempty(P_tmp)
    irf_log('proc',sprintf('Empty P%d.',cl_id))
    data = []; cd(old_pwd); return
  end
  
  % We always start at 0,4,8.. secs, so that we have
  % the same timelines an all SC at 2,6,10... sec
  t0 = fix(P_tmp(1,1)/4)*4 + 2;
  n = floor((P_tmp(end,1)-t0)/4) + 1;
  tvec = t0 + ( (1:n) -1)*4;
  
  if isfield(P_info,'useMax4hbiassa') && P_info.useMax4hbiassa==1
    P_tmp = irf_resamp(P_tmp,tvec','fsample',0.25,'max'); clear tvec %#ok<NASGU>
  else
    P_tmp = irf_resamp(P_tmp,tvec','fsample',0.25,'median'); clear tvec %#ok<NASGU>
  end
  c_eval('Ps?=P_tmp;save_list=[save_list ''Ps? '' ];',cl_id);
  
  if ~isempty(P_info)
    c_eval('Ps?_info=P_info;save_list=[save_list ''Ps?_info '' ];',cl_id);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % vce - E CIS PP [GSE+DSI]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vce')
  save_file='./mCIS.mat';
  
  if check_caa_sh_interval
    if ~exist('./.caa_sh_interval','file')
      irf_log('proc','Inside magnetosphere. No vce data fetched.')
      data = []; cd(old_pwd), return
    end
  end
  
  c_load('SAX?',cl_id)
  if ~exist('./mCISR.mat','file')
    irf_log('load','Please run ''vcis'' first (mCIS missing)')
    data = []; cd(old_pwd); return
  end
  
  CIS = load('mCISR'); %#ok<NASGU>
  % Load BPP. We use BPP for EDI as it must be a rather approximation
  [ok,B] = c_load('BPP?',cl_id);
  if ~ok
    [ok,B] = c_load('B?',cl_id);
    if ~ok
      irf_log('load',...
        irf_ssub('No B? and BPP?. Use getData(CDB,...,cl_id,''b'')',cl_id))
      data = []; cd(old_pwd); return
    end
  end
  
  vars = {'VCp', 'VCh'};
  varo = {'VCEp', 'VCEh'};
  for va=1:length(vars)
    eval(irf_ssub(['if isfield(CIS,''' vars{va} '?''); v=CIS.' vars{va} '?; else, v=[]; end; clear ' vars{va}], cl_id));
    if ~isempty(v)
      evxb=irf_tappl(irf_cross(v,B),'*(-1e-3)'); %#ok<NASGU>
      eval(irf_ssub([varo{va} '?=evxb;save_list=[save_list '' ' varo{va} '?'']; clear evxb'],cl_id));
      if ~exist(irf_ssub('SAX?',cl_id),'var')
        irf_log('load','must fetch spin axis orientation (option ''sax'')')
      else
        eval(irf_ssub(['di' varo{va} '?=c_gse2dsi(' varo{va} '?,SAX?);save_list=[save_list '' di' varo{va} '?''];'],cl_id));
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dibsc - despun B STAFF SC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dibsc')
  save_file='./mBSC.mat';
  
  [ok,wBSC,msg] = c_load('wBSC?',cl_id);
  if ~ok || isempty(wBSC)
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  [ok,pha, msg] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  aa = c_phase(wBSC(:,1),pha);
  diBSC = c_efw_despin(wBSC,aa);
  % DS-> DSI
  diBSC(:,3)=-diBSC(:,3);
  diBSC(:,4)=-diBSC(:,4); %#ok<NASGU>
  
  c_eval('diBSC?=diBSC;save_list=[save_list '' diBSC?''];',cl_id);
  
  [ok,sax,msg] = c_load('SAX?',cl_id);
  if ~ok || isempty(sax)
    irf_log('load',msg)
  else
    c_eval('BSC?=c_gse2dsi(diBSC?,sax,-1);save_list=[save_list '' BSC?''];',cl_id);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dibscburst - despun B STAFF SC burst
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dibscburst')
  save_file='./mEFWburst.mat';
  
  [ok,wBSC4kHz,msg] = c_load('wBSC4kHz?',cl_id);
  if ~ok || isempty(wBSC4kHz)
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  [ok,pha, msg] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  aa = c_phase(wBSC4kHz(:,1),pha);
  diBSC4kHz = c_efw_despin(wBSC4kHz,aa);
  % DS-> DSI
  diBSC4kHz(:,3)=-diBSC4kHz(:,3);
  if size(diBSC4kHz,2)>3
    diBSC4kHz(:,4)=-diBSC4kHz(:,4); %#ok<NASGU>
  end
  
  c_eval('diBSC4kHz?=diBSC4kHz;save_list=[save_list '' diBSC4kHz?''];',cl_id);
  
  [ok,sax,msg] = c_load('SAX?',cl_id);
  if ~ok || isempty(sax)
    irf_log('load',msg)
  else
    c_eval('BSC4kHz?=c_gse2dsi(diBSC4kHz?,sax,-1);save_list=[save_list '' BSC4kHz?''];',cl_id);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % manproblems: Reads manually-set problems from database.
  %              affects whip|sweep|bdump|badbias|probesa|hbiassa|wake
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'manproblems')
  save_file='./mEFW.mat';
  
  d=[c_ctl('get', 5, 'data_path') '/caa-control'];
  f_name = [d '/manual_problems_c' num2str(cl_id) '.dat'];
  if ~exist(f_name,'file')
    irf_log('load',['file ' f_name ' not found']);
    data = []; cd(old_pwd), return
  end
  fid = fopen(f_name);
  C = textscan(fid, '%s %n %1[+-] %s','commentStyle', '%');
  fclose(fid);
  
  [iso_t,dt] = caa_read_interval;
  int_start = iso2epoch(iso_t);
  int_end=int_start+dt;
  for i=1:length(C{1})
    st=iso2epoch(C{1}{i});
    dt=C{2}(i);
    if (st<int_end && (st+dt)>int_start)
      if ~exist(C{4}{i},'var')
        eval(['[ok,' C{4}{i} ',msg]=c_load(C{4}{i});'])
        if ~ok %#ok<NODEF>
          irf_log('load',['Load failed of ' C{4}{i}])
        else, irf_log('load',msg) %#ok<NODEF>
        end
        clear ok hbsa msg
      end
      if C{3}{i} == '+'
        irf_log('proc',['Setting manual problem:' C{4}{i}]);
        eval([C{4}{i} '=[' C{4}{i} ''' [st st+dt]'']'';']);
        eval('save_list=[save_list C{4}{i} '' ''];');
      elseif C{3}{i} == '-'
        irf_log('proc',['Removing manual problem: ' C{4}{i}]);
        eval(['prob=' C{4}{i} ';'])
        if any(prob)
          idx=find(prob(:,1)<st & prob(:,2)>(st+dt));
          if any(idx)
            prob=[prob' [prob(idx,1) st]' [st+dt prob(idx,2)]']';
            prob(idx,:)=0;
          end
          idx=find(prob(:,1)<st & prob(:,2)>st);
          if any(idx), prob(idx,2)=st;end
          idx=find(prob(:,1)<st+dt & prob(:,2)>st+dt);
          if any(idx), prob(idx,1)=st+dt;end
          idx=find(prob(:,1)>=st & prob(:,2)<=st+dt);
          if any(idx), prob(idx,1:2)=0;end
          idx=find(prob(:,1) ~= 0);
          if any(idx), prob=prob(idx,:);
          else, prob=[];
          end
          eval([C{4}{i} '=prob;'])
          eval('save_list=[save_list C{4}{i} '' ''];');
        end
      else
        error(['Corrupt manual problems file ' f_name])
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % hk - housekeeping: time bias1-4 puck(stub)1-4 guard1-4
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'hk')
  save_file='./mEFW.mat';
  
  global c_ct %#ok<TLEV> % includes hk calib values
  if isempty(c_ct)
    disp('c_ct undefined');
  end
  if c_ct{1,1}.ibias(2,1) == 0.
    %        irf_log('proc','c_ct hk cal not loaded: run c_ctl(''load_hk_cal'')');
    c_ctl('load_hk_cal');
  end
  
  [ok,HK,msg] = c_load('DSC?',cl_id);
  if ~ok || isempty(HK)
    irf_log('load',msg)
    data = []; cd(old_pwd); return
  end
  
  [ok,fdm] = c_load('FDM?',cl_id);
  if ~ok
    irf_log('load',...
      irf_ssub('No FDM?. Use getData(CDB,...,cl_id,''fdm'')',cl_id))
    data = []; cd(old_pwd); return
  end
  
  %    datestr(epoch2date(HK(1,1)))
  %    binval=HK(130:141,:);
  binval=HK(:,130:141);
  sbv = size(binval,1);
  sfdm=size(fdm,1);
  DSCindex=bitand(fdm(:,3),31);
  %    DSCindex
  
  % convert bin values like in isdat ./server/Wec/Efw/EfwCal.c
  ix = binval>127;
  binval(ix)=binval(ix)-256;
  binval=binval+128;
  ix = binval<0;
  binval(ix)=0;
  ix = binval>255;
  binval(ix)=255;
  
  binvalred=zeros(sbv,12);
  timest = zeros(sbv);
  
  k=0;
  lastfound=1;
  for j=1:sbv % check DSCindex 0-31 and romid as in isdat ./server/Wec/Index.c
    for i=lastfound:sfdm
      if abs(HK(j,1)-fdm(i,1))<.5
        cnt=0;
        ok=1;
        for l=i:sfdm
          if cnt~=DSCindex(l)
            ok=0;
            break;
          end
          if cnt>=31
            break;
          end
          cnt = cnt + 1;
        end
        romid=HK(j,80); % Does checking romid really work for a reset?
        if romid ~= hex2dec('b1') && romid ~= hex2dec('e1') && romid ~= hex2dec('f1') && romid ~= hex2dec('f2') &&...
            romid ~= hex2dec('f3') && romid ~= hex2dec('f4') && romid ~= hex2dec('f5') && romid ~= hex2dec('f6') &&...
            romid ~= hex2dec('f7') && romid ~= hex2dec('f8') && romid ~= hex2dec('f9')
          ok=0;
        end
        
        if ok==1
          k = k + 1;
          binvalred(k,:)=binval(j,:);
          timest(k)=HK(j,1);
        else
          irf_log('proc',sprintf('DSCindex and/or romid bad. i%d cnt%d ind%d rom%X',i,cnt,DSCindex(l),HK(j,80)));
        end
        lastfound=i+1;
      end
    end
  end
  redsize=k;
  
  calhk=zeros(12,redsize);
  % substitute to calib values from c_ct
  for i=1:redsize
    for k=1:4
      calhk(k,i)=c_ct{1,cl_id}.ibias(binvalred(i,k)+1,1+k);
      calhk(k+4,i)=c_ct{1,cl_id}.puck(binvalred(i,k+4)+1,1+k);
      calhk(k+8,i)=c_ct{1,cl_id}.guard(binvalred(i,k+8)+1,1+k);
    end
  end
  
  %    HK_info.probe = '1234'; % c_desc support for hk info needed
  %    c_eval('HK?_info=HK_info;save_list=[save_list ''HK?_info'' '' HK?''];',cl_id)
  % save time and only dsc 128->139: time,  BIAS_DAC_1 -> 4, STUB_DAC_1 -> 4, GUARD_DAC_1 -> 4
  c_eval('HK?=[timest(1:redsize)'' calhk''];',cl_id);
  c_eval('save_list=[save_list ''HK?''];',cl_id);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF DATA MANIPULATIONS
% saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If flag_save is set, save variables to specified file
if flag_save==1 && ~isempty(save_file) && ~isempty(save_list)
  irf_log('save',[save_list ' -> ' save_file])
  if exist(save_file,'file')
    eval(['save -append ' save_file ' ' save_list]);
  else
    eval(['save ' save_file ' ' save_list]);
  end
end

% prepare the output
if nargout > 0
  if isempty(save_list)
    data = [];
  else
    sl = tokenize(save_list);
    data = {sl};
    for k=1:length(sl)
      if exist(sl{k},'var'), eval(['data{k+1}=' sl{k} ';']); end
    end
  end
end

cd(old_pwd)
end

function s = probeS(p)
s = sprintf('p%d',p);
end
