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
%          correct_sw_wake - correct wakes in the Solar Wind
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
%   whip: WHIP{cl_id} -> mEFW	// Whisper pulses present +1 precceding sec
%   bdump: DBUMP{cl_id} -> mEFW	// Burst dump present
%   sweep: SWEEP{cl_id} -> mEFW	// Sweep + dump present
%   badbias: BADBIASRESET{cl_id}, BADBIAS{cl_id}p[1..4] -> mEFW	
%          // Bad bias settings
%   probesa: PROBESA{cl_id}p[1..4] -> mEFW	// Probe saturation
%   rawspec: RSPEC{cl_id}p{12/32,34} -> mEFW // Spectrum of raw signal (1,2..5 omega)
%   wake : PSWAKE{cl_id}p{12/34}, LOWAKE{cl_id}p{12/34} -> mEFW // wakes
%   edi : EDI{cl_id}, diEDI{cl_id} -> mEDI // EDI E in sc ref frame
%   br, brs : Br[s]{cl_id}, diBr[s]{cl_id} -> mBr // B resampled to E[s]
%   vedbs, vedb : VExB[s]{cl_id}, diVExB[s]{cl_id} -> mEdB // E.B=0 [DSI+GSE]
%   vce : VCE(p,h){cl_id},diVCE(p,h){cl_id} ->mCIS	// E CIS PP [GSE+DSI] 
%
% Example: 
%       getData(cp,4,'edbs','ang_fill','ang_limit',20,'probe_p',12)
%
% General options - one of the following:
%       nosave : do no save on disk
%       withwhip : do not remove time intervals with Whisper pulses
%       notusesavedoff : recalculating everything instead of using saved offsets
%
% See also C_GET, C_CTL
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(3,15,nargin))
if nargin > 3, have_options = 1; args = varargin;
else have_options = 0;
end

% default options
flag_save = 1;
flag_usesavedoff = 1;
flag_edb = 1;
sfit_ver = -1;
correct_sw_wake = 0;

CAA_MODE = c_ctl(0,'caa_mode');

flag_rmwhip = c_ctl(cl_id,'rm_whip');
if isempty(flag_rmwhip), flag_rmwhip = 1; end 
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
	case 'notusesavedoff'
		flag_usesavedoff = 0;
	case 'ang_limit'
		if length(args)>1
			if isnumeric(args{2})
				ang_limit = args{2};
				l = 2;
            else irf_log('fcal,','wrongArgType : ang_limit must be numeric')
			end
        else irf_log('fcal,','wrongArgType : ang_limit value is missing')
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
            else probe_p_tmp = str2double(args{2});
			end
			if (probe_p_tmp==12 || probe_p_tmp==34), probe_p = probe_p_tmp;
            else irf_log('fcal,','wrongArgType : probe_p must be 12 or 34')
			end
        else irf_log('fcal,','wrongArgType : ang_limit value is missing')
		end
	case 'sfit_ver'
		if length(args)>1
			if isnumeric(args{2})
				l = 2;
				if	args{2}>=0 && args{2}<2, sfit_ver = args{2};
                else irf_log('fcal,','wrongArgType : sfit_ver must be 0 or 1')
				end
            else irf_log('fcal,','wrongArgType : sfit_ver must be numeric')
			end
        else irf_log('fcal,','wrongArgType : sfit_ver value is missing')
		end
	case 'correct_sw_wake'
		correct_sw_wake = 1;
	otherwise
		irf_log('fcal,',['Option ''' args{1} '''not recognized'])
	end
	if length(args) > l, args = args(l+1:end);
	else break
	end
end


save_list = '';

old_pwd = pwd;
cd(cp.sp) %enter the storage directory
if cp.sp~='.', irf_log('save',['Storage directory is ' cp.sp]), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ec - correct raw Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'ec')
	save_file = './mERC.mat';
	
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
		
		problems = 'reset|bbias|probesa|probeld|sweep|bdump'; %#ok<NASGU>
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
			[da, n_corrected,wake_dsc] = c_efw_swwake(da,probe,pha,whip); %#ok<NASGU,ASGLU>
			
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
% dies - spin fiting of Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dies')
	save_file = './mEDSI.mat';
	
	% Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
	[ok,pha] = c_load('Atwo?',cl_id);
	if ~ok || isempty(pha)
		irf_log('load',...
			irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
		data = []; cd(old_pwd); return
	end
	
	p12 = 12; e12 = []; e34 =[];
	n_ok = 0;
	corrected_raw_data_p12 = 1;
	corrected_raw_data_p34 = 1;
	for probe = [12,32,34]
		[ok,da] = c_load(irf_ssub('wcE?p!',cl_id,probe));
		if ~ok || isempty(da)
			irf_log('load', irf_ssub('No/empty wcE?p!',cl_id,probe));
			[ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
			if ~ok || isempty(da)
				irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
				continue
			end
			irf_log('load','using raw (not corrected) data')
			if probe==34, corrected_raw_data_p34 = 0;
			else corrected_raw_data_p12 = 0;
			end
		end
		n_ok = n_ok + 1;
		
		if probe==32
			p12 = 32;
			e12 = da;
        else c_eval('e?=da;',probe)
		end
		clear ok da
	end
	if ~n_ok, data = []; cd(old_pwd), return, end
	
	% If we have different timelines for p1(3)2 and p34 we try to make 
	% a common timeline, so that the resulting spinfits will also have 
	% the common timeline
	not_same = 0;
	if n_ok==2
		if length(e12(:,1))==length(e34(:,1))
			not_same = 0;
			% Check fo same length but different timelines
			if length(find((e12(:,1)-e34(:,1))~=0)), not_same = 1; end
        else not_same = 1;
		end
	end
	if not_same && n_ok==2
		ts = e34(1,1); te = e34(end,1);
		ts1 = e12(1,1); te1 = e12(end,1);
		if ts>ts1, ts = ts1; end
		if te<te1, te = te1; end
		dt = double(te - ts); clear te ts1 te1
		
		% Guess the sampling frequency
		fsamp = c_efw_fsample(e34(:,1),'hx');
		if ~fsamp, error('no sampling frequency'),end
		
		irf_log('proc','Using new time line for spinfits')
		t_new = double(0):double(1/fsamp):dt+double(1/fsamp); 
		t_new = t_new';
		
		d_new = zeros(length(t_new),2);
		d_new(:,1) = ts + t_new;
		d_new(:,2) = NaN;
		
		ii = round((e34(:,1)-ts)*fsamp+1);
		if ii(end)>length(t_new), error('problemo with new time line'), end
		d_new(ii,2) = e34(:,2); e34 = d_new;
		d_new(:,2) = NaN;
		
		ii = round((e12(:,1)-ts)*fsamp+1);
		if ii(end)>length(t_new), error('problemo with new time line'), end
		d_new(ii,2) = e12(:,2); e12 = d_new;
		
		clear t_new ts dt d_new ii
	end
	clear not_same
	
	aa = [];
	n_sig = 0;
	for pr=[12,34]
		
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
			irf_log('load',sprintf('No spinfits for C%d p%d',cl_id,probe))
			continue
		end
		
		if corrected_raw_data, ss = 'c';
		else ss = '';
		end
		irf_log('proc',sprintf('Spin fit w%sE%dp%d -> diEs%dp%d',...
			ss,cl_id,probe,cl_id,probe))

		if isempty(aa), aa = c_phase(tt(:,1),pha); end
		if isempty(aa)
			irf_log('proc','Empty phase')
			continue
		end
		
		fsamp = c_efw_fsample(tt,'hx');
		if ~fsamp, error('no sampling frequency'),end
		
		problems = 'reset|bbias|probesa|probeld|sweep|bdump'; %#ok<NASGU>
		%if flag_rmwhip, problems = [problems '|whip']; end
		signal = tt; %#ok<NASGU>
		remove_problems
		tt = res; %#ok<NODEF>
		clear res signal problems
		
		% Check if we have at least 1 spin of data left
		if length(find(~isnan(tt(:,2)))) < 4*fsamp
			irf_log('proc',irf_ssub('No p? data after removals',probe))
			continue
		end

		
		if sfit_ver>=0
			irf_log('proc',['using SFIT_VER=' num2str(sfit_ver)])
			sp = c_efw_sfit(probe,3,10,20,tt(:,1),tt(:,2),aa(:,1),...
				aa(:,2),sfit_ver);
		else
			sp = c_efw_sfit(probe,3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
		end
		
		% Check if we have any data left
		if isempty(sp)
			irf_log('load',sprintf('No data left after spinfit for C%d p%d',...
				cl_id,probe))
			continue
		end
		
		% Remove point with zero time
		ind = find(sp(:,1)>0);
		if length(ind)<length(sp(:,1))
			irf_log('proc',...
				[num2str(length(sp(:,1))-length(ind)) ' spins removed (bad time)']);
			sp = sp(ind,:);
		end
		
		% ADC offsets
		adc_off = sp(:,[1 4]);
		% Warn about points with sdev>.8
		ii = find(sp(:,6)>.8);
		if length(ii)/size(sp,1)>.05,
			irf_log('proc',[sprintf('%.1f',100*length(ii)/size(sp,1)) '% of spins have SDEV>.8 (ADC offsets)']);
		end
		adc_off = irf_waverage(adc_off,1/4);
		adc_off(adc_off(:,2)==0,2) = mean(adc_off(abs(adc_off(:,2))>0,2)); %#ok<NASGU>
		
		sp = sp(:,[1:4 6]);
		sp(:,4) = 0*sp(:,4); % Z component
		
		% Remove spins with bad spin fit (obtained E > 10000 mV/m)
		ind = find(abs(sp(:,3))>1e4); sp(ind,:) = []; %#ok<NASGU>
		if ind, disp([num2str(length(ind)) ' spins removed due to E>10000 mV/m']);end
		eval(irf_ssub('diEs?p!=sp;Dadc?p!=adc_off;',cl_id,probe)); 
		eval(irf_ssub('save_list=[save_list ''diEs?p! Dadc?p! ''];',cl_id,probe));
		n_sig = n_sig + 1;
		clear tt sp adc_off
	end

	% Delta offsets (offsets between E DSI obtained from p12/32 and p34)
	if n_sig==2
		
		% To compute delta offsets we remove points which are > deltaof_sdev_max*sdev
		% as this must de a stable quantity
		eval(irf_ssub('[ii1,ii2] = irf_find_comm_idx(diEs?p!,diEs?p34);',cl_id,p12))
		eval(irf_ssub('df=diEs?p!(ii1,2:3)-diEs?p34(ii2,2:3);',cl_id,p12))
		clear ii1 ii2
		df(isnan(df)) = []; %#ok<NODEF>
		iia = [];
		if size(df,1)>2
			sdev = std(df);
			for comp = 1:2
				ii = find(abs(df(:,comp)-mean(df(:,comp))) > deltaof_sdev_max*sdev(comp)); 
				iia = [iia; ii];
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
			
		irf_log('calb',sprintf('delta offsets are: %.2f [x] %.2f [y]', ...
			Del(1), Del(2)))

		% Check for unreallistically large Del. 
		% If it is larger than deltaof_max, we set it to zero and 
		% NOT doing any corrections.
		if ( abs(Del(1)) > deltaof_max ) || ( abs(Del(2)) > deltaof_max )
			irf_log('calb',...
				irf_ssub('DELTA OFFSET TOO BIG >!. Setting D?p12p34=[0 0]',...
				cl_id,deltaof_max))
		else
			% Always correct p12/p32.
			% Deprecated behavior: real offset is applied to p12, imaginary to p34.
			irf_log('calb',irf_ssub('correcting p?',p12))
			eval(irf_ssub('diEs?p!(:,2:3)=diEs?p!(:,2:3)-ones(size(diEs?p!,1),1)*Del;',cl_id,p12));
		end
		
		eval(irf_ssub('D?p12p34=Del;',cl_id))
		clear m12 m34 Del

		eval(irf_ssub('save_list=[save_list ''D?p12p34 ''];',cl_id));
	elseif n_sig==1
		% If we have only p12, we should use it for sfits&co
		if exist(irf_ssub('diEs?p!',cl_id,p12), 'var')
			irf_log('calb',irf_ssub('Will use p? for spin res data',p12))
			caa_sfit_probe(cl_id,p12);
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die - despin of full resolution data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'die') || strcmp(quantity,'dief') || ...
	strcmp(quantity,'diespec') || strcmp(quantity,'dieburst')
	
	if strcmp(quantity,'dieburst'), do_burst = 1; else do_burst = 0; end
	if strcmp(quantity,'dief'), do_filter = 1; else do_filter = 0; end
	if do_burst
		save_file = './mEFWburst.mat';
		var1_name = 'dibE?p1234';
	else
		if strcmp(quantity,'diespec'), save_file = './mEDSI.mat';
        else save_file = './mEDSIf.mat';
		end
		if do_filter, var1_name = 'diEF?p1234';
        else var1_name = 'diE?p1234';
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
	
    % Make electric field for the burst
	if do_burst
		var_name = 'wbE?p!';
		p12 = 12; e12 = []; e34 =[];
		param={'180Hz','4kHz','32kHz'};
		p_ok = [];
		loaded = 0;
		for k=1:length(param)
			for probe=[2 4]
				[ok,p2] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,probe));
				if ok
					loaded = 1;
					p_sep = .088;
					if probe==2
						[ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,1));
						if ~ok
							[ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,3));
							p_sep = .066;
							p12 = 32;
						end
						if ~ok
							irf_log('load', irf_ssub(['No P' param{k} '?p1/3'],cl_id));
							continue
						end

						e12(:,1) = p2(:,1);
						e12(:,2) = ( p2(:,2) - p1(:,2) )/p_sep;
						p_ok = [p_ok 12]; %#ok<AGROW>
						eval(irf_ssub('wbE?p!=e12;save_list=[save_list ''wbE?p! ''];',cl_id, p12));
					else
						[ok,p1] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,3));
						if ~ok
							irf_log('load', irf_ssub(['No P' param{k} '?p3'],cl_id));
							continue
						end

						e34(:,1) = p2(:,1);
						e34(:,2) = ( p2(:,2) - p1(:,2) )/p_sep;
						p_ok = [p_ok 34]; %#ok<AGROW>
					end
				else
					irf_log('load', irf_ssub(['No P' param{k} '?p!'],cl_id, probe));
				end
			end
			if loaded, break, end
		end
		if ~loaded
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
	else % Not burst
		p12 = 12; e12 = []; e34 =[];
		p_ok = [];
		for probe = [12 32 34]
			[ok,da] = c_load(irf_ssub('wcE?p!',cl_id,probe));
			if ~ok || isempty(da)
				irf_log('load', irf_ssub('No/empty wcE?p!',cl_id,probe));
				[ok,da] = c_load(irf_ssub('wE?p!',cl_id,probe));
				if ~ok || isempty(da)
					irf_log('load', irf_ssub('No/empty wE?p!',cl_id,probe));
					continue
				end
				irf_log('load','using raw (not corrected) data')
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
	if ~length(p_ok), data = []; cd(old_pwd), return, end

	% Load ADC offsets
	for probe = p_ok
		if probe==12, ps=p12; else ps=probe; end
		[ok,dadc] = c_load(irf_ssub('Dadc?p!',cl_id,ps));
		if ~ok
			if CAA_MODE, error(irf_ssub('Cannot load Dadc?p!',cl_id,ps)), end
			irf_log('load',irf_ssub('Cannot load Dadc?p!',cl_id,ps))
		end
		if isempty(dadc)
			if CAA_MODE, error(irf_ssub('Empty Dadc?p!',cl_id,ps)), end
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
			if ~isempty(dadc)
				irf_log('calb','using saved ADC offsets')
				tmp_adc = irf_resamp(dadc,tt,'fsample',c_efw_fsample(tt,'ib'));
				tt(:,2) = tt(:,2) - tmp_adc(:,2); %#ok<NASGU>
				clear tmp_adc
			else irf_log('calb','saved ADC offset empty')
			end
			n_sig = n_sig + 1;
			c_eval('e?=tt;',p)
	
		else
			fsamp = c_efw_fsample(tt,'hx');
			if ~fsamp, error('no sampling frequency'),end
			
			problems = 'reset|bbias|probesa|probeld|sweep|bdump';
			if flag_rmwhip, problems = [problems '|whip']; end %#ok<NASGU,AGROW>
			signal = tt; %#ok<NASGU>
			probe = ps; %#ok<NASGU>
			remove_problems
			tt = res; %#ok<NODEF>
			clear res signal problems probe
		
			% Check if we have at least 1 sec of data left
			if length(find(~isnan(tt(:,2)))) < fsamp
				irf_log('proc',irf_ssub('No p? data after removals',ps))
				c_eval('e?=[];',p)
				continue
			end
			
			% Correct ADC offset
			if flag_usesavedoff
				% Correct ADC offset
				if ~isempty(dadc)
					irf_log('calb','using saved ADC offset')
					tmp_adc = irf_resamp(dadc,tt,'fsample',fsamp);
					tt(:,2) = tt(:,2) - tmp_adc(:,2);
					clear tmp_adc
				else
					irf_log('calb','saved ADC offset empty')
					flag_usesavedoff = 0;
				end
			end
			if ~flag_usesavedoff
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
					e34_tmp = e12;
					e34_tmp(~isnan(e12(:,2)),2) = 0;
					e34_tmp(ii12,2) = e34(ii34);
					e34 = e34_tmp;
					clear e34_tmp
				else
					irf_log('proc',['Setting Ep' num2str(p12)...
						' to 0, except for ' num2str(length(ii12)) ' data points'])
					e12_tmp = e34;
					e12_tmp(~isnan(e34(:,2)),2) = 0;
					e12_tmp(ii34,2) = e12(ii12);
					e12 = e12_tmp;
					clear e12_tmp
				end
				irf_log('proc','!!! REDUCED DATA QUALITY !!!')
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
        else fsamp = c_efw_fsample(e12,'hx');
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
			irf_log('proc','!!! REDUCED DATA QUALITY !!!')
			irf_log('proc','using one probe pair some part of the interval')
		end
		
		% Use WEC coordinate system E=[t,0,p34,p12]
		full_e = zeros(length(e12),4);
		full_e(:,[1,4]) = e12;
		full_e(:,3) = e34(:,2);
		clear e12 e34
		
		% Load Delta offsets D?p12p34
		[ok,Del] = c_load('D?p12p34',cl_id);
		if ~ok || isempty(Del)
			irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id,ps))
		else
			
			if isreal(Del)
				% Real Del means we must correct p12/p32.
				irf_log('calb',['correcting delta offset on p' num2str(p12)])
				i_c = 1;
			else
				% Correcting p34 is now DEPRECATED
				irf_log('calb','correcting delta offset on p34')
				Del = imag(Del);
				i_c = 2;
			end
			coef(i_c,3) = Del(1) -Del(2)*1j;
			clear Del
			
	
% 			irf_log('calb',['correcting delta offset on p' num2str(p12)])
% 
% 			% Correcting p34 is now DEPRECATED [~isreal(Del)]
% 			if ~isreal(Del), Del = -imag(Del); end
% 			coef(1,3) = Del(1) -Del(2)*1j;

		
		end
	else
	
		% We have one probe pair
		if ~isempty(e12)
			pp = 12;
			E_info.probe = num2str(p12);
			EE = e12;
			clear e12
		else
			pp = 34;
			E_info.probe = '34';
			EE = e34;
			clear e34
		end
		% Use WEC coordinate system E=[t,0,p34,p12]
		full_e = zeros(length(EE),4);
		full_e(:,1) = EE(:,1);
		if pp==12, full_e(:,4) = EE(:,2);
        else full_e(:,3) = EE(:,2);
		end
		clear EE pp
	end

	c_eval([var1_name '_info=E_info;save_list=[save_list ''' var1_name '_info ''];'],cl_id);

	% Do actual despin
	aa = c_phase(full_e(:,1),pha);
	if p12==32, full_e=c_efw_despin(full_e,aa,coef,'asym');
    else full_e=c_efw_despin(full_e,aa,coef);
	end
	
	if strcmp(quantity,'diespec')
		% Make a spectrum and save it.
		sfreq = c_efw_fsample(full_e(:,1),'hx');
		if ~sfreq, error('no sampling frequency'),end
		if sfreq == 25, nfft = 512; %#ok<NASGU>
        else nfft = 8192; %#ok<NASGU>
		end
		c_eval(...
		'diESPEC?p1234=caa_powerfft(full_e(:,1:3),nfft,sfreq);save_list=[save_list ''diESPEC?p1234 ''];',cl_id);
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
            else err_s = [err_s ', ' var_s{k}]; %#ok<AGROW>
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
    else inert = 0; 
	end
	
	if inert, save_file = './mEdBI.mat';
    else save_file = './mEdB.mat';
	end
	
	if strcmp(quantity,'edb') || strcmp(quantity,'iedb')
		var_s = irf_ssub('diE?p1234',cl_id); e_opt = 'die';
		varo_s = irf_ssub('E?',cl_id);
		var_b = 'diBr?'; b_opt ='br';
	else
		e_opt = 'dies';
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
		otherwise
			error(['Invalid probe pair ' num2str(probe_p)])
		end
		varo_s = irf_ssub('Es?',cl_id);
		var_b = 'diBrs?'; b_opt ='brs';
	end
	
	% Load resampled B
	[ok,diB] = c_load(var_b,cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No ' var_b ' in mBr. Use getData(CP,cl_id,''' b_opt ''')'],cl_id))
		data = []; cd(old_pwd); return
	end

	% Load V if we need to do SC->Inertial transformation
	if inert
		[ok,diV] = c_load('diV?',cl_id);
		if ~ok
			irf_log('load',...
				irf_ssub('No diV? in mR. Use getData(CDB,...,cl_id,''v'')',cl_id))
			data = []; cd(old_pwd); return
		end
	end

	[ok,diE] = c_load(var_s);
	if ~ok
		dsc = c_desc(var_s);
		irf_log('load',...
			irf_ssub(['No ' var_s ' in ' dsc.file '. Use getData(CP,cl_id,''' e_opt ''')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
	% Save stand-dev in the 6-th column
	if strcmp(quantity,'edbs') || strcmp(quantity,'iedbs'), diE(:,6) = diE(:,5); end
	
	dsiof = c_ctl(cl_id,'dsiof');
	if isempty(dsiof)
		st = diE(~isnan(diE(:,1)),1);
		if ~isempty(st), st = st(1); else st = 0; end
		[dsiof_def, dam_def] = c_efw_dsi_off(st,cl_id); clear st

		[ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
		[ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end

		if ok1 || ok2, irf_log('calb','Using saved DSI offsets')
		else irf_log('calb','Using default DSI offsets')
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
    else s = '';
	end
	
 	% DSI->GSE
	if c_load('SAX?',cl_id)
		c_eval([s varo_s '=c_gse2dsi(diE(:,1:4),SAX?,-1);' s varo_s '(:,5)=diE(:,5);save_list=[save_list ''' s varo_s ' ''];'],cl_id);
	else
		irf_log('load',irf_ssub('No SAX? in mEPH. Use getData(CDB,...,cl_id,''sax'')',cl_id))
	end

	eval([ s 'di' varo_s '=diE;']); clear diE
	eval(irf_ssub('ang_limit?=ang_limit;',cl_id)) 
	save_list=[save_list s 'di' varo_s ' ang_limit' num2str(cl_id) ' '];

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
	ii_px = find(fdm_px(:,1)==1 & fdm_px(:,2)==1);
	if ~isempty(ii) || ~isempty(ii_px)
		if isempty(ii)
			% We hit sweep dump directly
			irf_log('dsrc','found loonely sweep dump')
			if length(ii_px)>1, irf_log('proc','WARNING: too many loonely sweep dumps'), end
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
                else bdump(k,2) = t_e_px(jj(end)) +1;
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
	[ok,efwt] = c_load('EFWT?',cl_id);
	if ok
		c_eval(['BADBIASRESET?=[];'...
			'save_list=[save_list '' BADBIASRESET? ''];'],cl_id);
				
		ii = find(efwt(:,2)<DELTA_PLUS);
		if ~isempty(ii)
			t0 = efwt(ii(1),1) - efwt(ii(1),2);
			irf_log('proc', ['EFW reset at ' epoch2iso(t0,1)]);
			c_eval('BADBIASRESET?=[double(t0-DELTA_MINUS)'' double(t0+DELTA_PLUS)''];',cl_id);
		end

		% 2006-06-16 23:30 bias current was lowered to 100 nA
		if efwt(1,1)>iso2epoch('2006-06-16T00:00:00Z'), GOOD_BIAS = -95; end
		
    else irf_log('dsrc',irf_ssub('Cannot load EFWT?',cl_id))
	end
	clear ok t0 efwt ii
	
	% The reason we remove 300 seconds (DELTA_MINUS) of data before a bad
	% bias is that we get rid of all the EFW resets in a clean way.
	% The 64 seconds (DELTA_PLUS) afterwards is because we have trouble deciding
	% the timing of the bias current, so we takek 2 x DSC interval.
	DELTA_MINUS = 300;
	DELTA_PLUS = 64;
	
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
				if ibias(1,2)==0, ii = [1; ii]; end
				if ibias(end,2)==0, ii = [ii; length(ibias(:,2))]; end
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
								irf_log('proc',[epoch2iso(bb_st(ii),1) ' -- ' ...
									epoch2iso(bb_et(ii),1)])
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
	
	% Saturation level nA
	SA_LEVEL = 66; %#ok<NASGU>
	SA_LEV_POS = 0; % 0.0 V
	
	% N_CONST sets the minimum number of points of constant potential
	% which we consider bad
	N_CONST = 4;
	% DT_PLUMIN is the interval by which we extend saturation 
	% sintervals from each side
	DT_PLUMIN = 4;
	% Delta = .1 sec, half sampling interval for P.
	DELTA = .1;
	
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
	
	[iso_t,dt] = caa_read_interval;
	start_time = iso2epoch(iso_t);
	c_eval('sa_int_p?=[];')
	for pro=1:4
		[ok,p] = c_load(irf_ssub('P10Hz?p!',cl_id,pro));
		if pro==3 && ~isempty(start_time) && ...
			start_time>toepoch([2001 07 23 00 00 00]) && cl_id==2
			
			sa_int = []; %#ok<NASGU>
			if ~isempty(ns_ops)
				sa_int = caa_get_ns_ops_int(start_time,dt,ns_ops,'no_p3');
				if ~isempty(sa_int), irf_log('proc','Found no_p3 in NS_OPS'), end
			end
			
			irf_log('dsrc',...
				irf_ssub('Using fake PROBELD?p!',cl_id,pro));
			c_eval(['p?=[];sa_int_p?=sa_int;PROBELD' num2str(cl_id) ...
				'p?=[];save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
			continue
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
	
	% Points below SA_LEVEL should be excluded from E, but not from
	% P ans they atill contain valuable physical information.
	% This is not the case with points with positive and/or 
	% constant potential (latched probe).
	
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
			if p_tmp(1,2)==0, ii = [1; ii]; end
			if p_tmp(end,2)==0, ii = [ii; length(p_tmp(:,2))]; end
			ii = reshape(ii,2,length(ii)/2);
			res = [p_tmp(ii(1,:))-DELTA; p_tmp(ii(2,:))+DELTA]'; %#ok<NASGU>
			c_eval(['PROBELD' num2str(cl_id) 'p?=res;'...
			'save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
			
			clear res ii
		end
		c_eval('p?=p;',pro)
		clear ii_god ii_bad p
	end
	
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
		ii_bad = find( p(:,2) >=SA_LEV_POS );
		ii_god = find( p(:,2) < SA_LEV_POS );
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
                            else jj = ll(1);
							end
						end
					end
				end
			end
			p(ii_bad,2) = 0; %#ok<FNDSB>
			ii = irf_find_diff(p(:,2));
			if p(1,2)==0, ii = [1; ii]; end
			if p(end,2)==0, ii = [ii; length(p(:,2))]; end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rawspec - Spectrum of raw EFW signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'rawspec')
	save_file = './mEFW.mat';
	
	% Src quantities: Atwo?, wE?p12/wE?p32, wE?p34
	[ok,pha] = c_load('Atwo?',cl_id);
	if ~ok || isempty(pha)
		irf_log('load',...
			irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
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
	n_ok = 0;
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
			else corrected_raw_data_p12 = 0;
			end
		end
		n_ok = n_ok + 1;
		
		if size(da,1) > length(tpharef), tpharef = da(:,1); end
		
		if probe==32
			p12 = 32;
			e12 = da;
        else c_eval('e?=da;',probe)
		end
		clear ok da
	end
	if ~n_ok, data = []; cd(old_pwd), return, end
	
	%Compute spin period
	ph = c_phase(tpharef,pha); 
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
    else irf_log('proc',['spin period is ' num2str(spin_p) ' sec'])
	end
		
	for pr=[12,34]
		
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
		else ss = '';
		end
		irf_log('proc',sprintf('Raw spectrum w%sE%dp%d -> RSPEC%dp%d',...
			ss,cl_id,probe,cl_id,probe))
		
		fsamp = c_efw_fsample(tt,'hx');
		if ~fsamp, error('no sampling frequency'),end
		
		problems = 'reset|bbias|probesa|probeld|sweep|bdump';
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
				IFFT_Difference = ifft( tt(eind,2) - mean(tt(eind,2)) );
				
				AmpSin = 2*imag(IFFT_Difference);
				AmpCos = 2*real(IFFT_Difference);
				
				rspec(i,2)  =  AmpCos(1+1);
				rspec(i,3)  = -AmpSin(1+1);
				rspec(i,4)  =  AmpCos(2+1);
				rspec(i,5)  = -AmpSin(2+1);
				rspec(i,6)  =  AmpCos(3+1);
				rspec(i,7)  = -AmpSin(3+1);
				rspec(i,8)  =  AmpCos(4+1);
				rspec(i,9)  = -AmpSin(4+1);
				rspec(i,10) =  AmpCos(5+1);
				rspec(i,11) = -AmpSin(5+1);
            else n_gap = n_gap + 1;
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
	
	var_s = sprintf('diEs%dp%d',cl_id, probe_p);
	[ok,diEs,msg] = c_load(var_s);
	if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
	
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
	[Ddsi,Damp] = c_efw_dsi_off(diEs(1,1),cl_id);
	diEs = caa_corof_dsi(diEs,Ddsi,Damp);

	pswake = c_efw_corrot(cl_id,diEs,diBrs,Ps,R,SAX,diV);
	diEs = caa_rm_blankt(diEs,pswake);
	[lowake,dEx] = c_efw_lobewake(cl_id,diEs,diBrs,Ps,R,diEDI); %#ok<ASGLU>
	if ~isempty(dEx)
		cmd = 'DdsiX?=Ddsi+dEx;save mXTRA DdsiX?';
		if exist('./mXTRA.mat','file'), cmd = [cmd ' -append']; end
		c_eval(cmd,cl_id), clear cmd
	end
	eval(irf_ssub(...
		'PSWAKE?p!=pswake;LOWAKE?p!=lowake;save_list=[save_list ''PSWAKE?p! LOWAKE?p!''];',...
		cl_id,probe_p));
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edi (sc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edi')
	
	save_file = './mEDI.mat';
	
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
	
	if strcmp(quantity,'br')
		var_b = 'Br?';
		[ok,E_tmp,msg] = c_load('diE?p1234',cl_id);
		if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
	else
		var_b = 'Brs?'; var_e = {'diEs?p34', 'diEs?p12'};
		[ok,E_tmp] = c_load(var_e{1},cl_id);
		if ~ok
			[ok,E_tmp] = c_load(var_e{2},cl_id);
			if ~ok
				irf_log('load',sprintf('Canot load diEs%d(p12|p34). Please load it.',cl_id))
				data = []; cd(old_pwd); return
			end
		end
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
		B_tmp = irf_tlim(B_tmp,E_tmp(1,1) + [0 dt]);
		if isempty(B_tmp), bad_coverage = 1;
		else
			fgm_sf = 1/(B_tmp(2,1)-B_tmp(1,1));
			del_f = 1.5;
			if (fgm_sf > 22.5 - del_f) && (fgm_sf < 22.5 + del_f), fgm_sf = 22.5;
			elseif (fgm_sf > 67.5 - del_f) && (fgm_sf < 67.5 + del_f), fgm_sf = 67.5;
            else irf_log('proc','cannot guess sampling frequency for B')
			end
			cover = length(B_tmp(:,1))/(dt*fgm_sf);
			% We allow for 10% of data gaps. (should we??)
			if cover < .9, bad_coverage = 1; end
		end
    else bad_coverage = 1; cover = 0;
	end
	
	% Try to use BPP as a backup
	if bad_coverage
		BPP_tmp = c_load('BPP?',cl_id,'var');
		if isempty(BPP_tmp)
			% Use FR data if there is any (cover != 0)
			if cover==0
				irf_log('load','Canot load B. Please load B FGM or B PP.')
				data = []; cd(old_pwd); return
			end
		else
			BPP_tmp = irf_tlim(BPP_tmp,E_tmp(1,1) + [0 dt]);
			if isempty(BPP_tmp)
				irf_log('load','Canot find any usefull B data. Please load B FGM or B PP.')
				data = []; cd(old_pwd); return
			end
	
			fgm_sf = 1/(BPP_tmp(2,1)-BPP_tmp(1,1));
			del_f = .1;
			if (fgm_sf > .25 - del_f) && (fgm_sf < .25 + del_f), fgm_sf = .25;
            else irf_log('proc','cannot guess sampling frequency for B PP')
			end
			cover_pp = length(BPP_tmp(:,1))/(dt*fgm_sf);
			
			% If there is more PP data, then use it.
			% Take .99 to avoid marginal effects.
			if .99*cover_pp > cover
				B_tmp = BPP_tmp;
				Binfo = 'PP'; %#ok<NASGU>
				irf_log('proc','Using B PP to calculate Br')
            else irf_log('proc',sprintf('Use B has %2.2f%% coverage',cover*100))
			end
		end
	end
	
	% Resample the data
	Br = irf_resamp(B_tmp,E_tmp,'fsample',c_efw_fsample(E_tmp)); %#ok<NASGU>
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
        save_file = './mEFWburst.mat';
        do_burst = 1;
    else
        do_burst = 0;
        save_file = './mP.mat';
    end
    
	n_ok = 0;
    if do_burst, param={'180Hz','4kHz','32kHz'};
    else param={'10Hz'};
    end
        
    loaded = 0;
    for k=1:length(param)
        for probe=1:4
            [ok,p] = c_load(irf_ssub(['P' param{k} '?p!'],cl_id,probe));
            if ~loaded && ok, loaded=1; loaded_param = param{k}; end
            c_eval('p?=p;',probe)
            if ~ok || isempty(p)
                irf_log('load', irf_ssub(['No/empty P' param{k} '?p!'],cl_id,probe));
                continue
            end
            n_ok = n_ok + 1;
            clear ok p
        end
        if loaded, break, end
    end
	if ~n_ok, data = []; cd(old_pwd), return, end
	
	if ~do_burst
		% Remove probe saturation
		for probe = 1:4
			if eval(irf_ssub('~isempty(p?)',probe))
				[ok,sa] = c_load(irf_ssub('PROBESA?p!',cl_id,probe));
				if ~ok
					irf_log('load',irf_ssub('Cannot load PROBESA?p!',cl_id,probe))
					continue
				end
				if ~isempty(sa)
					irf_log('proc',['blanking saturated P' num2str(probe)])
					c_eval('if ~isempty(p?), p? = caa_rm_blankt(p?,sa); end',probe)
				end
				clear ok sa
			end
		end
		c_eval('if ~isempty(p?), p?=p?(find(~isnan(p?(:,2))),:); end')

		problems = 'reset|bbias|sweep'; %#ok<NASGU>
		n_ok = 0;
		for probe=1:4
			c_eval('signal=p?;',probe)
			if ~isempty(signal) %#ok<NODEF>
				remove_problems
				res(isnan(res(:,2)),:) = []; %#ok<AGROW>
				c_eval('p?=res;',probe)
				n_ok = n_ok + 1;
			end
		end
		clear res signal problems probe
		if ~n_ok, data = []; cd(old_pwd), return, end
	end
	
	if flag_rmwhip
		problems = 'whip'; %#ok<NASGU>
		for probe=1:4
			c_eval('signal=p?;',probe)
			if ~isempty(signal)
				remove_problems
				c_eval('p?=res;',probe)
			end
		end
	end

	% Check for problem with one probe pair
	MAX_CUT = .1;
	if isempty(p1), l1=0; else l1 = length(find(~isnan(p1(:,2)))); end %#ok<NODEF>
	if isempty(p2), l2=0; else l2 = length(find(~isnan(p2(:,2)))); end %#ok<NODEF>
	if ~isempty(p1) && ~isempty(p2)
		if abs(l1-l2) > MAX_CUT*( max([p1(end,1) p2(end,1)]) - max([p1(1,1) p2(1,1)]) )*5
			if l1>l2, p2=[]; irf_log('proc','throwing away p2')
			else p1=[]; irf_log('proc','throwing away p1')
			end
		end
	end
	if isempty(p3), l3=0; else l3 = length(find(~isnan(p3(:,2)))); end %#ok<NODEF>
	if isempty(p4), l4=0; else l4 = length(find(~isnan(p4(:,2)))); end %#ok<NODEF>
	if ~isempty(p3) && ~isempty(p4)
		if abs(l3-l4)> MAX_CUT*( max([p3(end,1) p4(end,1)]) - max([p3(1,1) p4(1,1)]) )*5
			if l3>l4, p4=[]; irf_log('proc','throwing away p4')
			else p3=[]; irf_log('proc','throwing away p3')
			end
		end
	end
	
	if ~isempty(p1) && all(size(p1)==size(p2)) && all(size(p1)==size(p3)) && ...
            all(size(p1)==size(p4) )
		p = [p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4]; %#ok<NASGU>
		Pinfo.probe = 1234;
		irf_log('proc','computing from p1234')
	elseif ~isempty(p3) && all(size(p3)==size(p4)) && ~(all(size(p1)==size(p2)) && l1>l3)
		p = [p3(:,1) (p3(:,2)+p4(:,2))/2]; %#ok<NASGU>
		Pinfo.probe = 34;
		irf_log('proc','computing from p34')
	elseif ~isempty(p1) && all(size(p1)==size(p2))
		p = [p1(:,1) (p1(:,2)+p2(:,2))/2]; %#ok<NASGU>
		Pinfo.probe = 12;
		irf_log('proc','computing from p12')
	elseif ~isempty(p4) && l4>=max([l1 l2 l3])
		p = p4; %#ok<NASGU>
		Pinfo.probe = 4;
		irf_log('proc','computing from p4')
	elseif ~isempty(p2) && l2>=max([l1 l3])
		p = p2; %#ok<NASGU>
		Pinfo.probe = 2;
		irf_log('proc','computing from p2')
	elseif ~isempty(p3) && l3>=l1
		p = p3; %#ok<NASGU>
		Pinfo.probe = 3;
		irf_log('proc','computing from p3')
	elseif ~isempty(p1)
		p = p1; %#ok<NASGU>
		Pinfo.probe = 1;
		irf_log('proc','computing from p1')
    else irf_log('dsrc','Cannot compute P'), cd(old_pwd); return
	end
	
	c_eval(['P' loaded_param '?=p;save_list=[save_list ''P' loaded_param '? ''];'],cl_id);
    if do_burst
        c_eval('bP?=p;bP?_info=Pinfo;bNVps?=c_efw_scp2ne(p);bNVps?(:,end+1)=p(:,2); save_list=[save_list '' bP? bP?_info bNVps?''];',cl_id)
    else
        c_eval('P?=p;P?_info=Pinfo;NVps?=c_efw_scp2ne(p);NVps?(:,end+1)=p(:,2); save_list=[save_list '' P? P?_info NVps?''];',cl_id)
    end
	clear p p1 p2 p3 p4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P spin resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'ps')
	save_file = './mP.mat';
	
	[ok,P_tmp,msg] = c_load('P?',cl_id);
	if ~ok, irf_log('load',msg), data = []; cd(old_pwd); return, end
	
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
	
	P_tmp = irf_resamp(P_tmp,tvec','fsample',0.25); clear tvec %#ok<NASGU>
	c_eval('Ps?=P_tmp;save_list=[save_list ''Ps? '' ];',cl_id);
	
	[ok,P_info] = c_load('P?_info',cl_id); %#ok<NASGU>
	if ok
		c_eval('Ps?_info=P_info;save_list=[save_list ''Ps?_info '' ];',cl_id);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vce - E CIS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vce')
	save_file='./mCIS.mat';

	c_load('SAX?',cl_id)
	if ~exist('./mCISR.mat','file')
		irf_log('load','Please run ''vcis'' first (mCIS missing)')
		data = []; cd(old_pwd); return
	end
	if ~exist('./mBPP.mat','file')
		irf_log('load','Please run ''b'' first (mBPP missing)')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF DATA MANIPULATIONS
% saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If flag_save is set, save variables to specified file
if flag_save==1 && length(save_file)>0 && ~isempty(save_list)
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
