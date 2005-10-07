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
%   dieburst : dibE{cl_id}p1234 -> mEFWburst // despun ib(8kHz) E [DSI]
%          ADC offsets are NOT corrected
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
%	whip: WHIP{cl_id} -> mFDM	// Whisper pulses present +1 precceding sec
%	bdump: DBUMP{cl_id} -> mFDM	// Burst dump present
%	sweep: SWEEP{cl_id} -> mFDM	// Sweep + dump present
%	badbias: BADBIASRESET{cl_id}, BADBIAS{cl_id}p[1..4] -> mFDM	
%          // Bad bias settings
%	probesa: PROBESA{cl_id}p[1..4] -> mFDM	// Probe saturation
%   edi : EDI{cl_id}, diEDI{cl_id} -> mEDI // EDI E in sc ref frame
%   br, brs : Br[s]{cl_id}, diBr[s]{cl_id} -> mBr // B resampled to E[s]
%   vedbs, vedb : VExB[s]{cl_id}, diVExB[s]{cl_id} -> mEdB // E.B=0 [DSI+GSE]
%	vce : VCE(p,h){cl_id},diVCE(p,h){cl_id} ->mCIS	// E CIS PP [GSE+DSI] 
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

% Copyright 2004,2005 Yuri Khotyaintsev
% Parts of the code are (c) Andris Vaivads

error(nargchk(3,15,nargin))
if nargin > 3, have_options = 1; args = varargin;
else, have_options = 0;
end

% default options
flag_save = 1;
flag_usesavedoff = 1;
flag_edb = 1;
sfit_ver = -1;

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
			else, probe_p_tmp = str2num(args{2});
			end
			if (probe_p_tmp==12 | probe_p_tmp==34), probe_p = probe_p_tmp;
			else, irf_log('fcal,','wrongArgType : probe_p must be 12 or 34')
			end
		else, irf_log('fcal,','wrongArgType : ang_limit value is missing')
		end
	case 'sfit_ver'
		if length(args)>1
			if isnumeric(args{2})
				l = 2;
				if	args{2}>=0 & args{2}<2
					sfit_ver = args{2};
				else, irf_log('fcal,','wrongArgType : sfit_ver must be 0 or 1')
				end
			else, irf_log('fcal,','wrongArgType : sfit_ver must be numeric')
			end
		else, irf_log('fcal,','wrongArgType : sfit_ver value is missing')
		end
	otherwise
		irf_log('fcal,',['Option ''' args{1} '''not recognized'])
	end
	if length(args) > l, args = args(l+1:end);
	else break
	end
end


save_file = '';
save_list = '';

old_pwd = pwd;
cd(cp.sp) %enter the storage directory
if cp.sp~='.', irf_log('save',['Storage directory is ' cp.sp]), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dies - spin fiting of Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'dies')
	save_file = './mEDSI.mat';
	
	if ~c_load('Atwo?',cl_id)
		irf_log('load',...
			irf_ssub('No Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
		data = []; cd(old_pwd); return
	end
	if ~(c_load('wE?p12',cl_id) | c_load('wE?p32',cl_id) | c_load('wE?p34',cl_id)) 
		irf_log('load',...
			irf_ssub('No wE?p12/32 and/or wE?p34 in mER. Use getData(CDB,...,cl_id,''e'')',cl_id))
		data = []; cd(old_pwd); return
	end
	
	if exist(irf_ssub('wE?p32',cl_id),'var'), p12 = 32;
	else, p12 = 12;
	end
	
	% If we have different timelines for p1(3)2 and p34 we try to make 
	% a common timeline, so thet the resulting spinfits will also have 
	% the common timeline
	if exist(irf_ssub(['wE?p' num2str(p12)],cl_id),'var') & exist(irf_ssub('wE?p34',cl_id),'var')
		if eval(irf_ssub('length(wE?p!(:,1))==length(wE?p34(:,1))',cl_id,p12))
			not_same = 0;
			if eval(irf_ssub('length(find((wE?p!(:,1)-wE?p34(:,1))~=0))',cl_id,p12))
				not_same = 1;
			end
		else, not_same = 1;
		end
		if not_same
			c_eval('ts = wE?p34(1,1); te = wE?p34(end,1);',cl_id);
			c_eval(['ts1=wE?p' num2str(p12) '(1,1); te1=wE?p' num2str(p12) '(end,1);'],cl_id);
			if ts>ts1, ts = ts1; end
			if te<te1, te = te1; end
			dt = double(te - ts); clear te ts1 te1
			
			% Guess the sampling frequency
			c_eval('fsamp=c_efw_fsample(wE?p34(:,1));',cl_id);
			if fsamp
				irf_log('proc','Using new time line')
				t_new = double(0):double(1/fsamp):dt+double(1/fsamp); 
				t_new = t_new';
				
				d_new = zeros(length(t_new),2);
				d_new(:,1) = ts + t_new;
				d_new(:,2) = NaN;
				
				c_eval('ii = round((wE?p34(:,1)-ts)*fsamp+1);',cl_id)
				if ii(end)>length(t_new), error('problemo with new time line'), end
				c_eval('d_new(ii,2) = wE?p34(:,2); wE?p34 = d_new;',cl_id)
				
				d_new(:,2) = NaN;
				
				c_eval(['ii = round((wE?p' num2str(p12) '(:,1)-ts)*fsamp+1);'],cl_id)
				if ii(end)>length(t_new), error('problemo with new time line'), end
				c_eval(['d_new(ii,2) = wE?p' num2str(p12) '(:,2); wE?p' num2str(p12) ' = d_new;'],cl_id)
				
				clear t_new ts dt d_new ii
			end
		end
		clear not_same
	end
	
	pl=[p12,34];
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(irf_ssub(['wE?p' ps],cl_id),'var')
			c_eval(['tt=wE?p' ps ';'],cl_id)
			irf_log('proc',sprintf('Spin fit wE%dp%d -> diEs%dp%d',cl_id,pl(k),cl_id,pl(k)))

			c_eval('aa=c_phase(tt(:,1),Atwo?);',cl_id)
			
			% Remove bad bias around EFW reset
			[ok,bbias] = c_load('BADBIASRESET?',cl_id);
			if ok
				irf_log('proc','blanking bad bias due to EFW reset')
				tt = caa_rm_blankt(tt,bbias);
			end
			clear ok bbias
			% Remove bad bias from bias current indication
			for kk = [num2str(ps(1)) num2str(ps(2))]
				if ~exist(irf_ssub('BADBIAS?p!',cl_id,kk),'var')
					c_load(irf_ssub('BADBIAS?p!',cl_id,kk))
				end
				if exist(irf_ssub('BADBIAS?p!',cl_id,kk),'var')
					eval(irf_ssub('bbias=BADBIAS?p!;',cl_id,kk))
					if ~isempty(bbias)
						irf_log('proc',['blanking bad bias on P' num2str(kk)])
						tt = caa_rm_blankt(tt,bbias);
					end
					clear bbias
				end
			end
			
			% Remove probe saturation electronics & low density
			for kk = [num2str(ps(1)) num2str(ps(2))]
				if ~exist(irf_ssub('PROBESA?p!',cl_id,kk),'var')
					c_load(irf_ssub('PROBESA?p!',cl_id,kk))
				end
				if ~exist(irf_ssub('PROBELD?p!',cl_id,kk),'var')
					c_load(irf_ssub('PROBELD?p!',cl_id,kk))
				end
				
				if exist(irf_ssub('PROBESA?p!',cl_id,kk),'var')
					eval(irf_ssub('sa=PROBESA?p!;',cl_id,kk))
					if ~isempty(sa)
						irf_log('proc',['blanking saturated P' num2str(kk)])
						tt = caa_rm_blankt(tt,sa);
					end
					clear sa
				end
				if exist(irf_ssub('PROBELD?p!',cl_id,kk),'var')
					eval(irf_ssub('sa=PROBELD?p!;',cl_id,kk))
					if ~isempty(sa)
						irf_log('proc',...
							['blanking low density saturation on P' num2str(kk)])
						tt = caa_rm_blankt(tt,sa);
					end
					clear sa
				end
			end
			
			% Remove whisper pulses
			if flag_rmwhip
				[ok,whip] = c_load('WHIP?',cl_id);
				if ok
					irf_log('proc','blanking Whisper pulses')
					tt = caa_rm_blankt(tt,whip);
					clear whip
				else
					irf_log('load',...
						irf_ssub('No WHIP? in mFDM. Use getData(CP,cl_id,''whip'')',cl_id))
				end
			end
			
			% Remove sweeps and burst dumps
			[ok,sweep] = c_load('SWEEP?',cl_id);
			if ok
				if ~isempty(sweep)
					irf_log('proc','blanking sweeps')
					tt = caa_rm_blankt(tt,sweep);
					clear sweep
				end
			else
				irf_log('load',...
					irf_ssub(['No SWEEP?. Use getData(CP,cl_id,''sweep'')'],cl_id))
			end
			[ok,bdump] = c_load('BDUMP?',cl_id);
			if ok
				if ~isempty(bdump)
					irf_log('proc','blanking burst dumps')
					tt = caa_rm_blankt(tt,bdump);
					clear bdump
				end
			else
				irf_log('load',...
					irf_ssub(['No BDUMP?. Use getData(CP,cl_id,''bdump'')'],cl_id))
			end
			
			if sfit_ver>=0
				irf_log('proc',['using SFIT_VER=' num2str(sfit_ver)])
				sp = c_efw_sfit(pl(k),3,10,20,tt(:,1),tt(:,2),aa(:,1),...
					aa(:,2),sfit_ver);
			else
				sp = c_efw_sfit(pl(k),3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
			end
			
			% Check if we have any data left
			if isempty(sp)
				irf_log('load',sprintf('No p%d data left for sc%d',pl(k),cl_id))
				continue
			end
			
			% remove point with zero time
			ind = find(sp(:,1)>0);
			if length(ind)<length(sp(:,1))
				irf_log('proc',[num2str(length(sp(:,1))-length(ind)) ' spins removed (bad time)']);
				sp = sp(ind,:);
			end
			
			adc_off = sp(:,[1 4]);
			% warn about points with sdev>.8
			ii = find(sp(:,6)>.8);
			if length(ii)/size(sp,1)>.05,
				irf_log('proc',[sprintf('%.1f',100*length(ii)/size(sp,1)) '% of spins have SDEV>.8 (ADC offsets)']);
			end
			%adc_off(ii,2) = 0;
			adc_off = irf_waverage(adc_off,1/4);
			ii = find(adc_off(:,2)==0);
			adc_off(ii,2) = mean(adc_off(find(abs(adc_off(:,2))>0),2));
			
			sp = sp(:,[1:4 6]);
			sp(:,4) = 0*sp(:,4); % Z component
			
			% remove spins with bad spin fit (obtained E > 10000 mV/m)
			ind = find(abs(sp(:,3))>1e4); sp(ind,:) = [];
			if ind, disp([num2str(length(ind)) ' spins removed due to E>10000 mV/m']);end
			eval(irf_ssub(['diEs?p' ps '=sp;Dadc?p' ps '=adc_off;'],cl_id)); 
			clear tt aa sp adc_off
			eval(irf_ssub(['save_list=[save_list ''diEs?p' ps ' Dadc?p' ps ' ''];'],cl_id));
		else
			irf_log('load',sprintf('No p%d data for sc%d',pl(k),cl_id))
		end
	end

	% Delta offsets
	if (exist(irf_ssub('diEs?p12',cl_id),'var') | ...
	exist(irf_ssub('diEs?p32',cl_id),'var')) & exist(irf_ssub('diEs?p34',cl_id),'var')
		
		% To compute delta offsets we remove points which are > deltaof_sdev_max*sdev
		% as this must de a stable quantity
		eval(irf_ssub('[ii1,ii2] = irf_find_comm_idx(diEs?p!,diEs?p34);',cl_id,p12))
		eval(irf_ssub('df=diEs?p!(ii1,2:3)-diEs?p34(ii2,2:3);',cl_id,p12))
		clear ii1 ii2
		sdev = std(df);
		iia = [];
		for comp = 1:2
			ii = find(abs(df(:,comp)-mean(df(:,comp))) > deltaof_sdev_max*sdev(comp)); 
			iia = [iia; ii];
		end
		iia = sortrows(iia(:));
		iia(find(diff(iia)==0)) = [];
		irf_log('calb',sprintf('%d points are removed for delta offsets',...
				length(iia)))
		for comp = 1:2
			ddd = df(:,comp); ddd(iia) = [];
			Del(comp) = mean(ddd);
		end
			
		irf_log('calb',sprintf('delta offsets are: %.2f [x] %.2f [y]', ...
			Del(1), Del(2)))

		% Check for unreallistically large Del. 
		% If it is larger than deltaof_max, we set it to zero and 
		% NOT doing any corrections.
		if (Del(1)>deltaof_max) | (Del(2)>deltaof_max)
			Del = [0 0];
			irf_log('calb',...
				irf_ssub('DELTA OFFSET TOO BIG >!. Setting D?p12p34=[0 0]',...
				cl_id,deltaof_max))
		else
			% We suppose that smaller field is more realistic
			% and will correct the largest signal.
			% If we have p32, we always correct it, not p34.
			% Real offset is applied to p12, imaginary to p34.
			if Del(1)>0 & p12==12, Del = -Del*j; end
			if real(Del)
				irf_log('calb',irf_ssub('correcting p?',p12))
				eval(irf_ssub('diEs?p!(:,2:3)=diEs?p!(:,2:3)-ones(size(diEs?p!,1),1)*Del;',cl_id,p12));
			else
				irf_log('calb','correcting p34')
				c_eval('diEs?p34(:,2:3)=diEs?p34(:,2:3)-ones(size(diEs?p34,1),1)*imag(Del);',cl_id);
			end
		end
		
		eval(irf_ssub('D?p12p34=Del;',cl_id))
		clear m12 m34 Del

		eval(irf_ssub(['save_list=[save_list ''D?p12p34 ''];'],cl_id));
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die - despin of full resolution data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'die') | strcmp(quantity,'dief') | ...
	strcmp(quantity,'diespec') | strcmp(quantity,'dieburst')
	
	do_burst = 0;
	if strcmp(quantity,'dieburst'), do_burst = 1; else do_burst = 0; end
	if strcmp(quantity,'dief'), do_filter = 1; else do_filter = 0; end
	if do_burst
		save_file = './mEFWburst.mat';
		var_name = 'wbE?p';
		var1_name = 'dibE?p1234';
	else
		if strcmp(quantity,'diespec'), save_file = './mEDSI.mat';
		else, save_file = './mEDSIf.mat';
		end
        var_name = 'wE?p';
		if do_filter, var1_name = 'diEF?p1234';
		else, var1_name = 'diE?p1234';
		end
	end

	if ~c_load('A?',cl_id)
		irf_log('load',...
			irf_ssub('No A? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
		data = []; cd(old_pwd); return
	end
	if ~exist('./mEDSI.mat','file') & ~do_filter & ~strcmp(quantity,'diespec')
		irf_log('load','Please compute spin averages (mEDSI)')
		data = []; cd(old_pwd); return
	end
	
	if do_burst
		if exist('./mEFWburst.mat','file')
			c_eval(['load mEFWburst ' var_name '12 ' var_name '34;'],cl_id);
		else
			irf_log('load','Cannot find mEFWburst.mat')
			data = []; cd(old_pwd); return
		end
	else
		if exist('./mER.mat','file')
			c_eval(['load mER ' var_name '12 ' var_name '32 ' var_name '34;'],cl_id);
		else
			irf_log('load','Cannot find mER.mat')
			data = []; cd(old_pwd); return
		end
	end
	
	% calibration coefficients // see c_efw_despin
	coef=[[1 0 0];[1 0 0]];

	pl=[12,32,34];
	full_e = [];
	n_sig = 0;
	p12 = 12;
	
	% Load Whisper pulses
	whip = [];
	if flag_rmwhip
		[ok,whip] = c_load('WHIP?',cl_id);
		if ~ok
			irf_log('load',...
				irf_ssub('No WHIP? in mFDM. Use getData(CP,cl_id,''whip'')',cl_id))
		end
	end
	
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(irf_ssub([var_name ps],cl_id),'var')
			if pl(k)==32, p12 = 32; end
			n_sig = n_sig + 1;
			if do_burst
				c_eval(['Ep' ps '=' var_name ps ';'],cl_id);
				% Correct ADC offset
				c_load(['Dadc?p' ps],cl_id)
				c_load(['Da?p' ps],cl_id)
				if exist(irf_ssub(['Dadc?p' ps],cl_id),'var')
					c_eval(['irf_log(''calb'',''using saved Dadc?p' ps ''')'],cl_id)
					c_eval(['tmp_adc = irf_resamp(Dadc?p' ps ',Ep' ps ');'],cl_id)
					c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-tmp_adc(:,2);'],cl_id)
					clear tmp_adc
				elseif exist(irf_ssub(['Da?p' ps],cl_id),'var')
					c_eval(['irf_log(''calb'',sprintf(''Da?dp' ps ' (using saved) : %.2f'',Da?p' ps '))'],cl_id)
					c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-Da?p' ps ';'],cl_id)
				else
					irf_log('calb','ADC offset not corrected')
				end
			else
				% Remove bad bias around EFW reset
				[ok,bbias] = c_load('BADBIASRESET?',cl_id);
				if ok
					irf_log('proc','blanking bad bias due to EFW reset')
					c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',bbias);'],cl_id)
				end
				clear ok bbias
				% Remove bad bias from bias current indication
				for k = [num2str(ps(1)) num2str(ps(2))]
					if ~exist(irf_ssub('BADBIAS?p!',cl_id,k),'var')
						c_load(irf_ssub('BADBIAS?p!',cl_id,k))
					end
					if exist(irf_ssub('BADBIAS?p!',cl_id,k),'var')
						eval(irf_ssub('bbias=BADBIAS?p!;',cl_id,k))
						if ~isempty(bbias)
							irf_log('proc',['blanking bad bias on P' num2str(k)])
							c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',bbias);'],cl_id)
						end
						clear bbias
					end
				end
				
				% Remove probe saturation electronics & low density
				for kk = [num2str(ps(1)) num2str(ps(2))]
					if ~exist(irf_ssub('PROBESA?p!',cl_id,kk),'var')
						c_load(irf_ssub('PROBESA?p!',cl_id,kk))
					end
					if ~exist(irf_ssub('PROBELD?p!',cl_id,kk),'var')
						c_load(irf_ssub('PROBELD?p!',cl_id,kk))
					end
					
					if exist(irf_ssub('PROBESA?p!',cl_id,kk),'var')
						eval(irf_ssub('sa=PROBESA?p!;',cl_id,kk))
						if ~isempty(sa)
							irf_log('proc',['blanking saturated P' num2str(kk)])
							c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',sa);'],cl_id)
						end
						clear sa
					end
					if exist(irf_ssub('PROBELD?p!',cl_id,kk),'var')
						eval(irf_ssub('sa=PROBELD?p!;',cl_id,kk))
						if ~isempty(sa)
							irf_log('proc',...
								['blanking low density saturation on P' num2str(kk)])
							c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',sa);'],cl_id)
						end
						clear sa
					end
				end
				
				% Remove sweeps and burst dumps
				[ok,sweep] = c_load('SWEEP?',cl_id);
				if ok
					if ~isempty(sweep)
						irf_log('proc','blanking sweeps')
						c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',sweep);'],cl_id)
						clear sweep
					end
				else
					irf_log('load',...
						irf_ssub(['No SWEEP?. Use getData(CP,cl_id,''sweep'')'],cl_id))
				end
				[ok,bdump] = c_load('BDUMP?',cl_id);
				if ok
					if ~isempty(bdump)
						irf_log('proc','blanking burst dumps')
						c_eval(['wE?p' ps '=caa_rm_blankt(wE?p' ps ',bdump);'],cl_id)
						clear bdump
					end
				else
					irf_log('load',...
						irf_ssub(['No BDUMP?. Use getData(CP,cl_id,''bdump'')'],cl_id))
				end
				
				% Remove Whisper pulses
				if flag_rmwhip & ~isempty(whip)
					irf_log('proc',['blanking Whisper pulses on p' ps])
					c_eval(['if ~isempty(wE?p' ps '),wE?p' ps '=caa_rm_blankt(wE?p' ps ',whip);end'],cl_id)
				end
				
				% Check if we have at least 1 sec of data left
				if eval(irf_ssub(['length(find(~isnan(wE?p' ps '(:,2))))'],cl_id))<25
					irf_log('proc','No data after removals')
					data = []; cd(old_pwd); return
				end
				
				% Correct ADC offset
				if flag_usesavedoff & ~do_filter
					if c_load(['Dadc?p' ps],cl_id)
						c_eval(['irf_log(''calb'',''using saved Dadc?p' ps ''')'],cl_id)
						c_eval(['Ep' ps '=wE?p' ps '; tmp_adc = irf_resamp(Dadc?p' ps ',Ep' ps ');'],cl_id)
						c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-tmp_adc(:,2);'],cl_id)
						clear tmp_adc
					elseif c_load(['Da?p' ps],cl_id)
						c_eval(['disp(sprintf(''Da?dp' ps ' (using saved) : %.2f'',Da?p' ps '))'],cl_id)
						c_eval(['Ep' ps '=wE?p' ps '; Ep' ps '(:,2)=Ep' ps '(:,2)-Da?p' ps ';'],cl_id)
					else, flag_usesavedoff = 0;
					end
				end
				if ~flag_usesavedoff | do_filter
					c_eval(['[Ep' ps ',Da?p' ps ']=caa_corof_adc(wE?p' ps ');'],cl_id)
					c_eval(['irf_log(''calb'',sprintf(''Da?dp' ps ' : %.2f'',Da?p' ps '))'],cl_id)
					c_eval(['save_list=[save_list '' Da?p' ps ' ''];'],cl_id);
				end
			end
		end
	end
	if n_sig==0
		irf_log('load','No raw data found in mER')
		data = []; cd(old_pwd); return
	end
	if n_sig==2
		if p12==32 
			Ep12 = Ep32; clear Ep32
			E_info.probe = '3234';
		else
			E_info.probe = '1234';
		end
		if abs(length(Ep12)-length(Ep34))>0
			% different timelines. Need to correct
			irf_log('proc','using common timeline')
			[ii12,ii34] = irf_find_comm_idx(Ep12,Ep34);
			irf_log('proc',['Ep' num2str(p12) ' ' num2str(length(Ep12)) '->' num2str(length(ii12)) ' data points'])
			Ep12 = Ep12(ii12,:);
			irf_log('proc',['Ep34 ' num2str(length(Ep34)) '->' num2str(length(ii34)) ' data points'])
			Ep34 = Ep34(ii34,:);
		end
		% use WEC coordinate system E=[t,0,p34,p12]
		full_e = zeros(length(Ep12),4);
		full_e(:,[1,4]) = Ep12;
		full_e(:,3) = Ep34(:,2);
		clear Ep12 Ep34
	else
		if exist('Ep12','var')
			pp = 12;
			E_info.probe = '12';
			EE = Ep12;
			clear Ep12
		else
			pp = 34;
			E_info.probe = '34';
			EE = Ep34;
			clear Ep34
		end
		% use WEC coordinate system E=[t,0,p34,p12]
		full_e = zeros(length(EE),4);
		full_e(:,1) = EE(:,1);
		if pp==12, full_e(:,4) = EE(:,2);
		else, full_e(:,3) = EE(:,2);
		end
		clear EE pp
	end

	% load Delta offsets D?p12p34
	if c_load('D?p12p34',cl_id)
		eval(irf_ssub('Del=D?p12p34;',cl_id))
		if real(Del) % Real del means we must correct p12. real(Del)==imag(Del)
			irf_log('calb',['correcting p' num2str(p12)])
			i_c = 1;
		else
			irf_log('calb','correcting p34')
			Del = imag(Del);
			i_c = 2;
		end
		coef(i_c,3) = Del(1) -Del(2)*1j;
		clear Del

	else, irf_log('calb','no Delta offsets found in mEDSI, not doing correction...')
	end
	c_eval([var1_name '_info=E_info;save_list=[save_list ''' var1_name '_info ''];'],cl_id);

	% Do actual despin
	if p12==32, c_eval('full_e=c_efw_despin(full_e,A?,coef,''asym'');',cl_id);
	else, c_eval('full_e=c_efw_despin(full_e,A?,coef);',cl_id);
	end
	
	if strcmp(quantity,'diespec')
		%We make a spectrum and save it.
		sfreq = c_efw_fsample(full_e(:,1));
		if sfreq == 25, nfft = 512;
		else, nfft = 8192;
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
		full_e(:,3)=-full_e(:,3);
		
		c_eval([var1_name '=full_e; save_list=[save_list ''' var1_name '''];'],cl_id);
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idie, idies - DSI inertial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'idie') | strcmp(quantity,'idies')
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
			irf_ssub(['No diV? in mR. Use getData(CDB,...,cl_id,''v'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
	evxb = irf_tappl(irf_cross(diB,irf_resamp(diV,diB)),'*1e-3*(-1)');
	
	err_s = '';
	for k=1:length(var_s)
		[ok,diE] = c_load(var_s{k},cl_id);
		if ~ok
			if isempty(err_s), err_s = var_s{k};
			else, err_s = [err_s ', ' var_s{k}];
			end
			continue 
		end

		enew = diE;
		% We take only X and Y components. Z must remain zero.
		etmp = irf_resamp(evxb(:,1:3),enew(:,1));
		enew(:,2:3) = diE(:,2:3) - etmp(:,2:3);
		clear etmp
		
		c_eval(['i' var_s{k} '= enew; clear enew'],cl_id)
		save_list=[save_list 'i' irf_ssub(var_s{k},cl_id) ' '];
	end
	if ~isempty(err_s)
		irf_log('load',...
			irf_ssub(['No ' err_s ' in mEDSI(f). Use getData(CP,cl_id,''' e_opt ''')'],cl_id))
		data = []; cd(old_pwd); return
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edb,edbs,iedb,iedbs - E.B=0 (sc,inertial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edb') | strcmp(quantity,'edbs') | ...
	strcmp(quantity,'iedb') | strcmp(quantity,'iedbs')
	
	if strcmp(quantity,'iedb') | strcmp(quantity,'iedbs'), inert = 1; 
	else, inert = 0; 
	end
	
	if inert, save_file = './mEdBI.mat';
	else, save_file = './mEdB.mat';
	end
	
	if strcmp(quantity,'edb') | strcmp(quantity,'iedb')
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
				irf_ssub(['No diV? in mR. Use getData(CDB,...,cl_id,''v'')'],cl_id))
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
	if strcmp(quantity,'edbs') | strcmp(quantity,'iedbs'), diE(:,6) = diE(:,5); end
	
	dsiof = c_ctl(cl_id,'dsiof');
	if isempty(dsiof), dsiof = [1+0i 1]; end
	[ok,Dxy] = c_load('Ddsi?',cl_id);
	if ~ok, Dxy = dsiof(1); end
	[ok,Da] = c_load('Damp?',cl_id);
	if ~ok, Da = dsiof(2); end
	clear dsiof

	diE = caa_corof_dsi(diE,Dxy,Da);

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
		diE(:,2:4) = diE(:,2:4) - evxb(:,2:4); clear evxb
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
	eval(irf_ssub('ang_limit?=ang_limit;',cl_id)) 
	save_list=[save_list s 'di' varo_s ' ang_limit' num2str(cl_id) ' '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whip - Whisper present.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'whip')
	save_file = './mFDM.mat';
	
	[ok,fdm] = c_load('FDM?',cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No FDM?. Use getData(CDB,...,cl_id,''fdm'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
	[t_s,t_e,fdm_r] = caa_efw_mode_tab(fdm, 'r');
	ii = find(fdm_r==1);
	
	if ~isempty(ii)
		% add 1 sec before
		t_s = t_s(ii) - 1;
		t_e = t_e(ii);
		c_eval('WHIP?=[double(t_s)'' double(t_e)''];',cl_id); 
		c_eval('save_list=[save_list '' WHIP? ''];',cl_id);
	else
		irf_log('dsrc','No data')
		c_eval('WHIP?=[];save_list=[save_list '' WHIP? ''];',cl_id);
	end
	clear t_s t_e fdm_r ii fdm ok
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bdump - Burst dump
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bdump')
	save_file = './mFDM.mat';
	
	[ok,fdm] = c_load('FDM?',cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No FDM?. Use getData(CDB,...,cl_id,''fdm'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
	[t_s,t_e,fdm_px] = caa_efw_mode_tab(fdm, 'px');
	ii = find(fdm_px(:,1)==1 & fdm_px(:,2)==0);
	
	if ~isempty(ii)
		t_s = t_s(ii);
		% We add one second to the end of the interval for safety
		t_e = t_e(ii) +1;
		c_eval('BDUMP?=[double(t_s)'' double(t_e)''];',cl_id); 
		c_eval('save_list=[save_list '' BDUMP? ''];',cl_id);
	else
		irf_log('dsrc','No data')
		c_eval('BDUMP?=[];save_list=[save_list '' BDUMP? ''];',cl_id);
	end
	clear t_s t_e fdm_px ii fdm ok
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sweep - Sweep present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sweep')
	save_file = './mFDM.mat';
	
	[ok,fdm] = c_load('FDM?',cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No FDM?. Use getData(CDB,...,cl_id,''fdm'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
	[t_s,t_e,fdm_w] = caa_efw_mode_tab(fdm, 'w');
	ii = find(fdm_w==1);

	[t_s_px,t_e_px,fdm_px] = caa_efw_mode_tab(fdm, 'px');
	ii_px = find(fdm_px(:,1)==1 & fdm_px(:,2)==1);
	if ~isempty(ii) | ~isempty(ii_px)
		bdump = zeros(length(ii),2);
		for k=1:length(ii)
			bdump(k,1) = t_s(ii(k));
			% We look for dump of the sweep in the FDM which follows the wseep 
			% or in the next one
			jj = find(t_s_px>=t_e(ii(k)) & t_s_px<t_e(ii(k))+1.1);
			% We add one second to the end of the interval for safety
			if isempty(jj)
				bdump(k,2) = t_e(ii(k)) +1;
				irf_log('dsrc','no dump after sweep')
			else, bdump(k,2) = t_e_px(jj(end)) +1;
			end
		end 
		c_eval('SWEEP?=bdump;save_list=[save_list '' SWEEP? ''];',cl_id);
	else
		irf_log('dsrc','No data')
		c_eval('SWEEP?=[];save_list=[save_list '' SWEEP? ''];',cl_id);
	end
	clear t_s t_e fdm_px ii fdm ok
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% badbias - Bad bias settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'badbias')
	save_file = './mFDM.mat';
	
	DELTA_MINUS = 60;
	DELTA_PLUS = 3*60;
	
	% First we check for EFW resets which usually mean bas bias settings
	% in 2-3 minutes following the reset
	[ok,efwt] = c_load('EFWT?',cl_id);
	if ok
		ii = find(efwt(:,2)<DELTA_PLUS);
		if ~isempty(ii)
			t0 = efwt(ii(1),1) - efwt(ii(1),2);
			irf_log('proc', ['EFW reset at ' epoch2iso(t0,1)]);
			c_eval(['BADBIASRESET?=[double(t0-DELTA_MINUS)'' double(t0+DELTA_PLUS)''];'...
				'save_list=[save_list '' BADBIASRESET? ''];'],cl_id);
		end
	else, irf_log('dsrc','No EFWT')
	end
	clear ok efwt t0 ii
	
	DELTA_MINUS = 300;
	DELTA_PLUS = 64;
	GOOD_BIAS = -130;
	
	c_eval(['i_p? = c_load(''IBIAS' num2str(cl_id) 'p?'',''var'');']);
	if isempty(i_p1) & isempty(i_p2) & isempty(i_p3) & isempty(i_p4)
		irf_log('load',...
			irf_ssub(['No IBIAS?p[1..4]. Use getData(CDB,...,cl_id,''ibias'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	for pro=1:4
		c_eval('ibias=i_p?;',pro)
		% Good & bad points
		if isempty(ibias), continue, end
		ii_bad = find(ibias(:,2)>GOOD_BIAS);
		if isempty(ii_bad)
			irf_log('dsrc','No data')
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
				res = [ibias(ii(1,:))-DELTA_MINUS; ibias(ii(2,:))+DELTA_PLUS]';
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
	save_file = './mFDM.mat';
	
	% Saturation level nA
	SA_LEVEL = 66;
	% N_CONST sets the minimum number of points of constant potential
	% which we consider bad
	N_CONST = 4;
	% DT_PLUMIN is the interval by which we extend saturation 
	% sintervals from each side
	DT_PLUMIN = 4;
	% Delta = .1 sec, half sampling interval for P.
	DELTA = .1;
	
	c_eval(['p? = c_load(''P10Hz' num2str(cl_id) 'p?'',''var'');']);
	if isempty(p1) & isempty(p2) & isempty(p3) & isempty(p4)
		irf_log('load',...
			irf_ssub(['No P10Hz?p[1..4]. Use getData(CDB,...,cl_id,''p'')'],cl_id))
		data = []; cd(old_pwd); return
	end
	for pro=1:4
		c_eval('p=p?;',pro)
		if isempty(p), continue, end
		
		% Points below SA_LEVEL should be excluded from E, but not from
		% P ans they atill contain valuable physical information.
		% This is not the case with points with positive and/or 
		% constant (stuck probe) potential.
		
		% Bad points are points below SA_LEVEL
		ii_bad = find(p(:,2)<-SA_LEVEL);
		ii_god = find(p(:,2)>=-SA_LEVEL);

		if isempty(ii_bad)
			c_eval(['PROBELD' num2str(cl_id) ...
				'p?=[];save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
		elseif isempty(ii_god)
			c_eval(['PROBELD' num2str(cl_id) ...
				'p?=[double(p(1,1))'' double(p(end,1))''];'...
				'PROBESA' num2str(cl_id) 'p?=[];'...
				'save_list=[save_list '' PROBELD' num2str(cl_id) ...
				'p? PROBESA' num2str(cl_id) 'p? ''];'],pro);
			p = [];
		else
			p_tmp = p;
			p_tmp(ii_god,2) = 1;
			p_tmp(ii_bad,2) = 0;
			ii = irf_find_diff(p_tmp(:,2));
			if p_tmp(1,2)==0, ii = [1; ii]; end
			if p_tmp(end,2)==0, ii = [ii; length(p_tmp(:,2))]; end
			ii = reshape(ii,2,length(ii)/2);
			res = [p_tmp(ii(1,:))-DELTA; p_tmp(ii(2,:))+DELTA]';
			c_eval(['PROBELD' num2str(cl_id) 'p?=res;'...
			'save_list=[save_list '' PROBELD' num2str(cl_id) 'p? ''];'],pro);
			clear res ii
			
			% Leave only good points for further exploration
			p = p(ii_god,:);
		end
		clear ii_god ii_bad
		
		if isempty(p), continue, end
		
		% Bad points are points with positive and/or constant potential
		ii_bad = find(p(:,2)>=0);
		ii_god = find(p(:,2)<0);
		if isempty(ii_god)
			c_eval(['PROBESA' num2str(cl_id) ...
				'p?=[double(p(1,1))'' double(p(end,1))''];'...
				'save_list=[save_list '' PROBESA' num2str(cl_id) 'p? ''];'],pro);
		else
			dd = diff(p(ii_god,2));
			p(ii_god,2) = 1;
			
			% check for constant P, means probe is in a strange state
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
			p(ii_bad,2) = 0;
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
			c_eval(['PROBESA' num2str(cl_id) 'p?=res;'...
			'save_list=[save_list '' PROBESA' num2str(cl_id) 'p? ''];'],pro);
			clear ii res
		end
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edi (sc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edi')
	
	save_file = './mEDI.mat';
	
	var_s = 'iEDI?'; e_opt = 'edi';
	varo_s = 'EDI?';
	
	% Load BPP. We use BPP for EDI as it must be a rather approximation
	[ok,B] = c_load('BPP?',cl_id);
	if ~ok
		[ok,B] = c_load('B?',cl_id);
		if ~ok
			irf_log('load',...
				irf_ssub(['No B? and BPP?. Use getData(CDB,...,cl_id,''b'')'],cl_id))
			data = []; cd(old_pwd); return
		end
	end

	% Load V if we need to do SC->Inertial transformation
	[ok,V] = c_load('V?',cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No diV? in mR. Use getData(CDB,...,cl_id,''v'')'],cl_id))
		data = []; cd(old_pwd); return
	end

	% Load E EDI (inertial)
	[ok,E] = c_load(var_s,cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No ' var_s ' in mEDI. Use getData(CP,cl_id,''' e_opt ''')'],cl_id))
		data = []; cd(old_pwd); return
	end

	% SC -> Inertial
	B = irf_resamp(B,E);
	evxb = irf_tappl(irf_cross(B,irf_resamp(V,B)),'*1e-3*(-1)');
	E(:,2:4) = E(:,2:4) + evxb(:,2:4); clear evxb
	
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
elseif strcmp(quantity,'vedb') | strcmp(quantity,'vedbs')
	save_file = './mEdB.mat';

	if strcmp(quantity,'vedb')
		var_s = 'diE?'; e_opt = 'edb';
		varo_s = 'VExB?';
		var_b = 'diBr?'; b_opt ='br';
	else
		var_s = 'diEs?'; e_opt = 'edbs';
		varo_s = 'VExBs?';
		var_b = 'diBrs?'; b_opt ='brs';
	end
	
	% Load resampled B
	[ok,diB] = c_load(var_b,cl_id);
	if ~ok
		irf_log('load',...
			irf_ssub(['No ' var_b ' in mBr. Use getData(CP,cl_id,''' b_opt ''')'],cl_id))
		data = []; cd(old_pwd); return
	end
	
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
elseif strcmp(quantity,'br') | strcmp(quantity,'brs')
	save_file = './mBr.mat';
	
	if strcmp(quantity,'br')
		var_b = 'Br?';
		[ok,E_tmp] = c_load('diE?p1234',cl_id);
		if ~ok
			irf_log('load',sprintf('Canot load diE%dp1234. Please load it.',cl_id))
			data = []; cd(old_pwd); return
		end
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
		Binfo = 'FR';
		bad_coverage = 0;
		cover = 0;
		B_tmp = irf_tlim(B_tmp,E_tmp(1,1) + [0 dt]);
		if isempty(B_tmp), bad_coverage = 1;
		else
			fgm_sf = 1/(B_tmp(2,1)-B_tmp(1,1));
			del_f = 1.5;
			if (fgm_sf > 22.5 - del_f) & (fgm_sf < 22.5 + del_f), fgm_sf = 22.5;
			elseif (fgm_sf > 67.5 - del_f) & (fgm_sf < 67.5 + del_f), fgm_sf = 67.5;
			else, irf_log('proc','cannot guess sampling frequency for B')
			end
			cover = length(B_tmp(:,1))/(dt*fgm_sf);
			% We allow for 10% of data gaps. (should we??)
			if cover < .9, bad_coverage = 1; end
		end
	else, bad_coverage = 1; cover = 0;
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
			if (fgm_sf > .25 - del_f) & (fgm_sf < .25 + del_f), fgm_sf = .25;
			else, irf_log('proc','cannot guess sampling frequency for B PP')
			end
			cover_pp = length(BPP_tmp(:,1))/(dt*fgm_sf);
			
			% If there is more PP data, then use it.
			% Take .99 to avoid marginal effects.
			if .99*cover_pp > cover
				B_tmp = BPP_tmp;
				Binfo = 'PP';
				irf_log('proc','Using B PP to calculate Br')
			else, irf_log('proc',sprintf('Use B has %2.2f%% coverage',cover*100))
			end
		end
	end
	
	% Resample the data
	Br = irf_resamp(B_tmp,E_tmp);
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
elseif strcmp(quantity,'p')
	save_file = './mP.mat';
	
	c_eval(['p?=c_load(''P10Hz' num2str(cl_id) 'p?'',''var'');'])
	
	% Remove bad bias around EFW reset
	[ok,bbias] = c_load('BADBIASRESET?',cl_id);
	if ok
		irf_log('proc','blanking bad bias due to EFW reset')
		c_eval('p? = caa_rm_blankt(p?,bbias);')
	end
	clear ok bbias
	% Remove bad bias from bias current indication
	for kk = 1:4
		if ~exist(irf_ssub('BADBIAS?p!',cl_id,kk),'var')
			c_load(irf_ssub('BADBIAS?p!',cl_id,kk))
		end
		if exist(irf_ssub('BADBIAS?p!',cl_id,kk),'var')
			eval(irf_ssub('bbias=BADBIAS?p!;',cl_id,kk))
			if ~isempty(bbias)
				irf_log('proc',['blanking bad bias on P' num2str(kk)])
				c_eval('p? = caa_rm_blankt(p?,bbias);',kk)
			end
			clear bbias
		end
	end
	clear BADBIAS*
	
	% Remove probe saturation
	for k = 1:4
		if ~isempty(eval(irf_ssub('p?',k)))
			c_load(irf_ssub('PROBESA?p!',cl_id,k))
			if exist(irf_ssub('PROBESA?p!',cl_id,k),'var')
				eval(irf_ssub('sa=PROBESA?p!;',cl_id,k))
				if ~isempty(sa)
					irf_log('proc',['blanking saturated P' num2str(k)])
					% We remove saturations simultaneously on probes 1&2, 3&4
					% so that the resulting timelines are similar
					if k==1|k==2
						if ~isempty(p1), p1 = caa_rm_blankt(p1,sa); end
						if ~isempty(p2), p2 = caa_rm_blankt(p2,sa); end
					else
						if ~isempty(p3), p3 = caa_rm_blankt(p3,sa); end
						if ~isempty(p4), p4 = caa_rm_blankt(p4,sa); end
					end
				end
				clear sa
			end
		end
	end
	c_eval('if ~isempty(p?), p?=p?(find(~isnan(p?(:,2))),:); end')
	
	% Remove sweeps and burst dumps
	[ok,sweep] = c_load('SWEEP?',cl_id);
	if ok
		if ~isempty(sweep)
			irf_log('proc','blanking sweeps')
			c_eval('p?=caa_rm_blankt(p?,sweep);',cl_id)
			clear sweep
		end
	else
		irf_log('load',...
			irf_ssub(['No SWEEP?. Use getData(CP,cl_id,''sweep'')'],cl_id))
	end
	[ok,bdump] = c_load('BDUMP?',cl_id);
	if ok
		if ~isempty(bdump)
			irf_log('proc','blanking burst dumps')
			c_eval('p?=caa_rm_blankt(p?,bdump);',cl_id)
			clear bdump
		end
	else
		irf_log('load',...
			irf_ssub(['No BDUMP?. Use getData(CP,cl_id,''bdump'')'],cl_id))
	end
	
	% Remove Whisper pulses
	if flag_rmwhip
		[ok,whip] = c_load('WHIP?',cl_id);
		if ok & ~isempty(whip)
			irf_log('proc','blanking Whisper pulses')
			c_eval('if ~isempty(p?),p?=caa_rm_blankt(p?,whip);end')
			clear whip
		else
			irf_log('load',...
				irf_ssub('No WHIP? in mFDM. Use getData(CP,cl_id,''whip'')',cl_id))
		end
	end
	
	if size(p1)==size(p2)&size(p1)==size(p3)&size(p1)==size(p4) & size(p1)~=[0 0]
		p = [p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
		Pinfo.probe = 1234;
	elseif size(p3)==size(p4) & size(p3)~=[0 0] & cl_id~=2
		p = [p3(:,1) (p3(:,2)+p4(:,2))/2];
		Pinfo.probe = 34;
	elseif size(p1)==size(p2) & size(p1)~=[0 0]
		p = [p1(:,1) (p1(:,2)+p2(:,2))/2];
		Pinfo.probe = 12;
	elseif size(p4)~=[0 0]
		p = p4;
		Pinfo.probe = 4;
	else, irf_log('dsrc','No data'), cd(old_pwd); return
	end
	
	c_eval('P10Hz?=p;save_list=[save_list ''P10Hz? ''];',cl_id);
	c_eval('P?=p;P?_info=Pinfo;NVps?=c_efw_scp2ne(p);NVps?(:,end+1)=p(:,2); save_list=[save_list '' P? P?_info NVps?''];',cl_id)
	clear p p1 p2 p3 p4
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P spinn resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'ps')
	save_file = './mP.mat';
	
	[ok,P_tmp] = c_load('P?',cl_id);
	if ~ok
		irf_log('load',sprintf('No P? in mP. Use getData(CP,...,cl_id,''p'')',cl_id))
		data = []; cd(old_pwd); return
	end
	
	P_tmp = P_tmp(find(~isnan(P_tmp(:,2))),:);
	
	t0 = '';
	% Try to use time from spin fit
	% TODO: This code can be made smarter.
	[ok,Es_tmp] = c_load('diEs?p34',cl_id);
	if ok
		ii = find(abs(Es_tmp(:,1)-P_tmp(1,1))<2.1);
		if ~isempty(ii)
			irf_log('proc',irf_ssub('using timeline of diEs?p34',cl_id))
			t0 = Es_tmp(ii,1);
		end
	end
	clear Es_tmp
	
	if isempty(t0)
		irf_log('proc','using new timeline')
		if isempty(P_tmp)
		  t0=[];
		else
		  t0 = P_tmp(1,1) + 2; 
		end
	end
	
	if ~isempty(P_tmp),
	  n = floor((P_tmp(end,1)-t0)/4) + 1;
	  tvec = t0 + ( (1:n) -1)*4;
	
      P_tmp = irf_resamp(P_tmp,tvec'); clear tvec
	  c_eval('Ps?=P_tmp;save_list=[save_list ''Ps? '' ];',cl_id);
	
	  [ok,P_info] = c_load('P?_info',cl_id);
	  if ok
	    c_eval('Ps?_info=P_info;save_list=[save_list ''Ps?_info '' ];',cl_id);
	  end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vce - E CIS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vce')
	save_file='./mCIS.mat';

	if exist('./mEPH.mat','file'), eval(irf_ssub('load mEPH SAX?',cl_id)), end
	if ~exist('./mCIS.mat','file')
		irf_log('load','Please run ''vcis'' first (mCIS missing)')
		data = []; cd(old_pwd); return
	end
	if ~exist('./mBPP.mat','file')
		irf_log('load','Please run ''b'' first (mBPP missing)')
		data = []; cd(old_pwd); return
	end

	CIS = load('mCIS');
	% Load BPP. We use BPP for EDI as it must be a rather approximation
	[ok,B] = c_load('BPP?',cl_id);
	if ~ok
		[ok,B] = c_load('B?',cl_id);
		if ~ok
			irf_log('load',...
				irf_ssub(['No B? and BPP?. Use getData(CDB,...,cl_id,''b'')'],cl_id))
			data = []; cd(old_pwd); return
		end
	end

	vars = {'VCp', 'VCh'};
	varo = {'VCEp', 'VCEh'};
	for va=1:length(vars)
		eval(irf_ssub(['if isfield(CIS,''' vars{va} '?''); v=CIS.' vars{va} '?; else, v=[]; end; clear ' vars{va}], cl_id));
		if ~isempty(v)
			evxb=irf_tappl(irf_cross(v,B),'*(-1e-3)');
			eval(irf_ssub([varo{va} '?=evxb;save_list=[save_list '' ' varo{va} '?'']; clear evxb'],cl_id));
			if ~exist(irf_ssub('SAX?',cl_id),'var')
				irf_log('load','must fetch spin axis orientation (option ''sax'')')
			else
				eval(irf_ssub(['di' varo{va} '?=c_gse2dsi(' varo{va} '?,SAX?);save_list=[save_list '' di' varo{va} '?''];'],cl_id));
			end
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF DATA MANIPULATIONS
% saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If flag_save is set, save variables to specified file
if flag_save==1 & length(save_file)>0 & ~isempty(save_list)
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
			if exist(sl{k}), eval(['data{k+1}=' sl{k} ';']); end
		end
	end
end

cd(old_pwd)
