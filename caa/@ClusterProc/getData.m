function data = getData(cp,cl_id,quantity,varargin)
%GETDATA(cp) produce Cluster level 0 data from the raw data
% data = getData(cp,cl_id,quantity,options)
%
% Input:
%	cp - ClusterProc object
%	cl_id - SC#
%	quantity - one of the following:
%
%	dies : diEs{cl_id}p12, diEs{cl_id}p34 -> mEDSI // spin fits [DSI]
%		also creates delta offsets D{cl_id}p12p34.
%		If the offset is real then it must be applied to p12,
%		if imaginary - to p34
%	die : diE{cl_id}p1234 -> mEDSI // despun full res E [DSI]
%		also created ADC offsets Da{cl_id}p12 and Da{cl_id}p34
%	dieburst : dibE{cl_id}p1234 -> mEFWburst // despun ib(8kHz) E [DSI]
%		ADC offsets are NOT corrected
%	edbs, edb : E[s]{cl_id}, diE[s]{cl_id} -> mEdB // Ez from E.B=0 [DSI+GSE]
%		has the following options:
%		ang_limit - minimum angle(B,spin plane) [default 10 deg]
%		ang_blank - put Ez to NaN for points below ang_limit [default]
%		ang_fill - fill points below ang_limit with 1e27
%		ang_ez0 - use Ez=0 for points below ang_limit
%	vedbs, vedb : VExB[s]{cl_id}, diVExB[s]{cl_id} -> mEdB // E.B=0 [DSI+GSE]
%
%	Example usage: getData(cp,4,'edbs','ang_fill','ang_limit',20)
%
%	options - one of the following:
%	nosave : do no save on disk
%	leavewhip : do not remove time intervals with Whisper pulses
%	notusesavedoff : recalculating everything instead of using saved offsets
%
%
% $Id$
%
% see also C_GET

% Copyright 2004 Yuri Khotyaintsev
% Parts of the code are (c) Andris Vaivads

error(nargchk(3,15,nargin))
if nargin > 3, have_options = 1; args = varargin;
else, have_options = 0;
end

% default options
flag_save = 1;
flag_usesavedoff = 0;
flag_edb = 1;
flag_rmwhip = 1; 
ang_limit = 10;

while have_options
	l = 1;
    switch(args{1})
    case 'nosave'
        flag_save = 0;
	case 'leavewhip'
		flag_rmwhip = 0;
	case 'notusesavedoff'
		flag_usesavedoff = 0;
	case 'ang_limit'
		if length(args)>1
			if isnumeric(args{2})
				ang_limit = args{2};
				l = 2;
			else
				c_log('fcal,','wrongArgType : ang_limit must be numeric')
			end
		else
			c_log('fcal,','wrongArgType : ang_limit value is missing')
		end
	case 'ang_blank'
	   flag_edb = 1;	% [default]
	case 'ang_fill'
	   flag_edb = 2;	% fill points below ang_limit with 1e27
	   fill_val = 1e27;
	case 'ang_ez0'
	   flag_edb = 0;	% use Ez=0 for points below ang_limit
    otherwise
        c_log('fcal,',['Option ''' args{i} '''not recognized'])
    end
	if length(args) > l, args = args(l+1:end);
	else break
	end
end


save_file = '';
save_list = '';

old_pwd = pwd;
cd(cp.sp) %enter the storage directory
c_log('save',['Storage directory is ' cp.sp])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dies - spin fiting of Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'dies')
	save_file = './mEDSI.mat';

	if ~(c_load(av_ssub('wE?p12',cl_id)) & c_load(av_ssub('wE?p12',cl_id)) & ...
	c_load(av_ssub('A?',cl_id)))
		c_log('load','Please load raw data (mER) and phase (mA)')
		data = [];
		return
	end

	pl=[12,34];
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(av_ssub(['wE?p' ps],cl_id),'var')
			c_eval(['tt=wE?p' ps ';aa=A?;'],cl_id)
			c_log('proc',sprintf('Spin fit wE%dp%d -> diEs%dp%d',cl_id,pl(k),cl_id,pl(k)))

			if flag_rmwhip
				if exist('./mFDM.mat','file')
					c_eval('load mFDM WHIP?',cl_id)
				end
				if exist(av_ssub('WHIP?',cl_id),'var')
					c_log('proc','not using times with Whisper pulses')
					c_eval('tt=blankTimes(tt,WHIP? );clear WHIP?',cl_id)
				end
			end
			
			sp = EfwDoSpinFit(pl(k),3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2),0);
			
			% remove point with zero time
			ind = find(sp(:,1)>0);
			if length(ind)<length(sp(:,1))
				c_log('proc',[num2str(length(sp(:,1))-length(ind)) ' spins removed (bad time)']);
				sp = sp(ind,:);
			end
			
			adc_off = sp(:,[1 4]);
			% warn about points with sdev>.8
			ii = find(sp(:,6)>.8);
			if length(ii)/size(sp,1)>.05,
				c_log('proc',[sprintf('%.1f',100*length(ii)/size(sp,1)) '% of spins have SDEV>.8 (ADC offsets)']);
			end
			%adc_off(ii,2) = 0;
			adc_off = wAverage(adc_off,1/4);
			ii = find(adc_off(:,2)==0);
			adc_off(ii,2) = mean(adc_off(find(abs(adc_off(:,2))>0),2));
			
			sp = sp(:,1:4);
			sp(:,4) = 0*sp(:,4); % Z component
			
			% remove spins with bad spin fit (obtained E > 10000 mV/m)
			ind = find(abs(sp(:,3))>1e4); sp(ind,:) = [];
			if ind, disp([num2str(length(ind)) ' spins removed due to E>10000 mV/m']);end
			eval(av_ssub(['diEs?p' ps '=sp;Dadc?p' ps '=adc_off;'],cl_id)); 
			clear tt aa sp adc_off
			eval(av_ssub(['save_list=[save_list ''diEs?p' ps ' Dadc?p' ps ' ''];'],cl_id));
		else
			c_log('load',sprintf('No p%d data for sc%d',pl(k),cl_id))
		end
	end

	% delta offsets
	if exist(av_ssub('diEs?p12',cl_id),'var') & exist(av_ssub('diEs?p34',cl_id),'var')
		eval(av_ssub(['m12=mean(diEs?p12);m34=mean(diEs?p34);'],cl_id))
		Del = m12(2:3) - m34(2:3);

		c_log('calb',sprintf('delta offsets are: %.2f [x] %.2f [y]', ...
		abs(Del(1)), abs(Del(2))))

		% we suppose that smaller field is more realistic
		% and will correct the largest signal
		% real offset is applied to p12, imaginary to p34
		if abs(m34)>abs(m12), Del = -Del*j; end
		eval(av_ssub('D?p12p34=Del;',cl_id))

		if real(Del)
			c_log('calb','correcting p12')
			eval(av_ssub('diEs?p12(:,2:3)=diEs?p12(:,2:3)-ones(length(diEs?p12),1)*Del;',cl_id));
		else
			c_log('calb','correcting p34')
			Del = imag(Del);
			eval(av_ssub('diEs?p34(:,2:3)=diEs?p34(:,2:3)-ones(length(diEs?p34),1)*Del;',cl_id));
		end
		clear m12 m34 Del

		eval(av_ssub(['save_list=[save_list ''D?p12p34 ''];'],cl_id));
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die - despin of full resolution data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'die') | strcmp(quantity,'dieburst')
	if strcmp(quantity,'dieburst'), do_burst = 1; else do_burst = 0; end
	if do_burst
		save_file = './mEFWburst.mat';
                var_name = 'wbE?p';
		var1_name = 'dibE?p1234';
	else
		save_file = './mEDSI.mat';
                var_name = 'wE?p';
		var1_name = 'diE?p1234';
	end

	if ~(exist('./mA.mat','file') & exist('./mEDSI.mat','file'))
		c_log('load','Please compute spin averages (mER) and load phase (mA)')
		data = [];
		return
	end
	if do_burst, c_eval(['load mEFWburst ' var_name '12 ' var_name '34;'],cl_id);
	else, c_eval(['load mER ' var_name '12 ' var_name '34;'],cl_id);
	end
	
	% calibration coefficients // see c_despin
	coef=[[1 0 0];[1 0 0]];

	pl=[12,34];
	full_e = [];
	n_sig = 0;
	
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(av_ssub([var_name ps],cl_id),'var')
			n_sig = n_sig + 1;
			if do_burst
				c_eval(['Ep' ps '=' var_name ps ';'],cl_id);
				% correct ADC offset
				if exist('./mEDSI.mat','file')
					eval(av_ssub(['load mEDSI Da?p' ps ' Dadc?p' ps],cl_id))
				end
				if exist(av_ssub(['Dadc?p' ps],cl_id),'var')
					c_eval(['c_log(''calb'',''using saved Dadc?p' ps ''')'],cl_id)
					c_eval(['tmp_adc = av_interp(Dadc?p' ps ',Ep' ps ');'],cl_id)
					c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-tmp_adc(:,2);'],cl_id)
					clear tmp_adc
				elseif exist(av_ssub(['Da?p' ps],cl_id),'var')
					c_eval(['c_log(''calb'',sprintf(''Da?dp' ps ' (using saved) : %.2f'',Da?p' ps '))'],cl_id)
					c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-Da?p' ps ';'],cl_id)
				else
					c_log('calb','ADC offset not corrected')
				end
			else
				% correct ADC offset
				if flag_usesavedoff & exist('./mEDSI.mat','file')
					eval(av_ssub(['load mEDSI Da?p' ps ' Dadc?p' ps],cl_id))
				end
				if exist(av_ssub(['Dadc?p' ps],cl_id),'var')
					c_eval(['c_log(''calb'',''using saved Dadc?p' ps ''')'],cl_id)
					c_eval(['Ep' ps '=wE?p' ps '; tmp_adc = av_interp(Dadc?p' ps ',Ep' ps ');'],cl_id)
					c_eval(['Ep' ps '(:,2)=Ep' ps '(:,2)-tmp_adc(:,2);'],cl_id)
					clear tmp_adc
				elseif exist(av_ssub(['Da?p' ps],cl_id),'var')
					c_eval(['disp(sprintf(''Da?dp' ps ' (using saved) : %.2f'',Da?p' ps '))'],cl_id)
					c_eval(['Ep' ps '=wE?p' ps '; Ep' ps '(:,2)=Ep' ps '(:,2)-Da?p' ps ';'],cl_id)
				else
					if flag_rmwhip & exist('./mFDM.mat','file')
						c_eval('load mFDM WHIP?',cl_id)
					end
					if flag_rmwhip & exist(av_ssub('WHIP?',cl_id),'var')
						%removing times with Whisper pulses
						c_eval(['[Ep' ps ',Da?p' ps ']=corrADCOffset(wE?p' ps ',WHIP?);clear WHIP?'],cl_id)
					else
						c_eval(['[Ep' ps ',Da?p' ps ']=corrADCOffset(wE?p' ps ');'],cl_id)
					end
					c_eval(['c_log(''calb'',sprintf(''Da?dp' ps ' : %.2f'',Da?p' ps '))'],cl_id)
					c_eval(['save_list=[save_list '' Da?p' ps ' ''];'],cl_id);
				end
			end
		end
	end
	if n_sig==0
		c_log('load','No raw data found in mER')
		data = [];
		return
	end
	if n_sig==2
		if abs(length(Ep12)-length(Ep34))>0
			% different timelines. Need to correct
			c_log('proc','using common timeline')
			[ii12,ii34] = findCommInd(Ep12,Ep34);
			c_log('proc',['Ep12 ' num2str(length(Ep12)) '->' num2str(length(ii12)) ' data points'])
			Ep12 = Ep12(ii12,:);
			c_log('proc',['Ep34 ' num2str(length(Ep34)) '->' num2str(length(ii34)) ' data points'])
			Ep34 = Ep34(ii34,:);
		end
		% use WEC coordinate system E=[t,0,p34,p12]
		full_e = zeros(length(Ep12),4);
		%keyboard
		full_e(:,[1,4]) = Ep12;
		full_e(:,3) = Ep34(:,2);
		clear Ep12 Ep34
	else
		if exist('Ep12','var')
			pp = 12;
			EE = Ep12;
			clear Ep12
		else
			pp = 34;
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

	if ~do_burst
		% load Delta offsets D?p12p34
		if exist('./mEDSI.mat','file')
			eval(av_ssub('load mEDSI D?p12p34;',cl_id));
		end
		if exist(av_ssub('D?p12p34',cl_id))
			eval(av_ssub('Del=D?p12p34;',cl_id))
			if real(Del)                               % ?????????????????????????
				c_log('calb','correcting p12')
				i_c = 1;
			else
				c_log('calb','correcting p34')
				Del = imag(Del);
				i_c = 2;
			end
			eval(av_ssub('coef(i_c,3)=Del(1)-Del(2)*1j;',cl_id));
			clear Del
	
		else, c_log('calb','no Delta offsets found in mEDSI, not doing correction...')
		end
	end

	% Do actual despin
	c_eval('load -mat mA.mat A?;',cl_id);
	c_eval([var1_name '=c_despin(full_e,A?,coef);'],cl_id);
	% DS-> DSI
	c_eval([var1_name '(:,3)=-' var1_name '(:,3);'],cl_id);
	c_eval(['save_list=[save_list ''' var1_name '''];'],cl_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edb,edbs - E.B=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edb') | strcmp(quantity,'edbs')
	save_file = './mEdB.mat';

	if ~(exist('./mEDSI.mat','file') & exist('./mBPP.mat','file'))
		c_log('load','Please despin E (mEDSI) and load B PP (mBPP)')
		data = [];
		return
	end

	eval(av_ssub('load mBPP diBPP?; diB=diBPP?;',cl_id));

	if strcmp(quantity,'edb')
		var_s = av_ssub('diE?p1234',cl_id);
		varo_s = av_ssub('E?',cl_id);
	else
		c_log('proc','using p34')
		var_s = av_ssub('diEs?p34',cl_id);
		varo_s = av_ssub('Es?',cl_id);
	end

	Dxy_s =  av_ssub('Ddsi?',cl_id);
	Dx_s =  av_ssub('real(Ddsi?)',cl_id);
	Dy_s =  av_ssub('imag(Ddsi?)',cl_id);
	Da_s =  av_ssub('Damp?',cl_id);

	eval(['load mEDSI ' var_s ' ' Dxy_s ' ' Da_s])
	if exist(var_s,'var'), eval(['diE=' var_s ';'])
	else
		c_log('load','Please despin E (no diE in mEDSI)')
		data = [];
		return
	end
	if exist(Dxy_s,'var'), eval(['Dx=real(' Dxy_s ');Dy=imag(' Dxy_s ');'])
	else, c_log('calb','using Dx,Dy=0'), Dx = 0; Dy=0;
	end
	if exist(Da_s,'var'), eval(['Da=' Da_s ';'])
	else, disp('using Da=1'), Da = 1;
	end

	diE = corrDSIOffsets(diE,Dx,Dy,Da);

	c_log('proc',['using angle limit of ' num2str(ang_limit) ' degrees'])
	[diE,angle]=av_ed(diE,diB,ang_limit);
	diE(:,5) = angle; clear angle

	ii = find(abs(diE(:,5)) < ang_limit);
	if length(ii) > 1
		switch(flag_edb)
		case 0 % Ez=0, do nothing
			c_log('proc','using Ez=0')
		case 1 % remove points
			c_log('proc','setting points < ang_limit to NaN')
			diE(ii,4) = diE(ii,4)*NaN;
		case 2 % fill with fill_val
			c_log('proc','setting points < ang_limit to 1e27')
			diE(ii,4) = ones(size(diE(ii,4)))*fill_val;
		end
	end

	% DSI->GSE
 	if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
	if ~exist(av_ssub('SAX?',cl_id),'var')
		c_log('load','must fetch spin axis orientation (option ''sax'')')
	else
		eval(av_ssub([varo_s '=c_gse2dsc(diE(:,1:4),SAX?,-1);' varo_s '(:,3:4)=-' varo_s '(:,3:4);' varo_s '(:,5)=diE(:,5);save_list=[save_list ''' varo_s ' ''];'],cl_id));
	end

	eval(['di' varo_s '=diE;']); clear diE
	eval(av_ssub('ang_limit?=ang_limit;',cl_id)) 
	save_list=[save_list 'di' varo_s ' ang_limit' num2str(cl_id) ' '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vedb,Vedbs = ExB with E.B=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vedb') | strcmp(quantity,'vedbs')
	save_file = './mEdB.mat';

	if ~(exist('./mEdB.mat','file') & exist('./mBPP.mat','file'))
		c_log('load','Please calculate Ez (mEdB)')
		data = [];
		return
	end

	eval(av_ssub('load mBPP diBPP?;',cl_id));

	if strcmp(quantity,'vedb')
		var_s = av_ssub('diE?',cl_id);
		varo_s = av_ssub('VExB?',cl_id);
	else
		var_s = av_ssub('diEs?',cl_id);
		varo_s = av_ssub('VExBs?',cl_id);
	end

	eval(['load mEdB ' var_s])
	if exist(var_s,'var')
		eval(av_ssub(['di' varo_s '=av_e_vxb(' var_s '(:,1:4),diBPP?,-1);'],cl_id))
		eval(['di' varo_s '(:,5)=' var_s '(:,5);'])
	else
		c_log('load','Please calculate Ez (no diE in mEDSI)')
		data = [];
		return
	end

	save_list=[save_list 'di' varo_s ' '];

	% DSI->GSE
 	if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
	if ~exist(av_ssub('SAX?',cl_id),'var')
		c_log('load','must fetch spin axis orientation (option ''sax'')')
	else
		eval(av_ssub([varo_s '=c_gse2dsc(di' varo_s '(:,1:4),SAX?,-1);' varo_s '(:,3:4)=-' varo_s '(:,3:4);' varo_s '(:,5)=di' varo_s '(:,5);save_list=[save_list ''' varo_s ' ''];'],cl_id));
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
	c_log('save',[save_list ' -> ' save_file])
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
