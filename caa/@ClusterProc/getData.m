function data = getData(cp,cl_id,quantity,varargin)
%GETDATA(cdb) produce Cluster level 0 data from the raw data
% data = getData(cdb,cl_id,quantity,options)
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
%	usesavedoff : use saved offsets instead of recalculate everything
%
% $Revision$  $Date$
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
ang_limit = 10;

while have_options
	l = 1;
    switch(args{1})
    case 'nosave'
        flag_save = 0;
	case 'usesavedoff'
		flag_usesavedoff = 1;
	case 'ang_limit'
		if length(args)>1
			if isnumeric(args{2})
				ang_limit = args{2};
				l = 2;
			else
				warning('caa:wrongArgType','ang_limit must be numeric')
			end
		else
			warning('caa:wrongArgType','ang_limit value is missing')
		end
	case 'ang_blank'
	   flag_edb = 1;	% [default]
	case 'ang_fill'
	   flag_edb = 2;	% fill points below ang_limit with 1e27
	   fill_val = 1e27;
	case 'ang_ez0'
	   flag_edb = 0;	% use Ez=0 for points below ang_limit
    otherwise
        disp(['Option ''' args{i} '''not recognized'])
    end
	if length(args) > l, args = args(l+1:end);
	else break
	end
end


save_file = '';
save_list = '';

old_pwd = pwd;
cd(cp.sp) %enter the storage directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dies - spin fiting of Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'dies')
	save_file = './mEDSI.mat';

	if ~(exist('./mA.mat','file') & exist('./mER.mat','file'))
		warning('Please load raw data (mER) and phase (mA)')
		data = [];
		return
	end

	eval(av_ssub('load mER wE?p12 wE?p34;',cl_id));
	eval(av_ssub('load mA A?;',cl_id));

	pl=[12,34];
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(av_ssub(['wE?p' ps],cl_id),'var')
			eval(av_ssub(['tt=wE?p' ps ';aa=A?;'],cl_id))
			disp(sprintf('Spin fit wE%dp%d -> diEs%dp%d',cl_id,pl(k),cl_id,pl(k)))
			sp = EfwDoSpinFit(pl(k),3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
			sp = sp(:,1:4);
			sp(:,4) = 0*sp(:,4); % Z component
			ind=find(abs(sp(:,3))>1e4);sp(ind,:)=[]; % remove spins with bad spin fit (obtained E > 10000 mV/m)
			if ind, disp([num2str(length(ind)) ' spinds removed due to E>10000 mV/m']);end
			eval(av_ssub(['diEs?p' ps '=sp;'],cl_id)); clear tt aa sp
			eval(av_ssub(['save_list=[save_list '' diEs?p' ps ' ''];'],cl_id));
		else
			disp(sprintf('No p%d data for sc%d',pl(k),cl_id))
		end
	end

	% delta offsets
	if exist(av_ssub('diEs?p12',cl_id),'var') & exist(av_ssub('diEs?p34',cl_id),'var')
		eval(av_ssub(['m12=mean(diEs?p12);m34=mean(diEs?p34);'],cl_id))
		Del = m12(2:3) - m34(2:3);

		disp(sprintf('delta offsets are: %.2f [x] %.2f [y]', ...
		abs(Del(1)), abs(Del(2))))

		% we suppose that smaller field is more realistic
		% and will correct the largest signal
		% real offset is applied to p12, imaginary to p34
		if abs(m34)>abs(m12), Del = -Del*j; end
		eval(av_ssub('D?p12p34=Del;',cl_id))

		if real(Del)
			disp('correcting p12')
			eval(av_ssub('diEs?p12(:,2:3)=diEs?p12(:,2:3)-ones(length(diEs?p12),1)*Del;',cl_id));
		else
			disp('correcting p34')
			Del = imag(Del);
			eval(av_ssub('diEs?p34(:,2:3)=diEs?p34(:,2:3)-ones(length(diEs?p34),1)*Del;',cl_id));
		end
		clear m12 m34 Del

		eval(av_ssub(['save_list=[save_list '' D?p12p34 ''];'],cl_id));
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% die - despin of full resolution data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'die')
	save_file = './mEDSI.mat';

	if ~(exist('./mA.mat','file') & exist('./mA.mat','file') & exist('./mEDSI.mat','file')) % ?????????????
		warning('Please compute spin averages (mER) and load phase (mA)')
		data = [];
		return
	end
	eval(av_ssub('load mER wE?p12 wE?p34;',cl_id));

	% calibration coefficients // see c_despin
	coef=[[1 0 0];[1 0 0]];

	pl=[12,34];
	full_e = [];
	n_sig = 0;
	for k=1:length(pl)
		ps = num2str(pl(k));
		if exist(av_ssub(['wE?p' ps],cl_id),'var')
			n_sig = n_sig + 1;
			% correct ADC offset
			if flag_usesavedoff & exist('./mEDSI.mat','file')
			   	eval(av_ssub(['load mEDSI Da?p' ps ],cl_id))
			end
			if exist(av_ssub(['Da?p' ps],cl_id),'var')
				eval(av_ssub(['disp(sprintf(''Da?dp' ps ' (using saved) : %.2f'',Da?p' ps '))'],cl_id))
				eval(av_ssub(['Ep' ps '=wE?p' ps '; Ep' ps '(:,2)=Ep' ps '(:,2)-Da?p' ps ';'],cl_id))    % ????????????????????????
			else
				eval(av_ssub(['[Ep' ps ',Da?p' ps ']=corrADCOffset(cp,wE?p' ps ');'],cl_id))
				eval(av_ssub(['disp(sprintf(''Da?dp' ps ' : %.2f'',Da?p' ps '))'],cl_id))
				eval(av_ssub(['save_list=[save_list '' Da?p' ps ' ''];'],cl_id));
			end
		end
	end
	if n_sig==0
		warning('No raw data found in mER')
		data = [];
		return
	end
	if n_sig==2
		if abs(length(Ep12)-length(Ep34))>0
			% different timelines. Need to correct
			[ii12,ii34] = findCommInd(Ep12,Ep34);
			Ep12 = Ep12(ii12,:);
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

	% load Delta offsets D?p12p34
	if exist('./mEDSI.mat','file')
		eval(av_ssub('load mEDSI D?p12p34;',cl_id));
	end
	if exist(av_ssub('D?p12p34',cl_id))
		eval(av_ssub('Del=D?p12p34;',cl_id))
		if real(Del)                               % ?????????????????????????
			disp('correcting p12')
			i_c = 1;
		else
			disp('correcting p34')
			Del = imag(Del);
			i_c = 2;
		end
		eval(av_ssub('coef(i_c,3)=Del(1)-Del(2)*1j;',cl_id));
		clear Del

	else, disp('no Delta offsets found in mEDSI, not doing correction...')
	end

	% Do actual despin
	eval(av_ssub('load mA A?;',cl_id));
	eval(av_ssub('diE?p1234=c_despin(full_e,A?,coef);',cl_id));
	% DS-> DSI
	eval(av_ssub('diE?p1234(:,3)=-diE?p1234(:,3);',cl_id));
	eval(av_ssub(['save_list=[save_list '' diE?p1234 ''];'],cl_id));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edb,edbs - E.B=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'edb') | strcmp(quantity,'edbs')
	save_file = './mEdB.mat';

	if ~(exist('./mEDSI.mat','file') & exist('./mBPP.mat','file'))
		warning('caa:noSuchFile','Please despin E (mEDSI) and load B PP (mBPP)')
		data = [];
		return
	end

	eval(av_ssub('load mBPP diBPP?; diB=diBPP?;',cl_id));

	if strcmp(quantity,'edb')
		var_s = av_ssub('diE?p1234',cl_id);
		varo_s = av_ssub('E?',cl_id);
	else
		disp('using p34')
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
		warning('caa:noData','Please despin E (no diE in mEDSI)')
		data = [];
		return
	end
	if exist(Dxy_s,'var'), eval(['Dx=real(' Dxy_s ');Dy=imag(' Dxy_s ');'])
	else, disp('using Dx,Dy=0'), Dx = 0; Dy=0;
	end
	if exist(Da_s,'var'), eval(['Da=' Da_s ';'])
	else, disp('using Da=1'), Da = 1;
	end

	diE = corrDSIOffsets(diE,Dx,Dy,Da);

	disp(['using angle limit of ' num2str(ang_limit) ' degrees'])
	[diE,angle]=av_ed(diE,diB,ang_limit);
	diE(:,5) = angle; clear angle

	ii = find(abs(diE(:,5)) < ang_limit);
	if length(ii) > 1
		switch(flag_edb)
		case 0 % Ez=0, do nothing
			disp('using Ez=0')
		case 1 % remove points
			disp('settiong points < ang_limit to NaN')
			diE(ii,4) = diE(ii,4)*NaN;
		case 2 % fill with fill_val
			disp('settiong points < ang_limit to 1e27')
			diE(ii,4) = ones(size(diE(ii,4)))*fill_val;
		end
	end

	% DSI->GSE
 	if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
	if ~exist(av_ssub('SAX?',cl_id),'var')
		warning('must fetch spin axis orientation (option ''sax'')')
	else
		eval(av_ssub([varo_s '=c_gse2dsc(diE(:,1:4),SAX?,-1);' varo_s '(:,3:4)=-' varo_s '(:,3:4);' varo_s '(:,5)=diE(:,5);save_list=[save_list ''' varo_s ' ''];'],cl_id));
	end

	eval(['di' varo_s '=diE;']); clear diE
	save_list=[save_list 'di' varo_s ' ang_limit '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vedb,Vedbs - E.B=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vedb') | strcmp(quantity,'vedbs')
	save_file = './mEdB.mat';

	if ~(exist('./mEdB.mat','file') & exist('./mBPP.mat','file'))
		warning('caa:noSuchFile','Please calculate Ez (mEdB)')
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
		warning('caa:noData','Please calculate Ez (no diE in mEDSI)')
		data = [];
		return
	end

	save_list=[save_list 'di' varo_s ' '];

	% DSI->GSE
 	if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
	if ~exist(av_ssub('SAX?',cl_id),'var')
		warning('must fetch spin axis orientation (option ''sax'')')
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
	disp([save_list ' -> ' save_file])
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
