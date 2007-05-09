function out=c_ctl(varargin)
%C_CTL Cluster control
%
% c_ctl('init',[dir])
% c_ctl('get',cl_id,ctl_name)
% c_ctl(cl_id,ctl_name) % equvalent to get
% c_ctl('list',[sc_list])
% c_ctl('load_ns_ops',[dir])
% c_ctl('save',[dir])
% c_ctl('set',[sc_list],ctl_name,ctl_val)
% c_ctl(sc_list,ctl_name,value,[ctl_name1,value1...]) % equvalent to set
% 
% $Id$

% Copyright 2005 Yuri Khotyaintsev

default_mcctl_path = '.';

if nargin<1, c_ctl_usage, return, end
args = varargin;

if ischar(args{1})
	if strcmp(args{1},'init')
		if nargin>=2
			d = args{2};
			if exist(d,'dir')
				if exist([d '/mcctl.mat'],'file')
					clear global c_ct; global c_ct
					eval(['load ' d '/mcctl.mat'])
				else
					error(['No mcctl.mat in ' d])
				end
			elseif exist(d,'file')
				clear global c_ct; global c_ct
				eval(['load ' d])
			else
				error(['Directory or file ' d ' does not exist'])
			end
		else
			clear global c_ct; global c_ct
			if exist([default_mcctl_path '/mcctl.mat'],'file')
				eval(['load ' default_mcctl_path '/mcctl.mat'])
			else
				% init 
				% NOTE: this is the place to initialize defaults
				def_ct.ns_ops = [];
				def_ct.ang_lim = 10;	% algle limit for E.B=0
				def_ct.rm_whip = 1;		% remove times with WHI pulses
				def_ct.probe_p = 34;	% default probe pair to use
				def_ct.deltaof_max = 1.5;	
										% maximum reasonable value of deltaof
				def_ct.deltaof_sdev_max = 2; 
										% delta offsets we remove points which 
										% are > deltaof_sdev_max*sdev
				
				% DSI offsets
				% Values used:
				% 2002 : .33 .69 .73 .92
				% 2003 : .15 .53 .47 .71
				% 2004 01-04 : .23 .72 .50 .92
				def_ct.dsiof = [.23+0i 1.1];
				c_ct{1} = def_ct;
				def_ct.dsiof = [.72+0i 1.1];
				c_ct{2} = def_ct;
				def_ct.dsiof = [.55+0i 1.1];
				c_ct{3} = def_ct;
				def_ct.dsiof = [.92+0i 1.1];
				c_ct{4} = def_ct;
				clear def_ct
				
				% cell number 5 has global settings
				% this cell must be accessed as SC # 0
				def_ct.isdat_db = 'db:10';
				def_ct.data_path = '/data/cluster';
				def_ct.caa_mode = 0;
				c_ct{5} = def_ct;
			end
		end
		
	elseif strcmp(args{1},'get')
		if nargin<3, error('get: must be c_ctl(''get'',cl_id,''ct_name''))'), end
		cl_id = args{2}; if cl_id==0, cl_id = 5; end
		c = args{3};
		
		global c_ct
		if isempty(c_ct)
			irf_log('fcal','CTL is not initialized. Initializing...') 
			c_ctl('init') 
			global c_ct
		end
		
		if isfield(c_ct{cl_id},c)
			if nargout>0
				eval(['out=c_ct{cl_id}.' c ';'])
			else
				if cl_id < 5, disp(['C' num2str(cl_id) '->' c ':'])
                else disp(['GLOBAL->' c ':'])
				end
				eval(['disp(c_ct{cl_id}.' c ');'])
			end
		else
			if nargout>0, out=[]; end
			irf_log('fcal',['unknown ctl: ' c])
		end

	elseif strcmp(args{1},'load_ns_ops')
		global c_ct
		if isempty(c_ct)
			irf_log('fcal','CTL is not initialized. Initializing...') 
			c_ctl('init') 
			global c_ct
		end
		
		if nargin>1, d = args{2};
		else d = '.';
		end
		
		for j=1:4
			try 
				f_name = [d '/ns_ops_c' num2str(j) '.dat'];
				if exist(f_name,'file')
					c_ct{j}.ns_ops = load(f_name,'-ascii');
					
					% remove lines with undefined dt
					c_ct{j}.ns_ops(find(c_ct{j}.ns_ops(:,2)==-157),:) = [];
				else irf_log('load',['file ' f_name ' not found'])
				end
			catch
				disp(lasterr)
			end
		end
	
	elseif strcmp(args{1},'list')
		if nargin>=2 
			sc_list = args{2};
			ii = find(sc_list==0);
			if ~isempty(ii), sc_list(ii) = 5; end
		else sc_list=1:5;
		end
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		for cl_id=sc_list
			c = fieldnames(c_ct{cl_id});
			if length(c)>0
				for j=1:length(c)
					if cl_id < 5, disp(['C' num2str(cl_id) '->' c{j} ':'])
					else disp(['GLOBAL->' c{j} ':'])
					end
					eval(['disp(c_ct{cl_id}.' c{j} ');'])
				end
			end
		end
		
	elseif strcmp(args{1},'save')
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		
		if nargin>1
			d = args{2};
			if exist(d,'dir')
				disp(['Saving ' d '/mcctl.mat'])
				eval(['save -MAT ' d '/mcctl.mat c_ct'])
			else
				disp(['Saving ' d])
				eval(['save -MAT ' d ' c_ct'])
			end
		else
			disp('Saving mcctl.mat')
			save -MAT mcctl.mat c_ct
		end
		
	elseif strcmp(args{1},'set')
		if nargin<3, error('set: must be c_ctl(''set'',''ctl_name'',value))'), end
		if isnumeric(args{2})
			sc_list = args{2};
			ii = find(sc_list==0);
			if ~isempty(ii), sc_list(ii) = 5; end
			c = args{3};
			c_val = args{4};
		else
			sc_list = 1:4;
			c = args{2};
			c_val = args{3};
		end
		if ~ischar(c), error('ctl_name must be a string'), end
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		for cl_id=sc_list
			try
				eval(['c_ct{cl_id}.' c '=c_val;'])
			catch
				disp(lasterr)
				error('bad option')
			end	
			if cl_id < 5, disp(['C' num2str(cl_id) '->' c ':'])
			else disp(['GLOBAL->' c ':'])
			end
			eval(['disp(c_ct{cl_id}.' c ');'])
		end
		
	else
		error('Invalid argument')
	end
elseif isnumeric(args{1})
	sc_list = args{1};
	ii = find(sc_list==5);
	if ~isempty(ii), sc_list(ii) = 5; end
	
	if nargin>2, have_options = 1; args = args(2:end);
	elseif nargin==2
		have_options = 0;
		if nargout>0, out=c_ctl('get',sc_list,args{2});
		else c_ctl('get',sc_list,args{2});
		end
	else have_options = 0;
	end
	
	while have_options
		if length(args)>1
			if ischar(args{1})
				c_ctl('set',sc_list,args{1},args{2})
			else
				error('option must be a string')
			end
			if length(args) >= 2
				args = args(3:end);
				if isempty(args), break, end
			else break
			end
		else
			disp('Usage: c_ctl(sc_list,''ctl_name'',value)')
			break
		end
	end
else
	error('Invalid argument')
end

function c_ctl_usage
	disp('Usage:')
	disp('  c_ctl(''init'')')
	disp('  c_ctl(''init'',''/path/to/mcctl.mat/'')')
	disp('  c_ctl(''init'',''/path/to/alternative_mcctl.mat'')')
	disp('  c_ctl(''load_ns_ops'')')
	disp('  c_ctl(''load_ns_ops'',''/path/to/ns_ops_cN.dat/'')')
	disp('  c_ctl(''save'')')
	disp('  c_ctl(''save'',''/path/to/mcctl.mat/'')')
	disp('  c_ctl(''save'',''/path/to/alternative_mcctl.mat'')')
	disp('  c_ctl(''sc_list'',''ctl'',value)')
