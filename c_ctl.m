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

default_mcctl_path = '.';

if nargin<1, c_ctl_usage, return, end
args = varargin;

if isstr(args{1})
	if strcmp(args{1},'init')
		if nargin>=2
			d = args{2};
			if exist(d,'dir')
				if exist([d '/mcctl.mat'])
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
			if exist([default_mcctl_path '/mcctl.mat'])
				eval(['load ' default_mcctl_path '/mcctl.mat'])
			else
				% init with defaults
				def_ct.ns_ops = [];
				def_ct.ang_lim = 15;
				def_ct.rm_whip = 1;
				c_ct{1} = def_ct;
				c_ct{2} = def_ct;
				c_ct{3} = def_ct;
				c_ct{4} = def_ct;
				clear def_ct
			end
		end
		
	elseif strcmp(args{1},'get')
		if nargin<3, error('get: must be c_ctl(''get'',cl_id,''ct_name''))'), end
		cl_id = args{2};
		c = args{3};
		
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		
		if isfield(c_ct{cl_id},c)
			if nargout>0
				eval(['out=c_ct{cl_id}.' c ';'])
			else
				disp(['C' num2str(cl_id) '->' c ':'])
				eval(['disp(c_ct{cl_id}.' c ');'])
			end
		else
			error(['unknown ctl: ' c])
		end

	elseif strcmp(args{1},'load_ns_ops')
		global c_ct
		if isempty(c_ct), c_ctl('init'), end
		
		if nargin>1, d = args{2};
		else, d = '.';
		end
		
		for j=1:4
			try 
				f_name = [d '/ns_ops_c' num2str(j) '.dat'];
				if exist(f_name,'file')
					eval(['c_ct{j}.ns_ops=load(''' f_name ''',''-ascii'');'])
					
					% remove lines with undefined dt
					c_ct{j}.ns_ops(find(c_ct{j}.ns_ops(:,2)==-157),:) = [];
				end
			catch
				disp(lasterr)
			end
		end
	
	elseif strcmp(args{1},'list')
		if nargin>=2, sc_list = args{2};
		else, sc_list=1:4;
		end
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		for cl_id=sc_list
			c = fieldnames(c_ct{cl_id});
			if length(c)>0
				for j=1:length(c)
					disp(['C' num2str(cl_id) '->' c{j} ':'])
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
			c = args{3};
			c_val = args{4};
		else
			sc_list = 1:4;
			c = args{2};
			c_val = args{3};
		end
		if ~isstr(c), error('ctl_name must be a string'), end
		global c_ct
		if isempty(c_ct), disp('CTL is not initialized.'), return, end
		for cl_id=sc_list
			try
				eval(['c_ct{cl_id}.' c '=c_val;'])
			catch
				disp(lasterr)
				error('bad option')
			end	
			disp(['C' num2str(cl_id) '->' c ':'])
			eval(['disp(c_ct{cl_id}.' c ');'])
		end
		
	else
		error('Invalid argument')
	end
elseif isnumeric(args{1})
	sc_list = args{1};
	
	if nargin>2, have_options = 1; args = args(2:end);
	elseif nargin==2
		have_options = 0;
		if nargout>0, out=c_ctl('get',sc_list,args{2});
		else, c_ctl('get',sc_list,args{2});
		end
	else, have_options = 0;
	end
	
	while have_options
		if length(args)>1
			if isstr(args{1})
				c_ctl('set',sc_list,args{1},args{2})
			else
				error('option must be a string')
			end
			if length(args) >= 2
				args = args(3:end);
				if length(args) == 0, break, end
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
