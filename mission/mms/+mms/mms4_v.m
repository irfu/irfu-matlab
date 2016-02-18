function v=mms4_v(t,coord_sys)
%MMS4_V   Calculate velocity from timing between 4 spacecraft
% NOT VALID FOR Cluster! (use c_v instead)
%
% v =mms4_v(t,[coord_sys]);
% dt=mms4_v(v,[coord_sys]);
%
% Calculate velocity from timing between 4 spacecraft
% or calculate timing from the velocity
% t=[t1 t2 t3 t4]; in EpochUnix
% v=[t vx vy vz]; t in EpochUnix, v in GSE ref frame,
% dt=[0 t2-t1 t3-t1 t4-t1];
%
% coord_sys='gsm' - when calculate in GSM reference frame
%

persistent R
persistent V

if ~exist('R','var') || isempty(R)
	R=struct('time',[],'gseR1',[],'gseR2',[],'gseR3',[],'gseR4',[],...
        'gsmR1',[],'gsmR2',[],'gsmR3',[],'gsmR4',[]);
       
	tt=[];temp=[]; % needed because we have nested functions
end

if ~exist('V','var') || isempty(R)
	V=struct('time',[],'gseV1',[],'gseV2',[],'gseV3',[],'gseV4',[],...
        'gsmV1',[],'gsmV2',[],'gsmV3',[],'gsmV4',[]);
       
	tt=[];temp=[]; % needed because we have nested functions
end


if nargin==1, coord_sys='gse'; end

if t(2) > 1e8,
	flag='v_from_t';
	tint = [min(t)-61 max(t)+61];
else flag='dt_from_v';
	v=t;
	t=v(1);
	tint = [t-61 t+61];
end






if ~is_R_ok && exist('./mmsR.mat','file'),
	load ./mmsR.mat
    load ./mmsV.mat
    irf.log('warning','mms position loaded from current folder (./mmsR.dat)')
end

if ~is_R_ok && exist('/data/mms/irfu/mmsR.mat','file'),
	load /data/mms/irfu/mmsR.mat
    load /data/mms/irfu/mmsV.mat
    irf.log('warning','mms position loaded from local repository folder (./data/mms/irfu/mmsR.mat)')
end

if ~is_R_ok && exist('mmsR.mat','file'),
	load mmsR.mat
    load mmsV.mat
        irf.log('warning','mms position loaded from somewhere in the path (mmsR.mat)')
end

% if ~is_R_ok
% 	irf.log('warning','Trying to read CAA files C?_CP_AUX_POSGSE_1M...')
% 	var = {'sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C2_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C3_CP_AUX_POSGSE_1M','sc_r_xyz_gse__C4_CP_AUX_POSGSE_1M',...
% 		'sc_v_xyz_gse__C1_CP_AUX_POSGSE_1M','sc_v_xyz_gse__C2_CP_AUX_POSGSE_1M','sc_v_xyz_gse__C3_CP_AUX_POSGSE_1M','sc_v_xyz_gse__C4_CP_AUX_POSGSE_1M'};
% 	ttt=c_caa_var_get(var,'mat','tint',tint);
% 	R.R1 = ttt{1}; R.R2 = ttt{2}; R.R3 = ttt{3}; R.R4 = ttt{4}; 
% 	R.V1 = ttt{5}; R.V2 = ttt{6}; R.V3 = ttt{7}; R.V4 = ttt{8};
% end

% if ~is_R_ok && exist('CAA/CL_SP_AUX','dir')==7,
% 	irf.log('warning','Trying to read CAA files CL_CP_AUX ...')
% 	R.R=irf_get_data('sc_r_xyz_gse__CL_SP_AUX','caa','mat');
% 	if ~isempty(R.R)
% 		c_eval('R.dR?=irf_get_data(''sc_dr?_xyz_gse__CL_SP_AUX'',''caa'',''mat'');');
% 		c_eval('R.R?=irf_add(1,R.R,1,R.dR?);');
% 		R.V=irf_get_data('sc_v_xyz_gse__CL_SP_AUX','caa','mat');
% 		irf.log('warning','!!!! Assumes all s/c move with the same velocity !!!');
% 		c_eval('R.V?=R.V;');
% 	end
% end

% if ~is_R_ok
% 	irf.log('warning','Streaming s/c position from CAA');
% 	read_RV_from_caa_stream;
% end

% if ~is_R_ok
% 	disp('loading position from isdat. Connecting ...');
% 	tmin=min(t);tmax=max(t);
% 	for ic='1234'
% 		[tt,temp] = irf_isdat_get(['Cluster/' ic '/ephemeris/position'], tmin-30, tmax-tmin+30);
% 		R.(['R' ic])=[tt temp];
% 		fprintf('%s',R.(['R' ic]));
% 		[tt,temp] = irf_isdat_get(['Cluster/' ic '/ephemeris/velocity'], tmin-30, tmax-tmin+30);
% 		R.(['V' ic])=[tt temp];
% 		fprintf('%s',R.(['V' ic]));
% 	end
% 	disp('');
% end

if ~is_R_ok
	irf.log('warning','!!! Could not obtain position data !!!');
	return
end

%remove NaN
ind1 = find(isnan(R.gseR1(:,1)));
ind2 = find(isnan(R.gseR2(:,1)));
ind3 = find(isnan(R.gseR3(:,1)));
ind4 = find(isnan(R.gseR4(:,1)));

s1=numel(ind1);
s2=numel(ind2);
s3=numel(ind3);
s4=numel(ind4);
[aaa bbb]=max([s1 s2 s3 s4]);%pick the most restrictive
switch bbb
    case 1
        ind=ind1;
    case 2
        ind=ind2;
    case 3 
        ind=ind3;
    case 4 
        ind=ind4;
end
R.gseR1(ind,:)=[];
R.gsmR1(ind,:)=[];
R.gseR2(ind,:)=[];
R.gsmR2(ind,:)=[];
R.gseR3(ind,:)=[];
R.gsmR3(ind,:)=[];
R.gseR4(ind,:)=[];
R.gsmR4(ind,:)=[];
R.time(ind)=[];

ind1 = find(isnan(V.gseV1(:,1)));
ind2 = find(isnan(V.gseV2(:,1)));
ind3 = find(isnan(V.gseV3(:,1)));
ind4 = find(isnan(V.gseV4(:,1)));
[aaa bbb]=max([s1 s2 s3 s4]);%pick the most restrictive
switch bbb
    case 1
        ind=ind1;
    case 2
        ind=ind2;
    case 3 
        ind=ind3;
    case 4 
        ind=ind4;
end
V.gseV1(ind,:)=[];
V.gseV2(ind,:)=[];
V.gseV3(ind,:)=[];
V.gseV4(ind,:)=[];
V.time(ind)=[];
V.gsmV1(ind,:)=[];
V.gsmV2(ind,:)=[];
V.gsmV3(ind,:)=[];
V.gsmV4(ind,:)=[];


if strcmp(flag,'v_from_t'),
	t_center=0.5*t(1)+0.5*t;
	for ic='1234',
		i=ic-'0';
R.(['vsc'  ic]) = irf_resamp([irf_time(V.time,'ttns>epoch') V.([coord_sys 'V' ic])],t_center','spline');
		R.(['drsc' ic]) = irf_resamp([irf_time(R.time,'ttns>epoch') R.([coord_sys 'R' ic])-R.([coord_sys 'R' num2str(1)])],t(i),'spline');
		R.(['dr'   ic]) = R.(['drsc' ic])+[0 (t(i)-t(1))*R.(['vsc' ic])(1,2:4)];
		R.dt(i)         = t(i)-t(1);
		R.(['sdt' ic])  = num2str(R.dt(i),3);
	end
	D=[R.dr2(2:4);R.dr3(2:4);R.dr4(2:4)];
	T=[R.dt(2),R.dt(3), R.dt(4)]';
	m=D\T;
	clear v
	v=m/norm(m)/norm(m);v=v';	% velocity vector of the boundary
	%if strcmpi(coord_sys,'gsm'), v = irf_gse2gsm([t(1) v]); v(1) = []; end
	disp([ datestr(datenum(fromepoch(t(1))))])
	strdt=['dt=[' , num2str(R.dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
	vn=irf_norm(v);
	strv=['V=' num2str(irf_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s ' coord_sys];
	disp(strdt);disp(strv);

elseif strcmp(flag,'dt_from_v'),
  vOrig = v;
  %if strcmpi(coord_sys,'gsm'), v = irf_gse2gsm([t(1) v], -1); v(1)=[]; end 
	t_center=0.5*t(1)+0.5*t;
	for ic='1234',
		i=ic-'0';
		R.(['vsc'  ic]) = irf_resamp([irf_time(V.time,'ttns>epoch') V.([coord_sys 'V' ic])],t_center','spline');
		R.(['v' ic])   = v(2:4)-dot(R.(['vsc' ic])(2:4),v(2:4)).*v(2:4)./norm(v(2:4))^2;
		R.(['dr' ic])=irf_resamp([irf_time(R.time,'ttns>epoch') R.([coord_sys 'R' ic])-R.([coord_sys 'R' num2str(1)])],t,'spline');
		R.dt(i)=irf_dot(R.(['dr' ic]),R.(['v' ic]),1)./norm(R.(['v' ic]))^2;
	end
	% print result
	disp([ datestr(datenum(fromepoch(t(1))))])
	vn=irf_norm(vOrig);
	strv=['V=' num2str(irf_abs(vOrig,1),3) '*[ ' num2str(vn(end-2:end),' %5.2f') '] km/s ' coord_sys];
	strdt=['dt=[' , num2str(R.dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
	disp(strv);  disp(strdt);
	v=R.dt; % output variable is v
end


%% %%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

	function answer=is_R_ok(sc)
		% check if position data are ok for spacecraft number 'sc'
		% if input argument not given check if ok for all spacecraft that needs
		% to be plotted.
		if nargin == 0,
			scList = 1:4;
		else
			scList = sc;
		end
		for iSc=scList
                strSc = [coord_sys 'R' num2str(iSc)];
            if numel(R.(strSc)) < 8 % less than 2 time points
				answer=false;
				return;
			else
				tintR=[irf_time(R.time(1),'ttns>epoch') irf_time(R.time(end),'ttns>epoch')];
				if(tintR(1)>min(t)) || (tintR(2)<max(t)),
					answer=false;
					return;
                end
			end
		end
		answer=true;
    end
    
% 	function read_RV_from_caa_stream
% 		currentDir = pwd;
% 		tempDir = tempname;
% 		mkdir(tempDir);
% 		cd(tempDir);
% 		caa_download([min(t)-60,max(t)+60],'CL_SP_AUX','stream');
% 		cd('CAA/CL_SP_AUX');
% 		d=dir('*.cef.gz');
% 		cefFile = d.name;
% 		R.R = c_caa_cef_var_get('sc_r_xyz_gse',cefFile);
% 		for sc='1234'
% 			tempR = c_caa_cef_var_get(['sc_dr' sc '_xyz_gse'],cefFile);
% 			R.(['R' sc])=R.R+[zeros(size(R.R,1),1) tempR(:,2:end)];
% 		end
% 		R.V = c_caa_cef_var_get('sc_v_xyz_gse',cefFile);
% 		irf.log('warning','!!!! Assumes all s/c move with the same velocity !!!');
% 		c_eval('R.V?=R.V;');
% 		cd(currentDir);
% 		rmdir(tempDir,'s');
% 	end
% 
% end
end