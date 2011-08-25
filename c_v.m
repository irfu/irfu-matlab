function v=c_v(t,coord_sys)
%C_V   Calculate velocity from timing between 4 spacecraft
%
% v=c_v(t,[coord_sys]);
% dt=c_v(v,[coord_sys]);
%
% Calculate velocity from timing between 4 spacecraft
% or calculate timing from the velocity
% t=[t1 t2 t3 t4]; in isdat_epoch units
% v=[t vx vy vz]; t in isdat_epoch, v in GSE ref frame,
% dt=[0 t2-t1 t3-t1 t4-t1];
%
% coord_sys='GSM' - when calculate in GSM reference frame 
%
% $Id$

sc_list=1:4;
if nargin==1, coord_sys='GSE'; end 

if t(2) > 1e8, flag='v_from_t'; else flag='dt_from_v';v=t;t=v(1);end

if exist('CAA/C1_CP_AUX_POSGSE_1M','dir')==7, % checks if exist CAA data (STUPID SOLUTION)
    irf_log('dsrc','Trying to read CAA files...')
    c_eval('[~,~,R?]=c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'');');
    c_eval('[~,~,V?]=c_caa_var_get(''sc_v_xyz_gse__C?_CP_AUX_POSGSE_1M'');');
elseif exist('./mR.mat','file'),
    load mR R1 R2 R3 R4 V1 V2 V3 V4;
else
    disp('loading position from isdat');
		DB_S = c_ctl(0,'isdat_db');
    db = Mat_DbOpen(DB_S);
    for ic=sc_list, disp(['...R' num2str(ic)]);
     [tr,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     eval(irf_ssub('R?=[double(tr) double(data)''];',ic));clear tr data;
    end
	for ic=sc_list, disp(['...V' num2str(ic)]);
       [tv,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
       eval(irf_ssub('V?=[double(tv) double(data)''];',ic));clear tv data;
    end
		Mat_DbClose(db);
end

switch coord_sys
    case 'GSE'
        % do nothing
    case 'GSM'
        c_eval('R?=irf_gse2gsm(R?);V?=irf_gse2gsm(V?);');
    otherwise
        % do nothing, i.e. assume GSE
end

if strcmp(flag,'v_from_t'),
  t_center=0.5*t(1)+0.5*t;
  c_eval('vsc?=irf_resamp(V?,t_center,''spline'');');
  c_eval('drsc?=irf_resamp(irf_add(1,R?,-1,R1),t(?),''spline'');');
  c_eval('dr?=drsc?+[0 (t(?)-t(1))*vsc?(1,2:4)];');
  c_eval('dt(?)=t(?)-t(1);sdt?=num2str(dt(?),3);');
  D=[dr2(2:4);dr3(2:4);dr4(2:4)];
  T=[dt(2),dt(3), dt(4)]';
  m=D\T;
  clear v
  v=m/norm(m)/norm(m);v=v';	% velocity vector of the boundary

  disp([ datestr(datenum(fromepoch(t(1))))])
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  vn=irf_norm(v);
  strv=['V=' num2str(irf_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s ' coord_sys];
  disp(strdt);disp(strv);
elseif strcmp(flag,'dt_from_v'),
  t_center=0.5*t(1)+0.5*t;
  for ic=1:4,eval(irf_ssub('vsc?=irf_resamp(V?,t_center,''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('v?=v(2:4)-dot(vsc?(2:4),v(2:4)).*v(2:4)./norm(v(2:4))^2;',ic));end
  for ic=1:4,eval(irf_ssub('dr?=irf_resamp(irf_add(1,R?,-1,R1),t,''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('dt(?)=irf_dot(dr?,v?,1)./norm(v?)^2;',ic));end
  
  % print result
  disp([ datestr(datenum(fromepoch(t(1))))])
  vn=irf_norm(v);
  strv=['V=' num2str(irf_abs(v,1),3) '*[ ' num2str(vn(end-2:end),' %5.2f') '] km/s ' coord_sys];
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  disp(strv);  disp(strdt);
  v=dt; % output variable is v
end

