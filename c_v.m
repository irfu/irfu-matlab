function v=c_v(t);
%C_V   Calculate velocity from timing between 4 spacecraft
%
% v=c_v(t);
% dt=c_v(v);
%
% Calculate velocity from timing between 4 spacecraft
% or calculate timing from the velocity
% t=[t1 t2 t3 t4]; in isdat_epoch units
% v=[t vx vy vz]; t in isdat_epoch, v in GSE ref frame,
% dt=[0 t2-t1 t3-t1 t4-t1];
%
% $Id$

sc_list=1:4;

if t(2) > 1e8, flag='v_from_t'; else, flag='dt_from_v';v=t;t=v(1);end

if exist('./mR.mat','file'),
    load mR R1 R2 R3 R4 V1 V2 V3 V4;
else,
    disp('loading position from isdat');
    DATABASE='disco:10';db = Mat_DbOpen(DATABASE);
    for ic=sc_list, disp(['...R' num2str(ic)]);
     [t,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     eval(irf_ssub('R?=[double(t) double(data)''];',ic));clear t data;
    end
	for ic=sc_list, disp(['...V' num2str(ic)]);
       [t,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
       eval(irf_ssub('V?=[double(t) double(data)''];',ic));clear t data;
    end
end

if strcmp(flag,'v_from_t'),
  t_center=0.5*t(1)+0.5*t;
  for ic=1:4,eval(irf_ssub('vsc?=av_interp(V?,t_center,''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('drsc?=av_interp(irf_add(1,R?,-1,R1),t(?),''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('dr?=drsc?+[0 (t(?)-t(1))*vsc?(1,2:4)];',ic));end
  for ic=1:4,eval(irf_ssub('dt(?)=t(?)-t(1);sdt?=num2str(dt(?),3);',ic));end
  D=[dr2(2:4);dr3(2:4);dr4(2:4)];
  T=[dt(2),dt(3), dt(4)]';
  m=D\T;
  clear v
  v=m/norm(m)/norm(m);v=v';	% velocity vector of the boundary

  disp([ datestr(datenum(fromepoch(t(1))))])
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  vn=irf_norm(v);
  strv=['V=' num2str(av_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s GSE'];
  disp(strdt);disp(strv);
elseif strcmp(flag,'dt_from_v'),
  t_center=0.5*t(1)+0.5*t;
  for ic=1:4,eval(irf_ssub('vsc?=av_interp(V?,t_center,''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('v?=v(2:4)-dot(vsc?(2:4),v(2:4)).*v(2:4)./norm(v(2:4))^2;',ic));end
  for ic=1:4,eval(irf_ssub('dr?=av_interp(irf_add(1,R?,-1,R1),t,''spline'');',ic));end
  for ic=1:4,eval(irf_ssub('dt(?)=irf_dot(dr?,v?,1)./norm(v?)^2;',ic));end
  
  % print result
  disp([ datestr(datenum(fromepoch(t(1))))])
  vn=irf_norm(v);
  strv=['V=' num2str(av_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s GSE'];
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  disp(strv);  disp(strdt);
  v=dt; % output variable is v
end

