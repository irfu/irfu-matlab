function v=c_v(t);
%function v=c_v(t);
%function dt=c_v(v);
%
% Calculate velocity from timing between 4 spacecraft
% or calculate timing from the velocity
% t=[t1 t2 t3 t4]; in isdat_epoch units
% v=[t vx vy vz]; t in isdat_epoch, v in GSE ref frame,
% dt=[0 t2-t1 t3-t1 t4-t1];
%

sc_list=1:4;

if t(2) > 1e8, flag='v_from_t'; else, flag='dt_from_v';v=t;t=v(1);end

if exist('mR.mat'),
    load mR R1 R2 R3 R4;
else,
    disp('loading position from isdat');
    DATABASE='disco:10';db = Mat_DbOpen(DATABASE);
    for ic=sc_list, disp(['...R' num2str(ic)]);
     [t,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     eval(av_ssub('R?=[double(t) double(data)''];',ic));clear t data;
    end
end

if exist('mV.mat'),
      load mV V1 V2 V3 V4;
  else,
      disp('loading position from isdat');
      DATABASE='disco:10';db = Mat_DbOpen(DATABASE);
      for ic=sc_list, disp(['...V' num2str(ic)]);
       [t,data] = isGetDataLite( db,min(t)-1, max(t)-min(t)+2,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
       eval(av_ssub('V?=[double(t) double(data)''];',ic));clear t data;
      end
end

if strcmp(flag,'v_from_t'),
  for ic=1:4,eval(av_ssub('vsc?=av_interp(V?,t,''spline'');',ic));end
  for ic=1:4,eval(av_ssub('drsc?=av_interp(av_add(1,R?,-1,R1),t(?),''spline'');',ic));end
  for ic=1:4,eval(av_ssub('dr?=drsc?+[0 (t(?)-t(1))*vsc?(:,2:4)];',ic));end
  for ic=1:4,eval(av_ssub('dt(?)=t(?)-t(1);sdt?=num2str(dt(?),3);',ic));end
  D=[dr2(2:4);dr3(2:4);dr4(2:4)];
  T=[dt(2),dt(3), dt(4)]';
  m=D\T;
  clear v
  v=m/norm(m)/norm(m);v=v';	% velocity vector of the boundary

  disp([ datestr(datenum(fromepoch(t(1))))])
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  vn=av_norm(v);
  strv=['V=' num2str(av_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s GSE'];
  disp(strdt);disp(strv);
elseif strcmp(flag,'dt_from_v'),
  for ic=1:4,eval(av_ssub('vsc?=av_interp(V?,t,''spline'');v?=av_add(1,v,-1,vsc?);',ic));end
  for ic=1:4,eval(av_ssub('dr?=av_interp(av_add(1,R?,-1,R1),t,''spline'');',ic));end
  for ic=1:4,eval(av_ssub('dt(?)=av_dot(v?,dr?,1)./av_abs(v?,1)^2;',ic));end

  % print result
  disp([ datestr(datenum(fromepoch(t(1))))])
  vn=av_norm(v);
  strv=['V=' num2str(av_abs(v,1),3) ' [ ' num2str(vn(end-2:end),' %5.2f') '] km/s GSE'];
  strdt=['dt=[' , num2str(dt,' %5.2f') '] s. dt=[t1-t1 t2-t1 ...]'];
  disp(strv);  disp(strdt);
  v=dt; % output variable is v
end

