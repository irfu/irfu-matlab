% get cluster sc data into matlab files,
% despinn if necessary, uses the neareset calibration values
% export to ASCII files
% uses ISGETDATALITE

if exist('sc_list') == 0, sc_list=1:4;end % default values
%eval('help c_get ');
mmm =  ...
 ['q  end, quit                          ';
	'0  this menu                          ';
	'1  time interval                      ';
	'2  sc list [1:4]                      ';
	'a  load phase A1..A4 into mA          ';
	'e  load wE1...wE4 -> mE               ';
	'eph load ephemeris,r,v,dr,dv          ';
	'b  load BPP1...BPP4 -> mBPP           ';
	'bf load Hres FGM B1...B4 -> mB        ';
%	'bfgm from DDS  FGM B1...B4 -> mB      ';
	'bs load BS1...BS4 into mBS            ';
	'p  potential P1..P4 -> mP             ';
%	's  with/ without saving               ';
	'vc CIS VCp?,VCh?,dVCp?,dVCh?-> mCIS   ';
	'vce CIS E VCE1..&dVCE1..-> mCIS       ';
	'x  free format load                   ';
	'---------------------------           ';
	'db despinned dBPP1... -> mBPP         ';
	'dbf despinned dB1... -> mB            ';
	'dbs despinned dBS1..dBS4 -> mBS       ';
	'dve despinned dvE1..dvE4 -> mE        ';
	'de (dve+vsxBPP) dE1..dE4 -> mE        ';
	'deo dEo1..dEo4 Eo.B=0                 ';
	'es  E ->ascii E?                      ';
	'ps  P ->ascii P?                      ';
];
disp('-------------- Get Cluster II data ----------------');
disp(mmm); % show onces menu
q='0';flag_save=1;
while(q ~= 'q') % ====== MAIN LOOP =========
 q=input('input>','s');if isempty(q),q='0';end
 save_list='';
 if q == 'q', return,
 elseif q == '0', disp(mmm);
 elseif q == 's',
    if flag_save==1,flag_save=0;disp('not saving variables');
    else, flag_save=1;disp('saving variables to mfiles');
    end
 elseif q == '1',
  variable='DATABASE';default='disco:10';question='Give database as string [%]>';av_ask
  db = Mat_DbOpen(DATABASE);
  variable='start_time_s';default='1999 01 01 00 00 00';question='Start time [%]>';av_ask
  start_time=eval(['[' start_time_s ']']);
  variable='Dt';default=60;question='How many seconds of data [%]>';av_ask
  tint_epoch=toepoch(start_time)+[0 Dt];
 elseif strcmp(q,'2'),
  % define sc_list
  variable='sc_list';default=[1:4];question='Spacecraft list [%]>'; av_ask;

 elseif strcmp(q,'a'),
    for ic=sc_list, disp(['...A' num2str(ic)]);
     [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');
     eval(av_ssub('A?=[double(t) double(data)];',ic));%clear t data;
     if flag_save==1, eval(av_ssub('if exist(''./mA.mat''),save mA A? -append; else, save mA A?;end',ic));end
    end

 elseif strcmp(q,'eph'),
    for ic=sc_list, disp(['...ephemeris' num2str(ic) '...LT,MLT,ILAT,L->mEPH...R->mR...V->mV']);
     [tlt,lt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'lt', ' ', ' ', ' ');
     [tmlt,mlt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'mlt', ' ', ' ', ' ');
     [tL,Lshell] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'l_shell', ' ', ' ', ' ');
     [tilat,ilat] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'inv_lat', ' ', ' ', ' ');
     [tlat, lat] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_lat__CL_SP_AUX'], ' ', ' ',' ');
     [tlong, long] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_long__CL_SP_AUX'], ' ', ' ',' ');
     [tr,r] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     [tv,v] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
     eval(av_ssub('',ic));clear t data;
     eval(av_ssub('LT?=[double(tlt) double(lt)];MLT?=[double(tmlt) double(mlt)];L?=[double(tL) double(Lshell)];ILAT?=[double(tilat) double(ilat)];R?=[double(tr) double(r)''];V?=[double(tv) double(v)''];spinaxis_latlong?=[double(tlat) double(lat) double(long)];',ic));clear tlt tmlt tL tilat lt mlt Lshell ilat tr r tv v tlat lat long;
     eval(av_ssub('if exist(''./mEPH.mat''),save mEPH LT? MLT? L? ILAT? spinaxis_latlong? -append; else, save mEPH LT? MLT? L? ILAT? spinaxis_latlong?;end',ic));
     eval(av_ssub('tt=R?(1,1);dR?=c_gse2dsc(R?,[tt ic]);',ic));  % despinned coordinates
     eval(av_ssub('if exist(''./mR.mat''),save mR R? dR? -append; else, save mR R? dR? ;end',ic));
     eval(av_ssub('tt=V?(1,1);dV?=c_gse2dsc(V?,[tt ic]);',ic));  % despinned coordinates
     eval(av_ssub('if exist(''mV.mat''),save mV V? dV? -append; else, save mV V? dV? ;end',ic));
    end

 elseif strcmp(q,'x'),
    var_name=input('matlab variable name =','s');
    disp('? in input is substituted by cluster number');
    qstring=input('example {''Cluster'' ''1'' ''efw'' ''E'' ''p12'' ''10Hz'' ''any''}=>','s');
    for ic=sc_list,
     string=eval(av_ssub(qstring,ic));
     varic = [var_name num2str(ic)];
     disp(['...free format s/c' num2str(ic) ' ' varic]);
     for jj=1:size(string,2),str{jj}=av_ssub(string{jj},ic);end, for jj=size(string,2)+1:7,str{jj}=' ';end
     [t,data] = isGetDataLite( db, start_time, Dt,str{1}, str{2}, str{3}, str{4},str{5},str{6},str{7});
     eval([varic '=[double(t) double(data)];']);
    end

 elseif q == 'b',
    yyyymmdd=[datestr(datenum(start_time),10) datestr(datenum(start_time),5) datestr(datenum(start_time),7) ];
    for ic=sc_list, disp(['CSDS...BPP' num2str(ic)]);
      %eval(av_ssub(['data=cdfread(''/data/cluster/CSDS/PP/FGM/C?/C?_PP_FGM_' yyyymmdd '_V01.CDF'',{''Epoch__C?_PP_FGM'' ''B_xyz_gse__C?_PP_FGM''});'],ic));
      %eval(av_ssub(['t=cdfread(''/data/cluster/CSDS/FGM/C?/C?_PP_FGM_' yyyymmdd '_V01.CDF'',{''Epoch__C?_PP_FGM''};'],ic));
      %tei=toepoch(datevec(cdf2date(double(data{1}))));ind=find(tei>0);t=tei(ind);
      %eval(av_ssub('b=double(data{2}); BPP?=[t b(ind,:)];',ic));clear t tei data b;
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'FGM', ['B_xyz_gse__C' num2str(ic) '_PP_FGM'], ' ', ' ',' ');
      eval(av_ssub('BPP?=[double(t) double(data)''];',ic));clear t,data;
     eval(av_ssub('if exist(''./mBPP.mat''),save mBPP BPP? -append; else, save mBPP BPP?;end',ic));
    end

 elseif strcmp(q,'db'),
  for ic=sc_list,
   eval(av_ssub('load mBPP BPP?;dBPP?=c_gse2dsc(BPP?,[BPP?(1,1) ic]);save -append mBPP dBPP?;',ic));
  end

 elseif strcmp(q,'bf'),
  for ic=sc_list,
    disp(['Choose FGM GSE data for spacecraft ' num2str(ic) ]);
    fvs=fgmvec_stream('/home/andris/data/cluster/fgm/');
    disp(tavail(fvs,[]));
    fgm_t_interval=[['T' datestr(datenum(start_time),13) 'Z'] ['T' datestr(datenum(start_time)+Dt/86400,13) 'Z']];
    disp(['fgm_t_interval=' fgm_t_interval]);
    dat=get(fvs,'data','b',fgm_t_interval);
    eval(av_ssub('B?=[rem(dat.time,1)*3600*24+toepoch(start_time.*[1 1 1 0 0 0]) dat.b];save_list=[save_list '' B? ''];',ic));
  end
  eval(['save mB ' save_list]);

 elseif strcmp(q,'bfgm'),
  disp('CONTACT STEPHAN BUCHERT!!!!!!!!!!!!!!!!!!!');
  for ic=sc_list,
    eval(av_ssub('B?=c_get_bfgm(tint_epoch,ic);',ic));
  end
  eval(['save mB ' save_list]);

 elseif strcmp(q,'dbf'),
  for ic=sc_list,
   eval(av_ssub('load mB B?;tt=B?(1,1);',ic));
   eval(av_ssub('dB?=c_gse2dsc(B?,[B?(1,1) ic]);',ic));
   eval(av_ssub('save -append mB dB?;',ic));
  end

 elseif strcmp(q,'bs'),
  mode=input('Model L=1/M=2? If different give as vector. [1]');if isempty(mode),mode=1;end;
  for ic=sc_list,
   	if (length(mode)>1), mm=mode(ic);else, mm=mode;end
		if (mm == 1), param='0-10Hz';end;
		if (mm == 2), param='0-180Hz';end;
    disp(['STAFF...wBS' num2str(ic) ' ' param ' filter' ]);
    [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'staff', 'B_SC', 'Bx_By_Bz', param, '');
    eval(av_ssub('wBS?=[double(t) double(data)''];',ic));clear t data;
    end
    save mBS wBS1 wBS2 wBS3 wBS4;

 elseif strcmp(q,'dbs'),
  for ic=sc_list,
   eval(av_ssub('load mBS wBS?;tt=wBS?(1,1);',ic));
   eval(av_ssub('load mA.mat A?;',ic));
   eval(av_ssub('dBS?=c_despin(wBS?,A?);save -append mBS dBS?;',ic));
  end

 elseif strcmp(q,'bfgm'),
    disp('High time resolution FGM from isdat');
    save_list='';
    for ic=sc_list,
     disp(['...Bprim' num2str(ic) ', Bsec' num2str(ic)]);
     [t, data] = isGetDataLite( db, start_time, Dt, 'Cluster', num2str(ic), 'fgm', 'Bprim' , ' ', ' ',' ');
     eval(av_ssub('Bprim?=[double(t) double(real(data))''];',ic));clear t,data;
     eval(av_ssub('if size(Bprim?), save_list=[save_list '' Bprim? '']; end; ',ic));
     [t, data] = isGetDataLite( db, start_time, Dt, 'Cluster', num2str(ic), 'fgm', 'Bsec' , ' ', ' ',' ');
     eval(av_ssub('Bsec?=[double(t) double(real(data))''];',ic));clear t,data;
     eval(av_ssub('if size(Bsec?), save_list=[save_list '' Bsec? '']; end; ',ic));
    end
    if exist('./mB.mat'), eval(['save -append mB ' save_list]); else, eval(['save mB ' save_list]);end

 elseif q == 'e',
  mode=av_q('Filter 1)10Hz_lx 2)180Hz 3)10Hz_hx? If different give as vector. [%]','mode',1);
  for ic=sc_list,
   	if (length(mode)>1), mm=mode(ic);else, mm=mode;end
		if (mm == 1), param='10Hz';tmmode='lx';end;
		if (mm == 2), param='180Hz';tmmode='hx';end;
		if (mm == 3), param='10Hz';tmmode='hx';end;
    disp(['EFW...sc' num2str(ic) '...E ' param ' filter']);
    if toepoch(start_time)>toepoch([2001 12 28 03 00 00])&ic==1, % probe 1 probelm on sc1
       [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
       data=double(real(data));
       data=[data(:,1)*0 data(:,1) data(:,1)*0]';
       disp('            !Only p34 exist for sc1');
    elseif toepoch(start_time)>toepoch([2002 07 29 09 06 59 ])&ic==3, % probe 1 probelm on sc3
       [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
       data=double(real(data));
       data=[data(:,1)*0 data(:,1) data(:,1)*0]';
       disp('            !Only p34 exists for sc3');
    else,
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p1234', param, tmmode);
    end
    t_e=double(t);
    eval(av_ssub('wE?=[t_e double(real(data))''];',ic));clear t t_e data;
    eval(av_ssub('save_list=[save_list '' wE? ''];',ic));
  end
  if exist('./mE.mat'), eval(['save -append mE  ' save_list]); else, eval(['save mE  ' save_list]);end

 elseif strcmp(q,'dve'),
  q_efw_offset=av_q('How to treat offsets? 1)do nothing, 2) nearest callibration, 3) subtract probe signal mean and use nearest sunward offsets. [%]>','q_efw_offset',2);
  for ic=sc_list,
   eval(av_ssub('load mE wE?;tt=wE?(1,1);',ic));
   eval(av_ssub('load mA.mat A?;',ic));
   [coef1,coef2,coef3,coef4]=c_efw_calib(tt); % the time stamp of the first sample defines calibration
   if q_efw_offset==1,
     eval(av_ssub('dvE?=c_despin(wE?,A?);',ic));
   elseif q_efw_offset==2,
     eval(av_ssub('dvE?=c_despin(wE?,A?,coef?);',ic));
   elseif q_efw_offset==3,
     eval(av_ssub('dvE?=c_despin(wE?,A?,coef?,''efw_offs'');',ic));
   else, disp('wrong offest option, using option 2.');
     eval(av_ssub('dvE?=c_despin(wE?,A?,coef?);',ic));
   end
   eval(av_ssub('save_list=[save_list '' dvE? ''];',ic));
  end
  eval(['save -append mE  ' save_list]);

 elseif strcmp(q,'de'),
  for ic=sc_list,
   eval(av_ssub('load mE dvE?;load mBPP dBPP?;load mV dV?; tt=dvE?(1,1);db=dBPP?; dv=dV?;clear dBPP? dV?;',ic));
   evxb=av_interp(av_t_appl(av_cross(db,dv),'*1e-3*(-1)'),tt);
   eval(av_ssub('dE?=av_add(1,dvE?,-1,evxb);dE?(:,4)=0;',ic));
   eval(av_ssub('save_list=[save_list '' dE? ''];',ic));
  end
  eval(['save -append mE  ' save_list]);

 elseif strcmp(q,'deo'),
  disp('Estimating dEo where dEo.B=0, Eo (GSE) and Vo=Eo/B')
  deg=input('B angle with respect to the spin plane should be at least x deg, x=');
  qb=input('To use FGM high res (1) or PP data (2) [1] >');if isempty(qb),qb=1;end;
  for ic=sc_list,
     eval(av_ssub('load mE dE?;tt=dE?(1,1);',ic));
   if qb ==1, eval(av_ssub('load mB dB?; db=av_interp(dB?,dE?);clear dB?;',ic));
   else, eval(av_ssub('load mBPP dBPP?; db=av_interp(dBPP?,dE?);clear dBPP?;',ic));
   end
   eval(av_ssub('[dEo?,d?]=av_ed(dE?,db,deg);Eo?=c_gse2dsc(dEo?,[tt ?],-1);indzero=find(abs(d?)<10);Eo?(indzero,4)=0;',ic));
   eval(av_ssub('dVo?=av_e_vxb(dEo?,db,-1);Vo?=c_gse2dsc(dVo?,[tt ?],-1);',ic));
   eval(av_ssub('save_list=[save_list '' dEo? d? Eo? Vo? ''];',ic));
  end
  eval(['save -append mE  ' save_list]);

 elseif strcmp(q,'es'), % create E ascii files
%  variable='qf';default=1;question='Create E1.dat ... E4.dat ascii files with 1) E_GSE B_angle 2) E_DS [%]>';av_ask
  for ic=sc_list,
     % E_GSE file creation
     eval(av_ssub('load mE Eo? d?;tt=Eo?(1,1);x=Eo?;x(:,end+1)=d?;',ic));
     t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
     file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_E_GSE_sc' num2str(ic) '.dat'];
     disp(['E' num2str(ic) ' --> ' file_name '  ' num2str(size(x,1)) ' samples']);
     fid = fopen(file_name,'w');
     fprintf(fid,'%% E-field in GSE reference frame, assuming E_GSE.B=0 \n');
     fprintf(fid,'%% Ez_GSE is not reliable when the ambient magnetic field is close to the spin plane\n');
     fprintf(fid,'%% The last column shows the angle of the ambient field with respect to the spin plane\n');
     fprintf(fid,'%% Ez_GSE is set to 0 when this angle is less than 10 degrees \n');
     fprintf(fid,'%% Cluster %1d\n',ic);
     fprintf(fid,['%% Time is in seconds from ' time_ref '\n'],ic);
     fprintf(fid,['%%  time      Ex_GSE   Ey_GSE   Ez_GSE angle B/spin_plane\n']);
     fprintf(fid,['%%  (s)       (mV/m)   (mV/m)   (mV/m)  (degrees)\n']);
     x(:,1)=x(:,1)-t_ref;x=x';
     fprintf(fid,'%10.4f %8.3f %8.3f %8.3f %5.1f\n',x);
     fclose(fid);
     % E_DS file creation
     eval(av_ssub('load mE dE?;tt=dE?(1,1);x=dE?(:,[1 2 3]);',ic));
     t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
     file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_E_DS_sc' num2str(ic) '.dat'];
     disp(['dE' num2str(ic) ' --> ' file_name '  ' num2str(size(x,1)) ' samples']);
     fid = fopen(file_name,'w');
     fprintf(fid,'%% E-field in the despinned spacecraft reference frame DS\n');
     fprintf(fid,'%% Approximately Ex_DS=Ex_GSE, Ey_DS=-Ey_GSE\n');
     fprintf(fid,'%% Cluster %1d\n',ic);
     fprintf(fid,['%% Time is in seconds from ' time_ref '\n'],ic);
     fprintf(fid,['%%  time      Ex_DS    Ey_DS \n']);
     fprintf(fid,['%%  (s)       (mV/m)   (mV/m)\n']);
     x(:,1)=x(:,1)-t_ref;x=x';
     fprintf(fid,'%10.4f %8.3f %8.3f\n',x);
     fclose(fid);
   end

 elseif q == 'p',
    save_list = '';
    variable='mode';default=1;question='Model 1)10Hz_lx, 2)180Hz_hx, 4)32kHz_any, 11)10Hz_hx? If different give as vector. [%]';av_ask;
  for ic=sc_list,
  	if (length(mode)>1), mm=mode(ic);else, mm=mode;end
		if (mm == 1), param='10Hz';tmmode='lx';end;
		if (mm == 11), param='10Hz';tmmode='hx';end;
		if (mm == 2), param='180Hz';tmmode='hx';end;
		if (mm == 4), param='32kHz';tmmode='any';end;
    for probe=1:4;
      disp(['EFW...sc' num2str(ic) '...probe' num2str(probe) '->P' param num2str(ic) 'p' num2str(probe)]);
      [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', ['p' num2str(probe)],param, tmmode);
      eval(av_ssub(['p!=[double(t) double(real(data))];save_list=[save_list '' P' param '?p!''];P' param '?p!=p!;'],ic,probe)); clear t data;
    end
    clear p;
    if size(p1)==size(p2)&size(p1)==size(p3)&size(p1)==size(p4)&size(p1)~=[0 0]&ic~=2,  % sc2 has often problems with p3
       p=[p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
    elseif size(p1)==size(p2)&size(p1)~=[0 0],
       p=[p1(:,1) (p1(:,2)+p2(:,2))/2];
    elseif size(p3)==size(p4)&ic~=2,
       p=[p3(:,1) (p3(:,2)+p4(:,2))/2];
    else,
         p=p4;
    end
    eval(av_ssub(['P' param '?=p;save_list=[save_list '' P' param '? ''];'],ic));
    if ((mm==1) | (mm==11)); eval(av_ssub('P?=p;NVps?=c_n_Vps(p);save_list=[save_list '' P? NVps?''];',ic));end
    if (mm == 4),
      dtburst=input(['s/c' num2str(ic) ', time shift to obtain correct time (get from Anders Tjulin) =']);
      dtburst=double(dtburst);
      for probe=1:4,eval(av_ssub('P32kHz?p!(:,1)=P32kHz?p!(:,1)+dtburst;',ic,probe));end
      eval(av_ssub('P32kHz?(:,1)=P32kHz?(:,1)+dtburst;',ic));
    end
  end
  if exist('./mP.mat'), eval(['save mP ' save_list ' -append']); else eval(['save mP ' save_list]); end

 elseif strcmp(q,'ps'), % create V_sc n ascii files
  for ic=sc_list,
     % Vsc_N file creation
     eval(av_ssub('tt=P?(1,1);x=P?;x(:,end+1)=c_n_Vps(P?(:,end));',ic));
     t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
     file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_Vps_N_sc' num2str(ic) '.dat'];
     disp(['P' num2str(ic) ' --> ' file_name '  ' num2str(size(x,1)) ' samples']);
     fid = fopen(file_name,'w');
     fprintf(fid,'%% Vps - probe to spacecraft potential which is approximately \n');
     fprintf(fid,'%%       the same as satellite potential with respect to plasma.\n');
     fprintf(fid,'%% N   - density derived from the satellite potential based on \n');
     fprintf(fid,'%%       empirical fit to Cluster data.\n');
     fprintf(fid,'%%       It is NOT true density\n');
     fprintf(fid,'%% Cluster %1d\n',ic);
     fprintf(fid,['%% Time is in seconds from ' time_ref '\n'],ic);
     fprintf(fid,['%%  time        Vps        N\n']);
     fprintf(fid,['%%  (s)         (V)      (cm-3)\n']);
     x(:,1)=x(:,1)-t_ref;x=x';
     fprintf(fid,'%10.4f %8.2f %8.2f\n',x);
     fclose(fid);
   end

 elseif q == 'r',
    for ic=sc_list, disp(['...R' num2str(ic)]);
     [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     eval(av_ssub('R?=[double(t) double(data)''];',ic));clear t data;
     eval(av_ssub('if exist(''./mR.mat''),save mR R? -append; else, save mR R?;end',ic));
    end

 elseif strcmp(q,'v'),
    for ic=sc_list, disp(['...V' num2str(ic)]);
     [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
     eval(av_ssub('V?=[double(t) double(data)''];',ic));clear t data;
     eval(av_ssub('if exist(''./mV.mat''),save mV V? -append; else, save mV V?;end',ic));
    end

 elseif strcmp(q,'vc'),
    save_list='';
    for ic=sc_list,
     disp(['...VCp' num2str(ic) ', dVCp' num2str(ic)]);
     [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_p_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' ');
     eval(av_ssub('VCp?=[double(t) double(real(data))''];',ic));clear t,data;
     eval(av_ssub('if size(VCp?), dVCp?=c_gse2dsc(VCp?,[VCp?(1,1) ic]); save_list=[save_list '' VCp? dVCp? '']; end; ',ic));
     disp(['...VCh' num2str(ic) ', dVCh' num2str(ic)]);
     [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_HIA_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' ');
     eval(av_ssub('VCh?=[double(t) double(real(data))''];',ic));clear t,data;
     eval(av_ssub('if size(VCh?), dVCh?=c_gse2dsc(VCh?,[VCh?(1,1) ic]); save_list=[save_list '' VCh? dVCh? '']; end;',ic));
    end
    eval(['save mCIS ' save_list]);

 elseif strcmp(q,'vce'),
  CIS=load('mCIS');
  for ic=sc_list,
   eval(av_ssub('if isfield(CIS,''VCp?''); vp=CIS.VCp?;else,vp=[];end;  if isfield(CIS,''VCh?'');vh=CIS.VCh?;else,vh=[];end; load mBPP BPP?;b=BPP?; clear BPP? VCp? VCh?;',ic));
   %
    if min(size(vp)) ~= 0,
      disp(['...VCEp' num2str(ic)]);
      evxb=av_t_appl(av_cross(vp,b),'*(-1e-3)');
      eval(av_ssub('VCEp?=evxb;save -append mCIS VCEp?;',ic));
      disp(['...dVCEp' num2str(ic)]);
      eval(av_ssub('dVCEp?=c_gse2dsc(VCEp?,[VCEp?(1,1) ic]);save -append mCIS dVCEp?;',ic));
    end
    if min(size(vh)) ~= 0,
      disp(['...VCEh' num2str(ic)]);
      evxb=av_t_appl(av_cross(vh,b),'*(-1e-3)');
      eval(av_ssub('VCEh?=evxb;save -append mCIS VCEh?;',ic));
      disp(['...dVCEh' num2str(ic)]);
      eval(av_ssub('dVCEh?=c_gse2dsc(VCEh?,[VCEh?(1,1) ic]);save -append mCIS dVCEh?;',ic));
    end
  end

 else
  eval(q,'');
 end

end
if exist('db'), Mat_DbClose(db); end


