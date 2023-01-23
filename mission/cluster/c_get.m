% get cluster sc data into matlab files,
% despinn if necessary, uses the neareset calibration values
% export to ASCII files
% uses isGetDataLite
%

DP_S = c_ctl(0,'data_path');
csds_dir = [DP_S '/CSDS/'];

if ~exist('sc_list','var'), sc_list=1:4;end % default values

mmm =  ...
  ['q  end, quit                          ';
  %	'0  this menu                          ';
  '1  time interval                      ';
  '2  sc list [1:4]                      ';
  'a  load phase A1..A4     -> mA        ';
  'e  load wE1...wE4        -> mE        ';
  'e1  load wE1p12...wE4p34 -> mE        ';
  'eph load ephemeris,r,v,dr,dv          ';
  'b  load BPP1...BPP4 -> mBPP           ';
  'bf load Hres FGM B1...B4 -> mB        ';
  %	'bfgm from DDS  FGM B1...B4 -> mB      ';
  'bs load BS1...BS4 ->   mBS            ';
  'p  potential P1..P4 -> mP             ';
  %	's  with/ without saving               ';
  'vc CIS VCp?,VCh?,dVCp?,dVCh?-> mCIS   ';
  'vce CIS E VCE1...dVCE1..-> mCIS       ';
  'edi EDI E EDI1...dEDI1..-> mEDI       ';
  'wbdwf WBD E/B wfWBD..-> mWBD          ';
  'x  free format load                   ';
  '---------------------------           ';
  'db despinned dBPP1... -> mBPP         ';
  'dbf despinned dB1... -> mB            ';
  'dibf despinned diB1... -> mB          ';
  'dbs despinned dBS1..dBS4 -> mBS       ';
  'dve despinned dvE1..dvE4 -> mE        ';
  'dve1 spin fits dvE1p12..dvE4p34 -> mE ';
  'de (dve+vsxBPP) dE1..dE4 -> mE        ';
  'deo dEo1..dEo4 Eo.B=0                 ';
  'ea,esa (s-fit)  E ->ascii E?          ';
  'va,vsa (s-fit)  ExB ->ascii VExB?     ';
  'pa  P ->ascii P?                      ';
  ];
disp('-------------- Get Cluster II data ----------------');
disp(mmm); % show onces menu
q='0';flag_save=1;
while(q ~= 'q') % ====== MAIN LOOP =========
  q=input('input>','s');if isempty(q),q='0';end
  save_list='';save_file='';
  if q == 'q', return,
  elseif q == '0', disp(mmm);
  elseif q == 's'
    if flag_save==1,flag_save=0;disp('not saving variables');
    else, flag_save=1;disp('saving variables to mfiles');
    end
  elseif q == '1'
    DB_S = c_ctl(0,'isdat_db');
    DATABASE = irf_ask('Give database as string [%]>','DATABASE',DB_S);
    db = Mat_DbOpen(DATABASE);
    start_time_s = irf_ask('Start time [%]>','start_time_s','2001 02 01 00 00 00');
    start_time=eval(['[' start_time_s ']']);
    start_date_str=strrep(datestr(start_time,29),'-','');
    Dt = irf_ask('How many seconds of data [%]>','Dt',60);
    tint_epoch=toepoch(start_time)+[0 Dt];
  elseif strcmp(q,'2')
    % define sc_list
    sc_list = irf_ask('Spacecraft list [%]>','sc_list',1:4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'a')
    save_file='./mA.mat';save_list='';
    for ic=sc_list, c_eval('disp(''...A?...Atwo?..'');',ic);
      [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');
      eval(irf_ssub('A?=[double(t) double(data)];',ic));%clear t data;
      [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase_2', ' ', ' ', ' ');
      eval(irf_ssub('Atwo?=[double(t) double(data)];',ic));%clear t data;
      c_eval('save_list=[save_list '' A? Atwo? ''];',ic);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ephemeris
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'eph')
    for ic=sc_list, disp(['...ephemeris' num2str(ic) '...LT,MLT,ILAT,L->mEPH...R->mR...V->mR']);
      [tlt,lt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'lt', ' ', ' ', ' ');
      [tmlt,mlt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'mlt', ' ', ' ', ' ');
      [tL,Lshell] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'l_shell', ' ', ' ', ' ');
      [tilat,ilat] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'inv_lat', ' ', ' ', ' ');
      [tlat, lat] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_lat__CL_SP_AUX'], ' ', ' ',' ');
      [tlong, long] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_long__CL_SP_AUX'], ' ', ' ',' ');
      [tr,r] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
      [tv,v] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
      eval(irf_ssub('',ic));clear t data;
      eval(irf_ssub('LT?=[double(tlt) double(lt)];MLT?=[double(tmlt) double(mlt)];L?=[double(tL) double(Lshell)];ILAT?=[double(tilat) double(ilat)];R?=[double(tr) double(r)''];V?=[double(tv) double(v)''];',ic));clear tlt tmlt tL tilat lt mlt Lshell ilat tr r tv v;
      eval(irf_ssub('spinaxis_latlong?=[double(tlat) double(lat) double(long)];',ic)); clear tlat lat long;
      eval(irf_ssub('if exist(''./mEPH.mat''),save mEPH LT? MLT? L? ILAT? spinaxis_latlong? -append; else, save mEPH LT? MLT? L? ILAT? spinaxis_latlong?;end',ic));
      eval(irf_ssub('tt=R?(1,1);dR?=c_gse2dsc(R?,[tt ic]);',ic));  % despinned coordinates
      eval(irf_ssub('if exist(''./mR.mat''),save mR R? dR? -append; else, save mR R? dR? ;end',ic));
      eval(irf_ssub('tt=V?(1,1);dV?=c_gse2dsc(V?,[tt ic]);',ic));  % despinned coordinates
      eval(irf_ssub('if exist(''mR.mat''),save mR V? dV? -append; else, save mR V? dV? ;end',ic));
    end
    save_list = '';
    
  elseif strcmp(q,'x')
    var_name=input('matlab variable name =','s');
    disp('? in input is substituted by cluster number');
    qstring=input('example {''Cluster'' ''1'' ''efw'' ''E'' ''p12'' ''10Hz'' ''any''}=>','s');
    for ic=sc_list
      string=eval(irf_ssub(qstring,ic));
      varic = [var_name num2str(ic)];
      disp(['...free format s/c' num2str(ic) ' ' varic]);
      for jj=1:size(string,2),str{jj}=irf_ssub(string{jj},ic);end
      for jj=size(string,2)+1:7,str{jj}=' ';end
      [t,data] = isGetDataLite( db, start_time, Dt,str{1}, str{2}, str{3}, str{4},str{5},str{6},str{7});
      eval([varic '=[double(t) double(data)];']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Magnetic fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif q == 'b'
    save_file='./mBPP.mat';
    for ic=sc_list, disp(['CSDS...BPP' num2str(ic)]);
      eval(irf_ssub('BPP?=irf_cdf_read([csds_dir ''PP/FGM/C?/C?_PP_FGM_'' start_date_str ''*''],''B_xyz_gse__C?_PP_FGM'');BPP?=irf_tlim(BPP?,tint_epoch);save_list=[save_list '' BPP?''];',ic));
      %      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'FGM', ['B_xyz_gse__C' num2str(ic) '_PP_FGM'], ' ', ' ',' ');
      %      eval(irf_ssub('BPP?=[double(t) double(data)''];',ic));clear t,data;
      %     eval(irf_ssub('if exist(''./mBPP.mat''),save mBPP BPP? -append; else, save mBPP BPP?;end',ic));
    end
    
  elseif strcmp(q,'db')
    for ic=sc_list
      eval(irf_ssub('load mBPP BPP?;dBPP?=c_gse2dsc(BPP?,[BPP?(1,1) ic]);save -append mBPP dBPP?;',ic));
    end
    
  elseif strcmp(q,'bf')
    for ic=sc_list
      disp(['Choose FGM GSE data for spacecraft ' num2str(ic) ]);
      fvs=fgmvec_stream('/home/andris/data/cluster/fgm/');
      disp(tavail(fvs,[]));
      fgm_t_interval=[['T' datestr(datenum(start_time),13) 'Z'] ['T' datestr(datenum(start_time)+Dt/86400,13) 'Z']];
      disp(['fgm_t_interval=' fgm_t_interval]);
      dat=get(fvs,'data','b',fgm_t_interval);
      eval(irf_ssub('B?=[rem(dat.time,1)*3600*24+toepoch(start_time.*[1 1 1 0 0 0]) dat.b];save_list=[save_list '' B? ''];',ic));
    end
    eval(['save mB ' save_list]);
    save_list = '';
    
  elseif strcmp(q,'bfgm')
    disp('CONTACT STEPHAN BUCHERT!!!!!!!!!!!!!!!!!!!');
    save_file='./mB.mat';
    save_list = '';
    c_eval('B?=c_get_bfgm(tint_epoch,?);save_list=[save_list '' B? ''];',sc_list);
    c_eval('diB?=c_gse2dsc(B?,[B?(1,1) ?],2);save_list=[save_list '' diB? ''];',sc_list);
    
  elseif strcmp(q,'dbf')
    for ic=sc_list
      eval(irf_ssub('load mB B?;tt=B?(1,1);',ic));
      eval(irf_ssub('dB?=c_gse2dsc(B?,[B?(1,1) ic]);',ic));
      eval(irf_ssub('save -append mB dB?;',ic));
    end
    
  elseif strcmp(q,'dibf')
    for ic=sc_list
      eval(irf_ssub('load mB B?;',ic));
      eval(irf_ssub('diB?=c_gse2dsc(B?,?,2);',ic));
      eval(irf_ssub('save -append mB diB?;',ic));
    end
    
  elseif strcmp(q,'bs')
    mode=input('Model L=1/M=2? If different give as vector. [1]');if isempty(mode),mode=1;end
    for ic=sc_list
      if (length(mode)>1), mm=mode(ic);else, mm=mode;end
      if (mm == 1), param='0-10Hz';end
      if (mm == 2), param='0-180Hz';end
      disp(['STAFF...wBS' num2str(ic) ' ' param ' filter' ]);
      [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'staff', 'B_SC', 'Bx_By_Bz', param, '');
      eval(irf_ssub('wBS?=[double(t) double(data)''];',ic));clear t data;
    end
    save mBS wBS1 wBS2 wBS3 wBS4;
    
  elseif strcmp(q,'dbs')
    c_eval('load mBS.mat wBS?; c_load A?; if ~isempty(wBS?), dBS?=c_efw_despin(wBS?,A?); save -append mBS dBS?;end; clear wBS? A?',sc_list);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E p1234 or p12 & p34
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'e') || strcmp(q,'e1')
    save_file = './mE.mat';
    mode=irf_ask('Sampling 1)hx 2)lx ? If different give as vector. [%]','mode',1);
    for ic=sc_list
      if (length(mode)>1), mm=mode(ic);else, mm=mode;end
      if (mm == 2), param='10Hz';tmmode='lx';
      elseif (mm == 1)
        %% Find TapeMode
        % We read FDM from isdat and 5-th column contains the HX mode (undocumented)
        % 0 - normal mode  (V12L,V34L)
        % 1 - tape mode 1  (V12M,V34M)
        % 2 - tape mode 2  (V12M,V34M,)
        % 3 - tape mode 3  (V1M,V2M,V3M,V4M)
        %
        clear tm mTMode1 mTMode2 mTMode3 mTMode4
        [ok,tm] = c_load('mTMode?',ic);
        if ~ok
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic),'efw','FDM');
          if ~isempty(data), tm=double(data(5,:)); else, error('Cannot fetch FDM'), end
          eval(irf_ssub('mTMode?=tm;',ic));
          if exist('./mTMode.mat','file'), eval(irf_ssub('save -append mEFWR mTMode?;',ic));
          else, eval(irf_ssub('save mEFWR mTMode?;',ic));
          end
        end
        if tm~=tm(1)*ones(size(tm))
          warning('tape mode changes during the selected tile inteval')
        end
        tm=tm(1);
        tmmode='hx';
        if tm<1e-30, param='10Hz';	else, param='180Hz'; end
        clear tm
        tst = toepoch(start_time);
        if tst>toepoch([2001 07 31 00 00 00])&&tst<toepoch([2001 09 01 00 00 00])
          % all sc run on 180Hz filter in august 2001
          param='180Hz';
        elseif tst>toepoch([2001 07 31 00 00 00])&&ic==2 % 10Hz filtef probelm on sc2
          param='180Hz';
        end
      end
      
      if (toepoch(start_time)>toepoch([2001 12 28 03 00 00])&&ic==1) || (toepoch(start_time)>toepoch([2002 07 29 09 06 59 ])&&ic==3)
        p34_only = 1;
        disp(sprintf('            !Only p34 exists for sc%d',ic));
      else
        p34_only = 0;
      end
      
      clear sensor;
      if strcmp(q,'e')
        disp(['EFW...sc' num2str(ic) '...E ' param ' filter']);
        if p34_only
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
          data = double(real(data));
          data = [data(:,1)*0 data(:,1) data(:,1)*0]';
        else
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p1234', param, tmmode);
          data = double(real(data));
        end
        t = double(t);
        eval(irf_ssub('wE?=[t data''];',ic)); clear t data;
        eval(irf_ssub('save_list=[save_list '' wE? ''];',ic));
      else
        % do separate probes
        if ~p34_only
          disp(['EFW...sc' num2str(ic) '...E p12 ' param ' filter']);
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p12', param, tmmode);
          data = double(real(data));
          t = double(t);
          eval(irf_ssub('wE?p12=[t data];',ic)); clear t data;
          eval(irf_ssub('save_list=[save_list '' wE?p12 ''];',ic));
        end
        disp(['EFW...sc' num2str(ic) '...E p34 ' param ' filter']);
        [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
        data = double(real(data));
        t = double(t);
        eval(irf_ssub('wE?p34=[t data];',ic)); clear t data;
        eval(irf_ssub('save_list=[save_list '' wE?p34 ''];',ic));
      end
    end
    %if exist('./mE.mat'), eval(['save -append mE  ' save_list]); else, eval(['save mE  ' save_list]);end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spinfits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'dve1')
    save_file = './mE.mat';
    for ic=sc_list
      eval(irf_ssub('load mE wE?p12 wE?p34;',ic));
      eval(irf_ssub('load mA.mat A?;',ic));
      
      if exist(irf_ssub('wE?p12',ic),'var')
        eval(irf_ssub('tt=wE?p12;aa=A?;',ic))
        disp(sprintf('Spin fit wE%dp12 -> dvE%dp12 mean:%.2f',ic,ic,mean(tt(:,2))))
        sp = c_efw_sfit(12,3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
        sp = sp(:,1:4);
        sp(:,3) = -sp(:,3); % DSI->DS
        sp(:,4) = 0*sp(:,4);
        o12(1) = mean(sp(:,2)); o12(2) = mean(sp(:,3));
        eval(irf_ssub('dvE?p12=sp;',ic))
        eval(irf_ssub('save_list=[save_list '' dvE?p12 ''];',ic));
      else
        disp(sprintf('No p12 data for sc%d',ic))
      end
      
      if exist(irf_ssub('wE?p34',ic),'var')
        eval(irf_ssub('tt=wE?p34;aa=A?;',ic))
        disp(sprintf('Spin fit wE%dp34 -> dvE%dp34 mean:%.2f',ic,ic,mean(tt(:,2))))
        sp = c_efw_sfit(34,3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
        sp = sp(:,1:4);
        sp(:,3) = -sp(:,3); % DSI->DS
        sp(:,4) = 0*sp(:,4);
        o34(1) = mean(sp(:,2)); o34(2) = mean(sp(:,3));
        eval(irf_ssub('dvE?p34=sp;',ic))
        eval(irf_ssub('save_list=[save_list '' dvE?p34 ''];',ic));
      else
        disp(sprintf('No p34 data for sc%d',ic))
      end
      
      %display offsets
      if exist(irf_ssub('wE?p12',ic),'var') && exist(irf_ssub('wE?p34',ic),'var')
        disp(sprintf('[X Y] offsets <p12>-<p34> : [ %.2f %.2f ]', o12(1)-o34(1), o12(2)-o34(2)))
      end
      
    end
    %eval(['save -append mE  ' save_list]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % despin E
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'dve')
    save_file = './mE.mat';
    q_efw_offset=irf_ask('How to treat offsets? \n  1) do nothing, \n  2) before despin subtract probe signal mean \n  3)  2 + use nearest sunward offsets,\n  4) nearest full hand-tuned callibration,\n[%]>','q_efw_offset',2);
    for ic=sc_list
      eval(irf_ssub('load mE wE?;tt=wE?(1,1);',ic));
      eval(irf_ssub('load mA.mat A?;',ic));
      switch q_efw_offset
        case 1,      eval(irf_ssub('dvE?=c_efw_despin(wE?,A?);',ic));
        case 2,      eval(irf_ssub('dvE?=c_efw_despin(wE?,A?,?,''efw_a'');',ic));
        case 3,      eval(irf_ssub('dvE?=c_efw_despin(wE?,A?,?,''efw_b'');',ic));
        case 4,      eval(irf_ssub('dvE?=c_efw_despin(wE?,A?,?);',ic));
        otherwise,   disp('wrong offest option, using option 2.');
          eval(irf_ssub('dvE?=c_efw_despin(wE?,A?,?);',ic));
      end
      eval(irf_ssub('save_list=[save_list '' dvE? ''];',ic));
    end
    %eval(['save -append mE  ' save_list]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'de')
    for ic=sc_list
      eval(irf_ssub('load mE dvE?;load mBPP dBPP?;load mV dV?; tt=dvE?(1,1);db=dBPP?; dv=dV?;clear dBPP? dV?;',ic));
      evxb=irf_resamp(irf_tappl(irf_cross(db,dv),'*1e-3*(-1)'),tt);
      eval(irf_ssub('dE?=irf_add(1,dvE?,-1,evxb);dE?(:,4)=0;',ic));
      eval(irf_ssub('save_list=[save_list '' dE? ''];',ic));
    end
    eval(['save -append mE  ' save_list]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % deo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'deo')
    disp('Estimating dEo where dEo.B=0, Eo (GSE) and Vo=Eo/B')
    deg=input('B angle with respect to the spin plane should be at least x deg, x=');
    qb=input('To use FGM high res (1) or PP data (2) [1] >');if isempty(qb),qb=1;end
    for ic=sc_list
      eval(irf_ssub('load mE dE?;tt=dE?(1,1);',ic));
      if qb ==1, eval(irf_ssub('load mB dB?; db=irf_resamp(dB?,dE?);clear dB?;',ic));
      else, eval(irf_ssub('load mBPP dBPP?; db=irf_resamp(dBPP?,dE?);clear dBPP?;',ic));
      end
      eval(irf_ssub('[dEo?,d?]=irf_edb(dE?,db,deg);Eo?=c_gse2dsc(dEo?,[tt ?],-1);indzero=find(abs(d?)<10);Eo?(indzero,4)=0;',ic));
      eval(irf_ssub('dVo?=irf_e_vxb(dEo?,db,-1);Vo?=c_gse2dsc(dVo?,[tt ?],-1);',ic));
      eval(irf_ssub('save_list=[save_list '' dEo? d? Eo? Vo? ''];',ic));
    end
    eval(['save -append mE  ' save_list]);
    save_list = '';
    
  elseif strcmp(q,'es') % create E ascii files
    for ic=sc_list
      % E_GSE file creation
      eval(irf_ssub('load mE Eo? d?;tt=Eo?(1,1);x=Eo?;x(:,end+1)=d?;',ic));
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
      eval(irf_ssub('load mE dE?;tt=dE?(1,1);x=dE?(:,[1 2 3]);',ic));
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E ASCII
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'ea')||strcmp(q,'esa') % create E ascii files
    if strcmp(q,'ea'), s = ''; else, s = 's'; end
    for ic=sc_list
      eval(irf_ssub(['if ~exist(''E' s '?'') & exist(''mEdB.mat'',''file''), load mEdB E' s '? ang_limit?;disp(''Loading E' s '?, ang_limit? from mEdB'');end'],ic));
      eval(irf_ssub('if ~exist(''D?p12p34'') & ~exist(''Ddsi?''), load mEDSI D?p12p34 Ddsi? Da?p12 Da?p34 Damp?;disp(''Loading offset values from mEDSI.mat'');end',ic));
      eval(irf_ssub(['if ~exist(''diE' s '?'')& exist(''mEdB.mat'',''file''), load mEdB diE' s '? ang_limit?;disp(''Loading diE' s '?, ang_limit? from mEdB'');end'],ic));
      offset_comment = 'Offsets => ';
      if exist(irf_ssub('Damp?',ic),'var')
        eval(irf_ssub('offset_comment=[offset_comment '' Damp='' num2str(Damp?)];',ic));
      end
      if exist(irf_ssub('Da?p12',ic),'var')
        eval(irf_ssub('offset_comment=[offset_comment '' Dap12='' num2str(Da?p12)];',ic));
      end
      if exist(irf_ssub('Da?p34',ic),'var')
        eval(irf_ssub('offset_comment=[offset_comment '' Dap34='' num2str(Da?p34)];',ic));
      end
      if exist(irf_ssub('D?p12p34',ic),'var')
        eval(irf_ssub('offset_comment=[offset_comment '' Dp12p34='' num2str(D?p12p34)];',ic));
      end
      if exist(irf_ssub('Ddsi?',ic),'var')
        eval(irf_ssub('offset_comment=[offset_comment '' Ddsi(x,y)=('' num2str(real(Ddsi?)) '','' num2str(imag(Ddsi?)) '')''];',ic));
      end
      offset_comment=[offset_comment '\n'];
      eval(irf_ssub('ang_limit_s=num2str(ang_limit?);',ic)),
      if eval(irf_ssub(['exist(''E' s '?'')'],ic))
        eval(irf_ssub(['number_of_points=size(E' s '?,1);'],ic));
        disp(['E' s num2str(ic) ' --> E' s num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
        E_add_comment=[offset_comment '\nang_limit=' ang_limit_s '\n'];
        E_add_comment=[E_add_comment 'E.B=0 used only for points in which magnetic field makes an angle \nwith respect to the spin plane that is larger than ang_limit'];
        eval(irf_ssub(['c_export_ascii(E' s '?,''E' s '?'',''' E_add_comment ''');'],ic));
      end
      if eval(irf_ssub(['exist(''diE' s '?'')'],ic))
        eval(irf_ssub(['number_of_points=size(diE' s '?,1);'],ic));
        disp(['diE' s num2str(ic) ' --> diE' s num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
        diE_add_comment=[offset_comment '\nang_limit=' ang_limit_s '\nE.B=0 used to estimate Ez for points in which magnetic field makes an \nangle with respect to the spin plane that is larger than ang_limit'];
        eval(irf_ssub(['c_export_ascii(diE' s '?,''diE' s '?'',''' diE_add_comment ''');'],ic));
      end
      clear E_add_comment diE_add_comment number_of_points;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % V=ExB ASCII
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'va')||strcmp(q,'vsa')
    if strcmp(q,'va'), s = ''; else, s = 's'; end
    for ic=sc_list
      % GSE
      eval(irf_ssub(['if ~exist(''VExB' s '?'') & exist(''mEdB.mat'',''file''), load mEdB VExB' s '? ang_limit?;disp(''Loading VExB' s '?, ang_limit? from mEdB'');end'],ic));
      eval(irf_ssub('ang_limit_s=num2str(ang_limit?);',ic)),
      if eval(irf_ssub(['exist(''VExB' s '?'')'],ic))
        eval(irf_ssub(['number_of_points=size(VExB' s '?,1);'],ic));
        disp(['VExB' s num2str(ic) ' --> VExB' s num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
        E_add_comment=['\nang_limit=' ang_limit_s '\n'];
        E_add_comment=[E_add_comment 'E.B=0 used only for points in which magnetic field makes an angle \nwith respect to the spin plane that is larger than ang_limit'];
        eval(irf_ssub(['c_export_ascii(VExB' s '?,''VExB' s '?'',''' E_add_comment ''');'],ic));
      end
      % DSI
      eval(irf_ssub(['if ~exist(''diVExB' s '?'',''var'') & exist(''mEdB.mat'',''file''), load mEdB diVExB' s '? ang_limit?;disp(''Loading diVExB' s '?, ang_limit? from mEdB'');end'],ic));
      eval(irf_ssub('ang_limit_s=num2str(ang_limit?);',ic)),
      if eval(irf_ssub(['exist(''diVExB' s '?'')'],ic))
        eval(irf_ssub(['number_of_points=size(diVExB' s '?,1);'],ic));
        disp(['diVExB' s num2str(ic) ' --> diVExB' s num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
        diE_add_comment=['\nang_limit=' ang_limit_s '\nE.B=0 used to estimate Ez for points in which magnetic field makes an \nangle with respect to the spin plane that is larger than ang_limit'];
        eval(irf_ssub(['c_export_ascii(diVExB' s '?,''diVExB' s '?'',''' diE_add_comment ''');'],ic));
      end
      clear E_add_comment diE_add_comment number_of_points;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EDI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'edi')
    save_file='./mEDI.mat';
    for ic=sc_list
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'EDI', ['E_xyz_gse__C' num2str(ic) '_PP_EDI'], ' ', ' ',' ');
      eval(irf_ssub('EDI?=[double(t) double(data)''];',ic));clear t,data;
      if eval(['min(size(EDI' num2str(ic) '))';])==0 % if there are no data
        disp(['CSDS...EDI' num2str(ic) '... no data']);
      else
        eval(irf_ssub('dEDI?=c_gse2dsc(EDI?,?);',ic));
        [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'EDI', ['V_ed_xyz_gse__C' num2str(ic) '_PP_EDI'], ' ', ' ',' ');
        eval(irf_ssub('VEDI?=[double(t) double(data)''];',ic));clear t,data;
        disp(irf_ssub(['CSDS....EDI?,dEDI?,VEDI? -> ' save_file],ic));
        save_list=[save_list irf_ssub(' EDI? dEDI? VEDI? ',ic)];
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WBD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'wbdwf')
    save_file='./mWBD.mat';
    for ic=sc_list
      data = getData(ClusterDB,toepoch(start_time),Dt,ic,'wbdwf','nosave');
      if ~isempty(data)
        c_eval('wfWBD?=data{2};',ic);
        disp(irf_ssub(['wfWBD? -> ' save_file],ic));
        save_list=[save_list irf_ssub(' wfWBD? ',ic)];
      end
      clear data
    end
    disp('!!! please check whether the data is E or B !!!')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif q == 'p'
    mode = irf_ask(...
      'Sampling 1)lx, 2)hx, 3)4kHz_any, 4)32kHz_any? If different give as vector. [%]','mode',1);
    for ic=sc_list
      if (length(mode)>1), mm=mode(ic);else, mm=mode;end
      if (mm == 1), param='10Hz'; tmmode='lx';
      elseif (mm == 2)
        %% Find TapeMode
        [ok,tm] = c_load('mTMode?',ic);
        if ~ok
          [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', ...
            num2str(ic),'efw','FDM');
          if ~isempty(data), tm=double(data(5,:)); else, error('Cannot fetch FDM'), end
          eval(irf_ssub('mTMode?=tm;',ic));
          if exist('./mTMode.mat','file')
            eval(irf_ssub('save -append mEFWR mTMode?;',ic));
          else, eval(irf_ssub('save mEFWR mTMode?;',ic));
          end
        end
        if tm~=tm(1)*ones(size(tm))
          warning('tape mode changes during the selected tile inteval')
        end
        tm=tm(1);
        if tm==3, param='180Hz'; tmmode='hx';
        else, param='10Hz'; tmmode='lx'; end
        clear tm
      elseif (mm == 3), param='4kHz';tmmode='any';
      elseif (mm == 4), param='32kHz';tmmode='any';
      end
      
      for probe=1:4
        disp(['EFW...sc' num2str(ic) '...probe' num2str(probe) '->P' param num2str(ic) 'p' num2str(probe)]);
        [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', ['p' num2str(probe)],param, tmmode);
        eval(irf_ssub(['p!=[double(t) double(real(data))];save_list=[save_list '' P' param '?p!''];P' param '?p!=p!;'],ic,probe)); clear t data;
      end
      clear p;
      if size(p1)==size(p2)&&size(p1)==size(p3)&&size(p1)==size(p4)&&size(p1)~=[0 0]&&ic~=2  % sc2 has often problems with p3
        p=[p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
        disp('Vps = (p1+p2+p3+p4)/4, satellite potential is average over 4 probes');
      elseif size(p1)==size(p2)&&size(p1)~=[0 0]
        p=[p1(:,1) (p1(:,2)+p2(:,2))/2];
        disp('Vps = (p1+p2)/2, satellite potential is average over probes 1 and 2');
      elseif size(p3)==size(p4)&&ic~=2
        p=[p3(:,1) (p3(:,2)+p4(:,2))/2];
        disp('Vps = (p3+p4)/2, satellite potential is average over probes 3 and 4');
      else
        p=p4;
        disp('Vps = p4, satellite potential is put to potential of probe 4');
      end
      eval(irf_ssub(['P' param '?=p;save_list=[save_list '' P' param '? ''];'],ic));
      if ((mm==1) || (mm==11)); eval(irf_ssub('P?=p;NVps?=c_efw_scp2ne(p);save_list=[save_list '' P? NVps?''];',ic));end
      if (mm == 3) || (mm == 4)
        dtburst=input(['s/c' num2str(ic) ', time shift to obtain correct time (get from Anders Tjulin) =']);
        dtburst=double(dtburst);
        for probe=1:4,eval(irf_ssub(['P' param '?p!(:,1)=P' param '?p!(:,1)+dtburst;'],ic,probe));end
        c_eval(['P' param '?(:,1)=P' param '?(:,1)+dtburst;'],ic);
      end
    end
    if exist('./mP.mat','file'), eval(['save mP ' save_list ' -append']); else, eval(['save mP ' save_list]); end
    save_list = '';
    
  elseif strcmp(q,'pa') % create V_sc n ascii files
    for ic=sc_list
      % Vsc_N file creation
      %t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
      %file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_Vps_N_sc' num2str(ic) '.dat'];
      eval(irf_ssub('if ~exist(''NVps?''), load mP NVps?;end',ic));
      if eval(irf_ssub('exist(''NVps?'')',ic))
        eval(irf_ssub('number_of_points=size(NVps?,1);',ic));
        disp(['NVps' num2str(ic) ' --> NVps' num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
        eval(irf_ssub('c_export_ascii(NVps?);',ic));
      end
    end
    
  elseif strcmp(q,'r') || strcmp(q,'v')
    save_file='./mR.mat';save_list=[];
    for ic=sc_list
      disp(['...R' num2str(ic) '--> mR.mat']);
      [tr,data_r] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
      eval(irf_ssub('R?=[double(tr) double(data_r)''];',ic));
      disp(['...V' num2str(ic) '--> mR.mat']);
      [tv,data_v] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
      eval(irf_ssub('V?=[double(tv) double(data_v)''];',ic));
      clear tr tv data_r data_v;
      save_list=[save_list irf_ssub(' R? V? ', ic)];
    end
    
  elseif strcmp(q,'vc')
    for ic=sc_list
      disp(['...VCp' num2str(ic) ', dVCp' num2str(ic)]);
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_p_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' ');
      eval(irf_ssub('VCp?=[double(t) double(real(data))''];',ic));clear t,data;
      eval(irf_ssub('if size(VCp?), dVCp?=c_gse2dsc(VCp?,[VCp?(1,1) ic]); save_list=[save_list '' VCp? dVCp? '']; end; ',ic));
      disp(['...VCh' num2str(ic) ', dVCh' num2str(ic)]);
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_HIA_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' ');
      eval(irf_ssub('VCh?=[double(t) double(real(data))''];',ic));clear t,data;
      eval(irf_ssub('if size(VCh?), dVCh?=c_gse2dsc(VCh?,[VCh?(1,1) ic]); save_list=[save_list '' VCh? dVCh? '']; end;',ic));
    end
    eval(['save mCIS ' save_list]);
    save_list = '';
    
  elseif strcmp(q,'vce')
    CIS=load('mCIS');
    for ic=sc_list
      eval(irf_ssub('if isfield(CIS,''VCp?''); vp=CIS.VCp?;else,vp=[];end;  if isfield(CIS,''VCh?'');vh=CIS.VCh?;else,vh=[];end; load mBPP BPP?;b=BPP?; clear BPP? VCp? VCh?;',ic));
      %
      if min(size(vp)) ~= 0
        disp(['...VCEp' num2str(ic)]);
        evxb=irf_tappl(irf_cross(vp,b),'*(-1e-3)');
        eval(irf_ssub('VCEp?=evxb;save -append mCIS VCEp?;',ic));
        disp(['...dVCEp' num2str(ic)]);
        eval(irf_ssub('dVCEp?=c_gse2dsc(VCEp?,[VCEp?(1,1) ic]);save -append mCIS dVCEp?;',ic));
      end
      if min(size(vh)) ~= 0
        disp(['...VCEh' num2str(ic)]);
        evxb=irf_tappl(irf_cross(vh,b),'*(-1e-3)');
        eval(irf_ssub('VCEh?=evxb;save -append mCIS VCEh?;',ic));
        disp(['...dVCEh' num2str(ic)]);
        eval(irf_ssub('dVCEh?=c_gse2dsc(VCEh?,[VCEh?(1,1) ic]);save -append mCIS dVCEh?;',ic));
      end
    end
    
    % If line is not recognized evaluate it in matlab
  else
    eval(q,'');
  end
  
  % If flag_save is set, save variables to specified file
  if flag_save==1 && length(save_file)>0 && ~isempty(save_list)
    if exist(save_file,'file')
      eval(['save -append ' save_file ' ' save_list]);
    else
      eval(['save ' save_file ' ' save_list]);
    end
  end
  
end

if exist('db'), Mat_DbClose(db); end


