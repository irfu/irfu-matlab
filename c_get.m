% get cluster sc data into matlab files,
% despinn if necessary, uses the neareset calibration values
% export to ASCII files
% uses ISGETDATALITE
%
% $Revision$  $Date$

csds_dir='/data/cluster/CSDS/';

if exist('sc_list') == 0, sc_list=1:4;end % default values

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
	'x  free format load                   ';
	'---------------------------           ';
	'db despinned dBPP1... -> mBPP         ';
	'dbf despinned dB1... -> mB            ';
	'dbs despinned dBS1..dBS4 -> mBS       ';
	'dve despinned dvE1..dvE4 -> mE        ';
        'dve1 spin fits dvE1p12..dvE4p34 -> mE ';
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
 save_list='';save_file='';
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
  start_date_str=strrep(datestr(start_time,29),'-','');
  variable='Dt';default=60;question='How many seconds of data [%]>';av_ask
  tint_epoch=toepoch(start_time)+[0 Dt];
 elseif strcmp(q,'2'),
  % define sc_list
  variable='sc_list';default=[1:4];question='Spacecraft list [%]>'; av_ask;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ephemeris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(q,'a'),
    for ic=sc_list, disp(['...A' num2str(ic)]);
     [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');
     eval(av_ssub('A?=[double(t) double(data)];',ic));%clear t data;
     if flag_save==1, eval(av_ssub('if exist(''./mA.mat''),save mA A? -append; else, save mA A?;end',ic));end
    end
    save_list = '';

 elseif strcmp(q,'eph'),
    for ic=sc_list, disp(['...ephemeris' num2str(ic) '...LT,MLT,ILAT,L->mEPH...R->mR...V->mV']);
     [tlt,lt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'lt', ' ', ' ', ' ');
     [tmlt,mlt] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'mlt', ' ', ' ', ' ');
     [tL,Lshell] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'l_shell', ' ', ' ', ' ');
     [tilat,ilat] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'inv_lat', ' ', ' ', ' ');
%      eval(av_ssub('[lat,long]=av_read_cdf([csds_dir ''SP/AUX/*'' start_date_str ''*''],{''sc_at?_lat__CL_SP_AUX'',''sc_at?_long__CL_SP_AUX''});spinaxis_latlong?=[lat long(:,2)];',ic));
     [tlat, lat] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_lat__CL_SP_AUX'], ' ', ' ',' ');
     [tlong, long] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_long__CL_SP_AUX'], ' ', ' ',' ');
     [tr,r] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
     [tv,v] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'velocity', ' ', ' ', ' ');
     eval(av_ssub('',ic));clear t data;
     eval(av_ssub('LT?=[double(tlt) double(lt)];MLT?=[double(tmlt) double(mlt)];L?=[double(tL) double(Lshell)];ILAT?=[double(tilat) double(ilat)];R?=[double(tr) double(r)''];V?=[double(tv) double(v)''];',ic));clear tlt tmlt tL tilat lt mlt Lshell ilat tr r tv v;
     eval(av_ssub('spinaxis_latlong?=[double(tlat) double(lat) double(long)];',ic)); clear tlat lat long;
     eval(av_ssub('if exist(''./mEPH.mat''),save mEPH LT? MLT? L? ILAT? spinaxis_latlong? -append; else, save mEPH LT? MLT? L? ILAT? spinaxis_latlong?;end',ic));
     eval(av_ssub('tt=R?(1,1);dR?=c_gse2dsc(R?,[tt ic]);',ic));  % despinned coordinates
     eval(av_ssub('if exist(''./mR.mat''),save mR R? dR? -append; else, save mR R? dR? ;end',ic));
     eval(av_ssub('tt=V?(1,1);dV?=c_gse2dsc(V?,[tt ic]);',ic));  % despinned coordinates
     eval(av_ssub('if exist(''mV.mat''),save mV V? dV? -append; else, save mV V? dV? ;end',ic));
    end
    save_list = '';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif q == 'b',
    save_file='./mBPP.mat';
    for ic=sc_list, disp(['CSDS...BPP' num2str(ic)]);
      eval(av_ssub('BPP?=av_read_cdf([csds_dir ''PP/FGM/C?/C?_PP_FGM_'' start_date_str ''*''],''B_xyz_gse__C?_PP_FGM'');BPP?=av_t_lim(BPP?,tint_epoch);save_list=[save_list '' BPP?''];',ic));
%      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'FGM', ['B_xyz_gse__C' num2str(ic) '_PP_FGM'], ' ', ' ',' ');
%      eval(av_ssub('BPP?=[double(t) double(data)''];',ic));clear t,data;
%     eval(av_ssub('if exist(''./mBPP.mat''),save mBPP BPP? -append; else, save mBPP BPP?;end',ic));
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
  save_list = '';

 elseif strcmp(q,'bfgm'),
  disp('CONTACT STEPHAN BUCHERT!!!!!!!!!!!!!!!!!!!');
  for ic=sc_list,
    eval(av_ssub('B?=c_get_bfgm(tint_epoch,ic);',ic));
  end
  eval(['save mB ' save_list]);
  save_list = '';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E p1234 or p12 & p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 elseif strcmp(q,'e') | strcmp(q,'e1'),
  save_file = './mE.mat';
  mode=av_q('Sampling 1)hx 2)lx ? If different give as vector. [%]','mode',1);
  for ic=sc_list,
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
  		if exist('./mTMode.mat','file'), load mTMode; end
		if exist(av_ssub('mTMode?',ic),'var'), eval(av_ssub('tm=mTMode?;',ic)), end
		if ~exist('tm','var')
			[t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic),'efw','FDM');
			if ~isempty(data), tm=data(5,:); else, error('Cannot fetch FDM'), end
			if tm~=tm(1)*ones(size(tm)),warning('tape mode changes during the selected tile inteval'), end
			tm=tm(1);
			eval(av_ssub('mTMode?=tm;',ic));
			if exist('./mTMode.mat','file'), eval(av_ssub('save -append mTMode mTMode?;',ic));
			else, eval(av_ssub('save mTMode mTMode?;',ic));	
                        end
		end
		tmmode='hx';
		if tm<1e-30, param='10Hz';	else, param='180Hz'; end
		clear tm
		tst = toepoch(start_time);
		if tst>toepoch([2001 07 31 00 00 00])&tst<toepoch([2001 09 01 00 00 00]), 
			% all sc run on 180Hz filter in august 2001
			param='180Hz';	
		elseif tst>toepoch([2001 07 31 00 00 00])&ic==2, % 10Hz filtef probelm on sc2
			param='180Hz';
		end
	end
    
    if (toepoch(start_time)>toepoch([2001 12 28 03 00 00])&ic==1) | (toepoch(start_time)>toepoch([2002 07 29 09 06 59 ])&ic==3),
	p34_only = 1;
       	disp(sprintf('            !Only p34 exists for sc%d',ic));
    else,
	p34_only = 0;
    end
    
    clear sensor;
    if strcmp(q,'e')
        disp(['EFW...sc' num2str(ic) '...E ' param ' filter']);
    	if p34_only
       	    [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
       	    data = double(real(data));
       	    data = [data(:,1)*0 data(:,1) data(:,1)*0]';
        else,
            [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p1234', param, tmmode);
        end
        t = double(t);
        eval(av_ssub('wE?=[t data''];',ic)); clear t data;
        eval(av_ssub('save_list=[save_list '' wE? ''];',ic));
    else
        % do separate probes
        if ~p34_only
            disp(['EFW...sc' num2str(ic) '...E p12 ' param ' filter']);
            [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p12', param, tmmode);
            data = double(real(data));
            t = double(t);
            eval(av_ssub('wE?p12=[t data];',ic)); clear t data;
            eval(av_ssub('save_list=[save_list '' wE?p12 ''];',ic));
        end
        disp(['EFW...sc' num2str(ic) '...E p34 ' param ' filter']);
        [t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'efw', 'E', 'p34', param, tmmode);
        data = double(real(data));
        t = double(t);
        eval(av_ssub('wE?p34=[t data];',ic)); clear t data;
        eval(av_ssub('save_list=[save_list '' wE?p34 ''];',ic));
    end
  end
  %if exist('./mE.mat'), eval(['save -append mE  ' save_list]); else, eval(['save mE  ' save_list]);end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % spinfits
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(q,'dve1'), 
  save_file = './mE.mat';
  for ic=sc_list,
   eval(av_ssub('load mE wE?p12 wE?p34;',ic));
   eval(av_ssub('load mA.mat A?;',ic));
   
   if exist(av_ssub('wE?p12',ic),'var')
       eval(av_ssub('tt=wE?p12;aa=A?;',ic))
       disp(sprintf('Spin fit wE%dp12 -> dvE%dp12 mean:%.2f',ic,ic,mean(tt(:,2))))
       sp = EfwDoSpinFit(12,3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
       sp = sp(:,1:4); 
       sp(:,3) = -sp(:,3); % DSI->DS
       sp(:,4) = 0*sp(:,4);
       o12(1) = mean(sp(:,2)); o12(2) = mean(sp(:,3));
       eval(av_ssub('dvE?p12=sp;',ic))
       eval(av_ssub('save_list=[save_list '' dvE?p12 ''];',ic));
   else
       disp(sprintf('No p12 data for sc%d',ic))
   end
   
   if exist(av_ssub('wE?p34',ic),'var')
       eval(av_ssub('tt=wE?p34;aa=A?;',ic))
       disp(sprintf('Spin fit wE%dp34 -> dvE%dp34 mean:%.2f',ic,ic,mean(tt(:,2))))
       sp = EfwDoSpinFit(34,3,10,20,tt(:,1),tt(:,2),aa(:,1),aa(:,2));
       sp = sp(:,1:4); 
       sp(:,3) = -sp(:,3); % DSI->DS 
       sp(:,4) = 0*sp(:,4);
       o34(1) = mean(sp(:,2)); o34(2) = mean(sp(:,3));
       eval(av_ssub('dvE?p34=sp;',ic))
       eval(av_ssub('save_list=[save_list '' dvE?p34 ''];',ic));
   else
       disp(sprintf('No p34 data for sc%d',ic))
   end
   
   %display offsets
   if exist(av_ssub('wE?p12',ic),'var') & exist(av_ssub('wE?p34',ic),'var')
       disp(sprintf('[X Y] offsets <p12>-<p34> : [ %.2f %.2f ]', o12(1)-o34(1), o12(2)-o34(2)))
   end
 
  end
  %eval(['save -append mE  ' save_list]);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despin E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(q,'dve'),
  save_file = './mE.mat';
  q_efw_offset=av_q('How to treat offsets? \n  1) do nothing, \n  2) before despin subtract probe signal mean \n  3)  2 + use nearest sunward offsets,\n  4) nearest full hand-tuned callibration,\n[%]>','q_efw_offset',2);
  for ic=sc_list,
   eval(av_ssub('load mE wE?;tt=wE?(1,1);',ic));
   eval(av_ssub('load mA.mat A?;',ic));
   switch q_efw_offset
   case 1,      eval(av_ssub('dvE?=c_despin(wE?,A?);',ic));
   case 2,      eval(av_ssub('dvE?=c_despin(wE?,A?,?,''efw_a'');',ic));
   case 3,      eval(av_ssub('dvE?=c_despin(wE?,A?,?,''efw_b'');',ic));
   case 4,      eval(av_ssub('dvE?=c_despin(wE?,A?,?);',ic));
   otherwise,   disp('wrong offest option, using option 2.');
                eval(av_ssub('dvE?=c_despin(wE?,A?,?);',ic));
   end
   eval(av_ssub('save_list=[save_list '' dvE? ''];',ic));
  end
  %eval(['save -append mE  ' save_list]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(q,'de'),
  for ic=sc_list,
   eval(av_ssub('load mE dvE?;load mBPP dBPP?;load mV dV?; tt=dvE?(1,1);db=dBPP?; dv=dV?;clear dBPP? dV?;',ic));
   evxb=av_interp(av_t_appl(av_cross(db,dv),'*1e-3*(-1)'),tt);
   eval(av_ssub('dE?=av_add(1,dvE?,-1,evxb);dE?(:,4)=0;',ic));
   eval(av_ssub('save_list=[save_list '' dE? ''];',ic));
  end
  eval(['save -append mE  ' save_list]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deo 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  save_list = '';

 elseif strcmp(q,'es'), % create E ascii files
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

 elseif strcmp(q,'ea'), % create E ascii files
  for ic=sc_list,
     % E_GSE file creation
     eval(av_ssub('load mEdB ang_limit E? diE?;number_of_points=size(E?,1); ',ic));
     %t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
     %file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_E_GSE_sc' num2str(ic) '.dat'];
     disp(['E' num2str(ic) ' --> E' num2str(ic) '.dat ' num2str(number_of_points) ' samples']);
     E_add_comment=['ang_limit=' num2str(ang_limit) '\n'];
     E_add_comment=[E_add_comment 'E.B=0 used only for points in which magnetic field makes an angle \n with respect to the spin plane that is larger than ang_limit'];
     eval(av_ssub(['exportAscii(E?,''E?'',''' E_add_comment ''');'],ic));
     % E_DS file creation
     disp(['diE' num2str(ic) ' --> diE' num2str(ic) '.dat ' num2str(number_of_points) ' samples']);
     diE_add_comment=['ang_limit=' num2str(ang_limit) '\nE.B=0 used to estimate Ez for points in which magnetic field makes an angle with respect to the spin plane that is larger than ang_limit'];
     eval(av_ssub(['exportAscii(diE?,''E?'',''' diE_add_comment ''');'],ic));
     clear E_add_comment diE_add_comment number_of_points;
   end
 elseif strcmp(q,'edi'),
  save_file='./mEDI.mat';
    for ic=sc_list, 
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'EDI', ['E_xyz_gse__C' num2str(ic) '_PP_EDI'], ' ', ' ',' ');
      eval(av_ssub('EDI?=[double(t) double(data)''];',ic));clear t,data;
      if eval(['min(size(EDI' num2str(ic) '))';])==0 % if there are no data
       disp(['CSDS...EDI' num2str(ic) '... no data']);
      else
       eval(av_ssub('dEDI?=c_gse2dsc(EDI?,?);',ic));
       [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'EDI', ['V_ed_xyz_gse__C' num2str(ic) '_PP_EDI'], ' ', ' ',' ');
       eval(av_ssub('VEDI?=[double(t) double(data)''];',ic));clear t,data;
       disp(av_ssub(['CSDS....EDI?,dEDI?,VEDI? -> ' save_file],ic));
       save_list=[save_list av_ssub(' EDI? dEDI? VEDI? ',ic)];
      end
    end

 elseif q == 'p',
  variable='mode';default=1;question='Sampling 1)lx, 2)hx, 3)32kHz_any? If different give as vector. [%]';av_ask;
  for ic=sc_list,
  	if (length(mode)>1), mm=mode(ic);else, mm=mode;end
		if (mm == 1), param='10Hz'; tmmode='lx';
		elseif (mm == 2),
			%% Find TapeMode
  			if exist('./mTMode.mat','file'), eval(av_ssub('load mTMode;',ic)); end
			if exist(av_ssub('mTMode?',ic),'var'), eval(av_ssub('tm=mTMode?;',ic)), end
			if ~exist('tm','var')
				[t,data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic),'efw','FDM');
				if ~isempty(data), tm=data(5,:);, else, error('Cannot fetch FDM'), end
				if tm~=tm(1)*ones(size(tm)),warning('tape mode changes during the selected tile inteval'), end
				tm=tm(1);
				eval(av_ssub('mTMode?=tm;',ic));
				if exist('./mTMode.mat','file'), eval(av_ssub('save -append mTMode mTMode?;',ic));
				else, eval(av_ssub('save mTMode mTMode?;',ic));	end
			end
			if tm==3, param='180Hz'; tmmode='hx';
			else, param='10Hz'; tmmode='lx'; end
			clear tm
		elseif (mm == 3), param='32kHz';tmmode='any';end;

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
  save_list = '';

 elseif strcmp(q,'ps'), % create V_sc n ascii files
  for ic=sc_list,
     % Vsc_N file creation
     %t_ref=toepoch(fromepoch(tt).*[1 1 1 0 0 0]);time_ref=datestr(datenum(fromepoch(t_ref)),0);
     %file_name=  [time_ref([8 9 10 11 3 4 5 6 3 1 2]) '_Vps_N_sc' num2str(ic) '.dat'];
     eval(av_ssub('if exist(''NVps?''), number_of_points=size(NVps?,1);else load mP NVps?;end',ic));
     if eval(av_ssub('exist(''NVps?'')',ic)),
       disp(['NVps' num2str(ic) ' --> NVps' num2str(ic) '.dat  ' num2str(number_of_points) ' samples']);
       eval(av_ssub('exportAscii(NVps?);',ic));
     end
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
    save_list = '';

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

% If line is not recognized evaluate it in matlab 
 else
  eval(q,'');
 end
 
% If flag_save is set, save variables to specified file
 if flag_save==1 & length(save_file)>0 & ~isempty(save_list),
  if exist(save_file,'file'), 
   eval(['save -append ' save_file ' ' save_list]); 
  else, 
   eval(['save ' save_file ' ' save_list]);
  end
 end
 
end

if exist('db'), Mat_DbClose(db); end


