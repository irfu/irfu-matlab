% get cluster sc data into matlab files,
% despinn if necessary, uses the neareset calibration values
% export to ASCII files
% uses isGetDataLite
%
% See also C_GET_BATCH
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


if exist('sc_list','var') == 0, sc_list=1:4;end % default values

mmm =  ...
  ['q        end, quit                              ';
  '1        time interval                          ';
  '2        sc list [1:4]                          ';
  'a        load phase A?                          ';
  'eph      load r,v,sax R?, V?, diV?, SAX?        ';
  'whip     load WHI active pulses WHIP?           ';
  'e        load raw E-field wE?p12,34             ';
  'p        potential P?, P10Hz?p1..4, NVps?       ';
  'ps       spin resolution P Ps?                  ';
  '---------------------------                     ';
  'dies     spin fit E-field diEs?p{12,34}         ';
  'die      despun full-res E-field diE?p1234      ';
  'idies    diEs?p{12,34} to iner fr idiEs?p{12,34}';
  'idie     diEp1234 to inert ref fr idiE?p1234    ';
  'edbs     spin res E with Ez(E.B=0) Es?,diEs?    ';
  'edb      full res E with Ez(E.B=0) E?,diE?      ';
  'iedbs    Es?,diEs? to inert ref fr iEs?,idiEs?  ';
  'iedb     E?,diE? to inert ref frame iE?,idiE?   ';
  '---------------------------                     ';
  'eburst   burst E-field wbE?p{12,34}             ';
  'pburst   P{4kHz,32kHz}?p1..4, wbE?p{12,34}      ';
  'dieburst despun burst E-field dibE?p1234        ';
  '---------------------------                     ';
  'b        load B PP BPP?, diBPP?                 ';
  'bfgm     high-res B FGM B?, diB?                ';
  'br       resample B to E Br?, diBr?             ';
  'brs      resample B to Es Brs?, diBrs?          ';
  'bsc      load B STAFF-SC wBSC?                  ';
  'dibsc    despun B STAFF-SC diBSC? BSC?          ';
  'ncis     CIS density NC{p,h}?                   ';
  'tcis     CIS temperature T{par,perp}C{p,h}?     ';
  'vcis     CIS vel VC{p,h}?,diVC{p,h}?            ';
  'vce      CIS VxB VCE{p,h}?, diVCE{p,h}?         ';
  'iedi     EDI E inertial fr iEDI?, idiEDI?       ';
  'edi      EDI E sc fr EDI?, diEDI?               ';
  'wbdwf    WBD E/B wfWBD?                         ';
  'whinat   WHISPER E WHINAT?                      ';
  %	'---------------------------                     ';
  %	'ea,esa   (s-fit)  E ->ascii E?                  ';
  %	'va,vsa   (s-fit)  ExB ->ascii VExB?             ';
  %	'pa       P ->ascii P?                           ';
  %	'x        free format load                       ';
  %	'bf load high-res FGM B1...B4          ';
  %   'bs load BS1...BS4 ->   mBS            ';
  %	'db despinned dBPP1... -> mBPP         ';
  %	'dbf despinned dB1... -> mB            ';
  %	'dibf despinned diB1... -> mB          ';
  ];

disp('-------------- Get Cluster II data ----------------');
disp(mmm); % show onces menu

q='0';flag_save=1;
while(q ~= 'q') % ====== MAIN LOOP =========
  
  q=input('input>','s');if isempty(q),q='0';end
  save_list='';save_file='';
  if strcmp(q,'q'), return,
  elseif strcmp(q,'0') || strcmp(q,'h'), disp(mmm);
  elseif q == '1'
    DB_S = c_ctl(0,'isdat_db');
    DP_S = c_ctl(0,'data_path');
    DATABASE = irf_ask('Database as string [%]>','DATABASE',DB_S);
    data_dir = irf_ask('Data path [%]>','data_dir',DP_S);
    cdb = ClusterDB(DATABASE,data_dir,pwd);
    if ~exist('start_time_s','var')
      [iso_t,Dt] = caa_read_interval; %#ok<NASGU>
      if ~isempty(iso_t)
        start_time_s = sprintf('%s %s %s %s %s %s', iso_t(1:4), ...
          iso_t(6:7),iso_t(9:10),iso_t(12:13),iso_t(15:16),iso_t(18:19));%#ok<NASGU>
      end
      clear iso_t
    end
    start_time_s = irf_ask('Start time [%]>','start_time_s','2001 02 01 00 00 00');
    start_time = eval(['[' start_time_s ']']);
    start_date_str = strrep(datestr(start_time,29),'-','');
    Dt = irf_ask('How many seconds of data [%]>','Dt',60);
    tint_epoch = toepoch(start_time) + [0 Dt];
  elseif strcmp(q,'2')
    % define sc_list
    sc_list = irf_ask('Spacecraft list [%]>','sc_list',1:4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ephemeris
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'eph')
    var_list = {'r','v','sax'};
    for j=1:length(var_list)
      for ic=sc_list, getData(cdb,tint_epoch(1),Dt,ic,var_list{j}); end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'p')
    for ic=sc_list
      getData(cdb,tint_epoch(1),Dt,ic,q);
      getData(ClusterProc(pwd),ic,q);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ClusterDB/getData quantities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'r') || strcmp(q,'v') || strcmp(q,'a') || strcmp(q,'whip') || ...
      strcmp(q,'e') || strcmp(q,'b') || ...
      strcmp(q,'bfgm') || strcmp(q,'bsc') || ...
      strcmp(q,'p') || strcmp(q,'pburst') || strcmp(q,'eburst') || ...
      strcmp(q,'ncis') || strcmp(q,'tcis') || strcmp(q,'vcis') || ...
      strcmp(q,'vce') || strcmp(q,'wbdwf') || strcmp(q,'sax')
    for ic=sc_list, getData(cdb,tint_epoch(1),Dt,ic,q); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EDI Inertial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'iedi')
    for ic=sc_list, getData(cdb,tint_epoch(1),Dt,ic,'edi'); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WHISPER Natural
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'whinat')
    for ic=sc_list, getData(cdb,tint_epoch(1),Dt,ic,'whinat'); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ClusterProc/getData quantities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'edi') || strcmp(q,'ps') || strcmp(q,'idies') || strcmp(q,'idie') || ...
      strcmp(q,'dieburst') || strcmp(q,'vedbs') || strcmp(q,'vedb') || ...
      strcmp(q,'br') || strcmp(q,'brs') || strcmp(q,'dibsc')
    for ic=sc_list, getData(ClusterProc(pwd),ic,q); end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spinfits
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'dies')
    es_rmwhip = irf_ask('Use times when Whisper pulses are present(y/n) [%]',...
      'es_rmwhip','n');
    if strcmpi(es_rmwhip,'y')
      for ic=sc_list, getData(ClusterProc(pwd),ic,'dies','withwhip'); end
    else
      es_rmwhip = 'n';
      for ic=sc_list, getData(ClusterProc(pwd),ic,'dies'); end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % despin E
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'die')
    e_rmwhip = irf_ask('Use times when Whisper pulses are present for ADC offset (y/n) [%]',...
      'es_rmwhip','n');
    if strcmpi(e_rmwhip,'y')
      for ic=sc_list, getData(ClusterProc(pwd),ic,'die','withwhip'); end
    else
      e_rmwhip = 'n';
      for ic=sc_list, getData(ClusterProc(pwd),ic,'die'); end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E.B=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'edbs') || strcmp(q,'edb') || strcmp(q,'iedb') || strcmp(q,'iedbs')
    ang_limit = irf_ask('Minimum angle(B,spin plane) [%]>','ang_limit',10);
    below_ang_limit = ...
      irf_ask('Points < min_ang (0-Ez to NaN, 1-Ez to 1e27, 2-use Ez=0)[%]>',...
      'below_ang_limit',0);
    switch below_ang_limit
      case 0
        flag_below_ang_limit = 'ang_blank';
      case 1
        flag_below_ang_limit = 'ang_fill';
      case 2
        flag_below_ang_limit = 'ang_ez0';
      otherwise
        error('Must be 0, 1 or 2')
    end
    if q(end)=='s'
      probe_p = irf_ask('Probe pair [%]>','probe_p',34);
      for ic=sc_list
        getData(ClusterProc(pwd),ic,q,flag_below_ang_limit,...
          'ang_limit',ang_limit,'probe_p',probe_p);
      end
    else
      for ic=sc_list
        getData(ClusterProc(pwd),ic,q,flag_below_ang_limit,...
          'ang_limit',ang_limit);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  elseif strcmp(q,'x')
    var_name=input('matlab variable name =','s');
    disp('? in input is substituted by cluster number');
    qstring=input('example {''Cluster'' ''1'' ''efw'' ''E'' ''p12'' ''10Hz'' ''any''}=>','s');
    for ic=sc_list
      string=eval(irf_ssub(qstring,ic));
      varic = [var_name num2str(ic)];
      disp(['...free format s/c' num2str(ic) ' ' varic]);
      for jj=1:size(string,2),str{jj}=irf_ssub(string{jj},ic);end, for jj=size(string,2)+1:7,str{jj}=' ';end
      [t,data] = isGetDataLite( DATABASE, start_time, Dt,str{1}, str{2}, str{3}, str{4},str{5},str{6},str{7});
      eval([varic '=[double(t) double(data)];']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Magnetic fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
  elseif strcmp(q,'dbf')
    for ic=sc_list
      eval(irf_ssub('load mB B?;tt=B?(1,1);',ic));
      eval(irf_ssub('dB?=c_gse2dsc(B?,[B?(1,1) ic]);',ic));
      eval(irf_ssub('save -append mB dB?;',ic));
    end
    
  elseif strcmp(q,'dibf')
    for ic=sc_list
      eval(irf_ssub('load mB B?;',ic));
      eval(irf_ssub('load mEPH SAX?;',ic));
      eval(irf_ssub('diB?=c_gse2dsi(B?,SAX?);',ic));
      eval(irf_ssub('save -append mB diB?;',ic));
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
  elseif strcmp(q,'va') || strcmp(q,'vsa')
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
    % P ASCII
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VCp ASCII
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(q,'vc')
    for ic=sc_list
      disp(['...VCp' num2str(ic) ', dVCp' num2str(ic)]);
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_p_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' '); %#ok<NASGU>
      eval(irf_ssub('VCp?=[double(t) double(real(data))''];',ic)); clear t data
      eval(irf_ssub('if size(VCp?), dVCp?=c_gse2dsc(VCp?,[VCp?(1,1) ic]); save_list=[save_list '' VCp? dVCp? '']; end; ',ic));
      disp(['...VCh' num2str(ic) ', dVCh' num2str(ic)]);
      [t, data] = isGetDataLite( db, start_time, Dt, 'CSDS_PP', ['C' num2str(ic)], 'CIS', ['V_HIA_xyz_gse__C' num2str(ic) '_PP_CIS'], ' ', ' ',' ');
      eval(irf_ssub('VCh?=[double(t) double(real(data))''];',ic)); clear t data
      eval(irf_ssub('if size(VCh?), dVCh?=c_gse2dsc(VCh?,[VCh?(1,1) ic]); save_list=[save_list '' VCh? dVCh? '']; end;',ic));
    end
    eval(['save mCIS ' save_list]);
    save_list = '';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If line is not recognized evaluate it in matlab
  else
    eval(q,'');
  end
  
  % If flag_save is set, save variables to specified file
  if flag_save==1 && ~isempty(save_file) && ~isempty(save_list)
    if exist(save_file,'file')
      eval(['save -append ' save_file ' ' save_list]);
    else
      eval(['save ' save_file ' ' save_list]);
    end
  end
  
end

if exist('db','var'), Mat_DbClose(db); end


