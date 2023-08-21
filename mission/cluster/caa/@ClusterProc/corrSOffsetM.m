function corrSOffsetM(cp,cl_id)
%corrSOffset correct the Sunward offset and amplitude factor
% do a manual correction of the offsets by comparison with
% CIS and EDI data (must be loaded)
%
% corrSOffset(cp,cl_id)
% cp - ClusterProc object
% cl_id - SC#
%
% Offset variable names
%  Ddsi#     - sunward offset real(Ddsi#) is X offset and imag(Ddsi#) is Y offset
%  Damp#     - factor with which to multiply amplitude
%  D#p12p34  - difference in sunward offset between different probe pairs
%              first value gives offset difference in X- and second in Y-direction
% Example:
% corrSOffsetM(ClusterProc('/home/yuri/caa-data/20020304'),1)
%
%  See also: C_CAL_GUI
%

% Copyright 2004,2007 Yuri Khotyaintsev (yuri@irfu.se)
%


old_pwd = pwd;
cd(cp.sp) %enter the storage directory


[ok, diE] = c_load('diE?p1234',cl_id);
if ~ok, irf_log('load','no diE{cl_id}p1234 data')
end
[ok, diEs] = c_load('diEs?p34',cl_id);
if ok
  % Load offsets
  dsiof = c_ctl(cl_id,'dsiof');
  if isempty(dsiof)
    [dsiof_def, dam_def] = c_efw_dsi_off(diEs(1,1),cl_id);

    [ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
    [ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end

    if ok1 || ok2, irf_log('calb',...
        ['Saved DSI offsets on C' num2str(cl_id)])
      %else irf_log('calb','Using default DSI offsets')
    end
    clear dsiof_def dam_def
  else
    Ddsi = dsiof(1); Damp = dsiof(2);
    irf_log('calb',['User DSI offsets on C' num2str(cl_id)])
  end
  clear dsiof

  offset(1) = Ddsi; clear Ddsi
  irf_log('proc',sprintf('Ddsi=%.2f', offset(1)))
  offset(2) = Damp; clear Damp
  irf_log('proc',sprintf('Damp=%.2f', offset(2)))

else
  error('caa:noData','no diEs{cl_id}p34 data in mEDSI')
end

if exist('diE','var')
  have_hres = 1;
  % remove points larger then 1 V/m
  for j=2:3, diE(abs(diE(:,j)) > 1000, j) = NaN; end
else, have_hres = 0;
end
% we load full res data, but plot only spin.
var_list = 'diEs_tmp';
var_list1 = 'diEsp34';

% load CIS
var = {'diVCEp', 'diVCEh'};
if exist('./mCIS.mat','file')
  CIS = load('mCIS');
  for i=1:length(var)
    c_eval(['if isfield(CIS,''' var{i} '?''); ' var{i} '=CIS.' var{i} '?; end; clear ' var{i} '?'], cl_id)
  end
  clear CIS
end
if ~exist('diVCEp','var') || ~exist('diVCEh','var')
  warning('caa:noData','no CIS data loaded')
else
  for i=1:length(var)
    if exist(var{i},'var')
      eval(['ll=length(find(~isnan(' var{i} '(:,2:end))));'])
      if ll>0
        var_list = [var_list ',' var{i}];
        var_list1 = [var_list1 ',' var{i}];
      end
      clear ll
    end
  end
end
clear var

% load EDI
if exist('./mEDI.mat','file')
  EDI = load('mEDI');
  var = 'diEDI';
  eval(irf_ssub(['if isfield(EDI,''' var '?''); ' var '=EDI.' var '?; end; clear ' var '?'], cl_id));
  clear EDI
end

if ~exist('diEDI','var')
  warning('caa:noData','no EDI data loaded')
else
  if any(~isnan(diEDI(:,2:end)))
    var_list = [var_list ',diEDI'];
    var_list1 = [var_list1 ',diEDI'];
  end
end

var_list0 = [];
if exist('./mP.mat','file')
  c_eval('load mP P?; P=P?; clear P?', cl_id)
  var_list0 = 'P';
end
if exist('./mBPP.mat','file')
  c_eval('load mBPP diBPP?; diB=irf_abs(diBPP?); clear diBPP?', cl_id)
  var_list0 = [var_list0 ',diB'];
end

diE_tmp = diE;
diE_tmp(:,2) = diE_tmp(:,2) - real(offset(1));
diE_tmp(:,3) = diE_tmp(:,3) - imag(offset(1));
diE_tmp(:,2:3) = diE_tmp(:,2:3)*real(offset(2));
diEs_tmp = diEs;
diEs_tmp(abs(diEs(:,2))>1e4,:) = []; % remove spinfits that has given large values
diEs_tmp(:,2) = diEs_tmp(:,2) - real(offset(1));
diEs_tmp(:,3) = diEs_tmp(:,3) - imag(offset(1));
diEs_tmp(:,2:3) = diEs_tmp(:,2:3)*real(offset(2));

figure(17)
clf
t = tokenize(var_list1,',');
leg = ['''' t{1} ''''];
for i=2:length(t), leg = [leg ',''' t{i} '''']; end

t0 = tokenize(var_list0,',');
t1 = tokenize(var_list,','); dvar = t1{1};

for k=1:length(t1)
  for j=2:3, eval([ t1{k} '(find(abs(' t1{k} '(:,j)) > 1000), j) = NaN;']),end
end

dummy = dvar;
for j=2:2+length(t0), dummy = [dummy ',' dvar]; end
eval(['h=irf_plot({' dummy '});'])

%Ex,Ey
for co=1:2
  dummy = [t1{1} '(:,1),' t1{1} '(:,' num2str(co+1) ')'];
  if length(t1)>1
    for j=2:length(t1)
      dummy = [dummy ',' t1{j} '(:,1),' t1{j} '(:,' num2str(co+1) ')'];
    end
  end

  axes(h(co))
  eval(['plot(' dummy ');'])
  irf_timeaxis
  grid
  set(gca,'XTickLabel',[])
  if co==1, ylabel('E_x DSI'), else, ylabel('E_y DSI'), end
  xlabel('')
  zoom on
end

axes(h(1))
title(sprintf(...
  'Cluster %d : offset X %.2f [mV/m], offset Y %.2f [mV/m], amplitude factor %.2f',...
  cl_id,real(offset(1)),imag(offset(1)),offset(2)))


% AUX information
for j=1:length(t0)
  axes(h(2+j))
  eval(['irf_plot(' t0{j} '); ylabel(''' t0{j} ''')'])
end

irf_pl_add_info
irf_figmenu
for j=1:length(t0)+2, axes(h(j)), end
eval(['legend(h(1),' leg ')'])
zoom on

q='0';
while(q ~= 'q')
  flag_replot=0;
  q_s = 'Delta Ex,Ey [mV/m], amplitude factor (s,r,';
  if have_hres, q_s = [q_s 'f,']; end
  q_s = [q_s 'q,h-help)[%]>'];
  q=irf_ask(q_s,'', ...
    num2str([real(offset(1)) imag(offset(1)) real(offset(2))],'%.2f '));
  switch(q)
    case 'h'
      disp('Usage:')
      disp('  give 3 numders for offsets in Ex, Ey and amplitude')
      disp('  correction or one of the following commands:')
      disp('  s - save calibration parameters')
      disp('  r - read from disk')
      if have_hres, disp('  f - show/hide full resolution data'), end
      disp('  q - quit')
      disp('  h - help (this message)')
    case 'f'
      if strcmp(var_list(1:7),'diE_tmp')
        %do not plot high resolution data
        var_list = var_list(9:end);
        leg = leg(12:end);
        flag_replot=1;
      else
        %plot high resolution data if we have its
        if have_hres
          var_list = ['diE_tmp,' var_list];
          leg = ['''diEp1234'',' leg];
          flag_replot=1;
        end
      end
    case 'r'
      disp(sprintf('Reading Ddsi%d, Damp%d <- ./mEDSI.mat',cl_id,cl_id))
      eval(irf_ssub('load mEDSI  Ddsi? Damp?; if exist(''Ddsi?''), offset(1)=Ddsi?; offset(2)=Damp?; else, disp(''Cannot find callibrations''); offset=[0+0i 1]; end;',cl_id))
      flag_replot=1;
    case 's'
      disp(sprintf('Ddsi%d, Damp%d -> ./mEDSI.mat',cl_id,cl_id))
      eval(irf_ssub('Ddsi?=offset(1); Damp?=offset(2);save -append mEDSI Ddsi? Damp?',cl_id))
    case 'q'
      cd(old_pwd)
      return
    otherwise
      [o_tmp,ok] = str2num(q);
      if ok
        if length(o_tmp) > 2
          offset(1) = o_tmp(1)+1i*o_tmp(2);
          offset(2) = o_tmp(3);
        elseif length(o_tmp) > 1, offset(1) = o_tmp(1)+1i*o_tmp(2);
        else, offset(1) = o_tmp(1)+1i*imag(offset(1));
        end
      else, disp('invalid command')
      end
      flag_replot=1;
  end
  if flag_replot
    diE_tmp = diE;
    diE_tmp(:,2) = diE_tmp(:,2) - real(offset(1));
    diE_tmp(:,3) = diE_tmp(:,3) - imag(offset(1));
    diE_tmp(:,2:3) = diE_tmp(:,2:3)*real(offset(2));
    diEs_tmp = diEs;
    diEs_tmp(:,2) = diEs_tmp(:,2) - real(offset(1));
    diEs_tmp(:,3) = diEs_tmp(:,3) - imag(offset(1));
    diEs_tmp(:,2:3) = diEs_tmp(:,2:3)*real(offset(2));

    figure(17)
    zoom off
    t1 = tokenize(var_list,',');
    %Ex,Ey
    for co=1:2
      dummy = [t1{1} '(:,1),' t1{1} '(:,' num2str(co+1) ')'];
      if length(t1)>1
        for j=2:length(t1)
          dummy = [dummy ',' t1{j} '(:,1),' t1{j} '(:,' num2str(co+1) ')'];
        end
      end

      axes(h(co)), cla(h(co))
      ax_pos=get(h(co),'position');
      eval(['plot(h(co),' dummy ');'])
      set(h(co),'position',ax_pos);
      irf_timeaxis(h(co))
      grid(h(co)), set(h(co),'XTickLabel',[]), xlabel('')
      if co==1, ylabel('E_x DSI'), else, ylabel('E_y DSI'), end
    end

    axes(h(1))
    title(sprintf('Cluster %d : offset X %.2f [mV/m], offset Y %.2f [mV/m], amplitude factor %.2f',cl_id,real(offset(1)),imag(offset(1)),offset(2)))
    irf_pl_add_info
    for j=1:length(t0)+2, axes(h(j)), end
    eval(['legend(h(1),' leg ')'])
    zoom on
    flag_replot=0;
  end
end

