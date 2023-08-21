function [int1,int2,int3,int4] = caa_ex_filter(iso_st,iso_stop,cl_id)
%CAA_EX_FILTER  filer Ex data to find saturation intervals
%
% [int1,int2,int3,int4] = CAA_EX_FILTER(iso_st,iso_stop,[cl_id])
%

int1 = []; int2 = []; int3 = []; int4 = [];

if nargin<3, cl_id = (1:4); end

old_pwd = pwd;

dt = iso2epoch(iso_stop) - iso2epoch(iso_st);

h = irf_plot(4);

%%load data
for ic = cl_id
  result = [];

  dirs = caa_get_subdirs(iso_st, dt, ic);
  for dd = 1:length(dirs)
    d = dirs{dd};
    cd(d);

    irf_log('proc',sprintf('Processing C%d : %s',ic,d(end-12:end)))

    [sfit_probe,flag_lx,probeS] = caa_sfit_probe(ic);
    if flag_lx, tS = '(LX)'; else, tS = ''; end
    irf_log('proc',sprintf('L3_E probe pair : %s%s',probeS(1:2),tS))

    probe_numeric = sfit_probe;
    switch sfit_probe
      case 120, sfit_probe = 12;
      case 320, sfit_probe = 32;
      case 340, sfit_probe = 34;
      case 420, sfit_probe = 42;
      otherwise
    end

    if flag_lx, vs = irf_ssub('diELXs?p!',ic,sfit_probe);
    else, vs = irf_ssub('diEs?p!',ic,sfit_probe);
    end

    [ok,data] = c_load(vs);
    if ~ok, irf_log('load',['Failed to load ' vs]); continue, end

    data = data(:,1:2);

    if 0 % this does double job, skip this
      % Remove saturation due to too high bias current
      [ok,hbias,msg] = c_load(irf_ssub('HBIASSA?p!',ic,sfit_probe));
      if ok
        % Special trick for C2 after April 2011
        if ~isempty(hbias) && ic==2
          hbias(hbias(:,1)==iso2epoch('2011-04-30T06:00:00.00Z'),:)=[];
        end
        if ~isempty(hbias)
          irf_log('proc','blanking HB saturation')
          data = caa_rm_blankt(data,hbias);
        end
      else, irf_log('load',msg)
      end
      clear ok hbias msg
    end

    % Remove saturation
    if probe_numeric<50, probepair_list=[mod(probe_numeric,10),fix(probe_numeric/10)];
    else, probepair_list=[1 2 3 4];end
    for probe=probepair_list
      [ok,hbias,msg] = c_load(irf_ssub('PROBESA?p!',ic,probe));
      if ok
        if ~isempty(hbias)
          irf_log('proc','blanking probe saturation')
          data = caa_rm_blankt(data,hbias);
        end
      else, irf_log('load',msg)
      end
      clear ok hbias msg
    end

    % Remove whisper pulses
    [ok,whip,msg] = c_load('WHIP?',ic);
    if ok
      if ~isempty(whip)
        %irf_log('proc','blanking Whisper pulses')
        data = caa_rm_blankt(data,whip);
      end
    else, irf_log('load',msg)
    end
    clear ok whip msg

    % Remove ns_ops intervals
    ns_ops = c_ctl('get', ic, 'ns_ops');
    if isempty(ns_ops)
      c_ctl('load_ns_ops', [c_ctl('get', 5, 'data_path') '/caa-control'])
      ns_ops = c_ctl('get', ic, 'ns_ops');
    end
    if ~isempty(ns_ops)
      ns_ops_intervals = [caa_get_ns_ops_int(data(1,1), data(end,1)-data(1,1), ns_ops, 'bad_data')' ...
        caa_get_ns_ops_int(data(1,1), data(end,1)-data(1,1), ns_ops, 'bad_tm')' ...
        caa_get_ns_ops_int(data(1,1), data(end,1)-data(1,1), ns_ops, 'high_bias')']';
      if ~isempty(ns_ops_intervals)
        ns_ops_intervals(:,1)=ns_ops_intervals(:,1)-4;
        ns_ops_intervals(:,2)=ns_ops_intervals(:,2)+4;
        irf_log('proc', 'blanking NS_OPS')
        data = caa_rm_blankt(data,ns_ops_intervals);
      end
      clear ns_ops ns_ops_intervals
    end

    if isempty(result), result = data;
    else
      if ~isempty(data)
        t = result(:,1);
        tapp = data(:,1);

        if tapp(1) <= t(end)
          if tapp(1) < t(end)
            irf_log('proc',sprintf('Last point in data is %f seconds before first point in appended data',t(end)-tapp(1)))
            data = irf_tlim(data,tapp(1),t(end),1);
          elseif tapp(1) == t(end)
            irf_log('proc','Overwriting last data point with appended data with identical time stamp')
            result(end,:) = [];
          end
          irf_log('proc',sprintf('   Last point in data: %s',epoch2iso(t(end))))
          irf_log('proc',sprintf('   Attempt to append interval %s to %s',epoch2iso(tapp(1)),epoch2iso(tapp(end))))
        end
      end

      if ~isempty(data)
        result = [result; data]; %#ok<AGROW>
      end
    end

  end % for dirs
  irf_log('proc',sprintf('Loaded C%d',ic))
  if isempty(result), continue, end
  %c_eval('Es? = result;',ic)

  %create evenly sampled data and fill gaps with moving median
  tNew = (result(1,1):4:result(end,1))';
  newData = NaN(length(tNew),2);
  [~,idxComm,~] = intersect(tNew,result(:,1));
  if length(idxComm) ~= length(result(:,1)), error('something wrong'), end
  newData(idxComm,:) = result;
  ttt = newData(:,2);
  newData(:,2) = movmedian(newData(:,2),5,'omitnan');
  ii = ~isnan(ttt); newData(ii,2) = ttt(ii); clear ttt
  newData(:,1) = tNew;
  resultSave = result;
  result = newData;

  %% find Ex<threhsold
  EX_MAX = 15; % repeated points above this value
  EX_MAXXX = 30; % all points above this value
  idx = find(abs(result(:,2))>EX_MAX); idxxx = find(abs(result(:,2))>EX_MAXXX);
  didx = diff(idx); didx(didx<5) = 1;
  idx1 = find(didx==1);
  idx1min = idx1-1; idx1min(idx1min<1) = [];
  idx1plu = idx1+1; idx1plu(idx1plu>length(idx)) = [];
  idxBad = unique([idxxx; idx(idx1min); idx(idx1); idx(idx1plu)]);
  irf_plot(h(ic),result,'k.'), hold(h(ic),'on')
  irf_plot(h(ic),[resultSave(:,1), resultSave(:,2)+2],'g.')
  irf_plot(h(ic),result(idxBad,:),'mo')
  ylabel(h(ic),sprintf('Ex C%d [mV/m]',ic))

  if ~isempty(idxBad)
    dbad = diff(idxBad); ii = find(dbad>1);

    a = sort([idxBad(1); idxBad(ii);  idxBad(ii+1); idxBad(end)]);
    a = result(a,1);
    a = reshape(a',2,length(a)/2)';
    iCurr = 2;
    %Merge intervals if no data between them
    while true
      if iCurr > size(a,1), break, end
      if any(~isnan(resultSave(resultSave(:,1)>a(iCurr-1,2) & resultSave(:,1)<a(iCurr,1),2)))
        iCurr = iCurr + 1; % next interval
      else, a(iCurr,1) = a(iCurr-1,1); a(iCurr-1,:) = []; % merge with previous
      end
    end
    for iInt = 1:size(a,1)
      [~,~,idxComm] = intersect(a(iInt,:)',result(:,1));
      irf_plot(h(ic),result(idxComm,:),'kx-');
    end
    a(:,2) = a(:,2) - a(:,1)+4;   a(:,1) = a(:,1) - 2; %#ok<NASGU>
    c_eval('int? = a;',ic), clear a
  end
end % for ic

irf_zoom(h,'x',[iso2epoch(iso_st) iso2epoch(iso_stop)])

cd(old_pwd);
