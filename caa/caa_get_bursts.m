function ret = caa_get_bursts(filename, plot_flag)
%caa_get_bursts(filename) produce and plot Cluster burst data from the raw data
%
% data = caa_get_bursts(filename, plotflag)
%
% Input:
%   filename   : burst data file name
%   burst_plot : generate plot 0=off 1=on (default off)
%
% Burst files are read from directory /data/cluster/burst/
% Plot files are saved in directory $HOME/figures/
% 
% Example: 
%       caa_get_bursts('070831101839we.03')
%
% $Id$

error(nargchk(1,2,nargin));
if nargin < 2
    plot_flag = 0;
end
old_pwd = pwd;

ret=1;
plot_save=1;
plotpath=[getenv('HOME') '/figures/'];
flag_save = 1;

DP = c_ctl(0,'data_path');
DB = c_ctl(0,'isdat_db');
B_DT = 300;
B_DELTA = 60;

cp = ClusterProc;

%Remove old files
fn={'mEFWburstTM.mat' 'mEFWburstR.mat' 'mEFWburst.mat'};
for i=1:size(fn,2)
    if exist(fn{i},'file')
        delete(fn{i});
    end
end

cl_id=str2double(filename(end)); %Get the satellite number
fname=irf_ssub([plotpath 'p?-c!'],filename(1:12),cl_id); %Sets the name that will be used to save the plots
fnshort=filename;
s=filename;
full_time = iso2epoch(['20' s(1:2) '-' s(3:4) '-' s(5:6) 'T' s(7:8) ':' s(9:10) ':' s(11:12) 'Z']);
start_time=full_time;
st=full_time;

dirs = caa_get_subdirs(st, 90, cl_id);
if isempty(dirs)
    irf_log('proc',['Can not find L1 data dir for ' s]);
    return;
end
found=false;
for i=size(dirs,2):-1:1 % find start time directory
    d=dirs{i}(end-12:end);
    dtime=iso2epoch([d(1:4) '-' d(5:6) '-' d(7:8) 'T' d(10:11) ':' d(12:13) ':00Z']);
    if dtime<=start_time
        found=true;
        break;
    end
end
if ~found
    irf_log('proc','iburst start time does not match any L1 data dir');
    return;
end
cd(dirs{i})

varsb = c_efw_burst_param([DP '/burst/' filename]);
varsbsize = length(varsb);

for out = 1:varsbsize;
    [field,sen,filter] = get_ib_props(varsb{out});
    instrument = 'efw';
    
    [t,data] = caa_is_get(DB,st-B_DELTA,B_DT,cl_id,instrument,field,...
        sen,filter,'burst','tm');
    
    if (out==1) % Create data matrix for t and all 8 possible variables
        if size(t,1)<3 || size(data,1)<3 % sanity check
            irf_log('proc','No usable burst data');
            return;
        end
        data8 = NaN(size(data,1),9);
        
        start_satt = c_efw_burst_chkt(DB,[DP '/burst/' filename]);
        if isempty(start_satt)
            irf_log('dsrc','burst start time was not corrected')
        elseif isempty(t)
            irf_log('proc','t is empty. no iburst data?!');
            cd(old_pwd);
            return;
        else
            err_t = t(1) - start_satt;
            irf_log('dsrc',['burst start time was corrected by ' ...
                num2str(err_t) ' sec'])
            t = t - err_t;
        end
        data8(:,1)=t;   % corrected time
    end
    data8(:,out+1)=data;
end
clear data
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Checking the order of the data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if varsbsize>=2
    can_check_order = 1;
    [ok,pha] = c_load('Atwo?',cl_id);
    if ~ok || isempty(pha)
        irf_log('load',...
            irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
        can_check_order = 0;
    end
    
    if can_check_order
        % See which single-ended probes we got in the burst, e.g. V1H, V2H...
        probes = [];
        for i=1:length(varsb)
            if varsb{i}(1)=='V' && ischar(varsb{i}(3))
                probes = [probes str2double(varsb{i}(2))]; %#ok<AGROW>
            end
        end
        if isempty(probes)
            can_check_order = 0;
            irf_log('proc','Cannot check burst order: no single-ended probes')
        end
    end
    
    if can_check_order
        % See which NM data we may use
        probep = [];
        if any(probes==3) && any(probes==4), probep = 34; end
        if any(probes==1) && any(probes==2), probep = [probep 12]; end
        if any(probes==3) && any(probes==2), probep = [probep 32]; end
        if isempty(probep)
            can_check_order = 0;
            irf_log('proc','Cannot check burst order: no useful probe pairs')
        end
    end
    
    if can_check_order
        % Load NM data
        for i=1:length(probep)
            [ok,nmdata] = c_load(irf_ssub('wE?p!',cl_id,probep(i)));
            if ok && ~isempty(nmdata)
                ref_probep = probep(i);
                irf_log('proc',sprintf('Using p%d NM data as reference',probep(i)))
                break
            end
        end
        if isempty(nmdata)
            can_check_order = 0;
            irf_log('proc','Cannot check burst order: no reference NM data')
        end
    end
    
    if can_check_order
        data8ord = c_efw_burst_order_signals(data8, varsb, nmdata, pha, ref_probep);
    else
        data8ord = data8;
    end
else
    data8ord=data8; % assume no order change
end

% save raw ordered data
save_file = './mEFWburstTM.mat';
save_list='';
burst_info=sprintf('%s ',varsb{:}); %#ok<NASGU>
eval(irf_ssub(['ib?_info=burst_info(1:end-1);' 'save_list=[save_list ''ib?_info ''];'],cl_id));

eval(irf_ssub(['iburst?=data8ord;' 'save_list=[save_list ''iburst? ''];'],cl_id));
if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
    irf_log('save',[save_list ' -> ' save_file])
    if exist(save_file,'file')
        eval(['save -append ' save_file ' ' save_list]);
    else
        eval(['save ' save_file ' ' save_list]);
    end
end
% filter data
% calib data

[data8ordfc, ix] = caa_identify_ib_spikes(data8ord);

BSCpos=zeros(1,3);
BSCcnt=0;
for i=1:varsbsize
    if varsb{i}(1)=='B' % BP data
        % BP12 factor is unkrown. These data do not have an easy physical
        % meaning.
        %data8ordfc(:,i+1)=data8ordfc(:,i+1)*BPFACTOR;
    elseif varsb{i}(1)=='S' % SC data       
        xyzord=double(varsb{i}(3))-87; % X=1
        BSCpos(xyzord)=i;
        data8ordfc(:,i+1) = -data8ordfc(:,i+1)/7000;
        c_efw_burst_bsc_tf(data8ordfc(:,[1 BSCpos+1]),cl_id,xyzord);
        BSCcnt=BSCcnt+1;
    elseif varsb{i}(1)=='V'
        if length(varsb{i})>3 % Differential signals are special
            data8ordfc(:,i+1) = data8ordfc(:,i+1)*0.00212;
            data8ordfc(:,[1 i+1]) = c_efw_invert_tf(data8ordfc(:,[1 i+1]),'BP');
        else
            data8ordfc(:,i+1) = data8ordfc(:,i+1)*0.00212;
            data8ordfc(:,[1 i+1]) = c_efw_invert_tf(data8ordfc(:,[1 i+1]),...
                varsb{i}(end));
        end
    else
        error(['Unknown ib data type: ' varsb{i}]);
    end
end
data8ordfc(ix,2:end) = NaN; % Set spikes to NaNs

% Save data
save_file = './mEFWburstR.mat';
save_list='';
if BSCcnt % BSC
    if BSCcnt<3, error('Less than 3 BSC components'), end % Sanity check
    BSCtemp = data8ordfc(:,[1 BSCpos+1]);  %#ok<NASGU>
	c_eval('wBSC4kHz?=BSCtemp;save_list=[save_list '' wBSC4kHz? ''];',cl_id);
    clear BSCtemp
end

for i=1:varsbsize % Single ended
    if varsb{i}~='V'
        continue
    end
    [~,~,filter,probe] = get_ib_props(varsb{i});
    data=data8ordfc(:,[1 i+1]); %#ok<NASGU>
    if length(probe)>1
        eval(irf_ssub(['wE?!p$=data;' 'save_list=[save_list ''wE?!p$ ''];'],filter,cl_id,probe));
    else
        eval(irf_ssub(['P?!p$=data;' 'save_list=[save_list ''P?!p$ ''];'],filter,cl_id,probe));
    end
end

% Make spike problem time vector
mem=1;
ii=1;
tinter=[];
sz=length(ix);
if sz>1
    for i=2:sz
        if ix(i)-ix(i-1)>1
            tinter(ii,1)=data8ordfc(ix(mem),1); %#ok<AGROW>
            tinter(ii,2)=data8ordfc(ix(i-1),1); %#ok<AGROW>
            ii=ii+1;
            mem=i;
        end
    end
    tinter(ii,1)=data8ordfc(ix(mem),1);
    tinter(ii,2)=data8ordfc(ix(sz),1);
end
if ~isempty(tinter)
    eval(irf_ssub(['SPIKE?=tinter;' 'save_list=[save_list ''SPIKE? ''];'],cl_id));
end
clear tinter ii mem sz

if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
    irf_log('save',[save_list ' -> ' save_file])
    
    if exist(save_file,'file')
        eval(['save -append ' save_file ' ' save_list]);
    else
        eval(['save ' save_file ' ' save_list]);
    end
end

% make L2
getData(cp,cl_id,'whip');
getData(cp,cl_id,'pburst');
getData(cp,cl_id,'dieburst');
getData(cp,cl_id,'dibscburst');

if plot_flag  
    st=data8ordfc(1,1);
    sp=data8ordfc(end,1);
    if st==-Inf || isnan(st)
        return;
    end
    getData(ClusterDB(DB,DP,'.'),st-B_DELTA,B_DT,cl_id,'bfgm');
    
    clf;
    dt2=5;
    st_int=st-dt2;
    st_int2=sp-st+2*dt2;
    
    summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',st_int2,'vars',char([fnshort varsb]));
    
    if plot_save
        orient landscape
        print('-dpdf', fname);
    end
end

ret=0;
cd(old_pwd);
end

function [field,sen,filter,probe] = get_ib_props(ibvarstr)
% get filter from ib string
field='E';
switch ibvarstr(1)
    case 'B'
        probe=ibvarstr(3:4);
        sen = irf_ssub('p?',probe);
        filter='bp';
    case 'S'
        field='dB';
        probe=lower(ibvarstr(3));
        sen=probe;
        filter='4kHz';
    case 'V'
        if length(ibvarstr) > 3
            if ibvarstr(2)=='4'  % 43 check
                probe = ibvarstr(3:-1:2);
            else
                probe = ibvarstr(2:3);
            end
        else
            probe = ibvarstr(2);
        end
        sen = irf_ssub('p?',probe);
        switch ibvarstr(end)
            case 'U'
                filter='32kHz';
            case 'H'
                if length(probe)==1, filter='4kHz';
                else filter='8kHz';
                end
            case 'M'
                filter='180Hz';
            case 'L'
                filter='10Hz';
            otherwise
                error(['Unknown filter char for V: ' ibvarstr(vtlen)]);
        end
    otherwise
        error(['Unknown leading char in V: ' ibvarstr(vtlen)]);
end
end