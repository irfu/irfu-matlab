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
save_list = '';

DP = c_ctl(0,'data_path');
DB = c_ctl(0,'isdat_db');

B_DT = 300;
B_DELTA = 60;

cp = ClusterProc;

delete('mEFWburs*.mat') %Remove old files

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

filtv=zeros(1,4);   % remember V filter usage
for out = 1:varsbsize;
    vt=varsb{out};
    vtlen=length(vt);
    field='E';
    
    if vt(1)=='B'
        probe=vt(3:4);
        sen = irf_ssub('p?',probe);
        filter='bp';
    elseif vt(1)=='S'
        field='dB';
        probe=lower(vt(3));
        sen=probe;
        filter='4kHz';
    elseif vt(1)=='V'
        filt=vt(vtlen);
        switch filt
            case 'U'
                filter='32kHz';
                filtv(1)=filtv(1)+1;
            case 'H'
                filter='4kHz';
                filtv(2)=filtv(2)+1;
            case 'M'
                filter='180Hz';
                filtv(3)=filtv(3)+1;
            case 'L'
                filter='10Hz';
                filtv(4)=filtv(4)+1;
            otherwise
                error(['Unknown filter char for V: ' vt(vtlen)]);
        end
        if vtlen>3
            if vt(2)=='4'  % 43 check
                probe = vt(3:-1:2);
            else
                probe = vt(2:3);
            end
        else
            probe = vt(2);
        end
        sen = irf_ssub('p?',probe);
    end
    instrument = 'efw';
    
    [t,data] = caa_is_get(DB,st-B_DELTA,B_DT,cl_id,instrument,field,...
        sen,filter,'burst','tm');
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
    d_phys=data*0.00212;
    data_phys = [t d_phys];
    
    if (out==1) % Create data matrix for t and all 8 possible variables
        if size(t,1)<3 || size(data,1)<3 % sanity check
            irf_log('proc','No usable burst data');
            return;
        end
        dlen=size(data,1);
        data8=NaN(dlen,9);
        data8(:,1)=t;   % corrected time
    end
    data8(:,out+1)=data;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
    %%%%%    Saving the data to mEFWburstTM/mEFWburstR1    %%%%%%%%%%%%%%%%%%%%%
    %%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save_file = './mEFWburstTM.mat';
    if ~isempty(data)
        data = [t data]; %#ok<NASGU,AGROW>
        if vt(1)=='S'
            eval(irf_ssub(['tm?!b$=data;' 'save_list=[save_list ''tm?!b$ ''];'],filter,cl_id,probe));
        else
            eval(irf_ssub(['tm?!p$=data;' 'save_list=[save_list ''tm?!p$ ''];'],filter,cl_id,probe));
        end
    else
    end
    if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
        irf_log('save',[save_list ' -> ' save_file])
        if exist(save_file,'file')
            eval(['save -append ' save_file ' ' save_list]);
        else
            eval(['save ' save_file ' ' save_list]);
        end
    end
    
    
    if field=='E'
        data = data_phys; %#ok<NASGU>
    elseif field=='B'
        continue;
    elseif strcmp(field,'dB')
        continue;
    else
        disp(['Info: Unknown field ' field]);
    end
    
    save_file = './mEFWburstR1.mat';
    data = data_phys; %#ok<NASGU>
    eval(irf_ssub(['PP?!p$=data;' 'save_list=[save_list ''PP?!p$ ''];'],filter,cl_id,probe));
    
    if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
        irf_log('save',[save_list ' -> ' save_file])
        if exist(save_file,'file')
            eval(['save -append ' save_file ' ' save_list]);
        else
            eval(['save ' save_file ' ' save_list]);
        end
    end
    save_list = '';
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Checking the order of the data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
[ok,pha] = c_load('Atwo?',cl_id);
if ~ok || isempty(pha)
    irf_log('load',...
        irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
end

%test = load('mEFWburstTM.mat');
%load('mEFWburstTM.mat')
%load('mER.mat')
%fn = fieldnames(test);
%bla=char(fn);
%cc=size(bla);
test2 = load('mEFWburstR1.mat');
%load('mEFWburstR1.mat')
load('mER.mat')
fn2 = fieldnames(test2);
bla2=char(fn2);
cc2=size(bla2);
ff=0;

if cc2(1)>=4
    if exist(irf_ssub('wE?p!',cl_id,12),'var')
        %        name11 = irf_ssub('PP?!p$',filter,cl_id,1);
        %        name22 = irf_ssub('PP?!p$',filter,cl_id,2);
        tt1=eval(irf_ssub('wE?p!',cl_id,12));   %(t1(:,2)-t2(:,2))/88;
        probev=1;
    elseif exist(irf_ssub('wE?p!',cl_id,34),'var')
        %        name11 = irf_ssub('PP?!p$',filter,cl_id,3);
        %        name22 = irf_ssub('PP?!p$',filter,cl_id,4);
        tt1=eval(irf_ssub('wE?p!',cl_id,34));   %(t1(:,2)-t2(:,2))/88;
        probev=3;
    else
        tt1=eval(irf_ssub('wE?p!',cl_id,32));   %(t1(:,2)-t2(:,2))/88;
        probev=2;
    end
    
    % ff gg di
    bestguess=[-1 -1 -1];
    for di=1:2:varsbsize
        t1(:,1:2)=data8(:,[1 di+1]);
        t1m=mean(t1(:,2));
        t2(:,1:2)=data8(:,[1 di+2]);
        t2m=mean(t2(:,2));
        %di
        divf=t2m/t1m;
        %    if divf<0.67 || divf>1.5    % different data types skip
        if divf<0.5 || divf>2    % different data types skip
            %divf
            continue;
        end
        aa1=c_phase(tt1(:,1),pha);
        if isempty(aa1)
            continue;
        end
        sp1=c_efw_sfit(12,3,10,20,tt1(:,1),tt1(:,2),aa1(:,1),aa1(:,2),1,'hx'); % org sfit2
        distance=88;
        tt2=(t1(:,2)-t2(:,2))/distance;
        aa2=c_phase(t2(:,1),pha);
        if isempty(aa2)
            continue;
        end
        sp2=c_efw_sfit(12,3,10,20,t2(:,1),tt2,aa2(:,1),aa2(:,2),1,'ib');       % org sfit2
        [a,~] = size(sp2);
        gg=0;
        i=1;
        
        while i<a+1
            [c,]=find(sp1==sp2(i,1));
            ex1=sp1(c,2);
            ex2=sp2(i,2);
            ey1=sp1(c,3);
            ey2=sp2(i,3);
            
            %timevec=fromepoch(sp2(i,1));
            z1=atan2(ey1,ex1);
            z2=atan2(ey2,ex2);
            y=round(abs(((z1-z2)/pi)*180));
            
            if ~isempty(y)
                if isnan(y)==1
                    %                irf_log('proc','NaN');
                elseif  25<y && y<155
                    gg=gg+1;
                elseif  205<y && y<335
                    gg=gg+1;
                else
                    ff=ff+1;
                    %                irf_log('proc','the data match');
                end
            end
            
            i=i+1;
            
        end
        
        % Remember best guess so far
        %    if (ff-gg>bestguess(1)-bestguess(2)) || (ff-gg==bestguess(1)-bestguess(2) && gg<bestguess(2))
        if (ff>bestguess(1) && gg==0)
            bestguess=[ff gg di];
        end
    end
    %bestguess
    xy=1:varsbsize;
    svar=irf_ssub('V?',probev);
    svar1=irf_ssub('V?',probev+1);
    pfound=false;
    % find variable position
    for pos=1:varsbsize-1
        if strcmp(varsb{pos}(1:2),svar) && length(varsb{pos})==3 && strcmp(varsb{pos+1}(1:2),svar1) && length(varsb{pos+1})==3
            pfound=true;
            break;
        end
    end
    data8ord=data8; % assume no order change
    %bestguess(3)=3;
    %pos
    if ~pfound || bestguess(1)==-1
        irf_log('proc','Can not find burst order. standard order 1-n used');
    elseif pos~=bestguess(3)
        % make order vector
        pcnt=pos;
        for j=bestguess(3):varsbsize
            xy(j)=pcnt;
            pcnt=pcnt+1;
            if pcnt>varsbsize
                pcnt=1;
            end
        end
        if bestguess(3)>1
            pcnt=pos-1;
            for j=bestguess(3)-1:-1:1
                if pcnt<1
                    pcnt=varsbsize;
                end
                xy(j)=pcnt;
                pcnt=pcnt-1;
            end
        end
        % order data
        for j=1:varsbsize
            data8ord(:,xy(j)+1)=data8(:,j+1);
        end
    end
else
    data8ord=data8; % assume no order change
end

% save raw ordered data
save_file = './mEFWburstTM1.mat';
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
        data8ordfc(:,i+1) = data8ordfc(:,i+1)*0.00212;
        data8ordfc(:,[1 i+1]) = c_efw_invert_tf(data8ordfc(:,[1 i+1]),...
            varsb{i}(end));
    else
        error(['Unknown ib data type: ' varsb{i}]);
    end
end

data8ordfc(ix,2:end) = NaN; % Set spikes to NaNs

% Save BSC data
if BSCcnt
    if BSCcnt<3, error('Less than 3 BSC components'), end % Sanity check
    BSCtemp = data8ordfc(:,[1 BSCpos+1]);  %#ok<NASGU>
    
    save_file = './mBSCBurst.mat';
    save_list='';
	c_eval('wBSC4kHz?=BSCtemp;save_list=[save_list '' wBSC4kHz? ''];',cl_id);
    if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
        irf_log('save',[save_list ' -> ' save_file])
        if exist(save_file,'file')
            eval(['save -append ' save_file ' ' save_list]);
        else
            eval(['save ' save_file ' ' save_list]);
        end
    end
end

save_file = './mEFWburstR.mat';
save_list='';
for i=1:varsbsize
    if varsb{i}~='V'
        continue
    end
    if length(varsb{i})>3
        if vt(2)=='4'  % 43 check
            probe = vt(3:-1:2);
        else
            probe = vt(2:3);
        end

    else
        probe=varsb{i}(2);
    end
    filter=get_filter(varsb{i});
    data=data8ordfc(:,[1 i+1]); %#ok<NASGU>
    eval(irf_ssub(['P?!p$=data;' 'save_list=[save_list ''P?!p$ ''];'],filter,cl_id,probe));
end
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
delete('mEFWburstR1.mat'); %Remove some files
delete('mEFWburstTM.mat');
cd(old_pwd);

end

function filt = get_filter(ibstr)
% get filter from ib string

    if ibstr(1)=='B'
        filt='bp';
    elseif ibstr(1)=='S'
            filt='4kHz';
    elseif ibstr(1)=='V'
        switch ibstr(end)
            case 'U'
                filt='32kHz';
            case 'H'
                filt='4kHz';
            case 'M'
                filt='180Hz';
            case 'L'
                filt='10Hz';
            otherwise
                error(['Unknown filter char for V: ' ibstr(end)]);
        end
    end
end