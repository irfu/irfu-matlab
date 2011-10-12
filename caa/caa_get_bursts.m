function ret = caa_get_bursts(filename, burst_plot)
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

error(nargchk(1,2,nargin));
if nargin < 2
    burst_plot = 0;
end
old_pwd = pwd;

ret=1;

% WARNING: Do NOT use 1 for this on server L1 data
fetch_efw_data=0; % must be 1 if no efw data in current directory

plot_save=1;
plotpath=[getenv('HOME') '/figures/'];

flag_save = 1;
save_list = '';
DBNO='db:9';
dt=30; %Sets the time interval for the data retrieval
no_data=0;
B_DT = 300;
B_DELTA = 60;
cp = ClusterProc;

burstreadbytes=44;

    fns=size(filename,2);
    DP = c_ctl(0,'data_path');
    fid=fopen([DP '/burst/' filename],'rb'); % opens the binary file for reading
    if fid==-1
        error(['Can not find burst file ' filename ' in ' DP '/burst']);
    end
        
    fseek(fid,128,0); % skip first 128 bytes
    % run on intel little endian read
    data(fns+1:fns+burstreadbytes,1) = fread(fid, burstreadbytes, 'uint8=>uint8', 'ieee-be'); % start at 18 due to filname byte length (17)
    fclose(fid);

    delete('mEFWburs*.mat') %Remove old files
%    delete('*.mat') %Removes old files. THIS IS DANGEROUS!!! FIX ME!

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
    cd(dirs{i});

    date = fromepoch(full_time);
    sp = [pwd() '/' irf_fname(full_time)];  
    cdb=ClusterDB(DBNO,[ DP '/burst'],'.');
    %Sets the variables for gathering normal mode data   
    vars0 = {'tmode','fdm','efwt','ibias','p','e','a','sax','r','v','bfgm','bsc'};
	%Sets the variables thats needed for burst data.
    vars11 = {'whip','sweep','bdump','probesa','p','ps' 'dies','die','pburst','dieburst','dibsc','dibscburst'};
	%Scans the correct line from the text file to be used later.
%	data=sscanf(tline,'%s %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x');%,fname,&bcode,&bfreq,&btrig,&bchirp,&bpages,&bthrsh0,&bp0,&bp1,&bp2,&bp3,&st[0],&st[1],&st[2],&st[3],&st[4],&vt[0], &vt[1], &vt[2], &vt[3], &vt[4],&et[0],&et[1],&et[2],&et[3],&et[4],&sa[0],&sa[1],&sa[2],&ea[0],&ea[1],&ea[2],&lr[0],&lr[1],&lr[2],&spare1,&spare2,&f[0],&f[1],&f[2],&f[3],&f[4],&f[5],&f[6],&f[7]))
    filename=irf_ssub([ DP '/burst/?'],filename)
    dec2bin(data(19),8);
    OUT = c_efw_burst_geth(filename);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting frequency and number of parameters    %%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %getting the information from the data and makes the corresponding number smaller.
    fff = bitand(data(19),7);
    ss = bitand(data(19),48);
    SS = bitshift(ss,-4);

    switch SS
        case 0,
            switch fff
                case 0, 
                    output=[450  8];
                case 1,
                    output=[900 8];
                case 2, 
                    output=[2250 8];
                case 3,
                    output=[4500 8];
                case 4, 
                    output=[9000 4];
                case 7,
                    output=[25000 2];
                otherwise
                    output=[0];
            end
            
        case 1
            if fff==5
                output=[18000 2];
            else
                error('bad fff');
            end    
            
        case 2
            switch fff
                case 4
                    output=[9000 8];
                case 5
                    output=[18000 2];
                otherwise
                    output=[0];
            end
            
        case 3
            switch fff
                case 0
                    output=[450 16];
                case 1
                    output=[900 16];
                case 2
                    output=[2250 16];
                case 3
                    output=[4500 16];
                case 4
                    output=[9000 8];
                case 5
                    output=[18000 4];
                case 6
                    output=[36000 2];
                otherwise
                    output=[0];
            end
    end
    
    i=1;
    iii=2;
    ii=54;
%    vars=zeros(output(1,2),4);
    varsb={};
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting information about the burst data    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    madc0={ 'V1L' 'V1M' 'V1H' 'V1U' 'V3L' 'V3M' 'V3H' 'V3U'...
            'V12M' 'V12M' 'SCX' 'SCZ' 'BAD' };
    madc1={ 'V2L' 'V2M' 'V2H' 'V2U' 'V4L' 'V4M' 'V4H' 'V4U'...
            'V43M' 'V12H' 'SCY' 'BP12' 'BAD' };

    while length(data)>=ii && data(ii)~=63
        
        adc0 = bitand(data(ii),15);
        adc1 = bitand(data(ii),240);
        adc1 = bitshift(adc1,-4);
        adc1 = bitor(bitand(adc0,8),bitand(adc1,7));
        
        if adc0>11 || adc0<0
            adc0=11;
        end
        varsb{i} = madc0{adc0+1};

        if adc1>11 || adc1<0
            adc1=11;
        end
        varsb{iii} = madc1{adc1+1};
        
        i=i+2;
        ii=ii+1;
        iii=iii+2;

    end
%    x=char(vars(2,2));
%    vars1 = char(vars);
    
%    [s,l] = size(vars);
    varsbsize=size(varsb,2);
    if varsbsize==2
        varsbsize
        pause;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting the burst data from ISDAT database    %%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    probe_list=1:4;
    do_burst=1;
%%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
	switch cl_id
		case 1
			if start_time>toepoch([2009 10 14 07 00 00]) || ...
					(start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]))
				% p1 and p4 failure
				probe_list = 2:3;
				irf_log('dsrc',sprintf('p1 and p4 are BAD on sc%d',cl_id))
			elseif start_time>toepoch([2001 12 28 03 00 00]) 
				% p1 failure
				probe_list = 2:4;
				irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
			elseif ( (start_time>=toepoch([2001 04 12 03 00 00]) && start_time<toepoch([2001 04 12 06 00 00])) || ...
					(  start_time>=toepoch([2001 04 14 06 00 00]) && start_time<toepoch([2001 04 16 15 00 00])) || ...
					(  start_time>=toepoch([2001 04 18 03 00 00]) && start_time<toepoch([2001 04 20 09 00 00])) || ...
					(  start_time>=toepoch([2001 04 21 21 00 00]) && start_time<toepoch([2001 04 22 03 00 00])) || ...
					(  start_time>=toepoch([2001 04 23 09 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) )
				% The bias current is a bit too large
				% on p3 and p4 on C1&2 in April 2001.
				% Ignore p3, p4 and p34 and only use p1, p2 and p12.
				% Use only complete 3-hour intervals to keep it simple.
				probe_list = [1 2];
				irf_log('dsrc',sprintf('Too high bias current on p3&p4 sc%d',cl_id));
			end
		case 2
			if start_time>=toepoch([2007 06 01 17 20 00])
				% We use 180 Hz filter
				if ~do_burst, param={'180Hz'}; end
				irf_log('dsrc',sprintf('using 180Hz filter on sc%d',cl_id))
				probe_list = [2 4];
				irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
			elseif start_time>=toepoch([2007 05 13 03 23 48])
				probe_list = [2 4];
				irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
			elseif start_time+dt>toepoch([2001 07 23 13 54 18]) && ~do_burst
				% 10Hz filter problem on C2 p3
				% Any changes should also go to ClusterProc/getData/probesa
				probe_list = [1 2 4];
				irf_log('dsrc',sprintf('10Hz filter problem on p3 sc%d',cl_id))
			elseif ( (start_time>=toepoch([2001 04 09 21 00 00]) && start_time<toepoch([2001 04 10 06 00 00])) || ...
					(  start_time>=toepoch([2001 04 10 09 00 00]) && start_time<toepoch([2001 04 19 15 00 00])) || ...
					(  start_time>=toepoch([2001 04 20 03 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) || ...
					(  start_time>=toepoch([2001 04 24 00 00 00]) && start_time<toepoch([2001 04 24 15 00 00])) )
				% The bias current is a bit too large
				% on p3 and p4 on C1&2 in April 2001.
				% Ignore p3, p4 and p34 and only use p1, p2 and p12.
				% Use only complete 3-hour intervals to keep it simple.
				probe_list = [1 2];
				irf_log('dsrc',sprintf('Too high bias current on p3&p4 sc%d',cl_id));
			end
		case 3
			if start_time>toepoch([2002 07 29 09 06 59])
				% p1 failure
				probe_list = 2:4;
				irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
			end
    end
	pl = [12,34];
	switch cl_id
		case 1
			if start_time>toepoch([2009 10 14 07 00 00]) ||  ...
					(start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]))
				pl = 32;
				irf_log('dsrc',sprintf('  !Only p32 exists on sc%d',cl_id));
			elseif (start_time>toepoch([2003 9 29 00 27 0]) || ...
					(start_time>toepoch([2003 3 27 03 50 0]) && start_time<toepoch([2003 3 28 04 55 0])) ||...
					(start_time>toepoch([2003 4 08 01 25 0]) && start_time<toepoch([2003 4 09 02 25 0])) ||...
					(start_time>toepoch([2003 5 25 15 25 0]) && start_time<toepoch([2003 6 08 22 10 0])) )
				pl = [32, 34];
				irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
			elseif start_time>toepoch([2001 12 28 03 00 00])
				pl = 34;
				irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
			elseif  (start_time>=toepoch([2001 04 12 03 00 00]) && start_time<toepoch([2001 04 12 06 00 00])) || ...
					(  start_time>=toepoch([2001 04 14 06 00 00]) && start_time<toepoch([2001 04 16 15 00 00])) || ...
				    (  start_time>=toepoch([2001 04 18 03 00 00]) && start_time<toepoch([2001 04 20 09 00 00])) || ...
					(  start_time>=toepoch([2001 04 21 21 00 00]) && start_time<toepoch([2001 04 22 03 00 00])) || ...
					(  start_time>=toepoch([2001 04 23 09 00 00]) && start_time<toepoch([2001 04 23 15 00 00]))
				% The bias current is a bit too large
				% on p3 and p4 on C1&2 in April 2001.
				% Ignore p3, p4 and p34 and only use p1, p2 and p12.
				% Use only complete 3-hour intervals to keep it simple.
				pl = 12;
				irf_log('dsrc',sprintf('  !Too high bias current on p34 for sc%d',cl_id));
			end
		case 2
			if start_time>toepoch([2007 11 24 15 40 0])
				pl = [32, 34];
				irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
			elseif start_time>toepoch([2007 05 13 03 23 48])
				pl = 34;
				irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
			elseif (start_time>=toepoch([2001 04 09 21 00 00]) && start_time<toepoch([2001 04 10 06 00 00])) || ...
					(  start_time>=toepoch([2001 04 10 09 00 00]) && start_time<toepoch([2001 04 19 15 00 00])) || ...
					(  start_time>=toepoch([2001 04 20 03 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) || ...
					(  start_time>=toepoch([2001 04 24 00 00 00]) && start_time<toepoch([2001 04 24 15 00 00]))
				pl = 12;
				irf_log('dsrc',sprintf('  !Too high bias current on p34 for sc%d',cl_id));
			end
		case 3
			if start_time>toepoch([2003 9 29 00 27 0]) || ...
					(start_time>toepoch([2003 3 27 03 50 0]) && start_time<toepoch([2003 3 28 04 55 0])) ||...
					(start_time>toepoch([2003 4 08 01 25 0]) && start_time<toepoch([2003 4 09 02 25 0])) ||...
					(start_time>toepoch([2003 5 25 15 25 0]) && start_time<toepoch([2003 6 08 22 10 0])) 
				pl = [32, 34];
				irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
			elseif start_time>toepoch([2002 07 29 09 06 59])
				pl = 34;
				irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
			end
	end
%%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%
probes=[];
%pl
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

        [t,data] = caa_is_get(DBNO,st-B_DELTA,B_DT,cl_id,instrument,field,sen,filter,'burst','tm');
        start_satt = c_efw_burst_chkt(DBNO,filename);
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
            data = [t data];
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
        % prepare the output
        if nargout > 0      
            if ~isempty(save_list)    
                sl = tokenize(save_list);  
                out_data = {sl};    
                for i=1:length(sl)
                    eval(['out_data{i+1}=' sl{i} ';'])
                end   
            end  
        else    
            clear out_data 
        end
    save_list='';

    if field=='E'
            data = data_phys;
    elseif field=='B'
         continue;
    elseif strcmp(field,'dB')
         continue;
    else
            disp(['Info: Unknown field ' field]);
    end     
        
        save_file = './mEFWburstR1.mat';
        data = data_phys;
        eval(irf_ssub(['PP?!p$=data;' 'save_list=[save_list ''PP?!p$ ''];'],filter,cl_id,probe));

        if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
            irf_log('save',[save_list ' -> ' save_file])
            if exist(save_file,'file')
                eval(['save -append ' save_file ' ' save_list]);
            else
                eval(['save ' save_file ' ' save_list]);
            end
        end
        
        % prepare the output
         
        if nargout > 0 
            if ~isempty(save_list)
                sl = tokenize(save_list);
                out_data = {sl};
                for i=1:length(sl)
                    eval(['out_data{i+1}=' sl{i} ';'])
                end
            end
        else
            clear out_data
        end
        save_list='';
    end
    
data8(1:10,2:end)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting normal data from ISDAT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    st=data(1,1);
    sp=data(end,1);
    if st==-Inf || isnan(st)
        return;
    end
    if fetch_efw_data % must be 1 if no efw data in directory
        for v=1:length(vars0)
            data2 = getData(cdb,st-B_DELTA,B_DT,cl_id,vars0{v});
            if isempty(data2) && (strcmp(vars0{v},'tmode') || strcmp(vars0{v},'fdm'))
                irf_log('load','No EFW data')
                no_data = 1;
                break
            end
        end
    else
       data2 = getData(cdb,st-B_DELTA,B_DT,cl_id,'bfgm');
       data2 = getData(cdb,st-B_DELTA,B_DT,cl_id,'bsc');
    end
    clear data2;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Checking the order of the data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
[ok,pha] = c_load('Atwo?',cl_id);

test = load('mEFWburstTM.mat');
%load('mEFWburstTM.mat')
%load('mER.mat')
fn = fieldnames(test);
bla=char(fn);
cc=size(bla);
test2 = load('mEFWburstR1.mat');
%load('mEFWburstR1.mat')
load('mER.mat')
fn2 = fieldnames(test2);
bla2=char(fn2);
cc2=size(bla2);
ff=0;

if cc2(1)>=4
    if exist(irf_ssub('wE?p!',cl_id,12))~=0
%        name11 = irf_ssub('PP?!p$',filter,cl_id,1);
%        name22 = irf_ssub('PP?!p$',filter,cl_id,2);
        tt1=eval(irf_ssub('wE?p!',cl_id,12));   %(t1(:,2)-t2(:,2))/88;
        probev=1;
    elseif exist(irf_ssub('wE?p!',cl_id,34))~=0
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
di
    divf=t2m/t1m;
%    if divf<0.67 || divf>1.5    % different data types skip
    if divf<0.5 || divf>2    % different data types skip
divf
        continue;
    end
    aa1=c_phase(tt1(:,1),pha);
    if isempty(aa1)
        continue;
    end
    sp1=c_efw_sfit(12,3,10,20,tt1(:,1),tt1(:,2),aa1(:,1),aa1(:,2),1,output(1,1)); % org sfit2
    distance=88;
    tt2=(t1(:,2)-t2(:,2))/distance;
    aa2=c_phase(t2(:,1),pha);
    if isempty(aa2)
        continue;
    end
    sp2=c_efw_sfit(12,3,10,20,t2(:,1),tt2,aa2(:,1),aa2(:,2),1,output(1,1));       % org sfit2
    [a,b] = size(sp2);
    gg=0;                
    i=1;

    while i<a+1
        [c,d]=find(sp1==sp2(i,1));
        ex1=sp1(c,2);
        ex2=sp2(i,2);
        ey1=sp1(c,3);
        ey2=sp2(i,3);
                                
        timevec=fromepoch(sp2(i,1));
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
%t1(1:5,2)
%t2(1:5,2)
ff
gg
% Remember best guess so far
%    if (ff-gg>bestguess(1)-bestguess(2)) || (ff-gg==bestguess(1)-bestguess(2) && gg<bestguess(2))
    if (ff>bestguess(1) && gg==0)
        bestguess=[ff gg di];
    end
end
bestguess
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
%size(data8)
data8ord=data8; % assume no order change
%bestguess(3)=3;
pos
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
%    data8(1:5,2:end)
end

%data8ord(1:5,2:end)
xy
probe_list

    probe_l=probe_list;
                    
else
    data8ord=data8; % assume no order change
    probe_l=1:length(pl);
    if length(pl)==2
        xy=[1 2];
    else
        if pl==34
            xy=[1];
        else
            xy=[2];
        end
    end
end
save_time=[];

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
        % TODO: BP factor???
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
%    getData(cp,cl_id,'p');%
%    getData(cp,cl_id,'die');%
getData(cp,cl_id,'pburst');
getData(cp,cl_id,'dieburst');
getData(cp,cl_id,'dibscburst');

if burst_plot
    clf;
    %st=st+300
    %sp=sp+300
    dt2=5;
    st_int=st-dt2;
    st_int2=sp-st+2*dt2;
    %summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',st_int2,'vars',vars1);
    %        summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',st_int2,'vars',char(varsb));
    summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',st_int2,'vars',char([fnshort varsb]));
    %summaryPlot(cp,cl_id,'ib','st',st_int,'dt',st_int2);
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