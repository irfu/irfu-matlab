function output = read_LL(varStr,Tint)
% output = solo.read_solo_LL(varStr, Tint)
% 
% Read Solar Orbiter low latency data
% Warning: this data is not to be used for science
%
% varStr is one of the following 
% B_RTN B_SRF V_RTN V_SRF N
%
% Example
%   Tint = irf.tint('2020-07-29T10:20:00.000Z/2020-07-29T11:20:00.000Z');
%   V = solo.get_data('V_RTN',Tint) % Solar wind speed in RTN low latency
% 
% Script can be improved/optimised in the future

if ~isa(Tint,'GenericTimeArray')
    error('TINT must be of GenericTimeArray type');
elseif Tint.stop-Tint.start<=0
    error('TINT duration is zero or negative');
end

day_load = irf_time(Tint(1),'epochtt>date'):irf_time(Tint(2),'epochtt>date');

switch varStr
    case 'V_RTN', SWA_VAR  = 'SWA_PAS_VELOCITY_RTN'; conv=0;
    case 'V_SRF', SWA_VAR  = 'SWA_PAS_VELOCITY_RTN'; conv=1;
    case 'N', SWA_VAR  = 'SWA_PAS_DENSITY'; conv=nan;
    case 'B_RTN', coord = 'RTN'; conv=nan;
    case 'B_SRF', coord = 'SRF'; conv=nan;
end

% There are both I and C versions of the fiels that are conplete and
% incomplete. This script will load them both and combine them.
for k=1:length(day_load)
    switch varStr(1)
        case 'B', dir_file = ['/Volumes/solo/soar/mag/LL02/' datestr(day_load(k),'YYYY')  '/' datestr(day_load(k),'mm')  '/' datestr(day_load(k),'dd')];
            if exist(dir_file,'dir')
                cd(dir_file)
                file = dir('**/*.cdf'); % move to directory and search for file
                for kk=1:length(file)
                    if exist(file(kk).name,'file')
                        dat_time = spdfcdfread(file(kk).name, 'Variable', {'EPOCH'}); % read time
                        [dat, ~] = spdfcdfread(file(kk).name,'variable',varStr); % read the data
                        flag_rem = dat(:,1) == -1e31; % indices of flags
                        dat(flag_rem,:) = nan; % remove the flags
                        if kk==1
                            out = TSeries(irf_time(dat_time,'date>epochtt'),dat,'TensorOrder',1,'TensorBasis','xyz',...
                                'repres',{'x','y','z'});
                        else
                            out_ = TSeries(irf_time(dat_time,'date>epochtt'),dat,'TensorOrder',1,'TensorBasis','xyz',...
                                'repres',{'x','y','z'});
                            out = out.combine(out_);
                        end
                        out.coordinateSystem = coord;
                        out.name = varStr;
                        out.units = 'nT';
                        out.siConversion = '1.0E-9>T';
                    end
                end
            end

        case 'V', dir_file = ['/Volumes/solo/soar/swa/LL02/' datestr(day_load(k),'YYYY')  '/' datestr(day_load(k),'mm')  '/' datestr(day_load(k),'dd')];
            if exist(dir_file,'dir')
                cd(dir_file)
                file = dir('**/*.cdf'); % move to directory and search for file
                for kk=1:length(file) % when there is both I and C versions, read them all and take all data available
                    if exist(file(kk).name,'file')
                        dat_time = spdfcdfread(file(kk).name, 'Variable', {'EPOCH'});
                        [dat, ~] = spdfcdfread(file(kk).name,'variable',SWA_VAR);
                        flag_rem = dat(:,1) == -1e31;
                        dat(flag_rem,:) = nan;
                        if kk==1
                            out = TSeries(irf_time(dat_time,'date>epochtt'),dat,'TensorOrder',1,'TensorBasis','xyz',...
                                'repres',{'x','y','z'});
                        else
                            out_ = TSeries(irf_time(dat_time,'date>epochtt'),dat,'TensorOrder',1,'TensorBasis','xyz',...
                                'repres',{'x','y','z'});
                            out = out.combine(out_);
                        end
                        out.units = 'km s^-1';
                        out.siConversion = '1000.0 > m s^-1';
                    end
                end
            end

        case 'N', dir_file = ['/Volumes/solo/soar/swa/LL02/' datestr(day_load(k),'YYYY')  '/' datestr(day_load(k),'mm')  '/' datestr(day_load(k),'dd')];
            if exist(dir_file,'dir')
                cd(dir_file)
                file = dir('**/*.cdf'); % move to directory and search for file
                for kk=1:length(file) % when there is both I and C versions, read them all and take all data available
                    if exist(file(kk).name,'file')
                        dat_time = spdfcdfread(file(kk).name, 'Variable', {'EPOCH'});
                        [dat, ~] = spdfcdfread(file(kk).name,'variable',SWA_VAR);
                        flag_rem = dat(:,1) == -1e31;
                        dat(flag_rem,:) = nan;
                        if kk==1
                            out = TSeries(irf_time(dat_time,'date>epochtt'),dat);
                        else
                            out_ = TSeries(irf_time(dat_time,'date>epochtt'),dat);
                            out = out.combine(out_);
                        end
                        out.name = 'SWA_PAS_DENSITY';
                        out.units = 'cm^-3';
                        out.siConversion = '1.0E-6 > m^-3';
                    end
                end
            end
    end

    if k==1 || ~exist("output",'var')
        if exist('out','var');output = out; end
        clear out out_
    else
        if exist('out','var');output = output.combine(out);end
        clear out out_
    end

end

if exist('output','var')~=0
    output = output.tlim(Tint);
else
    output = [];
end

if conv==1
    output = solo.srf2rtn(output,-1);
    output.coordinateSystem = 'SRF';
    output.name = 'SWA_PAS_VELOCITY_SRF';
end

if conv==0
    output = solo.srf2rtn(output,-1);
    output.coordinateSystem = 'RTN';
    output.name = 'SWA_PAS_VELOCITY_RTN';
end

output.userData = 'Low Latency data, not for science!';