function [t,res]=irf_isdat_get(datasource,start_time,dt_interval)
% getting data from isdat (using ISDAT.jar)
% FUNCTIONS STILL IN DEVELOPMENT, do not use
% [t,res]=irf_isdat_get('Cluster/1/ephemeris/position','2002-03-04T10:00:00.000Z',60)
% irf_stdt to get st,dt

%% Fix path
if ~any(cellfun(@(x) ~isempty(x), regexpi(javaclasspath,'ISDAT.jar')))
  wdir = fileparts(which('isGetDataLite.m'));
  javaclasspath([wdir, filesep, 'ISDAT.jar']);
end

%% Get data from isdat
%Open connection
isdatConnection = ISDAT.NetConnection('db.irfu.se',0);
dataSrcName=ISDAT.DataSourceName(datasource);
%dataSrcName=ISDAT.DataSourceName('Cluster/1/ephemeris/position');
%dataSrcName=ISDAT.DataSourceName('Cluster/1/efw/E/p4/10Hz/lx');

isdatDbSession=ISDAT.DbSession(isdatConnection);
dataSrc=isdatDbSession.SourceName2ID(dataSrcName);
%isdatDbSession.nameThisSourceFully(dataSrc)

st = ISDAT.Epoch(start_time);
%st = ISDAT.Epoch('2002-03-04T10:00:00.000Z');
dt = ISDAT.Epoch(dt_interval);
%dt = ISDAT.Epoch(2.0);

dataReq = ISDAT.DbDataRequest(dataSrc,st,dt);
                        
sciData = isdatDbSession.getSciData(dataReq);
%disp(sprintf('totNrSamples : %d',sciData.totNrSamples))

isdatConnection.close()

%% Prepare data
res = []; t = [];
for seg = 1:sciData.segments
    data = reshape(sciData.data(seg,:),...
        length(sciData.data(seg,:))/sciData.samples(seg),...
        sciData.samples(seg))';
    res = [res; double(data)];
    time = ISDAT.Epoch.toDoubleArray(sciData.time(seg));
    t = [t; time];
end
