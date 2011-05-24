function [t,res]=irf_isdat_get(datasource,st,dt)
% getting data from isdat (using ISDAT.rar)
% FUNCTIONS STILL IN DEVELOPMENT, do not use
% [t,res]=irf_isdat_get('Cluster/1/ephemeris/position','2002-03-04T10:00:00.000Z',60)


wd=which('irf_plot.m');
ii=strfind(wd,filesep);
wdir=wd(1:ii(end));
javaclasspath([wdir 'ISDAT.jar']);

%Open connection
isdatConnection = ISDAT.NetConnection('db.irfu.se',0);
dataSrcName=ISDAT.DataSourceName('Cluster/1/ephemeris/position');
%dataSrcName=ISDAT.DataSourceName('Cluster/1/efw/E/p4/10Hz/lx');

isdatDbSession=ISDAT.DbSession(isdatConnection);
dataSrc=isdatDbSession.SourceName2ID(dataSrcName);
isdatDbSession.nameThisSourceFully(dataSrc)

st = ISDAT.Epoch('2002-03-04T10:00:00.000Z');
dt = ISDAT.Epoch(2.0);

dataReq = ISDAT.DbDataRequest(dataSrc,st,dt);
                        
res = isdatDbSession.getSciData(dataReq);
disp(sprintf('totNrSamples : %d',res.totNrSamples))

isdatConnection.close()
