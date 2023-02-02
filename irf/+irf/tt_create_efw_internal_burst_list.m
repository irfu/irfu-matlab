% IRF.TT_CREATE_EFW_INTERNAL_BURST_LIST
%	script used to generate efw internal burst list from CAA inventory
%   and save it at the end on IRF server.
%
% 	See also IRF.TT,IRF.TImeTable
%


tint=[irf_time([2001,1,1,0,0,0]) irf_time(now,'date>epoch')];
irf_log('dsrc',['Checking interval: ' irf_time(tint,'tint>utc')]);
c_eval('disp(''Downloading inventory for Cluster ?'');s?=caa_download(tint,''list:C?_CP_EFW_L2_PB'');');
disp('You can choose to save on server the final lists.')
disp('If you choose not to save then locally will be created')
disp('time tables with names TT_C1_EFW_internal_bursts,TT_C2...');

questionSaveOnServer = irf_ask('Save on server the obtained lists y/n? [%]','questionSaveOnServer','n');
if strcmpi(questionSaveOnServer,'y')
  saveOnServer=true;
else
  saveOnServer=false;
end

for ic=1:4
  disp(['Creating time table for Cluster ' num2str(ic) '..']);
  c_eval('TT=s?;',ic);
  TT.Header=[sprintf('Intervals of EFW internal bursts. Cluster %d\n\n',ic) TT.Header];
  ii=find(TT.UserData.number);
  TT=select(TT,ii);
  label=['C' num2str(ic) '_EFW_internal_bursts'];
  if saveOnServer
    irf.tt(TT,'write_IRF',label)
    disp(['Time table ''' label ''' saved on IRF server.']);
  else
    varName=['TT_' label];
    assignin('base',varName,TT);
    disp(['Variable ' varName ' created']);
  end
end


