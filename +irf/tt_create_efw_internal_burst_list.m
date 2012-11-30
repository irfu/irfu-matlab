% IRF.CREATE_EFW_INTERNAL_BURST_LIST
%	script used to generate efw internal burst list from CAA inventory
%   and save it at the end on IRF server. 
%
% 	See also IRF.TT,IRF.TImeTable
%
% $Id$


tint=[irf_time([2001,1,1,0,0,0]) irf_time(now,'date2epoch')];
irf_log('dsrc',['Checking interval: ' irf_time(tint,'tint2iso')]);
c_eval('s?=caa_download(tint,''list:C?_CP_EFW_L2_PB'');');
disp('You can choose to save on server the final lists.')
disp('If you choose not to save then locally will be created')
disp('time tables with names TT_C1_EFW_internal_bursts,TT_C2...');

questionSaveOnServer = irf_ask('Save on server the obtained lists y/n? [%]','questionSaveOnServer','n');
if strcmpi(questionSaveOnServer,'y'),
	saveOnServer=true;
else
	saveOnServer=false;
end

for ic=1:4,
	disp(['Creating time table for Cluster ' num2str(ic) '..']); 
	c_eval('s=s?;',ic);
	tint=zeros(numel(s),2);
	toRemove= false(numel(s),1);
	textLine=regexp(s,'^C\d?_CP_EFW_L2_PB\s*(?<start>[\d-]*\s[\d:]*)\s*(?<end>[\d-]*\s[\d:]*)\s*(?<number>\d*)\s*(?<version>\d*)','names');
%	textLine=textscan(s,'^C\d?_CP_EFW_L2_PB\s*(?<start>[\d-]*\s[\d:]*)\s*(?<end>[\d-]*\s[\d:]*)\s*(?<number>\d*)\s*(?<version>\d*)','names');
	for j=1:numel(textLine),
		if isempty(textLine{j})
			toRemove(j)=true;
		elseif strcmp(textLine{j}.number,'0')
			toRemove(j)=true;
		else
			tint(j,:)=irf_time([textLine{j}.start '/' textLine{j}.end],'iso2tint');
			tint(j,1)=tint(j,1)-1;
			tint(j,2)=tint(j,2)+1;
		end
	end
	tint(toRemove,:)=[];
	TT=irf.TimeTable(tint);
	c_eval('TT.Header={''Intervals of EFW internal bursts. Cluster ?''};',ic);
	label=['C' num2str(ic) '_EFW_internal_bursts'];
	if saveOnServer,
		irf.tt(TT,'write_IRF',label)
		disp(['Time table ''' label ''' saved on IRF server.']);
	else
		varName=['TT_' label];
		assignin('base',varName,TT);
		disp(['Variable ' varName ' created']);
	end
end


