%read THEMIS position data and save into mRth.mat file

dataDir = '/data/themis';
thIds = 'abcde';

for thId=thIds
  R = [];
  for year=2007:2020
    fullPath = sprintf('%s%sth%s%sssc%s%02d',...
      dataDir,filesep,thId,filesep,filesep,year);
    if ~exist(fullPath,'dir'), continue, end
    files = dir(sprintf('%s%sth%s_or_ssc_*_v*.cdf',fullPath,filesep,thId));

    if ~isempty(files)
      for iFile=1:length(files)
        fileToRead = [fullPath,filesep,files(iFile).name];
        fprintf('reading %s\n',fileToRead)
        tmpData = spdfcdfread(fileToRead,'CombineRecords',true,'Variable','XYZ_GSE');
        tmpData = tmpData*6371.2; % comvert to kilometers
        tmpEpoch = spdfcdfread(fileToRead,'CombineRecords',true,...
          'KeepEpochAsIs',true,'Variable','Epoch');
        tmpEpoch = irf_time(tmpEpoch,'cdfepoch>epoch');
        R = [R; tmpEpoch tmpData]; %#ok<AGROW>
        clear tmpData tmpEpoch
      end
    end
  end
  % remove repeating points at month boundary
  ii = find(diff(R(:,1))==0); R(ii,:) = [];
  eval(['Rth' thId '=R;'])
  fprintf('Rth%s >> mRth.mat\n',thId);
  if exist('./mRth.mat','file')
    eval(['save(''mRth'',''Rth' thId ''',''-append'')'])
  else, eval(['save(''mRth'',''Rth' thId ''')'])
  end
end



