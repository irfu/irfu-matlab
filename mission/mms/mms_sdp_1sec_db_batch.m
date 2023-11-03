% Script for processing one month of data

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%yymmList = {'201511','201512','201601','201602','201610','201611','201612','201701'};
yymmList = {'201706','201707','201708'};
mmsId = 'mms2'; dataPath = '/data/mms';
for idxMo=yymmList
  yymm = idxMo{:};
  %%
  %yymm = '201511';
  db.files = {}; db.tint = {}; db.data = []; db.mmsId = mmsId;
  for dd=1:31
    day = sprintf('%02d',dd)
    dataDir = [dataPath '/' mmsId '/edp/fast/l2a/dce2d/' yymm(1:4) '/' yymm(5:6) '/'];
    fName = mms_find_latest_version_cdf(...
      [dataDir '*' yymm day '*.cdf'])
    if isempty(fName), continue, end

    disp(['Processing: ' fName.name])
    [out,Tint] = mms_sdp_1sec_db(fName.name,dataDir);
    if isempty(out), continue, end
    if isempty(db.data), db.data = out;
    else
      % Append data
      fi = fields(out);
      for idx = 1:length(fi)
        db.data.(fi{idx}) = [db.data.(fi{idx}); out.(fi{idx})];
      end
    end
    % Save file names for Bookkeeping and possible reprocessing
    db.files = [db.files fName.name]; db.tint = [db.tint Tint];
  end
  save([mmsId '__' yymm], 'db');
  %%
end

