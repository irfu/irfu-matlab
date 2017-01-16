% Script for processing one moth of data

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

yymm = '2015/12'; mmsId = 'mms1'; dataPath = '/Users/yuri/Data/mms';
%%
db.files = {}; db.tint = {};db.data = {};
for dd=1:31
  day = sprintf('%02d',dd)
  dataDir = [dataPath '/' mmsId '/edp/fast/l2a/dce2d/' yymm '/'];
  fName = mms_find_latest_version_cdf(...
    [dataDir '*' yymm(1:4) yymm(6:7) day '*.cdf'])
  if isempty(fName), continue, end

  disp(['Processing: ' fName.name])
  [out,Tint] = mms_sdp_1sec_db(fName.name,dataDir);
  
  % Append data
  fi = fields(out);
  for idx = 1:length(fi)
    db.data.(fi{idx}) = [db.data.(fi{idx}) out.(fi{idx})];
  end
  % Save file names for Bookkeeping and possible reprocessing
  db.files = [db.files fName]; db.tint = [db.tint fName];
end