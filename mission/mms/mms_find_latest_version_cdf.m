function fName = mms_find_latest_version_cdf(dir_input)
%MMS_FIND_LATEST_VERSION_CDF  find the highest version file
%
% FNAME = MMS_FIND_LATEST_VERSION_CDF(DIR_INPUT)
%
% Find the highest version file from the list returned by dir(DIR_INPUT)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

fName = '';

fList=dir(dir_input); if isempty(fList), return, end
  
vMax = [0 0 0]; fName = '';
for i=1:length(fList)
  filenameData = mms_fields_file_info(fList(i).name);
  v = tokenize(filenameData.vXYZ,'.');
  v = [num2str(v{1}) num2str(v{2}) num2str(v{3})];
  if isempty(fName)
    fName = fList(i).name; vMax = v; continue
  end
  if v(1)>vMax(1) || (v(1)==vMax(1) && v(2)>vMax(2)) || ...
      (v(1)==vMax(1) && v(2)==vMax(2) && v(3)>vMax(3))
    fName = fList(i).name; vMax = v;
  end
end