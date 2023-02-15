function fName = mms_find_latest_version_cdf(dir_input)
%MMS_FIND_LATEST_VERSION_CDF  find the highest version file(-s)
%
% FNAME = MMS_FIND_LATEST_VERSION_CDF(DIR_INPUT)
%
% Find the highest version file(-s) from the list returned by dir(DIR_INPUT)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

fName = '';
fList = dir(dir_input); if isempty(fList), fName=fList; return, end

fNameTmp = '';
for i=length(fList):-1:1 % Reversed order of "dir"
  % Store only the superseeded versions of each file.
  filenameData = mms_fields_file_info(fList(i).name);
  if isempty(fNameTmp)
    fNameTmp = fList(i);
    prevVer = filenameData.vXYZ;
    prevDate = filenameData.date;
    continue
  end
  if prevDate == filenameData.date
    if is_version_larger(filenameData.vXYZ, prevVer)
      fNameTmp(end) = fList(i);
      prevVer = filenameData.vXYZ;
    end
  else
    % New date/time stamps (burst?)
    fNameTmp(end+1) = fList(i); %#ok<AGROW>
    prevVer = filenameData.vXYZ;
    prevDate = filenameData.date;
  end
end

% Reverse order back to "dir" default, increasing with names of file.
for i=1:length(fNameTmp)
  if i==1, fName = fNameTmp(end); continue; end
  fName(end+1) = fNameTmp(end-i+1); %#ok<AGROW>
end

end
