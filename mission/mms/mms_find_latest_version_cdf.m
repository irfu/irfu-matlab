function fName = mms_find_latest_version_cdf(dir_input)
%MMS_FIND_LATEST_VERSION_CDF  find the highest version file(-s)
%
% FNAME = MMS_FIND_LATEST_VERSION_CDF(DIR_INPUT)
%
% Find the highest version file(-s) from the list returned by dir(DIR_INPUT)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

fName = '';

fList = dir(dir_input); if isempty(fList), fName=fList; return, end
  
for i=length(fList):-1:1
  filenameData = mms_fields_file_info(fList(i).name);
  if isempty(fName)
    fName = fList(i);
    prevVer = filenameData.vXYZ;
    prevDate = filenameData.date;
    continue
  end
  if prevDate == filenameData.date
    if is_version_larger(filenameData.vXYZ, prevVer)
      fName = fList(i);
      prevVer = filenameData.vXYZ;
    end
  else
    % New date/time stamps (burst?)
    fName(end+1) = fList(i); %#ok<AGROW>
    prevVer = filenameData.vXYZ;
    prevDate = filenameData.date;
  end
end

end
