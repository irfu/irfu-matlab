function res = is_version_larger(newVersionString,oldVersionString)
% res = is_version_larger(newVersionString,oldVersionString)
% compares version in format X.Y.Z (or an arbitrary number of elements 
% each separated by a single ".").
%
% res is true if newVersionString is larger than oldVersionString. With
% highest priority given to first element X, if X is the same in both then
% it compares values of Y, etc.
%

verNew = get_ver(newVersionString);
verOld = get_ver(oldVersionString);
diffLen = length(verOld) - length(verNew);
if(diffLen < 0)
  % Append zeros to the end of verOld.
  warnStr=['Different length strings, appending zeros to the end of old str: ' oldVersionString];
  irf.log('warning', warnStr);
  verOld(end+1:end+abs(diffLen)) = zeros(1, abs(diffLen));
elseif(diffLen > 0)
  % Append zeros to the end of verNew.
  warnStr=['Different length strings, appending zeros to the end of new str: ' newVersionString];
  irf.log('warning', warnStr);
  verNew(end+1:end+abs(diffLen)) = zeros(1, abs(diffLen));
end
res = false; % Assume it is not newer.
if(verNew(1) > verOld(1)), res = true; return; end
for ii=2:length(verOld) % Loop through each version segment.
  if(verNew(ii) > verOld(ii) && all(verOld(1:ii-1) == verNew(1:ii-1)))
    % Newer subversion and same X. etc.
    res = true;
    return
  end
end

  function ver = get_ver(verS)
    ver = sscanf(verS,'%d.',inf); % Split it using "." as delimiter
  end
end