function res = is_version_geq(newVersionString,oldVersionString)
%IS_VERSION_GEQ  Is version greater of equal
%
% res = IS_VERSION_GEQ(newVersionString,oldVersionString)
% compares version in format X.Y.Z (or an arbitrary number of elements 
% each separated by a single ".").
%
% res is true if newVersionString is larger of equal than oldVersionString. With
% highest priority given to first element X, if X is the same in both then
% it compares values of Y, etc.
%
% See also: IS_VERSION_LARGER

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

if verNew(1) < verOld(1), res = false; return
elseif verNew(1) > verOld(1), res = true; return
end
for ii=2:length(verOld) % Loop through each version segment.
  if verNew(ii) < verOld(ii), res = false; return
  elseif verNew(ii) > verOld(ii), res = true; return
  end
end
res = true;

  function ver = get_ver(verS)
    ver = sscanf(verS,'%d.',inf); % Split it using "." as delimiter
  end
end
