function res = is_version_larger(newVersionString,oldVersionString)
% res = is_version_larger(newVersionString,oldVersionString)
% compares version in format X.Y.Z (or an arbitrary number of elements as
% long as both strings contain the same number of elements separated by a
% single ".").
%
% res is true if newVersionString is larger than oldVersionString. With
% highest priority given to first element X, if X is the same in both then
% it compares values of Y, etc.
%
verNew = get_ver(newVersionString);
verOld = get_ver(oldVersionString);
if(~all(size(verOld) == size(verNew)))
  errStr='Unexpected format of version strings';
  irf.log('critical',errStr); error(errStr);
end
res = false; % Assume it is not newer.
for ii=1:length(verOld) % Loop through each version segment.
  if(verNew{ii} > verOld{ii})
    res = true; % It is newer.
    return % Exit function with result.
  end
end

	function ver = get_ver(verS)
		vT = strsplit(verS,'.');
		for iTmp=1:length(vT), vT{iTmp} = str2double(vT{iTmp}); end
        ver = vT;
	end
end