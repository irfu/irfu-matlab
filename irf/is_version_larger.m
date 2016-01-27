function res = is_version_larger(newVersionString,oldVersionString)
% res = is_version_larger(newVersionString,oldVersionString)
% compares version in format X.Y.Z
%
% res is true if newVersionString is larger than oldVersionString
%
verNew = get_ver(newVersionString);
verOld = get_ver(oldVersionString);
if(verNew.maj>verOld.maj) || ... % Newer major version
		(verNew.maj==verOld.maj && verNew.min>verOld.min) || ... % Same major version, newer calibration
		(verNew.maj==verOld.maj && verNew.min==verOld.min && verNew.rev>verOld.rev) % Same major and calib. but newer revision, replace file
	res = true;
else
	res = false;
end
	function ver = get_ver(verS)
		vT = strsplit(verS,'.');
		for iTmp=1:length(vT), vT{iTmp} = str2double(vT{iTmp}); end
		ver = struct('maj',vT{1},'min',0,'rev',0);
		if isempty(vT{2}), return, end, ver.min = vT{2};
		if isempty(vT{3}), return, end, ver.rev = vT{3};
	end
end
