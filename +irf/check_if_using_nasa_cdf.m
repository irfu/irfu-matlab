function ok=check_if_using_nasa_cdf
% IRF.CHECK_IF_USING_NASA_CDF returns which cdfread routines are used
%
%   OK=IRF.CHECK_IF_USING_NASA_CDF
%		ok=true  - NASA cdfread used
%		ok=false - Matlab original cdfread used

% $Id$


fid=fopen(which('cdfread'));
while 1
	tline = fgetl(fid);
	if ~ischar(tline), break, end
	if strfind(tline,'Mike Liu')
		fprintf('\n\n\n');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		disp(' You are using NASA cdfread patch which we have bad support!')
		disp(' This may give errors reading in multidimensional data sets!')
		disp(' Also option ''tint'' in routine databoj is disabled.');
		disp(' We suggest you to use the MATLAB cdfread!');
		disp(' To use MATLAB cdfread please remove path to NASA cdfread patch.');
		disp(' You can execute and then continue:');
		a=which('cdfread');
		ai=strfind(a,'/');
		disp(['> rmpath ' a(1:ai(end))]);
		disp('> clear databoj');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		fprintf('\n\n\n');
		irf_log('fcal','Using NASA cdf is not supported!');
		ok=true;
		break;
	else
		ok=false;
	end
end
fclose(fid);
end
