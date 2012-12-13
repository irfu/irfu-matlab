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
		disp(' You are using NASA cdf patch which within ')
		disp(' irfu-matlab is still experimental. Please report')
		disp(' if you encounter some problems.');
		disp(' If you want to use matlab cdf execute:');
		a=which('cdfread');
		ai=strfind(a,'/');
		disp(['> rmpath ' a(1:ai(end))]);
		disp('> clear databoj');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		fprintf('\n\n\n');
		irf_log('fcal','Using NASA cdf is still experimental!');
		ok=true;
		break;
	else
		ok=false;
	end
end
fclose(fid);
end
