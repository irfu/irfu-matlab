function ok=check_if_using_nasa_cdf
% IRF.CHECK_IF_USING_NASA_CDF returns which cdfread routines are used
%
%   OK=IRF.CHECK_IF_USING_NASA_CDF
%		ok=true  - NASA cdfread used
%		ok=false - Matlab original cdfread used
%
% It is encouraged to use NASA cdfread as it gives faster reading
% times and it directly works with EPOCH16 type of variables.
% 
% To obtain NASA cdf patch see their web page:
% http://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html


% $Id$


fid=fopen(which('cdfread'));
while 1
	tline = fgetl(fid);
	if ~ischar(tline), break, end
	if strfind(tline,'Mike Liu')
		fprintf('\n\n\n');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		disp(' You are using NASA cdf patch which within ')
		disp(' irfu-matlab is still experimental. It should')
		disp(' give you faster reading times but please report')
		disp(' if you encounter some problems.');
		disp(' If you want to use original matlab cdf, execute:');
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
