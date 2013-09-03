function ok=check_if_using_nasa_cdf
% IRF.CHECK_IF_USING_NASA_CDF returns which cdfread routines are used
%	If NASA cdf patch is not used then routine adds the path to
%	irfu-matlab/matlab_cdf341_patch directory so that NASA cdf patch is used. If
%	this does not succeed, routine returns false value.
%
%   OK=IRF.CHECK_IF_USING_NASA_CDF
%		ok=true  - NASA cdfread used
%		ok=false - Matlab original cdfread used
%
% It is encouraged to use NASA cdfread as it gives faster reading
% times and it directly works with EPOCH16 type of variables.
% 
% NASA cdf patch is located under 'irfu-matlab/matlab_cdf341_patch' 
% The NASA cdf patch web page:
% http://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html

cdfDir = fileparts(which('cdfread'));
irfDir = fileparts(which('irf.m'));
irfNasaCdfDir = [irfDir filesep 'matlab_cdf341_patch'];
ok = false; % default

if strfind(cdfDir,'patch')
	usingNasaCdfRead = true;
	if ~strcmp(cdfDir,irfNasaCdfDir),
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		disp('You are using NASA cdf patch but it is not located');
		disp('in irfu-matlab directory. Please, check that it is');
		disp('latest version NASA cdf patch!');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		return
	end
else
	usingNasaCdfRead = false;
end

if usingNasaCdfRead,
		fprintf('\n\n\n');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		disp(' You are using NASA cdf patch. It gives')
		disp(' you faster reading times but please report')
		disp(' if you encounter some problems.');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		fprintf('\n\n\n');
		ok=true;
else
	addpath(irfNasaCdfDir);
	ok = irf.check_if_using_nasa_cdf;
	if ~ok,
		disp(['ERROR: could not add to path: ' irfNasaCdfDir]);
	end
end

if nargout == 0, clear ok; end