function ok=check_if_using_nasa_cdf
% IRF.CHECK_IF_USING_NASA_CDF returns which cdfread routines are used
%	If NASA cdf patch is not used then routine adds the path to
%	irfu-matlab/contrib/nasa_cdf_patch directory so that NASA cdf patch is used. If
%	this does not succeed, routine returns false value.
%
%   OK=IRF.CHECK_IF_USING_NASA_CDF
%		ok=true  - NASA spdfcdfread used
%		ok=false - Matlab original cdfread used
%
% It is encouraged to use NASA spdfcdfread as it gives faster reading
% times and it directly works with EPOCH16 and TT2000 type of variables.
% 
% NASA cdf patch is located under 'irfu-matlab/contrib/nasa_cdf_patch'
% The NASA cdf patch web page:
% https://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html

irfDir = fileparts(which('irf.m'));
irfNasaCdfDir = fileparts(which('spdftt2000todatenum.m'));

if isempty(irfNasaCdfDir)
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  disp('NASA cdf patch is not on your path !!!');
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  ok = false;
  patchPath = [irfDir, filesep, 'contrib', filesep, 'nasa_cdf_patch'];
  addpath(patchPath);
  disp(['Added to path: ' patchPath]);
  ok = irf.check_if_using_nasa_cdf;
  if ~ok
    disp(['ERROR: could not add to path: ' patchPath]);
  end
else
  ok = true;
  if nargout == 0, disp('NASA cdf patch is being used'), end
end

if nargout == 0, clear ok; end
