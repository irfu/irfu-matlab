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
irfNasaCdfDir = fileparts(which('tt2000todatenum.m'));

if isempty(irfNasaCdfDir)
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  disp('NASA cdf patch is not on your path !!!');
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  ok = false;
else
  if strcmp(cdfDir,irfNasaCdfDir), 
    ok = true;
    if nargout == 0, disp('NASA cdf patch is being used'), end
  else
    addpath(irfNasaCdfDir);
    disp(['Added to path: ' irfNasaCdfDir]);
    ok = irf.check_if_using_nasa_cdf;
    if ~ok,
      disp(['ERROR: could not add to path: ' irfNasaCdfDir]);
    end
  end
end

if nargout == 0, clear ok; end