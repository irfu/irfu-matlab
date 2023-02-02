function out = irf_convert(in,flag)
%IRF_CONVERT  common but nonstandard unit conversions
%   out=IRF_CONVERT(in,flag); % Convert in to out
%       flag is in format 'in2out'
%
%  Example:
%    out=irf_convert(in,'eV2nm'); % convert from eV to nm
%
% Some common conversions: eV2nm, nm2eV, eV2K
%
% To see all possibilities execute
%   irf_convert
%

if nargin==0
  disp(' ');
  disp('out=irf_convert(in,flag); % flag is in format ''in2out''');
  disp('Possible flag values:');
  disp(' ');
  a=fscanf(fopen(which('irf_convert')),'%s');
  %    b=regexp(a,'case\s.*','match');
  b=regexp(a,'case''\w*''','match');
  for ib=1:numel(b)
    disp(b{ib}(6:end-1));
  end
  disp(' ');
  return
end

Units=irf_units;

switch flag
  case 'eV2nm'
    out=Units.h*Units.c/Units.e*1e9./in;
  case 'nm2eV'
    out=Units.h*Units.c/Units.e*1e9./in;
  case 'eV2K'
    out=Units.e/Units.kB.*in;
  case 'eV2Hz'
    out=Units.e/Units.h*in;
  case 'keV2K'
    out=Units.e*1e3/Units.kB.*in;
  case 'eV2MK'
    out=Units.e/Units.kB*1e-6.*in;
  case 'keV2MK'
    out=Units.e*1e3/Units.kB*1e-6.*in;
  case 'K2eV'
    out=Units.kB/Units.e.*in;
  case 'K2keV'
    out=Units.kB*1e-3/Units.e.*in;
  case 'MK2eV'
    out=Units.kB*1e6/Units.e.*in;
  case 'MK2keV'
    out=Units.kB*1e3/Units.e.*in;
  case 'F2C' % Fahrenheit to Celsius
    out=(in-32)*10/18;
  otherwise
    disp(['!!! irf_convert: unknown flag ''' lower(flag) ''', not converting.'])
    disp('CONSIDER ADDING THIS CONVERSION;)!');
    out=in;
end
