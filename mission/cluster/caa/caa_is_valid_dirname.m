function res = caa_is_valid_dirname(s)
%CAA_IS_VALID_DIRNAME  check if directory was created by the CAA s/w
%
% res = caa_is_valid_dirname(dir_name)
%

% Copyright 2006 Yuri Khotyaintsev

res = 0;

if length(s)~=13, return, end

try
  y = str2double(s(1:4));
  if y<2000 || y>2024
    irf.log('warning',['Invalid year of dataset ' num2str(y)])
    return
  end
  m = str2double(s(5:6));
  if m<1 || m > 12, return, end
  d = str2double(s(7:8));
  if d<1 || d > 31, return, end
  h = str2double(s(10:11));
  if h<0 || h > 23, return, end
  m = str2double(s(12:13));
  if m<0 || m > 59, return, end
  res = 1;
catch
  return
end
