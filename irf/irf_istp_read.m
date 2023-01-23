function [time,data] = irf_istp_read(filename)
%IRF_ISTP_READ reads ISTP CDF files (ACE) for IMF and SW data
%
%  [TIME,DATA] = irf_istp_read(FILENAME) reads ISTP key parameter
%  (1-hour) CDF files for ACE.
%  TIME is returned as serial date number (see DATENUM).
%  DATA contains:
%  IMF [nT] for MFI files downloaded from
%         ftp://cdaweb.gsfc.nasa.gov/pub/istp/ace/mfi_k2/ (GSE coordinates)
%         ftp://cdaweb.gsfc.nasa.gov/pub/istp/ace/mfi_h2/ (GSM coordinates)
%  SW pressure [nPa] computed as 1.6726e-6*Np*Vp^2
%  for SWE files downloaded from
%         ftp://cdaweb.gsfc.nasa.gov/pub/istp/ace/swe_k1/
%         ftp://cdaweb.gsfc.nasa.gov/pub/istp/ace/swe_h2/
%
%   K* data are preliminary, H* data are final.
%
%   See also SPDFCDFREAD, DATENUM, DATESTR
%

% Copyright Yuri Khotyaintsev, 2003

[dt,info]=spdfcdfread(filename);

ndata = length(dt(:,1));

dm = length(dt(1,:));

time = zeros(ndata,1);

if strcmp(info.Variables{5,1},'BGSEc') %IMF K2 data
  data = zeros(ndata,3);
  for i = 1:ndata
    time(i) = todatenum(dt{i,1});
    data(i,:) = dt{i,5}';
  end
elseif strcmp(info.Variables{6,1},'BGSM') %IMF H2 data
  data = zeros(ndata,3);
  for i = 1:ndata
    time(i) = todatenum(dt{i,1});
    data(i,:) = dt{i,6}';
  end
elseif strcmp(info.Variables{6,1},'Np') && strcmp(info.Variables{7,1},'Vp')
  %SW data
  data = zeros(ndata,1);
  for i = 1:ndata
    time(i) = todatenum(dt{i,1});
    %convert to SW pressure in nPa
    data(i,:) = 1.6726e-6*double(dt{i,6})*double(dt{i,7})^2;
  end
else
  error('CDF file format not recognized.')
end
