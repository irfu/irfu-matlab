function res = c_peace_read_qjas_cdf(fname)
%C_PEACE_READ_QJAS_CDF  read QJAS CDF files 
%
% data = c_peace_read_qjas_cdf(fname)
%
% Read PEACE data exported from QJas, http://www.space-plasma.qmul.ac.uk
%
% Output: Data structure with the following fields
%              t: Time stamps (epoch)
%             dt: 
%            psd: Phase Space Density
%          level: Peace Energy (eV)
%    level_delta: 
%            phi: Peace tzero Az Flow Angle
%      phi_delta: 
%          theta: Peace Polar Flow Angle
%    theta_delta: 
%
%    See also IRF_CDF_READ, CDFREAD
%
% $Id$


% Copyright 2007 Yuri Khotyaintsev

if ~exist(fname,'file'), error(['file ' fname ' does not exist']), end

t = cdfread(fname,'Variable','timetags');

% Convert time to isdat epoch
temp = struct([t{:,1}]);
t = double([temp.date]);
t = (t(:)-62167219200000)/1000;
clear temp

% Time deltas
dt = cdfread('qjas_data.cdf','Variable','timetags_delta');
dt = double(cell2mat(dt));

% PSD
psdcell = cdfread(fname,'Variable','psd');
ndata = length(psdcell);
[n,m] = size(psdcell{1});
psd = zeros(ndata,n,m);
for i=1:ndata, psd(i,:,:) = double(psdcell{i}); end
clear psdcell

% Energy levels
level = cdfread(fname,'Variable','level');
level = cell2mat(level')';
if any(any(diff(level)))
	error('energy channels change withing the time interval')
end
level = double(level(1,:)');

level_delta = cdfread(fname,'Variable','level_delta');
level_delta = cell2mat(level_delta')';
if any(any(diff(level_delta)))
	error('energy channels (delta) change withing the time interval')
end
level_delta = double(level_delta(1,:)');

% Phi (??)
phi = cdfread(fname,'Variable','phi');
phi = double(cell2mat(phi));

phi_delta = cdfread(fname,'Variable','phi_delta');
phi_delta = double(cell2mat(phi_delta));

% Theta
theta = cdfread(fname,'Variable','theta');
theta = cell2mat(theta')';
if any(any(diff(theta)))
	error('thetas change withing the time interval')
end
theta = double(theta(1,:)');

theta_delta = cdfread(fname,'Variable','theta_delta');
theta_delta = cell2mat(theta_delta')';
if any(any(diff(theta_delta)))
	error('thetas (delta) change withing the time interval')
end
theta_delta = double(theta_delta(1,:)');

res = struct('t', t, 'dt', dt, 'psd', psd,...
	'phi', phi, 'phi_delta', phi_delta,...
	'level', level, 'level_delta', level_delta,...
	'theta', theta, 'theta_delta', theta_delta);
