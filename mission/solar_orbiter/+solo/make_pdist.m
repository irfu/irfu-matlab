function PDout = make_pdist(varargin)
% SOLO.MAKE_PDIST Make a PDist skymap variable from Solar Orbiter Particle data
%
% PDist = mms.make_pdist(filename)
%
% Input:
%   file - filename and directory of cdf file, or loaded cdf file
%
% Output:
%   PDout - PDist skymap variable
%
% Examples:
%   PDout = make_pdist('solo_L2_swa-pas-vdf_20200729_V01.cdf')
%
%   tmpD = dataobj('solo_L2_swa-pas-vdf_20200729_V01.cdf');
%   PDout = solo.make_pdist(tmpD)
%
% Written by D. B. Graham
%
% TO DO: Add electron distributions

if (nargin > 1)
  error('Too many inputs. Need to add electron data.')
end

file = varargin{1};
flipenergy = 0;

if isa(file,'dataobj')
  tmpDataObj = file;
elseif isa(file,'char')
  tmpDataObj = dataobj(file);
else
  error('Input not recognized.')
end

if isempty(tmpDataObj); return; end

% Identify variable type
fname = tmpDataObj.GlobalAttributes.Logical_file_id{1};
if ~isempty(regexp(fname,'swa-pas-vdf','once'))
  % Load ion data
  Energyarr = get_variable(tmpDataObj,'Energy');
  Energyarr = Energyarr.data;
  if Energyarr(1) > Energyarr(end)
    flipenergy = 1;
    Energyarr = flip(Energyarr);
  end
  %dEnergyplus = get_variable(tmpDataObj,'delta_p_Energy');
  %dEnergyminus = get_variable(tmpDataObj,'delta_m_Energy');
  Azimuthangle = get_variable(tmpDataObj,'Azimuth');
  deltaAzimuth = get_variable(tmpDataObj,'delta_Azimuth');
  vdf = get_variable(tmpDataObj,'vdf');
  vdfp  = permute(vdf.data,[1 4 2 3]);
  if flipenergy
    vdfp = flip(vdfp,2);
  end
  Polarangle = get_variable(tmpDataObj,'Elevation');
  deltaPolar = get_variable(tmpDataObj,'delta_Elevation');
  time = vdf.DEPEND_0.data;
  time = EpochTT(time);
  dt = median(diff(time.epochUnix));

  energymat = ones(size(time))*Energyarr';

  PDout = PDist(time,vdfp,'skymap',energymat,-Azimuthangle.data,90-Polarangle.data); % Negative of azimuthal angle is used so Velocity moments agree
  PDout.species = 'ions';
  PDout.units = 's^3/m^6';
  PDout.siConversion = '1e12';
  PDout.ancillary.energy0 = Energyarr;
  PDout.ancillary.energy1 = Energyarr;
  %PDout.ancillary.delta_energy_plus = Energyarr.data-dEnergyplus.data; % Need to check this
  %PDout.ancillary.delta_energy_minus = Energyarr.data-dEnergyminus.data;
  PDout.ancillary.esteptable = zeros(size(time));
  PDout.ancillary.delta_theta_plus = deltaPolar.data;
  PDout.ancillary.delta_theta_minus = deltaPolar.data;
  PDout.ancillary.delta_phi_plus = deltaAzimuth.data;
  PDout.ancillary.delta_phi_minus = deltaAzimuth.data;

  PDout.ancillary.dt_minus = dt/2; % Need to check if time tags are the beginning or middle of distributions
  PDout.ancillary.dt_plus = dt/2;



  PDout.name = 'solo_L2_swa-pas-vdf'; % Temp fix

  ud = [];
  ud.GlobalAttributes = tmpDataObj.GlobalAttributes;
  ud.VALIDMIN = tmpDataObj.VariableAttributes.VALIDMIN{13,2};
  ud.VALIDMAX = tmpDataObj.VariableAttributes.VALIDMAX{13,2};
  PDout.userData = ud;


else
  PDout = NaN;
  error('Data Object not recognized.')
end


end

