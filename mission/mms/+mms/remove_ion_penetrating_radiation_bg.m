function [enflux_new, enflux_BG, idist_new, idist_BG, Ni_new, gseVi_new, gsePi_new, Ni_bg, EnergySpectr_bg, Pres_bg, EnergySpectr_bg_self] = remove_idist_background(varargin)
%MMS.REMOVE_IDIST_BACKGROUND  remove penetrating radiation background from ion distribution function
%
% [idist_new, Ni_bg] = MMS.REMOVE_IDIST_BACKGROUND(idist, 'tint', tint)
%
% Input:
%       idist - Ion distribution as a TSeries
%
% Options:
%       'tint' - time interval; DEFAULT: time interval of input eDist;
%
% Output:
%       1. enflux_new      - corrected omni-directional energy flux
%       2. enflux_BG       - background omni-directional energy flux
%       3. idist_new       - corrected ion distribution functions
%       4. idist_BG        - background ion distribution functions
%       5. Ni_new          - corrected ion number density: Ni - Ni_bg
%       6. gseVi_new       - corrected ion bulk velocity:  Ni * gseVi / Ni_new
%       7. gsePi_new       - corrected ion pressure tensor: see the following link
%       8. Ni_bg           - background ion number density [loaded from dis-moms]
%       9. EnergySpectr_bg - background ion energy flux [loaded from dis-moms]
%       10. Pres_bg        - background ion pressure scaler [ntime X 1]
%       11. EnergySpectr_bg_self - self calculated energy flux background from the lowest five energy channels.
%
% References:
%   1. MMS DIS Penetrating radiation correction methods:
%       https://lasp.colorado.edu/galaxy/display/MFDPG/Penetrating+Radiation+in+DIS+Data
%
%   2020-06-09, wyli;

%%  1. basic
Ni_new = 0;
Vi_new = 0;
Presi_new = 0;
[~, args, nargs] = axescheck(varargin{:});
idist = args{1};
Tint = irf.tint(idist.time(1), idist.time(end));
args_tmp = args(2: end);
if nargs > 1, have_options = 1; else, have_options = 0; end
while have_options
  l = 1;
  switch(lower(args_tmp{1}))
    case 'tint'
      l = 2;
      tint = args_tmp{2};
  end
  args_tmp = args_tmp(l+1:end);
  if isempty(args_tmp), break, end
end
idist_name = idist.name; %linkprop % previously linkprop was not commented out
ic = idist_name(4);     ic = str2num(ic);
idist_BG = idist;
idist_new = idist;

%%  2. load BACKGROUND from 'DIS-MOMS' data
% 2.1. background data (was 1:17 - now 1:18 - and was 10:13 - now 15:18)
c_eval('Ni_bg = mms.db_get_ts([idist_name(1:17), ''dis-moms''], [''mms?_dis_numberdensity_bg_'', idist_name(10:13)], Tint);', ic);
c_eval('EnergySpectr_bg = mms.db_get_ts([idist_name(1:17), ''dis-moms''], [''mms?_dis_spectr_bg_'', idist_name(10:13)], Tint);', ic);
c_eval('Pres_bg = mms.db_get_ts([idist_name(1:17), ''dis-moms''], [''mms?_dis_pres_bg_'', idist_name(10:13)], Tint);', ic);
c_eval('omniEnergySpectr_bg = mms.db_get_ts([idist_name(1:17), ''dis-moms''], [''mms?_dis_energyspectr_omni_'', idist_name(10:13)], Tint);', ic);

% 2.2. official moments data (was 10:13 - now 15:18)
Ni = mms.get_data(['Ni_fpi_', idist_name(10:13), '_l2'], Tint, ic);
gseVi = mms.get_data(['Vi_gse_fpi_', idist_name(10:13), '_l2'], Tint, ic);
gsePi = mms.get_data(['Pi_gse_fpi_', idist_name(10:13), '_l2'], Tint, ic);

%%  3. compute
Units = irf_units;
% 3.1. create idist_BG & idist_new
Wi_table = idist.depend{1};         Wi_table_size = size(Wi_table);             % ntime X 32;
EnergySpectr_bg_tmp = EnergySpectr_bg.data;
if not(length(EnergySpectr_bg_tmp) == Wi_table_size(1))
  error('idist & background dimension do not match! ');
end
EnergySpectr_bg_tmp = repmat(EnergySpectr_bg_tmp, 1, Wi_table_size(2));
idist_bg = Units.mp^2 * double(EnergySpectr_bg_tmp) * 1e4 ./ (2 * Units.e^2 * double(Wi_table.^2));   % [s^3/m^6]
idist_bg = idist_bg / 1e12;
idist_size = size(idist.data);
idist_bg = repmat(idist_bg, 1, 1, idist_size(3), idist_size(4));
idist_BG.data = idist_bg;
idist_new.data = idist.data - idist_bg;
idist_new.data(idist_new.data < 0) = 0;

% 3.2. create omni specrec for BACKGROUND & NEW
enflux_original = idist.deflux.omni.specrec;
enflux_BG = enflux_original;        enflux_BG.p = EnergySpectr_bg_tmp;
enflux_new = enflux_original;
enflux_new.p = enflux_new.p - enflux_BG.p;
enflux_new.p(enflux_new.p<0) = NaN;
tmp = mean(enflux_original.p(:, 1:5), 2);
EnergySpectr_bg_self = EnergySpectr_bg;
EnergySpectr_bg_self.data = tmp;

% 3.3. compute new moments
Ni_new = Ni - Ni_bg;
gseVi_new = Ni * gseVi / Ni_new;
gsePi_new = gsePi;      gsePi_new_data = gsePi_new.data;
% PiXX
gsePi_new_data(:, 1, 1) = gsePi_new_data(:, 1, 1) - Pres_bg.data(:) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 1) .* gseVi.data(:, 1) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 1) .* gseVi_new.data(:, 1);
% PiYY
gsePi_new_data(:, 2, 2) = gsePi_new_data(:, 2, 2) - Pres_bg.data(:) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 2) .* gseVi.data(:, 2) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 2) .* gseVi_new.data(:, 2);
% PiZZ
gsePi_new_data(:, 3, 3) = gsePi_new_data(:, 3, 3) - Pres_bg.data(:) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 3) .* gseVi.data(:, 3) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 3) .* gseVi_new.data(:, 3);
% PiXY & PiYX
gsePi_new_data(:, 1, 2) = gsePi_new_data(:, 1, 2) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 1) .* gseVi.data(:, 2) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 1) .* gseVi_new.data(:, 2);
gsePi_new_data(:, 2, 1) = gsePi_new_data(:, 1, 2);
% PiXZ & PiZX
gsePi_new_data(:, 1, 3) = gsePi_new_data(:, 1, 3) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 1) .* gseVi.data(:, 3) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 1) .* gseVi_new.data(:, 3);
gsePi_new_data(:, 3, 1) = gsePi_new_data(:, 1, 3);
% PiYZ & PiZY
gsePi_new_data(:, 2, 3) = gsePi_new_data(:, 2, 3) ...
  + Units.mp * Ni.data(:) .* gseVi.data(:, 2) .* gseVi.data(:, 3) ...
  - Units.mp * Ni_new.data(:) .* gseVi_new.data(:, 2) .* gseVi_new.data(:, 3);
gsePi_new_data(:, 3, 2) = gsePi_new_data(:, 2, 3);
% corrected gse ion pressure
gsePi_new.data = gsePi_new_data;

end
%%