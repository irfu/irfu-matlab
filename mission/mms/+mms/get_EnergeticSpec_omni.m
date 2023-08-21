function [efeeps_omni_spec, ifeeps_omni_spec] = get_EnergeticSpec_omni(varargin)
%   mms.get_EnergeticSpec_omni: get particle flux of FEEPS/EIS energetic particle data
%       - load brst/svry mode FEEPS data; 20190923
%       - future upgrade: EIS electron and ion ion data.
%   INPUT:
%       ic:             spacecraft number
%       tint:           time interval
%       'instrument':   FEEPS or EIS
%       'svry_brst':    'svry' or 'brst' mode data;
%       'iandore';      'ie'/'i'/'e'; DEFAULT: 'ie';
%   Example:
%   [efeeps, ifeeps] = mms.get_EnergeticSpec_omni(ic, Tint, 'instrument', 'FEEPS', 'srvy_brst', 'srvy', 'iandore', 'i');

%ic, tint, 'instrument', instrument, 'svry_brst', svry_brst

% 1. get input variables;
ic = varargin{1};
tint = varargin{2};
args=varargin(3:end);       % rest variables
efeeps_omni_spec = '';      eflag = 1;
ifeeps_omni_spec = '';      iflag = 1;

if numel(args)>0
  options=1;
else
  options=0;
end

while options
  l = 2;
  switch(lower(args{1}))
    case 'instrument'
      if numel(args)>1 && ~isempty(args{2})
        instrument = args{2};
        irf.log('notice', ['Energetic particle instrument ', upper(instrument), ' loaded.'])
      end
    case 'srvy_brst'
      if numel(args)>1 && ~isempty(args{2})
        srvy_brst = args{2};
        irf.log('notice', ['Data resolution ', upper(srvy_brst), '.']);
      end
    case 'iandore'
      if numel(args)>1 && ~isempty(args{2})
        iandore = args{2};
        if ~contains(iandore, 'i'); iflag = 0; end
        if ~contains(iandore, 'e'); eflag = 0; end
      end
    otherwise
      irf.log('critical',['Unknown flag: ' args{1}]);
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), options=0; end
end

% 2. load data
if strcmp(upper(instrument), 'FEEPS')
  % 2.1. FEEPS electron
  if eflag
    c_eval('efile? = mms.get_filepath([''mms?_feeps_'' srvy_brst ''_l2_electron''], tint);', ic);
    c_eval('eobj? = dataobj(efile?);', ic);
    c_eval('eElow = get_variable(eobj?, ''electron_energy_lower_bound'');', ic);
    c_eval('eEupp = get_variable(eobj?, ''electron_energy_upper_bound'');', ic);
    c_eval('eTit1 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_1'']);', ic);
    c_eval('eTit2 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_2'']);', ic);
    c_eval('eTit3 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_3'']);', ic);
    c_eval('eTit4 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_4'']);', ic);
    c_eval('eTit5 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_5'']);', ic);
    c_eval('eTit9 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_9'']);', ic);
    c_eval('eTit10 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_10'']);', ic);
    c_eval('eTit11 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_11'']);', ic);
    c_eval('eTit12 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_top_intensity_sensorid_12'']);', ic);
    c_eval('eBit1 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_1'']);', ic);
    c_eval('eBit2 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_2'']);', ic);
    c_eval('eBit3 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_3'']);', ic);
    c_eval('eBit4 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_4'']);', ic);
    c_eval('eBit5 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_5'']);', ic);
    c_eval('eBit9 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_9'']);', ic);
    c_eval('eBit10 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_10'']);', ic);
    c_eval('eBit11 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_11'']);', ic);
    c_eval('eBit12 = get_ts(eobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_electron_bottom_intensity_sensorid_12'']);', ic);
    % 2.2. make electron omni flux
    eenergy = (eElow.data + eEupp.data)/2.;
    c_eval('eTit?.data(eTit?.data<0) = 0;', [1:5, 9:12]);      % needs to be updated --> NaN
    c_eval('eBit?.data(eBit?.data<0) = 0;', [1:5, 9:12]);
    etmp = (eTit1 + eTit2 + eTit3 + eTit4 + eTit5 + eTit9 + eTit10 + eTit11 + eTit12 ...
      + eBit1 + eBit2 + eBit3 + eBit4 + eBit5 + eBit9 + eBit10 + eBit11 + eBit12)/18;
    efeeps_omni_spec = struct('t', eTit1.time.epochUnix);
    efeeps_omni_spec.p = double(etmp.data);
    efeeps_omni_spec.p_label = {'intensity'};
    efeeps_omni_spec.f_label = {'Energy'};
    efeeps_omni_spec.f = double(eenergy);
  end

  if iflag
    % 2.3. FEEPS ions
    c_eval('ifile? = mms.get_filepath([''mms?_feeps_'' srvy_brst ''_l2_ion''], tint);', ic);
    c_eval('iobj? = dataobj(ifile?);', ic);
    c_eval('iElow = get_variable(iobj?, ''ion_energy_lower_bound'');', ic);
    c_eval('iEupp = get_variable(iobj?, ''ion_energy_upper_bound'');', ic);
    c_eval('iTit6 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_top_intensity_sensorid_6'']);', ic);
    c_eval('iTit7 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_top_intensity_sensorid_7'']);', ic);
    c_eval('iTit8 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_top_intensity_sensorid_8'']);', ic);
    c_eval('iBit6 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_bottom_intensity_sensorid_6'']);', ic);
    c_eval('iBit7 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_bottom_intensity_sensorid_7'']);', ic);
    c_eval('iBit8 = get_ts(iobj?, [''mms?_epd_feeps_'' srvy_brst ''_l2_ion_bottom_intensity_sensorid_8'']);', ic);

    % 2.4. make ion omni flux
    ienergy = (iElow.data + iEupp.data)/2.;
    c_eval('iTit?.data(iTit?.data<0) = 0;', 6:8);          % needs to be updated --> NaN
    c_eval('iBit?.data(iBit?.data<0) = 0;', 6:8);
    itmp = (iTit6 + iTit7 + iTit8 + iBit6 + iBit7 + iBit8)/6;
    ifeeps_omni_spec = struct('t', iTit6.time.epochUnix);
    ifeeps_omni_spec.p = double(itmp.data);
    ifeeps_omni_spec.p_label = {'intensity'};
    ifeeps_omni_spec.f_label = {'Energy'};
    ifeeps_omni_spec.f = double(ienergy);
  end
end

end