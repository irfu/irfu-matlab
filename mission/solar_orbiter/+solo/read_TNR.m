function out = read_TNR(tint)
%     Read L2 data from TNR.
%
%     @author: Louis Richard
%     Updated by: Jordi Boldu
%
%     NOTE: Relies on hardcoded file system path for finding CDFs and then
%     reading them itself. Does not use irfu-matlab's solo.db_get_ts() for
%     locating and obtaining CDF data due to one zVariable not being readable
%     with it.
%
%     Parameters
%     ----------
%     path : str
%         String of the filename in .cdf containing the L2 data
%
%     tint : EpochTT (2x1)
%         Time interval
%         NOTE: Will only load the TNR file (solo_L2_rpw-tnr-surv-cdag)
%         corresponding to the initial timestamp.
%
%     sensor : int
%         TNR sensor to be read:
%             * 1: V1
%             * 2: V2
%             * 3: V2,
%             * 4: V1 - V2
%             * 5: V2 - V3
%             * 6: V3 - V1
%             * 7: B
%
%     Returns
%     -------
%     out : struct
%         Spectrum of the measured signals.
%
%     Notes
%     -----
%     The script check if there are data from the two channel and put them
%     together.
%

%%

%tint = irf.tint('2021-05-18T00:00:00Z','2021-05-19T00:00:00Z');
date_vector = irf_time(tint,'epochTT>vector');
yyyy = num2str(date_vector(1,1));

mm = num2str(date_vector(1,2));
if numel(mm)==1;mm=['0' mm];end

dd = num2str(date_vector(1,3));
if numel(dd)==1;dd=['0' dd];end



sensor = 5;
sensor2 = 4;

% IMPLEMENTATION NOTE: solo.get_db_ts() does not work to get the zVariable
% TNR_BAND_FREQ from the TNR cdf file, therefore the dataobj(x) function
% is used instead, which requires giving the full path of the file.
% The solo.get_db_ts() function seems to fail to create the TSeries object
% because the DEPEND_0 field is of different size from the data.
% NOTE: Path with wildcard.
path = ['/data/solo/remote/data/L2/thr/' yyyy '/' mm '/solo_L2_rpw-tnr-surv-cdag_' yyyy mm dd '_V*.cdf'];

data_l2 = read_TNR_CDFs(path, tint);

% =============================================================================
% Information on selected solo_L2_rpw-tnr-surv zVariables used by this function
% (=fields in struct "data_l2")
% =============================================================================
%
% zVariable Information
% ---------------------
% AUTO1                     CDF_DOUBLE/1   1:[32]    T/T
% AUTO2                     CDF_DOUBLE/1   1:[32]    T/T
% FRONT_END                 CDF_UINT1/1    0:[]      T/
% MAGNETIC_SPECTRAL_POWER1  CDF_DOUBLE/1   1:[32]    T/T
% MAGNETIC_SPECTRAL_POWER2  CDF_DOUBLE/1   1:[32]    T/T
% SWEEP_NUM                 CDF_UINT4/1    0:[]      T/
% TNR_BAND                  CDF_UINT1/1    0:[]      T/
% TNR_BAND_FREQ             CDF_UINT4/1    2:[4,32]  F/TT
% SENSOR_CONFIG             CDF_UINT1/1    1:[2]     T/T
%
% Descriptions in zVariable attribute VAR_NOTES
% ---------------------------------------------
% AUTO1                     "Power spectral density at receiver+PA for channel 1 before applying antenna gain"
% AUTO2                     "Power spectral density at receiver+PA for channel 2 before applying antenna gain"
% FRONT_END                 "Indicates the TNR front end setting (GND=0, PREAMP=1, CAL=2)"
% MAGNETIC_SPECTRAL_POWER1  "Magnetic power spectral density from 1 search coil axis in channel 1"
% MAGNETIC_SPECTRAL_POWER2  "Magnetic power spectral density from 1 search coil axis in channel 2"
% SWEEP_NUM                 "TNR sweep index number in the current file"
% TNR_BAND                  "TNR band of the current record. Possible values are: 1=A, 2=B, 3=C, 4=D"
% TNR_BAND_FREQ             "Frequencies of analysis of the 4 TNR bands in Hz"
% SENSOR_CONFIG             "Indicates the THR sensor configuration
%                            (V1=1, V2=2, V3=3, V1-V2=4, V2-V3=5, V3-V1=6, B_MF=7,
%                            HF_V1-V2=9, HF_V2-V3=10, HF_V3-V1=11)"

n_freqs  = size(   data_l2.tnr_band_freq.data, 2) * 4;
freq_tnr = reshape(data_l2.tnr_band_freq.data', n_freqs, 1);

puntical_ = find(data_l2.front_end.data == 1);   % 1 = PREAMP

epoch_ = data_l2.auto1.time.epoch(puntical_, :);
auto1_ = data_l2.auto1.data(puntical_, :);
auto2_ = data_l2.auto2.data(puntical_, :);
sweep_ = data_l2.sweep_num.data(puntical_, :);
bande_ = data_l2.tnr_band.data(puntical_, :);
confg_ = data_l2.sensor_config.data(puntical_, :);

if sensor == 7
  % NOTE: Overwrite previously assigned variables.
  auto1_ = data_l2.magnetic_spectral_power1.data(puntical_, :);
  auto2_ = data_l2.magnetic_spectral_power2.data(puntical_, :);
end

sweep_num = sweep_;

delta_sw = abs(sweep_(2:end) - sweep_(1:end - 1));

xdelta_sw = find(delta_sw > 100)';

if ~isempty(xdelta_sw)
  xdelta_sw = [xdelta_sw, numel(sweep_)];

  for inswn = 1:numel(xdelta_sw) - 1
    idx_l = xdelta_sw(inswn) + 1;
    idx_r = xdelta_sw(inswn + 1);
    sweep_num(idx_l:idx_r) = sweep_num(idx_l:idx_r) ...
      + sweep_num(xdelta_sw(inswn));
  end
end


timet_ = epoch_';

sens0_1 = find(confg_(:, 1) == sensor)';
sens0_2 = find(confg_(:, 1) == sensor2)';
sens0_ = sort([sens0_1 sens0_2]);

sens1_1 = find(confg_(:, 2) == sensor)';
sens1_2 = find(confg_(:, 2) == sensor2)';
sens1_ = sort([sens1_1 sens1_2]);



if ~isempty(sens0_) && ~isempty(sens1_)
  auto_calib = [auto1_(sens0_, :); auto2_(sens1_, :)];
  sens_ = [sens0_, sens1_];
  timet_ici = [timet_(sens0_), timet_(sens1_)];
elseif ~isempty(sens0_)
  auto_calib = auto1_(sens0_, :);
  sens_ = sens0_;
  timet_ici = timet_(sens0_);
elseif ~isempty(sens1_)
  auto_calib = auto1_(sens1_, :);
  sens_ = sens1_;
  timet_ici = timet_(sens1_);
else
  out = 0;
  return;
  %irf.log('critical', 'no data at all ?!?')
end

[~, ord_time] = sort(timet_ici);
time_rr = timet_ici(ord_time);
sens_ = sens_(ord_time);
auto_calib = auto_calib(ord_time, :);

bande_e = bande_(sens_);
max_sweep = max(sweep_num(sens_));
min_sweep = min(sweep_num(sens_));
sweep_num = sweep_num(sens_);

v_        = zeros(128, 1);
sweep_tnr = zeros(  1, 1);
time_     = zeros(  1, 1);

for ind_sweep = min_sweep:max_sweep
  v1_ = zeros(128, 1);
  p_punt = find(sweep_num == ind_sweep);
  if ~isempty(p_punt)
    for indband = 1:numel(p_punt)
      idx_l = floor(32 * bande_e(p_punt(indband))) + 1;
      idx_r = floor(32 * (bande_e(p_punt(indband)) + 1));
      v1_(idx_l:idx_r) = auto_calib(p_punt(indband), :);
    end
  end

  if sum(v1_) > 0.0    % Might be a way of detecting CDF fill values (large, negative).
    punt0_ = find(v1_ == 0.0);
    if ~isempty(punt0_)
      v1_(punt0_) = NaN;
    end
    v_ = [v_, v1_];
    sweepn_tnr = [sweep_tnr, sweep_num(p_punt(1))];
  end

  if ~isempty(p_punt)
    time_ = [time_, time_rr(min(p_punt))];
  end
end

time_ = EpochTT(time_(2:end));



%select frequencies lower than 100 kHz
freq_tnr = freq_tnr(freq_tnr<100000);
f100_ind = length(freq_tnr);
vp = v_(1:f100_ind, 2:end)';

%==============Integration %needs more testing
%   for ii = 1:f100_ind
%         itg(ii) = trapz(vp(:,ii))/length(vp(:,1));
%   end
%
%    vp = vp-itg;
%   %===============
%
%
%   %==============Moving Window
%    for i = 1:f100_ind
%        movm(:,i) = movmean(vp(:,i),5);
%    end
%    vp = vp-movm;
%    vp(vp<0)=0;
%=============

out = struct('t', time_.epochUnix, 'f', freq_tnr, 'p',vp.^10);
out.p_label={'dB'};


%For plotting
%         h(10)=irf_panel('tnr');
%         irf_spectrogram(h(10),out,'log','donotfitcolorbarlabel')
%         %fpe_sc.units = 'kHz';
%         %fpe_sc.name = 'f [kHz]';
%         hold(h(10),'on');
%         %irf_plot(h(10),fpe_sc,'r','linewidth',lwidth);
%         hold(h(10),'off');
%         text(h(10),0.01,0.3,'f_{pe,RPW}','units','normalized','fontsize',18,'Color','r');
%         irf_legend(h(10),'(a)',[0.99 0.98],'color','k','fontsize',12)
%         irf_legend(h(10),'f_{pe}',[0.15 0.60],'color','k','fontsize',12)
%         set(h(10), 'YScale', 'log');
%         %set(h(10),'ColorScale','log')
%         %caxis(h(10),[.01 1]*10^-12)
%         ylabel(h(10),{'f';'(kHz)'},'interpreter','tex');
%         colormap(h(10),jet)
%         yticks(h(10),[10^1 10^2]);
%         irf_zoom(h(10),'y',[10^1 10^2])
end

%%
function out_struct = read_TNR_CDFs(path, tint)
%     Reads required field from TNR .cdf file.
%
%     @author: Louis Richard
%
%     Parameters
%     ----------
%     path : str
%       Path to L2 TNR CDF(s). May refer to multiple files by using wildcard
%       which dataobj() understands.
%
%     tint : EpochTT (2x1)
%       Time interval
%
%     Returns
%     -------
%     out_struct : struct
%       L2 data structure.

try
  data_obj = dataobj(path);
catch CauseExc
  % IMPLEMENTATION NOTE: Assumes that error is due to File-not-found.
  % "path" contains "*" i.e. it is not a real path, but dataobj() can
  % handle that. Can therefore not (easily) manually check for path
  % existence outside of dataobj().
  % SolO IRFU quicklooks (solo.qli.generate_quicklooks_all_types()) fails for
  % 2022-08-08 if not for some way of handling non-existent path here.
  % /Erik P G Johansson 2022-09-12
  Exc = MException('read_TNR:FileNotFound', 'Can not find/open file path="%s".', path);
  Exc = addCause(Exc, CauseExc);
  throw(Exc)
end

% Key of the TNR frequencies
freqs_key = 'TNR_BAND_FREQ';

% Names of the field in time series format
tseries_keys = {'AUTO1', 'AUTO2', 'SWEEP_NUM', 'TNR_BAND', ...
  'SENSOR_CONFIG', 'FRONT_END', ...
  'MAGNETIC_SPECTRAL_POWER1', 'MAGNETIC_SPECTRAL_POWER2'};

% Construct the data structure starting with the TNR frequencies
% (only not time series field)
out_struct = struct(lower(freqs_key), ...
  get_variable(data_obj, freqs_key));

for i_key = 1:numel(tseries_keys)
  key_upper = tseries_keys{i_key};
  out_struct.(lower(key_upper)) = get_ts(data_obj, key_upper).tlim(tint);
end
end

