function mms_edi_epd_comparison( scId, tint )
%MMS_EDI_EPD_COMPARISON Summary of this function goes here
%   Detailed explanation goes here

% Verify input
narginchk(2,2);
if(~isa(tint,'GenericTimeArray'))
  errStr='Unexpected tint. Should be created with irf.tint("start","stop")';
  irf.log('critical',errStr);  error(errStr);
end
if((tint.stop - tint.start) > 86400)
  % Script is currently not at all optimized...
  errStr='Large tint, this will take a while... Get coffee?';
  irf.log('warning',errStr);
end

% Some limits and settings used...
limit.outlier = 3; % +/-5 mV/m is somewhat good outlier limit..
limit.edp_edi_zero1 = 2.5; % EPD > 2.5 of the outliers
limit.edp_edi_zero2 = 2.0; % abs(EDI) < 2.0 of the outliers.
limit.aspocNAN = true; % remove intervals with ASPOC on. (Set to NAN).

% Setup
global MMS_CONST
if(isempty(MMS_CONST)), MMS_CONST=mms_constants; end
data_path_root=getenv('DATA_PATH_ROOT');
if(exist(data_path_root, 'dir'))
  mms.db_init('local_file_db',data_path_root); % Setup MMS db
else
  error('Setup failed');
end

% Get EDI e-field (dsl) data
edi_e_dsl_ts = mms.db_get_ts(['mms',scId,'_edi_srvy_l2_efield'], ...
  ['mms',scId,'_edi_e_dsl_srvy_l2'],tint);

% Get EDP E-field (dsl2d) data
edp_l2pre_ts = mms.db_get_ts(['mms',scId,'_edp_fast_l2pre_dce2d'], ...
  ['mms',scId,'_edp_dce_dsl_fast_l2pre'],tint);

% Get spinfits
% dce_l2a_fit_e12_ts = mms.db_get_ts('mms1_edp_fast_l2a_dce2d', ...
%   'mms1_edp_espin_p12_fast_l2a', tint2);
% irf.log('warning','EDP spinfits loaded');

if(limit.aspocNAN)
  bitmask_ts = mms.db_get_ts(['mms',scId,'_edp_fast_l2pre_dce2d'], ...
    ['mms',scId,'_edp_bitmask_fast_l2pre'], tint);
  edp_l2pre_ts.data(logical(bitand(bitmask_ts.data(:,1), MMS_CONST.Bitmask.ASPOC_RUNNING)),:) = NaN;
  irf.log('warning','ASPOC ON has been NaN:ed in EDP data only.');
end

%% FIXME: This should not required, but to avoid extremes simply blank it here while testing!
% Ignore extreme values in either EDI or EDP
edp_l2pre_ts.data(abs(edp_l2pre_ts.data)>20) = NaN;
edi_e_dsl_ts.data(abs(edi_e_dsl_ts.data)>20) = NaN;

if(false) %% Very quickly examine data, compare EDI with closest, in time, EDP datapoint.
  % Simply use the EDP datapoint closest in time to the EDI datapoint.
  edp_t_norm = double(edp_l2pre_ts.time.ttns - edp_l2pre_ts.time.start.ttns);%#ok<UNRCH>
  edi_t_norm = double(edi_e_dsl_ts.time.ttns - edp_l2pre_ts.time.start.ttns);
  edp_l2pre_tmp = interp1(edp_t_norm, edp_l2pre_ts.data, edi_t_norm, 'nearest');

  % Scatter plot to locate differances
  figure;
  % Note: Only DSL X for now...
  h = scatter(edp_l2pre_tmp(:,1), edi_e_dsl_ts.data(:,1));
  % Add axis labels
  ylabel('EDI E-field DSL X');
  xlabel('EDP E-field DSL X (resamp nearest)');
  hold on; plot([-10 10],[-10 10],'--black'); % What would be perfect aligned data
  plot( [-10 10]-limit.outlier, [-10 10], '--red', ...
    [-10 10]+limit.outlier, [-10 10], '--red'); % Outlier limits (for abs(EDP-EDI)>limit).

  if(false) %% Plot heat map.
    figure;
    xb = linspace(-15, 15, 50);
    yb = linspace(-10, 10, 50);
    n = hist3([edp_l2pre_tmp(:,1), edi_e_dsl_ts.data(:,1)], {xb, yb});
    h = pcolor(xb, yb, n');
    colormap(flipud(colormap(hot)));
    % "reverse hot", ie change so that zero = white, higher values red towards black.
    % Add axis labels
    ylabel('EDI E-field DSL X');
    xlabel('EDP E-field DSL X (resamp nearest)');
    hold on; plot([-10 10],[-10 10],'--black'); % What would be perfect aligned data
    plot( [-10 10]-limit.outlier, [-10 10], '--red', ...
      [-10 10]+limit.outlier, [-10 10], '--red'); % Outlier limits (for abs(EDP-EDI)>limit).
  end

  % Find index of outliers
  % ind = abs(edp_l2pre_tmp(:,1) - edi_e_dsl_ts.data(:,1)) > limit.outlier;
  % Plot only outliers (making sure it was as expected).
  %figure;
  %h2 = scatter(edp_l2pre_ts.data(ind,1),edi_resamp.data(ind,1));

  % Timestamps of outliers
  %time = edi_e_dsl_ts.time(ind)

  % Plot timeseries of outliers (to identify regions with a lot almost
  % continuous outliers)
  %figure; irf_plot(irf.ts_scalar(time, ones(size(time))),'*');
  %
  % outlier_ts = edi_e_dsl_ts(ind);

  % Plot both data and manually zoom in to the intervals with a lot of
  % outliers identified before..
  % Note: this takes a while...
  %figure; irf_plot({edp_l2pre_ts, edi_e_dsl_ts, outlier_ts}, 'comp', 'linestyle',{'-b','-r*','go'});

  % ind2 = ind & (edp_l2pre_tmp(:,1) > limit.edp_edi_zero1);
  % ind2 = ind2 & (abs(edi_e_dsl_ts.data(:,1)) < limit.edp_edi_zero2);
  %
  % outliers_zero_ts = edi_e_dsl_ts(ind2);
  %
  % time = outliers_zero_ts.time;
  %
  % figure;
  % h=irf_plot({edp_l2pre_ts, edi_e_dsl_ts, outlier_ts, outliers_zero_ts}, 'comp', 'linestyle',{'-b','-r*','go','blacko'});

end

% Compute EDP mean values
%[edp_mean, edp_raw] = compute_edp_mean(edi_e_dsl_ts, edp_l2pre_ts);
edp_mean = compute_edp_mean(edi_e_dsl_ts, edp_l2pre_ts); % Somewhat quicker.

if(false) %% Plot Timeseries of EDI, EDP mean values and EDP measurements used
  figure; %#ok<UNRCH>
  irf_plot({edi_e_dsl_ts.x, ...
    irf.ts_scalar(edi_e_dsl_ts.time, edp_mean), ...
    edp_raw}, ...
    'comp', ...
    'linestyle',{'r*','go','-b*'});
  legend({'EDI (measured over 5s)', 'EDP (mean over 5s)', 'EDP (measurement)'});
end

if(true) %% Scatterplot similar to the ones Ivan produced.
  figure;
  scatter(edp_mean, edi_e_dsl_ts.x.data);
  ylabel('EDI E-field DSL X');
  xlabel('EDP E-field DSL X (mean 5 sec interval)');
  hold on;
  plot([-10 10],[-10 10],'--black'); % What would be perfect aligned data
  plot( [-10 10]-limit.outlier, [-10 10], '--red', ...
    [-10 10]+limit.outlier, [-10 10], '--red');
end

if(false) %% Plot heat map.
  figure;%#ok<UNRCH>
  xb = linspace(-15, 15, 50);
  yb = linspace(-10, 10, 50);
  n = hist3([edp_mean, edi_e_dsl_ts.x.data], {xb, yb});
  h = pcolor(xb, yb, n');
  colormap(flipud(colormap(hot))); % heat map, but zero = white, higher values red towards black.
  % Add axis labels
  ylabel('EDI E-field DSL X');
  xlabel('EDP E-field DSL X (mean 5 sec interval)');
  hold on; plot([-10 10],[-10 10],'--black'); % What would be perfect aligned data
  plot( [-10 10]-limit.outlier, [-10 10], '--red', ...
    [-10 10]+limit.outlier, [-10 10], '--red');
end

%% Give manual control..
keyboard

% Example on further analysis
%ind = abs(edp_mean-edi_e_dsl_ts.x.data) > limit.outlier;
% Locate times to look more at:
%time = edi_e_dsl_ts.time(ind);



end


function [edp_mean, edp_raw] = compute_edp_mean(edi_ts, edp_ts, interval)
% Compute mean of EDP values from middle of edi_ts.time and with interval
% for instance [-2.5, 2.5] seconds.
%% THIS FUNCTION MUST BE IMPROVED!!

% Verify inputs
narginchk(2,3);
if(~isa(edi_ts,'TSeries') || ~isa(edp_ts,'TSeries')), error('TSeries input expected.'); end
% Default 5 seconds interval.
if(nargin<3), interval = [-2.5, 2.5]; end

%% FIXME: Using ONLY X for now
edp_ts = edp_ts.x;    edi_ts = edi_ts.x;

% Pre allocate output
edp_mean = NaN(length(edi_ts), 1);

timeStamp = edi_ts.time;
interval = int64(interval*10^9); % Convert second to ns.
if(nargout==2), edp_raw = edp_ts(1); end

for ii=1:length(edi_ts)
  % Check if EDI is valid, otherwise NaN in, NaN out.
  if(~isnan(edi_ts(ii).data))
    edp_segment = edp_ts.tlim( irf.tint( timeStamp(ii) + interval) );
    % Compute the mean, ignoring NaN.
    edp_mean(ii) = irf.nanmean(edp_segment.data);
    if(nargout==2)% Keep the data used, if TSeries is to be plotted...
      edp_raw = combine(edp_raw, edp_segment);
    end
  end
end

end