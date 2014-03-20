function phase = mms_sdc_sdp_phase( DataTimeStamp, SunpulseTimeStampIn )
% MMS_SDC_SDP_PHASE calculates the phase of each point in time in the DataTimeStamp.
%	phase = MMS_SDC_SDP_PHASE(DataTimeStamp, SunpulseTimeStamp) return the
%	phase calculated with linear extrapolation for each time in
%	DataTimeStamp assuming each SunpulseTimeStamp corresponds to times with 
%   phase 0 deg. 
%
%	Example:
%		phase = mms_sdc_sdp_phase(DataTimeStamp, SunpulseTimeStamp);
%
% 	See also MMS_SDC_SDP_DATAMANAGER.

narginchk(2,2);

global MMS_CONST;

% Nominal spinrate is expected to be 3 rpm +/- 0.2 rpm.
 % One pulse every 18'750'000'000 ns, = 60/(3+0.2).
sunpulse_limit_low  = int64(60*10^9/MMS_CONST.Spinrate.max);
 % One pulse every 21'428'571'429 ns, = 60/(3-0.2).
sunpulse_limit_high = int64(60*10^9/MMS_CONST.Spinrate.min);

% Get each new sunpulse. (SunpulseTimeStampIn comes from CDF with HK, the
% rate of HK is different from Sunpulse rate, some pulses are repeated in
% next HK packet).
sunpulseTimeStamp = unique(SunpulseTimeStampIn);

if( (all( diff(sunpulseTimeStamp) < sunpulse_limit_high ) ) && ...
        (all( diff(sunpulseTimeStamp) > sunpulse_limit_low) ) )
    % Ok, sunpulse times are as expected.
    % Reduce the size by subtracting the time corresponding to the first
    % datapacket, then transform it to double.
    % (interp() - can't handle int64). 
    sunpulseTimeStamp = double( sunpulseTimeStamp - sunpulseTimeStamp(1) );
    DataTimeStamp = double( DataTimeStamp -sunpulseTimeStamp(1) );

    % For each new sunpulse, add 360 degrees to the total phase. or should
    % it be radian 0 -> 2*pi?
    for i=1:length(sunpulseTimeStamp)
         % Each new pulse equals to an extra 360 deg phase.
        sunpulseTimeStamp(i,2) = (i-1)*360;
    end

    % Make a phase assumption as linear progression between each sunpulse
    % recieved. Sunpulse may be a real sunpulse or a CIDP psuedo pulse but
    % for SITL / Quicklook this should be just fine anyhow.

    phase = interp1( sunpulseTimeStamp(:,1), ...
        sunpulseTimeStamp(:,2), DataTimeStamp, ...
        'linear', 'extrap' );

    % Wrap it to the default 0-360 degrees or should it be radian 0-2*pi?
    % Or should it be continous?
    phase = wrapTo360(phase);
    
else
    % Not ok. Some sunpulses are extra or missing...
    err_str = ['MMS_SDC_SDP_PHASE timesteps between each unique ', ...
        'sunpulse are too large or too small to be consistent with ', ...
        'spinrates between ', MMS_CONST.Spinrate.max, ' rpm and ', ...
        MMS_CONST.Spinrate.min, ' rpm.'];
    irf.log('critical', err_str);
    error('MATLAB:mms_sdc_sdp_phase', err_str);
end


function lon = wrapTo360(lon)
% Wrap positive angles to the [0 - 360] interval.
positiveInput = (lon > 0);
lon = mod(lon, 360);
lon((lon == 0) & positiveInput) = 360;
