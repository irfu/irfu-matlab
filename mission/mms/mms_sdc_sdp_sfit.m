function spinfit = mms_sdc_sdp_sfit(pair, fout, maxit, minpts, te, data, tp, ph, method, tmmode)
%MMS_SDC_SDP_SFIT  Produce spin fit data values (ex, ey)
%   of MMS data frome given probe pair at 20 second resolution.
%
% SPINFIT = MMS_SDC_SDP_SFIT(PAIR,FOUT,MAXIT,MINPTS,TE,DATA,TP,PH,METHOD,[TMMODE])
%
% Input:
%  PAIR - probe pair used (12, 32, 34)
%  FOUT - minimum fraction of fit standard deviation that defines an outlier 
%         (zero means no removal of outliers). Has no effect if METHOD=1
%  MAXIT - maximum number of iterations (zero means infinity)
%  MINPTS - minimum number of data points to perform fit
%      (set to 5 if smaller number if given)
%  TE - EFW time in seconds (isGetDataLite time)
%  DATA - EFW data from pair in mV/m, should correspond to te
%  TP - Ephemeris time in seconds (isGetDataLite time)
%  PH - Ephemeris phase_2 in degr (sun angle for s/c Y axis), should 
%      correspond to tp 
%  METHOD - 0: mms_sdc_sdp_onesfit, based on Matlab routine for Cluster, by AIE
%           1: mms_sdc_sdp_spinfit_mx (default), based on routine for Cluster, BHN Fortran source obtained from KTH
%  TMMODE - 'hx' (default) or 'ib', used for determination of sampling freq
%
% Output:
%  SPINFIT = [TS,EX,EY,OFFSET,SDEV0,SDEV,ITER,NOUT]
%  TS - time vector in seconds <- FIXME: NANOSECONDS. 
%  EX - E-field x-component in DSI coordinates (almost GSE)
%  EY - E-field y-component in DSI coordinates (almost GSE)
%  OFFSET - mean value of input data
%  SDEV0 - standard deviation in first fit. Has no meaning if METHOD=1
%  SDEV - standard deviation in final fit
%  ITER - number of iterations (one if OK at once)
%  NOUT - number of outliers removed
%
% This function chops up the time series in 20 second
% intervals, each of which are analysed by mms_sdc_sdp_onesfit.
% Spins with less then 90% of data are disregarded.
%
% See also MMS_SDC_SDP_PHASE
%
% $Id$

% Based on original version for Cluster (c_efw_sfit.m) by 
% Anders.Eriksson@irfu.se 13 December 2002

narginchk(8,10)

global MMS_CONST

if( pair~=12 && pair~=32 && pair~=34 && pair~=42 )
    err_str = 'PAIR must be one of: 12, 34, 32, 42';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_SFIT:INPUT', err_str);
end

if ~isequal(size(te), size(data))
    err_str = 'TE and DATA vectors must be the same size.';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_SFIT:INPUT', err_str);
end
if ~isequal(size(tp), size(ph))
    err_str = 'TP and PH vectors must be the same size.';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_SFIT:INPUT', err_str);
end

if( nargin<10 )
    tmmode='hx';
end

% Check if the phase is good
[ie, ip] = irf_find_comm_idx(te, tp);
if( isempty(ie) || isempty(ip) )
	% Try to compute phase
	irf.log('warning','MMS_SDC_SDP_SFIT phase appeared inconsistent with data. Recalculating.')
% 	aa = c_phase(te,[tp ph]);
%   % Check if new phase is good.
% 	[ie,ip] = irf_find_comm_idx(te, aa(:,1));
% 	if( isempty(ie) || isempty(ip) )
% 		irf.log('critical','Phase and/or E is empty.')
% 		spinfit = [];
% 		return
% 	else
% 		tp = aa(:,1);
% 		ph = aa(:,2);
% 	end
end

te = te(ie); data = data(ie);
tp = tp(ip); ph = ph(ip);
clear ie ip

% Convert from int64 to double and seconds.
te = double(te)/10^9;

% Set default method to BHN
if( nargin < 9) 
    method = 1;
end

if( method==1 || method==2 )
	if( ~exist(['mms_sdc_sdp_spinfit_mx.', mexext],'file') )
        % Provide warning to users logfile to ensure they know.
        irf.log('warning',['MMS_SDC_SDP_SPINFIT_MX.', mexext,...
            ' file does not exist. Building it.']);
        % Locate the path to source and output directory.
        pathToSrc = which('mms_sdc_sdp_spinfit_mx.F');
        [s, ~, ~] = fileparts(pathToSrc);

        try
            mex(['FFLAGS=$FFLAGS -O2 -fPIC -mtune=opteron -funroll-loops'],...
                '-outdir',[s,filesep], [s,filesep,'mms_sdc_sdp_spinfit_mx.F']);
        catch err
            % If running at SDC, there will be NO GCC installed and this 
            % will fail. MEXA64 to be built at IRFU or on dedicated machine
            % at SDC after filing a ticket for it at SDC.
            irf.log('critical', 'Mex failed. Are we running at SDC?');
            irf.log('critical', err.message);
            
            irf.log('critical', 'Returning to Matlab script (method 0) instead.');
            method = 0;
        end
	end
end


% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .9;
% Sane values should have periods between max and min.
MAX_SPIN_PERIOD = 60/MMS_CONST.Spinrate.min; % Approx 21.43 sec (nominally)
MIN_SPIN_PERIOD = 60/MMS_CONST.Spinrate.max; % Approx 18.75 sec (nominally)

% Guess the sampling frequency
%sf = c_efw_fsample(te, tmmode);
%THONI: FIXME: static sampling frequency for test purpose.
sf = 31.99963;

if( sf == 0 )
    error('cannot guess the sampling frequency')
end

tpha = tocolumn(tp);
pha = tocolumn(ph);
% Convert phase to rad
pha = unwrap(pi*pha/180);

% FIXME: THONI check which angle each probe pair actually has on MMS.
% Find phase of given pair:
% if( pair == 12)
%     pha = pha + 3*pi/4;
% elseif( pair == 32 )
%     pha = pha + pi/2;
% elseif( pair == 34 )
%     pha = pha + pi/4;
% elseif( pair == 42 )
%     pha = pha + pi;
% else
%     err_str = 'PAIR must be one of: 12, 34, 32, 42';
%     irf.log('critical', err_str);
%     error('MATLAB:MMS_SDC_SDP_SFIT:INPUT', err_str);
% end


% Do it:
if( method==1 || method==2 )
	te = torow(te);
	data = torow(data);
	ind = find(~isnan(data));
	
    % Check if we have at least one spin of data
    if( length(ind) < N_EMPTY * 20 * sf )
        irf.log('critical', 'Not enough data points for spinfit. Return empty.');
        spinfit = [];
        return
    end
    
    if( method==2 )
        NCOLS = 10;
        NTERM = 5;
    else
        % Method == 1.
        NCOLS = 8;
        NTERM = 3;
    end
    
	pha = torow( pha(ind) );
	
    [ts, sfit, sdev, iter, nout] = mms_sdc_sdp_spinfit_mx(maxit,...
        N_EMPTY*20*sf,NTERM,te(ind),data(ind),pha);
	
    n = length(sdev);
    
    % Locate any fillvalues
    ind = find( sdev~=-159e7 ); % FillVal set in compiled file to this value.
    n_gap = n - length(ind);
    
    % Preallocate output
    spinfit = zeros(n,NCOLS);
	
    spinfit(:, 1) = ts;		        % time
	spinfit(:, 2:end) = NaN;
    
    if ~isempty(ind)
		spinfit(ind, 2) =  sfit(2, ind);	% Ex
		spinfit(ind, 3) = -sfit(3, ind);% Ey, - Because s/c is spinning upside down
		spinfit(ind, 4) =  sfit(1, ind);
		spinfit(ind, 5) =  sdev(ind);
		spinfit(ind, 6) =  sdev(ind);
		spinfit(ind, 7) =  iter(ind);
		spinfit(ind, 8) =  nout(ind);
        if( method==2 )
			spinfit(ind, 9)  = sfit(4, ind); % 2 omega signals
			spinfit(ind, 10) = sfit(5, ind);
        end
    end
    
else
	% Turn off warnings for badly conditioned polynomial:
	warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
	
	% Chop up time interval
	% We always start at 0,4,8.. secs, so that we have
	% the same timelines an all SC at 2,6,10... sec
	tstart = fix(min(te)/20)*20;
	tend = max(te);

	n = floor((tend - tstart)/20) + 1;
	
    % Preallocate output
    spinfit = zeros(n,8);
	
	n_gap = 0;
    for i=1:n
		tsfit = tstart + (i-1)*20 +10;
		ind = find( ( te >= tsfit-10 ) & ( te < tsfit+10 ) );

		% Check for data gaps inside one spin.
		if( sf > 0 && length(ind) < N_EMPTY * MAX_SPIN_PERIOD * sf )
			irf.log('notice','Data gap inside one spin.');
            n_gap = n_gap + 1;
			continue
		end

		% Compute spin period
		pol = polyfit(tpha(ind), pha(ind), 1);
		spinp = 2*pi / pol(1);
		if( spinp > MAX_SPIN_PERIOD || spinp < MIN_SPIN_PERIOD )
			irf.log('warning',sprintf('Bad spin period %.4f s.',spinp));
			n_gap = n_gap + 1;
			continue
		end

		% Clear NaNs
		ind(isnan(data(ind))) = [];

		% Check for data gaps inside one spin.
        if( sf > 0 && length(ind) < N_EMPTY * sf * 20 )
            irf.log('notice','Data gap inside one spin, second check.');
            % SHOULD THIS NOT ALSO HAVE n_gap = n_gap + 1; ??? <- THONI,
            % 20140325.
            continue
        end

		% Use Matlab version by AIE
		spinfit(i - n_gap,:) = c_efw_onesfit( pair, fout, maxit, ...
            minpts, te(ind), data(ind), te(ind), ph(ind) );
		spinfit(i - n_gap, 1) = tsfit;
    end
    
    spinfit = spinfit(1:n - n_gap, :);
	
    % Turn warnings back on
	warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
end

irf.log('notice', sprintf('%d spins processed, %d gaps found',n,n_gap) );
