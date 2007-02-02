function [data, n_corrected,wakedesc] = c_efw_swwake(e,pair,phase_2,plotflag)
%C_EFW_SWWAKE  Correct raw EFW E data for wake in the solar wind
%
% [data, n_spins_corrected,w_dsc] = c_efw_swwake(e,pair,phase_2 [,plotflag])
%
% Input:
%   e        - raw EFW data (wE?p12/34)
%   pair     - probe pair (12/34/32)
%   phase_2  - Phase 2
%   plotflag - 0 no plotting of debug 
%
% Output:
%   data - corrected data
%   
% Wakes are identified by max derivative. First we find a narrow proxy
% wake, correct for it and find a proxy DC field (ground tone). Then
% we subtract the ground tone from the original data and find a final
% fit for the wake. The procedure is performed on five spins, with the
% resulting fit being applied to the spin in the middle.
%
% This program was written in order to improve quality of the EFW data
% in the solar wind for the CAA.
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev
% Original idea by Anders Eriksson.
% Many useful suggestion by Per-Arne Lindqvist.

error(nargchk(3,4,nargin))
if nargin==4
	if isnumeric(plotflag)
		if plotflag~=0, plotflag = 1; end
	else plotflag = 1;
	end
else plotflag = 0;
end

n_corrected = 0;

if pair~=12 && pair~=32 && pair~=34, error('PAIR must be one of: 12, 32, 34'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

if pair==32, error('PAIR 32 is not implemented yet'), end

% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .9; 
MAX_SPIN_PERIOD = 4.3;
WAKE_MAX_HALFWIDTH = 45; % degrees
WAKE_MIN_AMPLITUDE = 0.7; % mV/m

data = e;

i0 = find(phase_2(:,2)==0);
if isempty(i0), irf_log('proc','empty phase'), return, end
if length(i0)<5, irf_log('proc','not enough spins for correction'), return, end

% Hack to add phase=360 at the end and phase=0 at the start of the interval
% in order to try to correct these incomplete spins
if phase_2(end,2)~=360 && ...
		phase_2(i0(end),1)-phase_2(i0(end-1),1)<MAX_SPIN_PERIOD
	phase_2 = [phase_2; phase_2(i0(end),1)*2-phase_2(i0(end-1),1) 0];
end
if phase_2(1,2)~=0 && phase_2(i0(2),1)-phase_2(i0(1),1)<MAX_SPIN_PERIOD
	phase_2 = [phase_2(i0(1),1)*2-phase_2(i0(2),1) 0; phase_2];
	i0 = find(phase_2(:,2)==0);
end

n_spins = length(i0);
NPOINTS = 361;
tt = zeros(NPOINTS,n_spins);
wakedesc = zeros(n_spins*2, 4)*NaN;
ttime = tt;
sf = c_efw_fsample(e,'hx');
iok = [];

for in = 1:n_spins
	ts = phase_2(i0(in),1);
	i360 = find( phase_2(:,1)>ts & phase_2(:,1)<ts+MAX_SPIN_PERIOD & ...
		phase_2(:,2)==360);
	if isempty(i360)
		irf_log('proc',['gap in phase at ' epoch2iso(ts,1)])
		te = ts + 4.0;
		empty = 1;
	else
		if length(i360)~=1, irf_log('proc',['bogus phase at ' epoch2iso(ts,1)]), end
		te = phase_2(i360(end),1);
		empty = 0;
	end
	
	ttime(:,in) = (ts + (0:1:360) *(te-ts)/360.0)';
	if empty, tt(:,in) = NaN;
	else
		eind = find((e(:,1) > ts-0.15) & (e(:,1) < te+0.15));
		eind(isnan(e(eind,2))) = [];
		% Check for data gaps inside one spin.
		if sf>0 && length(eind)<N_EMPTY*MAX_SPIN_PERIOD*sf
			irf_log('proc',['data gap at ' epoch2iso(ts,1)])
			tt(:,in) = NaN;
		else
			if sf==450, dtmp = c_resamp(e(eind,:), ttime(:,in));
			else dtmp = c_resamp(e(eind,:), ttime(:,in), 'spline');
			end
			% Fill small gaps (at edges only?) with zeroes
			% This has a minor influence on the correction procedure
			dtmp(isnan(dtmp(:,2)),2) = 0;
			tt(:,in) = dtmp(:,2);
		end
		
		% Identify spins for which we attempt to correct wake
		if (in>=6 && sum(isnan(tt(1,in - (0:1:5) )))<=1)
			% We allow max one data gap within 6 spins
			iok = [iok in-3];
		end
		if (in==n_spins || n_spins==5) && ~any(isnan(tt(1,in - (0:1:4) )))
			% Special case when we have only 5 spins or
			% it is the last spin we are working with
			iok = [iok in-2];
		end
	end
end

if isempty(iok)
	irf_log('proc','not enough spins for correction')
	return
end

i1 = [];

for in = iok
	ts = ttime(1,in);
	if in==n_spins-2
		% Last interval
		idx = -2:1:2;
	else
		% Check for a data gap
		nans = isnan(tt(1,in + (-2:1:3) ));
		if any(nans)
			idx = 1:6;
			idx = idx(xor(idx,nans)) -3;
		else idx = -2:1:2;
		end
	end

	% Spin in the middle has maximum weigth
	av12 = sum(tt(:, in + (idx) ).*([.1 .25 .3 .25 .1]'*ones(1,NPOINTS))',2);

	% First find a proxy wake fit	
	% Identify wakes by max second derivative
	d12 = [av12(1)-av12(end); diff(av12)];
	d12 = [d12(1)-d12(end); diff(d12)];
	% Average with 7 points to minimize danger of detecting a wrong maximum
	d12 = w_ave(d12,7);
	
	% If i1 is coming not from the previous spin
	if isempty(i1) || (in~=iok(1) && in-1~=iok(find(iok==in)-1))
		i1 = 1:1:361;
	end
	
	% We expect for the second maximum to be 180 degrees from the first
	if max(d12(i1))>abs(min(d12(i1)))
		ind1 = find(d12 == max(d12(i1))) -1;
		i1 = mod( (ind1-WAKE_MAX_HALFWIDTH+180:ind1+WAKE_MAX_HALFWIDTH+180) -1,...
			NPOINTS) +1;
		ind2 = find(d12 == min(d12( i1 ))) -1;
	else
		ind1 = find(d12 == min(d12(i1))) -1;
		i1 = mod( (ind1-WAKE_MAX_HALFWIDTH+180:ind1+WAKE_MAX_HALFWIDTH+180) -1,...
			NPOINTS) +1;
		ind2 = find(d12 == max(d12( i1 ))) -1;
	end

	% The proxy wake is naroow (1/2 of the final fit)	
	wake_width = fix(WAKE_MAX_HALFWIDTH/2);
	i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
	i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;

	% The proxy wake is symmetric	
	dav = (d12(i1)-d12(i2))/2;
	cdav = cumsum(dav);
	cdav = cdav - mean(cdav);
	ccdav = cumsum(cdav);
	
	% Save wake description
	fw = (mod(ind2,NPOINTS)+1<mod(ind1,NPOINTS)+1);
	wakedesc(in*2-1+fw,1) = ttime(mod(ind1,NPOINTS)+1,in);
	wakedesc(in*2-1+fw,2) = ind1;
	wakedesc(in*2-fw,1) = ttime(mod(ind2,NPOINTS)+1,in);
	wakedesc(in*2-fw,2) = ind2;
	wakedesc([in*2-1 in*2],3) = max(abs(ccdav));
	% Wake half-width
	ii = find(abs(ccdav)<max(abs(ccdav))/2);
	wakedesc([in*2-1 in*2],4) = min(ii(ii>23))-max(ii(ii<23));
	if wakedesc(in*2,3)< WAKE_MIN_AMPLITUDE
		%irf_log('proc',['wake is too small at ' epoch2iso(ts,1)])
		i1 = [];
		continue
	end
	
	wake = zeros(NPOINTS,1);
	wake( i1 ) = ccdav;
	wake( i2 ) = -ccdav;

	% Correct for the proxy wake	
	av12_corr = av12 - wake;
	% Find the ground tone and remove it from the data
	x = fft(av12_corr);
	x(3:359) = 0;
	av12_corr = av12 -ifft(x,'symmetric');

	% Now find the final fit	
	d12 = [av12_corr(1)-av12_corr(end); diff(av12_corr)];
	d12 = [d12(1)-d12(end); diff(d12)];
	if plotflag
		d12_tmp = d12; % save for plotting
	end
	% Average with only 5 points to get a more fne fit
	d12 = w_ave(d12,5);
	
	wake_width = WAKE_MAX_HALFWIDTH;
	i1 = mod( (ind1-wake_width:ind1+wake_width) -1, NPOINTS) +1;
	i2 = mod( (ind2-wake_width:ind2+wake_width) -1, NPOINTS) +1;
	
	% Allow the final fit to be asymmetric
	cdav = cumsum(d12(i1));
	cdav = cdav - mean(cdav);
	ccdav1 = cumsum(cdav);
	cdav = cumsum(d12(i2));
	cdav = cdav - mean(cdav);
	ccdav2 = cumsum(cdav);
	
	% Save wake description
	wakedesc(in*2-1+fw,3) = max(abs(ccdav1));
	% Wake half-width
	ii = find(abs(ccdav1)<max(abs(ccdav1))/2);
	wakedesc(in*2-1+fw,4) = min(ii(ii>45))-max(ii(ii<45));
	wakedesc(in*2-fw,3) = max(abs(ccdav2));
	% Wake half-width
	ii = find(abs(ccdav2)<max(abs(ccdav2))/2);
	%disp(ii)
	%if in==244,keyboard,end
	wakedesc(in*2-fw,4) = min(ii(ii>45))-max(ii(ii<45));
	clear ii
	
	if max(max(abs(ccdav1)),max(abs(ccdav2)))< WAKE_MIN_AMPLITUDE
		%irf_log('proc',['wake is too small at ' epoch2iso(ts,1)])
		i1 = [];
		continue
	end
	
	wake = zeros(NPOINTS,1);
	wake( i1 ) = ccdav1;
	wake( i2 ) = ccdav2;
	
	if plotflag
		clf
		subplot(4,1,1)
		ts = ttime(1,in);
		te = ttime(end,in);
		plot(ttime(:,in)-ts, tt(:, in), 'b',...
			ttime(:,in)-ts, tt(:, in + ([-2 -1 1 2]) ), 'g',...
			ttime(:,in)-ts, av12, 'k',...
			ttime(ind1,in)*[1 1]-ts, [-2 2], 'r',...
			ttime(ind2,in)*[1 1]-ts, [-2 2], 'r',...
			ttime(:,in)-ts, av12-wake,'r');
		ylabel('E12 [mV/m]');
		add_timeaxis(gca,ts); xlabel('');
		set(gca,'XLim',[0 te-ts])
		
		subplot(4,1,2)
		plot(ttime(:,in)-ts,d12_tmp,'g',ttime(:,in)-ts, d12,'b');
		ylabel(['D2(E' num2str(pair) ') [mV/m]']);
		add_timeaxis(gca,ts); xlabel('');
		set(gca,'XLim',[0 te-ts])

		subplot(4,1,3)
		plot(ttime(:,in)-ts, wake)
		ylabel('Wake [mV/m]');
		add_timeaxis(gca,ts);
		set(gca,'XLim',[0 te-ts])
	end

	% Correct the spin in the middle	
	ind = find(e(:,1)>=ttime(1,in) & e(:,1)<ttime(end,in));
	wake_e = c_resamp([ttime(:,in) wake], e(ind,1));
	data(ind,2) = data(ind,2) - wake_e(:,2);
	n_corrected = n_corrected + 1;
	
	if plotflag
		subplot(4,1,4)
		irf_plot({e(ind,:),data(ind,:)},'comp')
		ylabel(['E' num2str(pair) ' [mV/m]']);
		irf_zoom([ts te],'x',gca)
	end

	% Correct edge spins	
	cox = [];
	if in==iok(1) || (in~=iok(1) && in-1~=iok(find(iok==in)-1))
		% If the previous spin was not corrected
		% we correct it here
		if in==iok(1) && in>3
			% We try to correct the spin at the start of the entire 
			% intrval which usually has a data gap 
			cox = - (1:3);
			n_corrected = n_corrected + 3;
		else
			cox = - (1:2);
			n_corrected = n_corrected + 2;
		end
	end
	if in==iok(end) || (in~=iok(end) && in+1~=iok(find(iok==in)+1))
		if ~isempty(cox)
			irf_log('proc',['single interval at ' epoch2iso(ts,1)])
		else
			irf_log('proc',['stop   interval at ' epoch2iso(ts,1)])
		end
		if in==iok(end) && n_spins-in>=3
			% We try to correct the spin at the end of the entire 
			% intrval which usually has a data gap 
			cox = [cox 1:3];
			n_corrected = n_corrected + 3;
		else
			cox = [cox 1:2];
			n_corrected = n_corrected + 2;
		end
	elseif ~isempty(cox)
		irf_log('proc',['start  interval at ' epoch2iso(ts,1)])
	end
	if ~isempty(cox)
		for cx = cox
			ind = find(e(:,1)>=ttime(1,in+cx) & e(:,1)<ttime(end,in+cx));
			if ~isempty(ind)
				wake_e = c_resamp([ttime(:,in+cx) wake], e(ind,1));
				data(ind,2) = data(ind,2) - wake_e(:,2);
				%irf_log('proc',['correcting spin: ' num2str(cx)])
				if plotflag
					hold on
					irf_plot({e(ind,:),data(ind,:)},'comp')
					ylabel(['E' num2str(pair) ' [mV/m]']);
					hold off
					if ttime(1,in+cx)<ts, ts = ttime(1,in+cx); end
					if ttime(end,in+cx)>te, te = ttime(end,in+cx); end
					irf_zoom([ts te],'x',gca)
				end
			end
		end
	end
	
end

% Save wake position only inside 0-180 degrees
wakedesc(wakedesc(:,2)>180,2) = wakedesc(wakedesc(:,2)>180,2)-180;

irf_log('proc',['corrected ' num2str(n_corrected) ' out of ' ...
	num2str(n_spins) ' spins'])

function av = w_ave(x,np)
% Weighted average
NPOINTS = 361;
if nargin<2, np=5; end
av = zeros(size(x));
if np==7
	m = [.07 .15 .18 .2 .18 .15 .07]';
	idx = -3:1:3;
else
	m = [.1 .25 .3 .25 .1]';
	idx = -2:1:2;
end
MIDX = max(idx);

for j=1:length(x)
	ii = j + (idx);
	if j<=MIDX, ii(ii<1) = ii(ii<1) +NPOINTS; end
	if j>length(x)-MIDX, ii(ii>NPOINTS) = ii(ii>NPOINTS) -NPOINTS; end
	av(j) = sum(x(ii).*m);
end
