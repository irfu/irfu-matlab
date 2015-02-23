function [phase_out,fits] = c_phase(t,phase_2)
%C_PHASE  Find spacecraft phase for given time vector
%
% [phase_out,fits] = c_phase(t,phase_2)
%
% Find spacecraft phase for give time vector t
%
% Input:
%       t - column vector with time in isdat epoch
%       phase_2 - column vector [time phase_2]
% Output:
%       phase_out - column vector [t phase]
%       fits - structure providing details of linear fits

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if size(t,1)>1 && size(t,2)>1, error('t must be a vector'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

MAX_SPIN_PERIOD = 4.3; % Sane values for Cluster
MIN_SPIN_PERIOD = 3.6;
MAX_ERR_PHA = 0.5;     %  Max error in phase [deg] we tolerate
MAX_SPINS_EXTRAP = 100;% Max number of spins we extrapolate the phase

t=t(:); % t should be column vector
tRef = t(1);
phase_out = [];

% Sanity check
badt=find(phase_2(:,1)<iso2epoch('2000-01-01T00:00:00Z'));
if ~isempty(badt)
	irf_log('proc',['Bad time ignored ' epoch2iso(phase_2(badt(1),1))])
	phase_2(badt,:) = [];
	clear badt;
end

% Check for gaps
pos = 1; idxS = {};
while pos < size(phase_2,1)
  % Search for gaps
  idx = pos:size(phase_2,1);
  ii = find( diff(phase_2(idx,1)) > MAX_SPIN_PERIOD );
  if isempty(ii), idxS = [idxS {pos:size(phase_2,1)}]; break, end %#ok<AGROW>
  irf_log('proc',...
      ['Gap in phase at ' epoch2iso(phase_2(idx(ii(1)),1),1)])
  if idx(ii(1))==pos
    irf_log('proc','throwing away 1 bad point'), phase_2(1,:) = [];
    continue
  end
  idxS = [idxS {pos:idx(ii(1))}]; pos = idx(ii(1))+1; %#ok<AGROW>
end

nChunks = length(idxS); fits = cell(nChunks,1); iPrevGoodChunk = 0;
for iChunk=1:nChunks
  idx = idxS{iChunk};
  if length(idx)==2 && abs(diff(phase_2(idx,2)))==360
    irf_log('proc',['Bad point at ' epoch2iso(phase_2(idx(1),1))])
    continue
  end
  lFit = fit_spin(phase_2(idx,:)); fits{iChunk} = lFit;
  if isempty(lFit.spinPeriod) || ...
      lFit.errAngle>MAX_ERR_PHA || lFit.errAngle>MAX_ERR_PHA
    irf_log('proc',...
      ['Spin rate changes at ' irf_disp_iso_range(phase_2(idx([1 end]),1)')])
    if ~isempty(lFit.errAngle)
      irf_log('proc',...
        sprintf('Phase errors > %.1f deg. err=%.1f deg, meanErr=%.1f deg',...
        MAX_ERR_PHA,lFit.errAngle,lFit.errAngleMean))
    end
    irf_log('proc','Using simple interpolation')
    tt = t(t>=phase_2(idx(1),1) & t<=phase_2(idx(end),1));
    if ~isempty(tt)
      phase_out = [phase_out; tt interp1q(phase_2(idx,1),phase_2(idx,2),tt)]; %#ok<AGROW>
    end
    continue
  end
  tStart = phase_2(idx(1),1);
  if iPrevGoodChunk==0
    tStart = tStart-MAX_SPINS_EXTRAP*lFit.spinPeriod;
  else
    % If the phase did not change, we extrapolate via the data gap
    if abs(lFit.spinPeriod-fits{iPrevGoodChunk}.spinPeriod)/4*360<MAX_ERR_PHA && ...
        abs(mod(polyval(lFit.phcCoef,phase_out(end,1)-tRef)*180/pi,360)-...
        phase_out(end,2))<MAX_ERR_PHA
      tStart = phase_2(idxS{iPrevGoodChunk}(end),1);
    end
  end
  tEnd = phase_2(idx(end),1);
  if iChunk==nChunks, tEnd = tEnd + MAX_SPINS_EXTRAP*lFit.spinPeriod; end 
  tt = t(t>=tStart & t<tEnd);
  phase_out = [phase_out; tt mod(polyval(lFit.phcCoef,tt-tRef)*180/pi,360)];  %#ok<AGROW>
  iPrevGoodChunk = iChunk;
end

% Sanity check
if ~isempty(phase_out())
  phase_unwrapped = unwrap(phase_out(:,2)/180*pi);
  SpinRate = diff(phase_unwrapped)./diff(phase_out(:,1));
  ii = find( diff(phase_out(:,1))< 0.95*2*pi/median(SpinRate) &...
    (SpinRate<2*pi/MAX_SPIN_PERIOD | SpinRate>2*pi/MIN_SPIN_PERIOD) );
  if ~isempty(ii)
    ii_jump = find(diff(ii')>1); ii_jump = [1 ii_jump];
    for i=1:length(ii_jump)
      if i==length(ii_jump), kk = [ii(ii_jump(i)) ii(end)];
      else kk = [ii(ii_jump(i)) ii(ii_jump(i+1))-1];
      end
      if kk(1,2)==0, kk(1,2)=1; end;
      irf_log('proc',['Bad phase at ' irf_disp_iso_range(phase_out(kk,1)',1)])
    end
    phase_out(ii,:) = [];
  end
end

n_out = length(t) - length(phase_out);
if n_out>0, irf_log('proc',['throwing away ' num2str(n_out) ' points']), end

  function res = fit_spin(phase)
    % res = fit_spin(PHASE_2)
    %
    % Compute Cluster spin period from PHASE_2 using linear fit
    res = struct('t',phase([1 end],1),'spinPeriod',[],...
      'errAngle',[],'errAngleMean',[],'phcCoef',[],'tRef',tRef);   
    ph = phase(:,2); tp = phase(:,1) - tRef;
    
    tTmp = (ceil(phase(1,1)*2):fix(phase(end,1)*2))'/2 - tRef; % new time with 0.5 sec step
    while 1
      tTmp(tTmp<tp(1)) = []; % Avoid NaNs
      if isempty(tTmp), return, end
      pTmp = interp1q(tp,ph,tTmp);
      phc = unwrap(pTmp/180*pi); res.phcCoef = polyfit(tTmp,phc,1);
      if isnan(res.phcCoef(1))
        irf_log('proc','Cannot determine spin period!'), return
      end
      diffangle = mod(ph - polyval(res.phcCoef,tp)*180/pi,360);
      diffangle = abs(diffangle);
      diffangle = min([diffangle';360-diffangle']);
      % Throw away erroneous points, we suppose these correspond to a single
      % error in sun pulse
      iJump = find(diffangle>5);
      if ~isempty(iJump) && length(iJump)<10
        if all(diff(iJump)==1),
          irf_log('proc', sprintf('Throwing away %d erroneous points at %s',...
            length(iJump),epoch2iso(tp(iJump(1))+tRef)))
        else
          for iOut=1:length(iJump)
            irf_log('proc', sprintf('Throwing away 1 erroneous points at %s',...
            epoch2iso(tp(iJump(iOut))+tRef)))
          end
        end
        ph(iJump) = []; tp(iJump) = [];
      else break % exit loop
      end
    end
    res.errAngleMean = mean(diffangle); res.errAngle = std(diffangle);
    res.spinPeriod = 2*pi/res.phcCoef(1);
  end % fit_spin
end



