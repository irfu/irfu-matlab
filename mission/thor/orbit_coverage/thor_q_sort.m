function tsQnew = thor_q_sort(tsQ,dT)
% THOR_Q_SORT assign the maximum Q from tsQ for all points in each time step dT
% tsQnew = THOR_Q_SORT(tsQ,dT)
%		tsQ     - Timeseries with new quality factors where in 
%		dT      - time step bins used for taking max(Q), tsQ.time.start:dT:tsQ.time.stop;
%   tsQnew  - output Timeseries

% Define new time axis and to whic point in new axis each of the input time
% points correspond to
tNew       = tsQ.time.start:dT:tsQ.time.stop;
QNewValues = nan(size(tNew));
indBin     = floor((tsQ.time-tsQ.time.start)/dT)+1;

% find Q values for each point in new time series
Q = tsQ.data;
for j = 1:numel(Q)
	QNewValues(indBin(j)) = max( QNewValues(indBin(j)) , Q(j) );
end

% assign new Q values
tsQnew = tsQ;
tsQnew.data = QNewValues(indBin);
