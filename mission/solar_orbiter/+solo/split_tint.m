function times_out = split_tint(tint,split_times)
% SOLO.SPLIT_TINT  Splits an input time interval (tint) according to
% the times listed in split_times.
% times_out = solo.split_tint(tint,split_times)
%
% Inputs: tint - EpochTT object : 2 records (e.g. from irf.tint('',''))
%         split_times - EpochTT object : N>0 records (e.g. from
%         solo.ProbePotDiscontinuities)
%
% Outputs:
%   times_out - time in epoch (int64) format.
%
% Example:
%   tint=irf.tint('2020-01-01T00:00:00.00Z','2020-01-02T00:00:00.00Z');
%   split_times = EpochTT(['2020-01-01T10:00:00.00Z';'2020-01-01T13:00:00.00Z';...
%   '2020-01-02T01:00:00.00Z';]);
%   times_out=EpochTT(solo.split_tint(tint,split_times));
%
%   gives output:
%   times_out =
%   EpochTT object : 4 records
%   2020-01-01T00:00:00.000000000Z
%   2020-01-01T10:00:00.000000000Z
%   2020-01-01T13:00:00.000000000Z
%   2020-01-02T00:00:00.000000000Z

times_epoch = tint.epoch;
split_epoch = split_times.epoch;
split_flag=zeros(size(split_epoch));
for ii=1:length(split_epoch)
  split_flag(ii) = logical((split_epoch(ii)>times_epoch(1))*...
    (split_epoch(ii)<times_epoch(2)));
end
use_splits = 1:length(split_flag);
use_splits(~split_flag)=[];

t=[times_epoch;split_epoch(use_splits)];
times_out = sort(t,'ascend');

end