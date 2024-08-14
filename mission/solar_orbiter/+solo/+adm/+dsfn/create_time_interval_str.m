%
% Create time interval string. The function supports
% solo.adm.dsfn.create_dataset_filename().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function timeIntervalStr = create_time_interval_str(Dt1, Dt2, timeIntervalFormat)
irf.dt.assert_UTC(Dt1)
irf.dt.assert_UTC(Dt2)

if strcmp(timeIntervalFormat, 'DAY')
  irf.dt.assert_UTC_midnight(Dt1)
  irf.dt.assert_UTC_midnight(Dt2)
  expDt2 = datetime(Dt1) + caldays(1);
  assert(isequal(Dt2, expDt2))

  dateVec1 = datevec(Dt1);

  timeIntervalStr = sprintf(...
    '%04i%02i%02i', ...
    dateVec1(1:3));

elseif strcmp(timeIntervalFormat, 'DAY_TO_DAY')
  irf.dt.assert_UTC_midnight(Dt1)
  irf.dt.assert_UTC_midnight(Dt2)

  Dt2 = Dt2 - caldays(1);   % Decrease by one day, since end day is inclusive.

  dateVec1 = datevec(Dt1);
  dateVec2 = datevec(Dt2);

  timeIntervalStr = sprintf(...
    '%04i%02i%02i-%04i%02i%02i', ...
    dateVec1(1:3), dateVec2(1:3));

elseif strcmp(timeIntervalFormat, 'SECOND_TO_SECOND')
  timeIntervalStr = sprintf(...
    '%04i%02i%02iT%02i%02i%02i-%04i%02i%02iT%02i%02i%02i', ...
    datevec(Dt1), datevec(Dt2));

elseif strcmp(timeIntervalFormat, 'NO_TIME_INTERVAL')
  NAT = datetime('NaT', 'TimeZone', 'UTCLeapSeconds');
  assert(isequaln(Dt1, NAT))
  assert(isequaln(Dt2, NAT))

  timeIntervalStr = [];

else

  error('Can not interpret timeIntervalFormat="%s".', timeIntervalFormat)
end

end
