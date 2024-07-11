%
% Create time interval string. Function supports
% solo.adm.dsfn.create_dataset_filename().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function timeIntervalStr = create_time_interval_str(dateVec1, dateVec2, timeIntervalFormat)
    % PROPOSAL: Automatic test code.

if strcmp(timeIntervalFormat, 'DAY')
  assert(isequal(dateVec1(4:6), [0,0,0]))
  assert(isequal(dateVec2(4:6), [0,0,0]))
  expDateVec2 = datevec(datetime(dateVec1) + caldays(1));
  assert(isequal(dateVec2, expDateVec2))

  timeIntervalStr = sprintf(...
    '%04i%02i%02i', ...
    dateVec1(1:3));

elseif strcmp(timeIntervalFormat, 'DAY_TO_DAY')
  assert(isequal(dateVec1(4:6), [0,0,0]))
  assert(isequal(dateVec2(4:6), [0,0,0]))
  dateVec2 = datevec(datetime(dateVec2) - caldays(1));

  timeIntervalStr = sprintf(...
    '%04i%02i%02i-%04i%02i%02i', ...
    dateVec1(1:3), dateVec2(1:3));

elseif strcmp(timeIntervalFormat, 'SECOND_TO_SECOND')
  timeIntervalStr = sprintf(...
    '%04i%02i%02iT%02i%02i%02i-%04i%02i%02iT%02i%02i%02i', ...
    dateVec1, dateVec2);

elseif strcmp(timeIntervalFormat, 'NO_TIME_INTERVAL')
  assert(isequal(dateVec1, [0,0,0,0,0,0]))
  assert(isequal(dateVec2, [0,0,0,0,0,0]))
  timeIntervalStr = [];

else
  error('Can not interpret timeIntervalFormat="%s".', timeIntervalFormat)
end

end
