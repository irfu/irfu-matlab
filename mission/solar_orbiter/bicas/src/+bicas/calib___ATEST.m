function calib___ATEST

get_calibration_time___ATEST
end



function get_calibration_time___ATEST
new_test = @(Epoch, CalibEpochList, outputs) (EJ_library.atest.CompareFuncResult(@bicas.calib.get_calibration_time, {int64(Epoch)', int64(CalibEpochList)'}, outputs));
tl = {};

EV = zeros(1,0);

tl{end+1} = new_test(EV,     [3],   {EV'});
tl{end+1} = new_test(EV,     [3,7], {EV'});

tl{end+1} = new_test([1:10],      []    , 'MException');
tl{end+1} = new_test([1:10],      [3]   , 'MException');
tl{end+1} = new_test([1:10],      [3, 7], 'MException');
tl{end+1} = new_test([5:10],      [3]   , {[1,1,1,1,1,1]'});
tl{end+1} = new_test([5:10],      [3, 7], {[1,1,2,2,2,2]'});

EJ_library.atest.run_tests(tl)
end
