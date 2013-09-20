function tests = TimeTest
tests = functiontests(localfunctions);
end

function test_leap_second(testCase)
t1=irf.Time('2005-12-31T12:00:00Z');
t2=irf.Time('2006-01-01T12:00:00Z');
actSolution = t2-t1-(24*3600);
expSolution = 1;
verifyEqual(testCase,actSolution,expSolution);
end

function test_add_seconds(testCase)
t = irf.Time('2006-01-01T12:00:00Z')+5*3600+21*60+17;
actSolution = toUTC(t);
expSolution = '2006-01-01T17:21:17.000000000';
verifyEqual(testCase,actSolution,expSolution);
end

function test_time_array(testCase)
t = irf.Time('2005-12-31T12:00:00Z');
a = t + [1 2 3];
actSolution = length(a);
expSolution = 3;
verifyEqual(testCase,actSolution,expSolution);
end