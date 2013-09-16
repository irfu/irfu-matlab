function tests = TimeTest
tests = functiontests(localfunctions);
end

function testLeapSecond(testCase)
t1=irf.Time('2005-12-31T12:00:00Z');
t2=irf.Time('2006-01-01T12:00:00Z');
actSolution = t2-t1-int64((24*3600)*1e9);
expSolution = int64(1e9);
verifyEqual(testCase,actSolution,expSolution);
end