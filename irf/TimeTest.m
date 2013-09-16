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

function testAddSeconds(testCase)
t = irf.Time('2006-01-01T12:00:00Z')+5*3600+21*60+17;
actSolution = toUTC(t);
expSolution = '2006-01-01T17:21:17.000000000';
verifyEqual(testCase,actSolution,expSolution);
end