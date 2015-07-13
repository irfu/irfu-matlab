function tests = test_mms_dsl2gse
tests = functiontests(localfunctions);
end


function testNoRotation(testCase)
epoch = EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:3));
data4x3 = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
TsIn = irf.ts_vec_xyz(epoch,data4x3);
sax = [0 0 1];

actSolution = mms_dsl2gse(TsIn, sax);
expSolution = [0 1 0];
verifyEqual(testCase,actSolution.data(1,:,:),expSolution)
end

function test45Rotation(testCase)
epoch = EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:3));
data4x3 = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
TsIn = irf.ts_vec_xyz(epoch,data4x3);
sax = [0 cosd(45) cosd(45)];

actSolution = mms_dsl2gse(TsIn, sax);
expSolution = [0 cosd(45) -cosd(45)];
verifyEqual(testCase,actSolution.data(1,:,:),expSolution,'AbsTol',1e-12)

actSolution = mms_dsl2gse(TsIn, sax,-1);
expSolution = [0 cosd(45) cosd(45)];
verifyEqual(testCase,actSolution.data(1,:,:),expSolution,'AbsTol',1e-12)
end