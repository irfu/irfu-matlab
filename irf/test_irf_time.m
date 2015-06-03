classdef test_irf_time < matlab.unittest.TestCase
	%TEST_IRF_TIME
	
	properties
	end
	
	methods (Test)
		function test_UTC_to_TTns_to_UTC(testCase)
			% 1000 random utc>ttns>utc
			
			%generate vector with 10000 times during last 500 years
			tDateArray = now - 365*500*rand(1000,1);
			s1=irf_time(tDateArray,'date>utc');
			t=irf_time(s1,'utc>ttns');
			s2=irf_time(t,'ttns>utc');
			testCase.verifyEqual(s1,s2);
			ok = strcmp(s1,s2);
		end
		function test_tint_to_UTC_to_tint(testCase)
			% 1000 random tint>utc>tint
			tDateArray = now - 365*500*rand(1000,1);
			tEpoch=irf_time(tDateArray,'date>epoch');
			tint=[tEpoch tEpoch+rand(1000,1)*2*365*24*3600]; % random length intervals up to 2 year
			s1=irf_time(tint,'tint>utc');
			tt=irf_time(s1,'utc>tint');
			testCase.verifyEqual(tint,tt,'AbsTol',2e-9); % 2ns precision
		end
		function test_different_UTC_formats(testCase)
			tDateArray = now - 365*500*rand(1000,1);
			s1=irf_time(tDateArray,'date>utc');
			t=irf_time(s1,'utc>ttns');
			t1=irf_time(s1(:,1:end-1),'utc>ttns');
			t2=irf_time(reshape(strrep(s1(:)','Z',''),size(s1)-[0 1]),'utc>ttns');
			testCase.verifyEqual(t,t1,'AbsTol',2e-9); % 2ns precision
			testCase.verifyEqual(t,t2,'AbsTol',2e-9); % 2ns precision
		end
	end
end

