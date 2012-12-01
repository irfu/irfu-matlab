function ok=test(varargin)
% IRF.TEST test different irfu-matlab routines
%
%   IRF.TEST show list of possible tests
%   IRF.TEST('full')  - make full test
%	IRF.TEST(testname) - make test 'testname'
%
%	Example:
%   IRF.TEST('irf_time')  - test irf_time
%

% $Id$

testList={'irf_time'};

if nargin == 0,
	disp('Possible tests of irf.test');
	disp('--------------------------');
	for j=1:numel(testList),
		disp([num2str(j) '. ' testList{j}]);
	end
	disp('--------------------------');
	disp('Execute irf.test(''testname'')');
	return
else
	testName = varargin{1};
	if ~ischar(testName),
		irf_log('fcal','Unknown test');
		return;
	end
	isTest = strcmpi(testList,testName);
	if any(isTest)
		testNumber = find(isTest);
	elseif strcmpi('full',testName)
		testNumber=1:numel(testList);
	else
		disp(['!!! Test ''' testName ''' not in list!!!']);
		irf.test;
		return;
	end
end
nTests = numel(testNumber);
okArray = false(nTests);
for j=1:nTests
	disp('*********************');
	if numel(testNumber)>1,
		disp(['TEST 1/' num2str(nTests) ': ' testList{j}]);
	else
		disp(['TEST: ' testList{j}]);
	end
	disp('*********************');
	tic;
	functionName = ['test_' testList{testNumber(j)}];
	okArray(j)=eval(functionName);
	tElapsed = toc;
	disp(' ');
	disp(['Test ended. Elapsed time: ' num2str(tElapsed) 's.']);
end

if all(okArray),
	disp('ALL TESTS PASSED SUCCESFULL!');
else
	for j=find(~okArray)
		disp(' ');
		disp(['Test FAILED! ' testList{testNumber(j)}]);
	end
end

if nargout == 1, 
	ok = okArray;
end
end

function okTest=test_irf_time
try
	okTest = 1;
	% generate vector with 10000 times during last 20 years
	a=rand(1000,1);
	tDateArray = now - 365*20*a;
	s1=irf_time(tDateArray,'date2iso');
	t=irf_time(s1,'iso2epoch');
	s2=irf_time(t,'epoch2iso');
	ok = strcmp(s1,s2);
	plusminus(ok); disp('1000 random iso>epoch>iso ');	
	okTest = okTest * ok;
	% next test
	tint=[t t+a*1000];
	s1=irf_time(tint,'tint2iso');
	tt=irf_time(s1,'iso2tint');
	s2=irf_time(tt,'tint2iso');
	ok = strcmp(s1,s2);
	plusminus(ok); disp('1000 random iso>tint>iso ');	
	okTest = okTest * ok;
	% next test. different iso formats should be recognized
	s1=irf_time(tDateArray,'date2iso');
	t=irf_time(s1,'iso2epoch');
	t1=irf_time(s1(:,1:end-1),'iso2epoch');
	t2=irf_time(reshape(strrep(s1(:)','T',' '),size(s1)),'iso2epoch');
	t3=irf_time(reshape(strrep(s1(:)','Z',' '),size(s1)),'iso2epoch');
	if all(t==t1) && all(t==t2) && all(t==t3)
		ok=1;
	else
		ok=0;
	end
	plusminus(ok); disp('1000 random iso (4 formats) > tint');	
	okTest = okTest * ok;
	% next test
catch
	okTest=false;
end
end

function plusminus(ok)
if ok,
	fprintf('+ ');
else
	fprintf('- ');
end
end

