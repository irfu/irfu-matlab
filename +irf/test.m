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

testList={'irf_time','c_4'};

if nargin == 0,
	disp('Possible tests');
	disp('--------------------------');
	for j=1:numel(testList),
		disp([num2str(j) '. ' testList{j}]);
	end
	disp('--------------------------');
	disp('Execute: irf.test(''testname'') or irf.test(number)');
	return
else
	testName = varargin{1};
	if ischar(testName),
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
	elseif isnumeric(testName) %assume given test number
		testNumber=testName;
	else
		irf_log('fcal','Unknown test');
		return;		
	end
end
nTests = numel(testNumber);
okArray = false(nTests,1);
for j=1:nTests
	iTest=testNumber(j);
	disp(' ');
	disp('*********************');
	if numel(testNumber)>1,
		disp(['TEST ' num2str(j) '/' num2str(nTests) ': ' testList{iTest}]);
	else
		disp(['TEST: ' testList{iTest}]);
	end
	disp('*********************');
	tic;
	functionName = ['test_' testList{iTest}];
	okArray(j)=eval(functionName);
	tElapsed = toc;
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
function okTest=test_c_4
try
	okTest = true; % default for all test
	ok     = true; % default for subtest
	% position of tetrahedron, base plane in xy
	R1n=[1 0 0];
	R2n=[cos(1/3*2*pi) sin(1/3*2*pi) 0];
	R3n=[cos(2/3*2*pi) sin(2/3*2*pi) 0];
	R4n=[0 0 sqrt(norm(R2n-R1n)^2-1)];
	zMassCenter=R4n(3)/4;
	c_eval('R?n(3)=R?n(3)-zMassCenter;');
	c_eval('R?=R?n;');
	
	% construction of field
	Bconst	=@(x) [0 0 1]; %
	
	% Cylindrical field with symmetry axis along z and center at rc
	% such field gives curvature equal to 1/R where R is distance to
	% the cylinder symmetry axis
	Bcircle	=@(x,rc) [(x(:,2)-rc(2))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) ...
		-(x(:,1)-rc(1))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) 0];
	rc=[10 10 0]; % cylinder symmetry axis location, B field is clockwise
	%disp('Test curvature')
	%disp([' should be: [' num2str(1/norm(rc)^2*rc,'%5.2f') ']']);
	c_eval('B?=Bcircle(R?,rc);');
	curv=c_4_grad('R?','B?','curvature');
	%disp(['        is: [' num2str(curv,'%5.2f') ']'])
	if norm(1/norm(rc)^2*rc - curv) > 1e-1
		disp('!!!! Curvature failed!!!!')
		disp([' should be: [' num2str(1/norm(rc)^2*rc,'%5.2f') ']']);
		disp(['        is: [' num2str(curv,'%5.2f') ']'])
		disp(['     error: ' num2str(norm(1/norm(rc)^2*rc - curv))])
		ok = false;
	end
	plusminus(ok); disp('curvature');	
	okTest = okTest * ok;
	
	% Constant current in Z direction
	jz=10;
	Bjz = @(x) [0 jz*x(1) 0];
	c_eval('B?=Bjz(R?);');
	testCurl=c_4_grad('R?','B?','curl');
	if norm([0 0 jz] - testCurl) > 1e-10
		disp('!!!! Curl failed!!!!')
		disp([' should be: [0 0 ' num2str(jz) ']']);
		disp(['        is: [' num2str(testCurl,'%6.2f') ']'])
		disp(['     error: ' num2str(norm([0 0 jz] - testCurl))])
		ok = false;
	end
	plusminus(ok); disp('curl');	
	okTest = okTest * ok;
	
	% Test drift gradient
	% Fitzpatrick book page 30
	% V=mv^2/2/e B x gradB / B^3

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

