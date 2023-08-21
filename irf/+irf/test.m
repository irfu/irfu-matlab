function ok=test(varargin)
% IRF.TEST test different irfu-matlab routines
%
%   IRF.TEST show list of possible tests
%   IRF.TEST('all')  - make all tests
%   IRF.TEST(testname) - make test 'testname'
%
%	Example:
%   IRF.TEST('irf_time')  - test irf_time
%


testList={'irf_time','epoch','c_4','coord_sys'};

if nargin == 0
  disp('Possible tests');
  disp('--------------------------');
  for j=1:numel(testList)
    disp([num2str(j) '. ' testList{j}]);
  end
  disp('--------------------------');
  disp('Execute: irf.test(''testname'') or irf.test(number) or irf.test(''all'')');
  return
else
  testName = varargin{1};
  if ischar(testName)
    isTest = strcmpi(testList,testName);
    if any(isTest)
      testNumber = find(isTest);
    elseif strcmpi('all',testName)
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
  if numel(testNumber)>1
    disp(['TEST ' num2str(j) '/' num2str(nTests) ': ' testList{iTest}]);
  else
    disp(['TEST: ' testList{iTest}]);
  end
  disp('*********************');
  tic;
  functionName = ['test_of_' testList{iTest}];
  okArray(j)=eval(functionName);
  tElapsed = toc;
  disp(['Test ended. Elapsed time: ' num2str(tElapsed) 's.']);
end

if all(okArray)
  disp('ALL TESTS PASSED SUCCESFULL!');
else
  for j=find(~okArray)
    disp(' ');
    disp(['Test FAILED! ' testList{testNumber(j)}]);
  end
end

if nargout == 1
  ok = okArray;
end
end

% Tests
function okTest = test_of_irf_time
testResults = run(test_irf_time);
okTest = check_test_results(testResults);
end
function okTest = test_of_c_4
testResults = run(testC4);
okTest = check_test_results(testResults);
end
function okTest = test_of_epoch
try
  okTest		= true; % default for all test
  ntests=1e5;
  disp('Testing ISDAT epoch conversion tools.')
  disp([' Using ' num2str(ntests) ' times for random tests'])

  %% SUBTEST: roundoff error
  okSubtest	= true; % default for subtest
  tol=iso2epoch('2020-01-01T00:00:00Z')*eps;
  t=1272094620.1;  % Known to cause problems on earlier Matlab versions
  d=fromepoch(t);
  if ~isequal(d(1:5),[2010 04 24 07 37]) || abs(d(6)-0.1)>tol
    okSubtest = false;
  end
  plusminus(okSubtest); disp('known earlier roundoff error');
  okTest = okTest * okSubtest;
  %% SUBTEST: roundoff errors from random times
  okSubtest	= true; % default for subtest
  t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
  d=fromepoch(t);
  if any(d(:,6)>60)
    okSubtest = false;
  end
  plusminus(okSubtest); disp('random roundoff errors');
  okTest = okTest * okSubtest;

  %% SUBTEST: roundoff errors from times near second boundaries
  okSubtest	= true; % default for subtest
  t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
  t=fix(t)+rand(1,ntests)*0.0001;
  d=fromepoch(t);
  if any(d(:,6)>60)
    okSubtest = false;
  end
  plusminus(okSubtest); disp('roundoff errors near second boundaries');
  okTest = okTest * okSubtest;
  %% SUBTEST: epoch2iso/iso2epoch round trip
  okSubtest	= true; % default for subtest
  year=fix(rand(ntests,1)*(2020-1970)+1970);
  month=fix(rand(ntests,1)*12)+1;
  day=fix(rand(ntests,1).*(eomday(year,month)))+1;
  hour=fix(rand(ntests,1)*24);
  minute=fix(rand(ntests,1)*60);
  second=rand(ntests,1);
  t=[year month day hour minute second];
  ISO =sprintf('%04d-%02d-%02dT%02d:%02d:%09.6fZ',t');
  ISO =reshape(ISO,27,ntests)';
  ISO2=epoch2iso(iso2epoch(ISO));
  if strcmp(ISO,ISO2)~=1
    okSubtest = false;
  end
  plusminus(okSubtest); disp('epoch2iso/iso2epoch round trip');
  okTest = okTest * okSubtest;
  if ~okSubtest
    cnt=0;
    for i=1:ntests
      if strcmp(ISO(i,:),ISO2(i,:))~=1
        disp(['  Failed for time ' ISO(i,:) ' --> ' ISO2(i,:)])
        cnt=cnt+1;
      end
      if cnt > 5
        disp('  ... and possibly more times (not listing all failures).')
        break
      end
    end
  end

  %% SUBTEST: leap seconds
  okSubtest	= true; % default for subtest
  stepdates = [...
    'Jan 6 1980'
    'Jul 1 1981'
    'Jul 1 1982'
    'Jul 1 1983'
    'Jul 1 1985'
    'Jan 1 1988'
    'Jan 1 1990'
    'Jan 1 1991'
    'Jul 1 1992'
    'Jul 1 1993'
    'Jul 1 1994'
    'Jan 1 1996'
    'Jul 1 1997'
    'Jan 1 1999'
    'Jan 1 2006'
    'Jan 1 2009'];
  stepdates = datenum(stepdates)';
  ISOtime=epoch2iso(date2epoch(stepdates)-0.5);
  for i=1:length(stepdates)
    leap_s=0;
    if strcmp(ISOtime(1,18:19),'60')
      okSubtest = false;
      leap_s=1;

      break;
    end
  end
  plusminus(okSubtest); disp('leap seconds');
  okTest = okTest * okSubtest;
  if leap_s==0,disp('  Neither Matlab nor ISDAT uses leap seconds.'),end
  if leap_s==1 && ~okSubtest, disp('ISDAT does not use leap seconds.'); end


catch
  okTest = false;
end
end
function okTest = test_of_coord_sys
try
  okTest		= true; % default for all test
  %% SUBTEST: conversion between geo/gei/gse/gsm/sm/mag
  okSubtest	= true; % default for subtest
  %
  % Here comes the subtest1
  %
  coordsysList = {'gei','geo','gse','gsm','sm','mag'};
  % generate vector with 1000 times during last 50 years
  iTimes = 100;
  a=rand(iTimes,1);a=a(:);
  tDateArray = now - 365*50*a;
  t=irf_time(tDateArray,'datenum>tt');
  for iT = 1:numel(t)
    iCoord = [randi(numel(coordsysList),1,3) 0];
    iCoord(end)=iCoord(1);
    vecStart = rand(1,3)*rand(1)*1e4;
    vec=[t(iT) vecStart];
    for jCoord = 1:numel(iCoord)-1
      transformString = [coordsysList{iCoord(jCoord)} '>' coordsysList{iCoord(jCoord+1)}];
      vec = irf.geocentric_coordinate_transformation(vec,transformString);
    end
    if abs(vec(2:4)-vecStart)>1e-10
      okSubtest = false;
      disp(['failed time: ' irf_time(t(iT),'ttns>utc') ]);
      disp(['failed conversion: ' coordsysList(iCoord) ]);
      disp(['failed start vector: ' num2str(vecStart,'%9.2e')]);
      disp(['failed end vector: ' num2str(vec(2:4),'%9.2e')]);
      disp(['eps: ' num2str(abs(vec(2:4)-vecStart))]);
      break;
    end
  end
  plusminus(okSubtest); disp([num2str(iTimes) ' random cyclic transformations gei/geo/gse/gsm/sm/mag']);
  okTest = okTest * okSubtest;
catch
  okTest = false;
end
end
function okTest = template_test
try
  okTest		= true; % default for all test
  okSubtest	= true; % default for subtest

  %% SUBTEST: description of subtest1
  okSubtest	= true; % default for subtest
  %
  % Here comes the subtest1
  %
  plusminus(okSubtest); disp('subtest1 text');
  okTest = okTest * okSubtest;
  %% SUBTEST: description of test2
  okSubtest	= true; % default for subtest
  %
  % Here comes the subtest2
  %
  plusminus(okSubtest); disp('subtest2 text');
  okTest = okTest * okSubtest;

catch
  okTest = false;
end
end

% Subfunctions
function plusminus(ok)
if ok
  fprintf('+ ');
else
  fprintf('- ');
end
end
function okTest = check_test_results(testResults)
okTest = true;
for iTest = 1:length(testResults)
  if testResults(iTest).Failed
    okTest = false;
    testResults(iTest);
  elseif testResults(iTest).Incomplete
    okTest = false;
    testResults(iTest);
  end
end
end

