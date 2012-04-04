

% This is an automatic test script for checking whether the time routines
% for converting to/from ISADT epoch are working correctly. There have been
% a lot of problems over the years with these routines.


ntests=1e5;
disp('Starting tests of ISDAT epoch conversion tools.')
disp([' Using ' num2str(ntests) ' times for random tests'])
failed=0;

% Test 1: known roundoff errors
disp('Test 1: known roundoff errors...')
tol=iso2epoch('2020-01-01T00:00:00Z')*eps;
t=1272094620.1;  % Known to cause problems on earlier Matlab versions
d=fromepoch(t);
if ~isequal(d(1:5),[2010 04 24 07 37]) || abs(d(6)-0.1)>tol
	disp('FAILED test 1: known roundoff errors.')
	failed=1;
else disp('PASSED.')
end

% Test 2: roundoff errors from random times
disp('Test 2: roundoff errors from random times...')
t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
d=fromepoch(t);
if any(d(:,6)>60)
	disp('FAILED test 2: random check for roundoff errors.')
	failed=1;
else disp('PASSED.')
end

% Test 3: roundoff errors from times near second boundaries
disp('Test 3: roundoff errors from times near second boundaries...')
t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
t=fix(t)+rand(1,ntests)*0.0001;
d=fromepoch(t);
if any(d(:,6)>60)
	disp('FAILED test 3: targeted check for roundoff errors.')
	failed=1;
else disp('PASSED.')
end

% test 4: epoch2iso/iso2epoch round trip
disp('Test 4: toepoch/iso2epoch round trip...')
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
	disp('FAILED test 4: epoch2iso/iso2epoch round trip.')
	failed=1;
	cnt=0;
	for i=1:ntests
		if strcmp(ISO(i,:),ISO2(i,:))~=1
			disp(['Failed for time ' ISO(i,:) ' --> ' ISO2(i,:)])
			cnt=cnt+1;
		end
		if cnt > 5
			disp('... and possibly more times (not listing all failures).')
			break
		end
	end
else disp('PASSED.')
end


%test 5: leap seconds
disp('Test 5: leap seconds ...')
disp('Caution: this can cause Matlab to hang.')
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
		failed=1;
		leap_s=1;
		disp('FAILED test 5: ISDAT does not use leap seconds.')
		break;
	end
end
if leap_s==0,disp('PASSED: Neither Matlab nor ISDAT uses leap seconds.'),end

% End of test
disp(' ')
disp(' ')
if failed == 0, disp('Passed testing. OK.')
else disp('FAILED testing'),end
