function [out1,out2]=igrf(t,flag)
% MODEL.IGRF return IGRF model
%	[lambda,phi]=MODEL.IGRF(t,'dipole') return magnetic dipole latitude and longitude 

persistent hIGRF gIGRF yearsIGRF

if isempty(hIGRF)
	fileIGRF = [fileparts(which('irf.m')) '/+model/igrf11coeffs.txt'];
	irf.log('warning',['Reading IGRF coeficients from file:' fileIGRF]);
	%file reading
	fid=fopen(fileIGRF);
	out = textscan(fid, '%s', 'delimiter',sprintf('\n')); % cell array with lines
	out=out{1};
	fclose(fid);
	% construct igrf_coef
	% g(i,[n m years]) 0
	% h(i,[n,m,years]) 1
	numberCoefIGRF = numel(out) - 4;
	a=textscan(out{4},'%s');a=a{1};
	numYears = numel(a)-3;
	formatYears=repmat('%f',1,numYears);
	yearsIGRF = textscan(out{4},['%*s%*s%*s' formatYears]);
	yearsIGRF = horzcat(yearsIGRF{:});
	yearsIGRF(end)=yearsIGRF(end-1)+5;
	% read in all IGRF coefficients from file
	iIGRF=zeros(numberCoefIGRF,numYears+3);
	ii=1;
	for iLine=5:numel(out)
		d=textscan(out{iLine},['%s%f%f' formatYears]);
		iIGRF(ii,2:end)=horzcat(d{2:end});
		if d{1}{1}-'h'==0, iIGRF(ii,1)=1;end 
		ii=ii+1;
	end
	% the last column is derivative, make to value for +5 year
	iIGRF(:,end) = iIGRF(:,end-1)+5*iIGRF(:,end);
	hIGRF=iIGRF(iIGRF(:,1)==1,2:end);
	gIGRF=iIGRF(iIGRF(:,1)==0,2:end);
end

timeVec       = irf_time(t,'vector');
yearRef       = timeVec(:,1);
if min(yearRef) < min(yearsIGRF),
	irf_log('fcal','requested time before available IGRF');
	return;
end

year        = yearRef + ...
	(t - irf_time([yearRef (yearRef*0+1)*[1 1] yearRef*[0 0 0]],'vector2epoch'))/(365.25*86400);

switch flag
	case 'dipole'
		g01 = interp1(yearsIGRF,gIGRF(1,3:end),year,'linear','extrap');
		g11 = interp1(yearsIGRF,gIGRF(2,3:end),year,'linear','extrap');
		h11 = interp1(yearsIGRF,hIGRF(1,3:end),year,'linear','extrap');
		lambda = atand(h11./g11);
		phi = 90-asind((g11.*cosd(lambda)+h11.*sind(lambda))./g01);
		if nargout == 2, 
			out1 = lambda;
			out2 = phi;
		else
			disp(['lamda = ' num2str(lambda) ' deg']);
			disp(['  phi = ' num2str(phi) ' deg']);
		end
end

