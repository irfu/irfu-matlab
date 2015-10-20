function out=geocentric_coordinate_transformation(inp,flag)
% IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION 
%   coordinate transformation GEO/GEI/GSE/GSM/SM/MAG
%
%	[out]=IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION(inp,'coord1>coord2')
%		where coord1 and coord2 can be geo/gei/gse/gsm/sm/mag
%		inp can be matrix with 1st column is time, 2-4th columns are x,y,z.
%   inp can be TSeries of vector data.
%
% [out]=IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION(inp,'coord2')
%   if inp is TSeries object coord1 is obtained from inp.userData.COORDINATE_SYSTEM
%
%	[T]=IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION(t,'coord1>coord2')
%		if only column vector with time given, return transformation matrix
%
%	[out]=IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION(t,'dipoledirectiongse')
%		out - 1st column time, 2-4th column dipole direction in GSE 
%
%	Example:
%		xgsm=irf.geocentric_coordinate_transformation(xgse,'GSE>GSM');
%
% Ref: Hapgood 1997 (corrected version of Hapgood 1992)
% Planet. Space Sci.. Vol. 40, No. 5. pp. 71l-717, 1992

%% Alternative for T1, src: aa.usno.navy.mil/faq/docs/GAST.php (Last modified: 2011/06/14T14:04)
% Greenwich mean sidereal time (GMST), in hours (24.0h = 360.0 deg).
% GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T^2
% where;
% Let JD be the Julian date of the time of interest. Let JD0 be the Julian
% date of the previous midnight (0h) UT (the value of JD0 will end in .5
% exactly), and let H be the hours of UT elapsed since that time. Thus we
% have JD = JD0 + H/24. For both of these Julian dates, compute the number
% of days and fraction (+ or -) from 2000 January 1, 12h UT, Julian date
% 2451545.0:
% D = JD - 2451545.0;
% D0 = JD0 - 2451545.0;
% T = D/36525; % is the number of centuries since the year 2000; thus the
% last term can be omitted in most applications. It will be necessary to
% reduce GMST to the range 0h to 24h.
% ThoNi; In other words, GMST = rem(GMST,24);
% ThoNi; and theta = (360/24)*GMST; % Convert hours to degrees.
% The epoch designated "J2000.0" is specified as Julian date 2451545.0 TT,
% or 2000 January 1, 12h TT.
% This epoch can also be expressed as 2000 January 1, 11:59:27.816 TAI or
% 2000 January 1, 11:58:55.816 UTC.
% ThoNi; This is exactly what we have in EpochTT(int64(0)).
% ThoNi; therefor: D = double(EpochTT.ttns)/(10^9*86400);
% ThoNi; or simply: D = EpochTT.tts / 86400; % Julian date from J2000.
% ThoNi; and D0 = floor(EpochTT.tts/86400)-0.5;
% ThoNi; and H = 24*(JD-JD0) = 24*(D-D0);
%
%% Alterative for T2, src aa.usno.navy.mil/faq/docs/SunApprox.php (Last modified: 2012/11/06T14:12)
% compute D, the number of days and fraction (+/-) from the epoch referred
% to as "J2000.0", which is 2000 January 1.5, Julian date 2451545.0:
% D = JD - 2451545.0;
% ThoNi; Same as above: D = EpochTT.tts/86400;
% where JD is the Julian date of interest. Then compute
% Mean anomaly of the Sun:
% M = 357.529 + 0.98560028 * D;
% Mean longitude of the Sun:
% q = 280.459 + 0.98564736 * D;
% Geocentric apparent ecliptic longitude
% of the Sun (adjusted for aberration):
% L = q + 1.915*sind(M) + 0.020*sind(2*M);
% where all the constants (therefore M, q, and L) are in degrees
% The mean obliquity of the ecliptic, in degrees:
% e = 23.439 - 0.00000036 * D;
%
%% End of Alternative

persistent dipoleDirectionGSE

isInpTSeries = isa(inp,'TSeries');

if strfind(flag,'>') % if input and output reference frames are the same return input
	refSystIn  = flag(1:strfind(flag,'>')-1);
	refSystOut = flag(strfind(flag,'>')+1:end);
else
	refSystOut = lower(flag);
end
if isInpTSeries 
	refSystInVariable = lower(inp.coordinateSystem);
	refSystInVariable = refSystInVariable(1:strfind(refSystInVariable,'>')-1);
	if isempty(refSystInVariable) && isempty(refSystIn)
		errStr = 'input reference frame undefined';
		irf.log('critical',errStr);error(errStr);
	end
	if ~isempty(refSystInVariable) && ~isempty(refSystIn) && ~strcmpi(refSystInVariable,refSystIn)
		errStr = 'input reference frame as defined in variable and input flag differs';
		irf.log('critical',errStr);error(errStr);
	end
	if isempty(refSystIn)
		refSystIn = lower(refSystInVariable);
	end
	flag = [refSystIn '>' refSystOut];
end

if strcmpi(refSystIn,refSystOut),
	out=inp;
	return
end

if isInpTSeries
	t = inp.time.epochUnix;
	inpTs = inp;
	inp = double(inp.data);
else
	t=inp(:,1);
	inp(:,1) = [];
end

timeVec       = irf_time(t,'vector');
dayStartEpoch = irf_time([timeVec(:,[1 2 3]) timeVec(:,1)*[0 0 0]],'vector>epoch');
mjdRefEpoch   = irf_time([2000 1 1 12 0 0],'vector>epoch');
% Tzero is time measured in Julian centuries from 2000-01-01 12:00 UT to the previous midnight
Tzero         = (dayStartEpoch - mjdRefEpoch)/3600/24/36525.0; 
UT            = timeVec(:,4)+timeVec(:,5)/60+timeVec(:,6)/3600;

switch lower(flag)
	case 'gse>gsm', tInd = 3;
	case 'gsm>gse', tInd = -3;
	
	case 'gse>gei', tInd = -2;
	case 'gse>geo', tInd = [1 -2];
	case 'gse>sm',  tInd = [4 3];
	case 'gse>mag', tInd = [5 1 -2];
	
	case 'gsm>gei', tInd = [-2 -3];
	case 'gsm>geo', tInd = [1 -2 -3];
	case 'gsm>sm',  tInd = 4;
	case 'gsm>mag', tInd = [5 1 -2 -3];

	case 'sm>gei', tInd = [-2 -3 -4];
	case 'sm>geo', tInd = [1 -2 -3 -4];
	case 'sm>gse', tInd = [-3 -4];
	case 'sm>gsm', tInd = -4;
	case 'sm>mag', tInd = [5 1 -2 -3 -4];

	case 'mag>gei', tInd = [-1 -5];
	case 'mag>geo', tInd = -5;
	case 'mag>gse', tInd = [2 -1 -5];
	case 'mag>gsm', tInd = [3 2 -1 -5];
	case 'mag>sm',  tInd = [4 3 2 -1 -5];

	case 'geo>gei', tInd = -1;
	case 'geo>gse', tInd = [2 -1];
	case 'geo>gsm', tInd = [3 2 -1];
	case 'geo>sm',  tInd = [4 3 2 -1];
	case 'geo>mag', tInd = 5;
	
	case 'gei>geo', tInd = 1;
	case 'gei>gse', tInd = 2;
	case 'gei>gsm', tInd = [3 2];
	case 'gei>sm',  tInd = [4 3 2];
	case 'gei>mag', tInd = [5 1];
	
	case 'dipoledirectiongse'
		out = [t dipole_direction_gse];
		return
	otherwise
		irf.log('critical',['Transformation ''' flag ''' is unknown!']);
		error('Fix transformation!');
end
if size(inp,2)>=3, % input is time and 3 components
	out = mult(T(tInd),inp(:,1:3));
elseif size(inp,2)==0, % input is time, return only transformation matrix
	out=T(tInd);
end
if size(inp,2) > 3, % more columns than 3 components
	irf.log('warning','Input has more columns than 3 components! Replicating last columns in output');
	out(:,4:size(inp,2))=inp(:,4:end); % replicate last columns in output
end

	function Tout=T(tInd)		
		for j=numel(tInd):-1:1
			tNum=tInd(j);
			if     tNum == 1 || tNum == -1 % T1
				theta = 100.461 + 36000.770*Tzero + 15.04107*UT;
				T = triang(theta*sign(tNum),3); % invert if tInd=-1
			elseif tNum == 2 || tNum == -2 % T2
				eps = 23.439 - 0.013*Tzero;
				Ttemp1 = triang(eps,1);
				M = 357.528+35999.050*Tzero+0.04107*UT; % Suns mean anomaly
				L = 280.460+36000.772*Tzero+0.04107*UT; % Suns mean longitude
				lSun = L + (1.915-0.0048*Tzero).*sind(M)+0.020*sind(2*M);
				Ttemp2 = triang(lSun,3);
				T = mult(Ttemp2,Ttemp1);
				if tNum==-2, 
					T = inverse(T);
				end
			elseif tNum == 3 || tNum == -3 % T3
				if ~exist('dipoleDirectionGSE','var') || ...
						numel(t) ~= size(dipoleDirectionGSE,1) || ...
						~all(t == dipoleDirectionGSE(:,1))
					dipoleDirectionGSE=dipole_direction_gse;
				end
				ye  = dipoleDirectionGSE(:,3); % 1st col is time
				ze  = dipoleDirectionGSE(:,4);
				psi = atand(ye./ze);
%				disp(['gse2gsm angle = ' num2str(psi) ' deg']);
%				disp(['  dipole tilt = ' num2str(asind(dipoleDirectionGSE(1))) ' deg']);
				T  = triang(-psi*sign(tNum),1); % inverse if -3
			elseif tNum == 4 || tNum == -4 % T4
				if ~exist('dipoleDirectionGSE','var') || ...
						numel(t) ~= size(dipoleDirectionGSE,1) || ...
						any(t == dipoleDirectionGSE(:,1))
					dipoleDirectionGSE=dipole_direction_gse;
				end
				mu = atand(dipoleDirectionGSE(:,2)./...
					sqrt(dipoleDirectionGSE(:,3).^2 + dipoleDirectionGSE(:,4).^2));
				T = triang(-mu*sign(tNum),2);
			elseif tNum == 5 || tNum == -5 % T5
				[lambda,phi]=model.igrf(t,'dipole');
				T=mult(triang(phi-90,2),triang(lambda,3));
				if tNum == -5,
					T=inverse(T);
				end
			end			
			if j==numel(tInd),
				Tout=T;
			else
				Tout=mult(T,Tout);
			end
		end
	end
	function dipoleDirectionGSE=dipole_direction_gse
		[lambda,phi]=model.igrf(t,'dipole');
		cosPhi = cosd(phi);
		dipoleDirectionGEO = [cosPhi.*cosd(lambda) cosPhi.*sind(lambda) sind(phi)];
		dipoleDirectionGSE = ...
			irf.geocentric_coordinate_transformation([t dipoleDirectionGEO],'geo>gse'); % mult(T2,invT1);
	end

	function Tout=inverse(Tin)
		Tout=Tin;
		Tout(:,1,2)=Tin(:,2,1);
		Tout(:,1,3)=Tin(:,3,1);
		Tout(:,2,3)=Tin(:,3,2);
		Tout(:,2,1)=Tin(:,1,2);
		Tout(:,3,1)=Tin(:,1,3);
		Tout(:,3,2)=Tin(:,2,3);
	end

	function T=triang(angle,ax)
		% angle in degrees
		% ax=1 > X, ax=2 > y,ax=3 > z)
		
		cosAngle = cosd(angle);
		sinAngle = sind(angle);
		a = [1 2 3];
		axXX = a(ax~=a);
		
		T = zeros(numel(angle),3,3);
		T(:,axXX(1),axXX(1)) = cosAngle;
		T(:,axXX(2),axXX(2)) = cosAngle;
		T(:,ax,ax)=1;
		T(:,axXX(1),axXX(2)) = sinAngle;
		T(:,axXX(2),axXX(1)) = -sinAngle;
	end

	function out = mult(inp1,inp2)
		dimInp1 = numel(size(inp1));
		dimInp2 = numel(size(inp2));
		if (dimInp1==3) && (dimInp2==3)
			T = inp1;
			for ii=1:3,
				for jj=1:3,
					T(:,ii,jj)=...
						inp1(:,ii,1).*inp2(:,1,jj)+...
						inp1(:,ii,2).*inp2(:,2,jj)+...
						inp1(:,ii,3).*inp2(:,3,jj);
				end
			end
		elseif (dimInp1==3) && (dimInp2==2)
			T = inp2;
			for ii=1:3,
				T(:,ii)=inp1(:,ii,1).*inp2(:,1)+...
						inp1(:,ii,2).*inp2(:,2)+...
						inp1(:,ii,3).*inp2(:,3);
			end
		end
		out = T;
	end

if isInpTSeries,
	outData = out;
	out = inpTs;
	out.data = outData;
	if isfield(out.userData,'COORDINATE_SYSTEM')
		out.userData.COORDINATE_SYSTEM = upper(refSystOut);
	end
else
	out = [t out];
end

end
