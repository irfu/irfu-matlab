function f = irf_get_data_omni( tint, parameter , database)
%IRF_GET_DATA_OMNI download omni data
%
% f=IRF_GET_DATA_OMNI(tint,parameter) download parameters for specified time interval
%       f is column vector, first column time and parameters in next columns
%
% parameters - string, paramters separated by comma (case does not matter).
%               'B'     - <|B|>, magnetic field magnitude [nT]
%               'avgB'  - |<B>|, magnitude of average magnetic field [nT]
%               'Blat'  - latitude angle of average B field (GSE)
%               'Blong' - longitude angle of average B field (GSE)
%               'Bx'    - Bx GSE (the same as 'BxGSE')
%               'By'
%               'Bz'
%               'ByGSM' - By GSM
%               'BzGSM' - Bz GSM
%               'T'     - proton temperature (K)
%               'n'     - proton density (cc)
%               'NaNp'  - alpha/proton ratio
%               'v'     - bulk speed (km/s)
%               'E'     - electric field (mV/m)
%               'P'     - flow pressure (nPa)
%               'beta'  - plasma beta
%               'Ma'    - Alfven Mach number
%               'Ms'    - 1 AU IP Magnetosonic Mach number
%               'ssn'   - daily sunspot number
%               'dst'   - DST index
%               'f10.7' - F10.7 flux
%               'pc'    - PC index
%               'ae'    - AE index
%               'al'    - AL index
%               'au'    - AL index
%
% f=IRF_GET_DATA_OMNI(tint,parameter,database) download from specified database
%
% database:  'omni2'    - 1h resolution OMNI2 data (default)
%            'omni_min' - 1min resolution OMNI data
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
%   ff= irf_get_data_omni(tint,'b,bx,bygsm');
%   ff= irf_get_data_omni(tint,'f10.7');
%   ff= irf_get_data_omni(tint,'b','omni_min');

% http://omniweb.gsfc.nasa.gov/html/ow_data.html

% data description of high res 1min,5min omni data set
% http://omniweb.gsfc.nasa.gov/html/omni_min_data.html#4b
% %
% %
% Year			        I4	      1995 ... 2006
% Day			        I4	1 ... 365 or 366
% Hour			        I3	0 ... 23
% Minute			        I3	0 ... 59 at start of average
% ID for IMF spacecraft	        I3	See  footnote D below
% ID for SW Plasma spacecraft	I3	See  footnote D below
% # of points in IMF averages	I4
% # of points in Plasma averages	I4
% Percent interp		        I4	See  footnote A above
% Timeshift, sec		        I7
% RMS, Timeshift		        I7
% RMS, Phase front normal	        F6.2	See Footnotes E, F below
% Time btwn observations, sec	I7	DBOT1, See  footnote C above
% Field magnitude average, nT	F8.2
% Bx, nT (GSE, GSM)		F8.2
% By, nT (GSE)		        F8.2
% Bz, nT (GSE)		        F8.2
% By, nT (GSM)	                F8.2	Determined from post-shift GSE components
% Bz, nT (GSM)	                F8.2	Determined from post-shift GSE components
% RMS SD B scalar, nT	        F8.2
% RMS SD field vector, nT	        F8.2	See  footnote E below
% Flow speed, km/s		F8.1
% Vx Velocity, km/s, GSE	        F8.1
% Vy Velocity, km/s, GSE	        F8.1
% Vz Velocity, km/s, GSE	        F8.1
% Proton Density, n/cc		F7.2
% Temperature, K		        F9.0
% Flow pressure, nPa		F6.2	See  footnote G below
% Electric field, mV/m		F7.2	See  footnote G below
% Plasma beta		        F7.2	See  footnote G below
% Alfven mach number		F6.1	See  footnote G below
% X(s/c), GSE, Re		        F8.2
% Y(s/c), GSE, Re		        F8.2
% Z(s/c), GSE, Re		        F8.2
% BSN location, Xgse, Re	        F8.2	BSN = bow shock nose
% BSN location, Ygse, Re	        F8.2
% BSN location, Zgse, Re 	        F8.2
%
% AE-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% AL-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% AU-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% SYM/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% SYM/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% ASY/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% ASY/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% PC(N) index,                    F7.2    See World Data Center for Geomagnetism, Copenhagen
% Magnetosonic mach number        F5.1    See  footnote G below

if nargin < 3, % database not specified defaulting to omni2
	dataSource='omni2';
	dateFormat='utc_yyyymmdd';
elseif nargin == 3, % database specified
	if strcmpi(database,'omni2')
		dataSource='omni2';
		dateFormat='utc_yyyymmdd';
	elseif strcmpi(database,'omni_min') || strcmpi(database,'min')
		dataSource='omni_min';
		dateFormat='utc_yyyymmddHH';
	else
		irf_log('fcal','Unknown database, using omni2.');
		dataSource='omni2';
		dateFormat='utc_yyyymmdd';
	end
end
httpRequest=['http://omniweb.gsfc.nasa.gov/cgi/nx1.cgi?activity=retrieve&spacecraft=' dataSource '&'];
startDate=irf_time(tint(1),dateFormat);
endDate=irf_time(tint(2),dateFormat);

i=strfind(parameter,',');
iEnd=[i-1 length(parameter)];
iStart=[1 i+1];
vars='';nVar=0;
for jj=1:length(iStart)
	variable=parameter(iStart(jj):iEnd(jj));
	switch lower(variable)
		case 'b', varNumberOmni2=8;varNumberOmni1min=13;
		case 'avgb', varNumberOmni2=9;varNumberOmni1min=-1;
		case 'blat', varNumberOmni2=10;varNumberOmni1min=-1;
		case 'blong', varNumberOmni2=11;varNumberOmni1min=-1;
		case {'bx','bxgse','bxgsm'}, varNumberOmni2=12;varNumberOmni1min=14;
		case {'by','bygse'}, varNumberOmni2=13;varNumberOmni1min=15;
		case {'bz','bzgse'}, varNumberOmni2=14;varNumberOmni1min=16;
		case 'bygsm', varNumberOmni2=14;varNumberOmni1min=17;
		case 'bzgsm', varNumberOmni2=15;varNumberOmni1min=18;
		case 't', varNumberOmni2=22;varNumberOmni1min=26;
		case 'n', varNumberOmni2=23;varNumberOmni1min=25;
		case 'nanp', varNumberOmni2=27;varNumberOmni1min=-1;
		case 'v', varNumberOmni2=24;varNumberOmni1min=21;
		case 'p', varNumberOmni2=28;varNumberOmni1min=27;
		case 'e', varNumberOmni2=35;varNumberOmni1min=28;
		case 'beta', varNumberOmni2=36;varNumberOmni1min=29;
		case 'ma', varNumberOmni2=37;varNumberOmni1min=30;
        case 'ms', varNumberOmni2=56;varNumberOmni1min=45;
		case 'ssn', varNumberOmni2=39;varNumberOmni1min=-1;
		case 'dst', varNumberOmni2=40;varNumberOmni1min=-1;
		case 'ae', varNumberOmni2=41;varNumberOmni1min=37;
		case 'al', varNumberOmni2=52;varNumberOmni1min=38;
		case 'au', varNumberOmni2=53;varNumberOmni1min=39;
		case 'kp', varNumberOmni2=38;varNumberOmni1min=-1;
		case 'pc', varNumberOmni2=51;varNumberOmni1min=44;
		case 'f10.7', varNumberOmni2=50;varNumberOmni1min=-1;
		otherwise, varNumberOmni2=0;varNumberOmni1min=-1;
	end
	if strcmp(dataSource,'omni2'),
		if varNumberOmni2>0,
			vars=[vars '&vars=' num2str(varNumberOmni2)];
			nVar=nVar+1;
		end
	else % datasource omni_min
		if varNumberOmni1min>0,
			vars=[vars '&vars=' num2str(varNumberOmni1min)];
			nVar=nVar+1;
		end
	end
end

url=[httpRequest 'start_date=' startDate '&end_date=' endDate vars];
disp(['url:' url]);
[c,status]=urlread(url);

if status==1, % success in downloading from internet
	cstart=strfind(c,'YEAR'); % returned by omni2 databse
	if isempty(cstart),
		cstart=strfind(c,'YYYY'); % returned by omni_min database
	end
	if isempty(cstart), % no data returned
		irf.log('warning','Can not get OMNI data from internet!');
		f=[];
		return
	end
	cend=strfind(c,'</pre>')-1;
	if strcmp(dataSource,'omni2')
		fmt=['%f %f %f' repmat(' %f',1,nVar)];
		cc=textscan(c(cstart:cend),fmt,'headerlines',1);
		xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
		f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600;
		for jj=1:nVar,
			f(:,jj+1)=cc{jj+3};
		end
	else
		fmt=['%f %f %f %f' repmat(' %f',1,nVar)];
		cc=textscan(c(cstart:cend),fmt,'headerlines',1);
		xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
		f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600+cc{4}*60;
		for jj=1:nVar,
			f(:,jj+1)=cc{jj+4};
		end
	end
	% Remove FILLVAL, however not good solution. TODO: improve, so that
	% only FILLVAL for the corresponding variable is used
	f(f==9999999.)=NaN;
	f(f==999999.9)=NaN;
	f(f==99999.9)=NaN;
	f(f==9999.99)=NaN;
	f(f==999.99)=NaN;
	f(f==999.9)=NaN;
	f(f==999)=NaN;
	f(f==99.99)=NaN;
else % no success in getting data from internet
	irf.log('warning','Can not get OMNI data form internet!');
	f=[];
end

