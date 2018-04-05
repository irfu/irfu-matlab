function [Lstar,traces] = onera_desp_lib_grid_lstar(kext,options,sysaxes,matlabd,x1,x2,x3,maginput,varargin)
% [Lstar,traces] = onera_desp_lib_grid_lstar(kext,options,sysaxes,matlabd,x1,x2,x3,maginput,...)
% compute Lstar on a grid
% kext, options, sysaxes, matlabd, and maginput are as required by onera_desp_lib_trace_field_line
% x1, x2, x3 are [N1,N2] matrices of launch points
% ***Option Bm, a0, or K is required
% keyword options:
% 'G' - set Bunit to G (default Bunit is nT)
% 'Phi' - return Phi in nT-RE^2 rather than Lstar
% 'periodic1' - space is period in i1 (e.g., longitude, local time)
% 'periodic2' - space is period in i2 (e.g., longitude, local time)
% param/value options:
% 'Bm' - grid of Bmirror values in Bunit [N3x1]
% 'a0' - grid of equatorial pitch angles in degrees [N3x1] (slower)
% 'K' - grid of K values in RE-sqrt(Bunit) [N3x1]
% 'k0' - magnetic moment for Lstar = 2*pi*k0/Phi, Bunit-RE^3
%   default is 0.301153e5 nT-RE^3 for epoch 2000.0
% 'R0' - minimum tracing radius, RE, passed on to
%   onera_desp_lib_trace_field_line (default is 0.9)
% 'NptsI' - number of points for I integral (default 1000)
% 'dlat' - latitude spacing for Phi integral (degrees, default 0.1)
% 'verbose' - verbose setting for diagnostic/status output
%   (default is 0, set to 1 for about 1 report per second.
%    higher integers have more reports)
% 'maxLm' - stop tracing in dimension 1 if |Lm| exceeds maxLm
%   set maxLm negative to trace from N1 down to 1
% 'digits' - number of digits to round off to (a kludge for edge effects, default is inf)
%
% output:
% Lstar - L* [N1 x N2 x N3], in RE. NaN on open/shabansky orbits
% traces - cell array [N1 x N2] of structures containing trace information
%             Nt: is number of points in field line trace
%             Lm: McIlwain's L, RE
%         Blocal: [Ntx1] local field strength, Bunit
%           Bmin: minimum field strength on field line, Bunit, i.e., min(Blocal)
%           iBmin: index of Bmin: Blocal(iBmin) == Bmin
%              J: I for particle mirroring at launch point x1(i1,i2), x2(i1,i2), x3(i1,i2), RE
%            GEO: [Nx3] field points in GEO
%            RLL: [Ntx3] geographic radius, lat, lon (RE, deg, deg)
%     NorthPoint: [lat, lon] geographic location of R=1 in northern foot point, (deg)
%          nmins: 1
%          imins: [nmins x 1] index of minima
%          Bmins: [nmins x 1] B at minima, Bunit
%          imaxs: [nmaxs x 1] index of maxima
%          Bmaxs: [nmaxs x 1] B at maxima, Bunit
%             Bm: [N3x1] Bmirror, Bunit
%      alpha0deg: [N3x1] equatorial pitch angle, deg
%     alpha0rads: [N3x1] equatorial pitch angle, radians
%              I: [N3x1] I, RE
%              K: [N3x1] K, RE*sqrt(Bunit)
%              Phi: [N3x1] Phi, RE^2*Bunit
%              Lstar: [N3x1] Lstar, RE
%              foot_points: [N3x1] cell array of [Nx2] lat,lon (deg) of northern foot points of drift shell
%              NorthMirrorPoint: [N3x3] GEO coordinates of northern mirror point
%              SouthMirrorPoint: [N3x3] GEO coordinates of southern mirror point
%
% call with no inputs to run test/demo
%
% Initial implementation by Paul O'Brien
% Algorithm adapted from
% Min, K., J. Bortnik, and J. Lee (2013), A novel technique for rapid 
% L* calculation using UBK coordinates, J. Geophys. Res. Space Physics, 
% 118, doi:10.1029/2012JA018177.

if nargin == 0 % test
    [Lstar,traces] = test;
    return;
end

nT = 1; % assume nT
k0 = nan;
Bm_option = '';
output_option = 'Lstar';
R0 = 0.9;
NptsI = 1000; % number of points for I integral
dlat = 0.1; % latitude spacing for Phi integrals
verbose = 0;
maxLm = inf;
periodic = false(2,1);
digits = inf;

i = 1;
while i <= length(varargin)
    switch(lower(varargin{i}))
        case 'maxlm'
            i = i+1;
            maxLm = varargin{i};
        case 'r0'
            i = i+1;
            R0 = varargin{i};
        case 'k0'
            i = i+1;
            k0 = varargin{i};
        case 'bm'
            Bm_option = 'Bm';
            i = i+1;
            z = varargin{i};
        case 'k'
            Bm_option = 'K';
            i = i+1;
            z = varargin{i};
        case {'a0','alpha0','alpha0deg','a0deg'}
            Bm_option = 'a0';
            i = i+1;
            z = varargin{i};
        case 'periodic1'
            periodic(1) = true;
        case 'periodic2'
            periodic(2) = true;
        case 'phi'
            output_option = 'Phi';
        case 'g'
            nT = 1e-5; % multiply nT by this to get Gauss
        case 'nptsi'
            i = i+1;
            NptsI = varargin{i};
        case 'digits'
            i = i+1;
            digits = varargin{i};
        case 'verbose'
            i = i+1;
            verbose = varargin{i};
        otherwise
            error('Unknown option %s',varargin{i});
    end
    i = i+1;
end
if isempty(Bm_option)
    error('Option a0 or Bm is required');
end

z = z(:); % column vector

if ~isfinite(k0)
    k0 = 0.301153e5*nT; % k0 in nT-RE^3 for 2000.0
end

maginput = onera_desp_lib_maginputs(maginput);
[N1,N2] = size(x1);
N3 = length(z);
if verbose
    fprintf('grid size is [%d,%d,%d]\n',N1,N2,N3);
end

if maxLm<0
    I1 = N1:-1:1; % trace from high to low in dim 1
else
    I1 = 1:N1;
end

traces = cell([N1,N2]);
Phi = nan([N1,N2,N3]);

last_t = now;
for i2 = 1:N2
    for i1 = I1
        clear trace
        if verbose >= 2
            fprintf('Tracing i=(%d,%d), x=(%g,%g,%g)\n',i1,i2,x1(i1,i2),x2(i1,i2),x3(i1,i2));
        end
        [trace.Lm,trace.Blocal,trace.Bmin,trace.J,trace.GEO] = onera_desp_lib_trace_field_line(kext,options,sysaxes,matlabd,x1(i1,i2),x2(i1,i2),x3(i1,i2),maginput,R0);
        if isempty(trace.Blocal)
            break;
        end
        if (abs(trace.Lm)>abs(maxLm))
            break;
        end
        trace.Nt = length(trace.Blocal);
        % convert nT to Bunit
        trace.Blocal = trace.Blocal*nT;
        [Bmin,iBmin] = min(trace.Blocal);
        Bmin = Bmin(1); % just in case minimum occurs twice
        iBmin = iBmin(1);
        trace.Bmin = Bmin; % replace library value which can sometimes be smaller than min Blocal
        trace.iBmin = iBmin;
        % check for Shabansky
        dir = -1; % B decreasing
        imins = [];
        imaxs = [];
        for i = 2:length(trace.Blocal)
            newdir = sign(trace.Blocal(i)-trace.Blocal(i-1));
            if newdir == 0
                newdir = dir;
            end
            
            if newdir ~= dir
                if (dir == +1) && (newdir == -1)
                    % change from increasing to decreasing
                    imaxs(end+1) = i-1;
                elseif (dir == -1) && (newdir == +1)
                    % change from decreasing to increasing
                    imins(end+1) = i-1;
                end
                dir = newdir;
            end
        end
        trace.imins = imins;
        trace.imaxs = imaxs;
        trace.Bmins = trace.Blocal(imins);
        trace.Bmaxs = trace.Blocal(imaxs);
        trace.nmins = length(imins);
        Bm = unique(trace.Blocal); % try all Blocal as Bm
        Bm = Bm(Bm <= min(trace.Blocal([1,end]))); % remove Bm in bounce loss cone
        
        if trace.nmins > 1 % at least two local minima
            Bm = Bm(Bm>max(trace.Bmaxs)); % omit points in local magnetic bottle
        end
        
        I = nan(size(Bm)); % allocate space for I
        ds = sqrt(sum(diff(trace.GEO,1,1).^2,2)); % step
        s = [0;cumsum(ds)]; % distance along field line
        NorthMirrorPoint = nan(length(I),3);
        SouthMirrorPoint = nan(length(I),3);
        for i = 1:length(Bm)
            if Bm(i) <= Bmin % equatorial case
                I(i) = 0;
                NorthMirrorPoint(i,:) = trace.GEO(iBmin,:);
                SouthMirrorPoint(i,:) = trace.GEO(iBmin,:);
            else
                bm = Bm(i); % Bmirror
                s1 = interp1(trace.Blocal(1:imins(1)),s(1:imins(1)),bm,'linear'); % one mirror point
                s2 = interp1(trace.Blocal(imins(1):end),s(imins(1):end),bm,'linear'); % other mirror point
                if s1==s2
                    I(i) = 0;
                    NorthMirrorPoint(i,:) = trace.GEO(imins(1),:);
                    SouthMirrorPoint(i,:) = trace.GEO(imins(1),:);
                else
                    si = linspace(s1,s2,NptsI); % grid in s betwen mirror points
                    bi = interp1(s,trace.Blocal,si,'linear');
                    assert(all(bi./bm<=1.001)); % allow only a little fudge for numerics
                    I(i) = trapz(si,sqrt(1-min(bi./bm,1))); % integrate
                    NorthMirrorPoint(i,:) = interp1(trace.Blocal(1:imins(1)),trace.GEO(1:imins(1),:),bm,'linear'); % one mirror point
                    SouthMirrorPoint(i,:) = interp1(trace.Blocal(imins(1):end),trace.GEO(imins(1):end,:),bm,'linear'); % other mirror point
                end
            end
        end
        K = I.*sqrt(Bm); % RE*sqrt(Bunit)
        switch(Bm_option)
            case 'Bm'
                trace.Bm = z;
                trace.I = interp1(Bm,I,z,'linear');
                trace.K = trace.I.*sqrt(trace.Bm);
                trace.alpha0deg = asind(sqrt(trace.Bmin./trace.Bm));
                % cover two ways we could get not-quite-90 degrees due to
                % round off
                trace.alpha0deg(trace.Bm <= trace.Bmin) = 90;
                trace.alpha0deg(trace.I == 0) = 90;
            case 'a0'
                trace.alpha0deg = z;
                trace.Bm = trace.Bmin./sind(trace.alpha0deg).^2;
                trace.I = interp1(Bm,I,trace.Bm,'linear');
                trace.I(trace.alpha0deg==90) = 0; % force equatorial case
                trace.K = trace.I.*sqrt(trace.Bm);
            case 'K'
                trace.K = z;
                trace.Bm = interp1(K,Bm,trace.K,'linear');
                trace.I = trace.K./sqrt(trace.Bm);
                trace.alpha0deg = asind(sqrt(trace.Bmin./trace.Bm));
                trace.alpha0deg(trace.K==0) = 90; % force equatorial case
            otherwise
                error('Unknown Bm_option %s',Bm_option);
        end
        trace.alpha0rads = trace.alpha0deg*pi/180;
        trace.RLL = onera_desp_lib_coord_trans(trace.GEO,'GEO2RLL',matlabd);
        trace.NorthMirrorPoint = interp1(Bm,NorthMirrorPoint,trace.Bm,'linear');
        trace.SouthMirrorPoint = interp1(Bm,SouthMirrorPoint,trace.Bm,'linear');
        ifix = trace.Bm < Bmin; % should never happen, now that trace.Bmin = min(trace.Blocal)
        if any(ifix) % use the point with smallest B, which should be near trace.Bmin
            Nfix = sum(ifix);
            trace.NorthMirrorPoint(ifix,:) = repmat(trace.GEO(iBmin,:),Nfix,1);
            trace.SouthMirrorPoint(ifix,:) = repmat(trace.GEO(iBmin,:),Nfix,1);
        end
        % find NorthPoint geographic LAT,LON in deg of northern crossing of
        % R=1
        fNorth = find(trace.RLL(:,2)>=median(trace.RLL(:,2)));
        [Rmin,imin] = min(trace.RLL(fNorth,1));
        [Rmax,imax] = max(trace.RLL(fNorth,1));
        if Rmin>1
            trace.NorthPoint = trace.RLL(fNorth(imin),2:3); % project along radius to R=1
        elseif Rmax<1
            trace.NorthPoint = trace.RLL(fNorth(imax),2:3); % project along radius to R=1
        else
            trace.NorthPoint = nan(1,2);
            trace.NorthPoint(1) = interp1(trace.RLL(fNorth,1),trace.RLL(fNorth,2),1.0,'linear'); % interpolate to R=1
            clon = interp1(trace.RLL(fNorth,1),cosd(trace.RLL(fNorth,3)),1.0,'linear'); % interpolate to R=1
            slon = interp1(trace.RLL(fNorth,1),sind(trace.RLL(fNorth,3)),1.0,'linear'); % interpolate to R=1
            trace.NorthPoint(2) = atan2(slon,clon)*180/pi; % interpolate in periodic coordinate
        end
        if (verbose>=1) && (now-last_t>1/24/60/60)
            fprintf('Traced i=(%d,%d), x=(%g,%g,%g): Lm=%g (dipoleL=%g)\n',i1,i2,x1(i1,i2),x2(i1,i2),x3(i1,i2),trace.Lm,cosd(trace.NorthPoint(1))^-2);
            last_t = now;
        end
        traces{i1,i2} = trace;
    end
end

switch(Bm_option)
    case 'Bm' % contours of constant I on grid in Bm with I(i1,i2,Bm)
        traces = Phi_at_fixed_i3(traces,'I',periodic,digits,verbose);
    case 'a0' % contours of constant I at fixed Bm but with I(i1,i2,a0)
        traces = Phi_at_fixed_a0(traces,periodic,digits,verbose);
    case 'K' % contours of constant Bm at fixed K with Bm(i1,i2,K)
        traces = Phi_at_fixed_i3(traces,'Bm',periodic,digits,verbose);
    otherwise
        error('Unknown Bm_option %s',Bm_option);
end

% now, with foot_point populated, compute Phi and Lstar
last_t = now;
for i1 = 1:N1
    for i2 = 1:N2
        for i3 = 1:N3
            if ~isempty(traces{i1,i2}) && ~isempty(traces{i1,i2}.foot_points{i3})
                lat = traces{i1,i2}.foot_points{i3}(:,1);
                lon = traces{i1,i2}.foot_points{i3}(:,2);
                traces{i1,i2}.Phi(i3) = computePhi(kext,options,matlabd,lat,lon,maginput,dlat)*nT;
                Phi(i1,i2,i3) = traces{i1,i2}.Phi(i3);
                traces{i1,i2}.Lstar(i3) = 2*pi*k0/traces{i1,i2}.Phi(i3);
                if (verbose>=1) && (now-last_t>1/24/60/60)
                    fprintf('at i=(%d,%d,%d), x=(%g,%g,%g),%s=%g: Phi = %g (Lstar = %g,Ldip=%g)\n',i1,i2,i3,...
                        x1(i1,i2),x2(i1,i2),x3(i1,i2),Bm_option,z(i3),Phi(i1,i2,i3),traces{i1,i2}.Lstar(i3),...
                        cosd(traces{i1,i2}.NorthPoint(1))^-2);
                    last_t = now;
                end
            end % if ~isempty...
        end % for i3
    end % for i2
end % for i1


switch(output_option)
    case 'Phi'
        Lstar = Phi;
    case 'Lstar'
        Lstar = 2*pi*k0./Phi;
    otherwise
        error('Unknown output option %s',output_option);
end

function Phi = computePhi(kext,options,matlabd,lat,lon,maginput,dlat)
% compute northern poalr cap Phi integral inside trajectory lat,lon with
% latitude spacing dlat

% we'll precompute the latitude line integral on a grid
persistent settings grid
newsettings = {kext,options,matlabd,maginput,dlat}; % can we re-use the grid?
if ~isequal(settings,newsettings)
    grid.lat = linspace(-20,90,ceil((90--20)/dlat));
    grid.NLAT = length(grid.lat);
    grid.lon = linspace(0,360,2*length(lon)); % overkill longitude spacing, just in case
    grid.NLON = length(grid.lon);
    [grid.LAT,grid.LON] = meshgrid(grid.lat,grid.lon);
    grid.N = numel(grid.LAT);
    [Bgeo,B] = onera_desp_lib_get_field(kext,options,'RLL',repmat(matlabd,grid.N,1),ones(grid.N,1),grid.LAT(:),grid.LON(:),repmat(maginput,grid.N,1));
    rhat = [cosd(grid.LAT(:)).*cosd(grid.LON(:)), cosd(grid.LAT(:)).*sind(grid.LON(:)), sind(grid.LAT(:))];
    Bdotr = reshape(sum(rhat.*Bgeo,2),size(grid.LAT));
    grid.partialPhi = cumtrapz(sind(grid.lat),Bdotr,2);
    grid.partialPhi = grid.partialPhi-repmat(grid.partialPhi(:,end),1,grid.NLAT); % includes minus sign for northern hemi
    % grid.partialPhi(ilon,ilat) is the line integral from
    % [grid.lat(ilat),grid.lon(ilon)] to [90,grid.lon(ilon)] of {B dot rhat dsin(LAT)}
    settings = newsettings; % success, so we can re-use grid if same settings are requested.
end

% prepare for lon integral
dlon = [lon(2)-lon(end);lon(3:end)-lon(1:(end-2));lon(1)-lon(end-1)]; % centered longitude spacing, with wrap
dlon = rem(2*360+180+dlon,360)-180; % deal with mod 360
dlon = dlon/2; % dlon is half of distance between neighbors i+1 and i-1
dlon = pi/180*dlon; % convert to radians
lon = rem(lon+360*2,360); % force onto 0,360 interval
% interpolate partialPhi (lat line integral) onto lat/lon's
partialPhi = interp2(grid.LAT,grid.LON,grid.partialPhi,lat,lon,'linear');
% perform longitude integral
Phi = dlon'*partialPhi; % northern hemisphere integral is negative B.r<0

function foot_points = Phi_contours(LAT,LON,I,I0,periodic,digits)
% general purpose contouring and interpolating routine
% LAT,LON are N1xN2
% I is N1xN2xN3, gives I as a function of LAT,LON
% I0 is N3x1, a set of values of I0 at which to find contours
% periodic defines whether i1 or i2 is periodic, e.g., [false,true]
% digits is inf for no rounding, an integer>0 for rounding
% foot_points is a N3x1 cell array of {lat,lon}
% where lat and lon are Nx1, and provide the contours

[N1,N2] = size(LAT);
N3 = length(I0);
foot_points = cell(N3,1);
% facilitate periodic interpolation
CLON = cosd(LON);
SLON = sind(LON);

I = round_to_digits(I,digits);
I0 = round_to_digits(I0,digits);

for i3 = 1:N3
    I3 = I(:,:,i3); % matrix to contour
    if ~isfinite(I0(i3)) || ~any(isfinite(I3(:)))
        continue;
    end
    C = contourc(I3,I0(i3)+[0 0]);
    %         C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
    %              pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]
    if isempty(C)
        continue;
    end
    ipair = 1;
    lat = nan;
    lon = nan;
    while ipair < size(C,2)
        Npairs = C(2,ipair);
        xy = C(:,ipair+(1:Npairs))'; %[i2,i1]
        ipair = ipair+Npairs+1; % next
        if Npairs == 1
            continue;
        end
        % remove repeats
        ikeep = setdiff((1:Npairs)',find(all(diff(xy,1,1)==0)));
        xy = xy(ikeep,:);
        
        % check for closed contour
        closed = false;
        if isequal(xy(1,:),xy(end,:)) % 1st and last match, so first is redundant
            closed = true;
            xy = xy(2:end,:); % remove repeated point     
        elseif periodic(1) && all(ismember([1,N1],xy(:,2)))
            closed = true; % end-to-end in i1
        elseif periodic(2) && all(ismember([1,N2],xy(:,1)))
            closed = true; % end-to-end in i2
        end
        if closed
            lat = interp2(1:N2,1:N1,LAT,xy(:,1),xy(:,2),'linear');
            % interpolate periodic using sin/cos
            clon = interp2(1:N2,1:N1,CLON,xy(:,1),xy(:,2),'linear');
            slon = interp2(1:N2,1:N1,SLON,xy(:,1),xy(:,2),'linear');
            lon = atan2(slon,clon)*180/pi;
            break;
        end
    end % while ipair
    if all(isfinite(lat)) && all(isfinite(lon)) && (length(lon)>=2)
        foot_points{i3} = [lat,lon];
    end
end % for i3

function traces = Phi_at_fixed_i3(traces,yvar,periodic,digits,verbose)
% yvar is Bm for traces at fixed K, yvar is I for traces at fixed Bm
last_t = now;
[N1,N2] = size(traces);
% northern foot points at R=1
LAT = nan([N1,N2]); % geographic latitude, degrees
LON = nan([N1,N2]); % geographic longitude, degrees
% I vs i1, i2, z
y = [];
for i1 = 1:N1
    for i2 = 1:N2
        if ~isempty(traces{i1,i2})
            if isempty(y)
                N3 = length(traces{i1,i2}.(yvar));
                y = nan([N1,N2,N3]);
            end
            traces{i1,i2}.Phi = nan(N3,1);
            traces{i1,i2}.Lstar = nan(N3,1);
            traces{i1,i2}.foot_points = cell(N3,1);
            y(i1,i2,:) = traces{i1,i2}.(yvar);
            LAT(i1,i2) = traces{i1,i2}.NorthPoint(1);
            LON(i1,i2) = traces{i1,i2}.NorthPoint(2);
        end % if ~isempty
    end % for i2
end % for i1

% make contours
for i1 = 1:N1
    for i2 = 1:N2
        if ~isempty(traces{i1,i2}) && isfinite(LAT(i1,i2)) && isfinite(LON(i1,i2))
            y0 = squeeze(y(i1,i2,:));
            traces{i1,i2}.foot_points = Phi_contours(LAT,LON,y,y0,periodic,digits);
            if (verbose>=1) && (now-last_t>1/24/60/60)
                fprintf('Traced foot points for i=(%d,%d): lat = %g, lon = %g)\n',i1,i2,LAT(i1,i2),LON(i1,i2));
                last_t = now;
            end
        end % if ~isempty
    end % for i2
end % for i1

function traces = Phi_at_fixed_a0(traces,periodic,digits,verbose)
last_t = now;
[N1,N2] = size(traces);
% northern foot points at R=1
LAT = nan([N1,N2]); % geographic latitude, degrees
LON = nan([N1,N2]); % geographic longitude, degrees
% I,Bm vs i1, i2, z
I = [];
Bm = [];
for i1 = 1:N1
    for i2 = 1:N2
        if ~isempty(traces{i1,i2})
            if isempty(I)
                N3 = length(traces{i1,i2}.alpha0deg);
                I = nan([N1,N2,N3]);
                Bm = nan([N1,N2,N3]);
            end
            traces{i1,i2}.Phi = nan(N3,1);
            traces{i1,i2}.Lstar = nan(N3,1);
            traces{i1,i2}.foot_points = cell(N3,1);
            I(i1,i2,:) = traces{i1,i2}.I;
            Bm(i1,i2,:) = traces{i1,i2}.Bm;
            LAT(i1,i2) = traces{i1,i2}.NorthPoint(1);
            LON(i1,i2) = traces{i1,i2}.NorthPoint(2);
        end % if ~isempty
    end % for i2
end % for i1

% now have I(LAT,LON,a0) and Bm(LAT,LON,a0)

% make contours of constant Bm at fixed I
for i1 = 1:N1
    for i2 = 1:N2
        if ~isempty(traces{i1,i2}) && isfinite(LAT(i1,i2)) && isfinite(LON(i1,i2))
            I0 = squeeze(I(i1,i2,:));
            Bm0 = squeeze(Bm(i1,i2,:));
            ib = find(isfinite(Bm0) & isfinite(I0)); % logical index into finite Bm0
            if isempty(ib)
                continue;
            end
            I0 = I0(ib); % only interpolate to finite Bm0
            tmpBm = nan(size(Bm)); % I(LAT,LON,Bm0)
            for j1 = 1:N1
                for j2 = 1:N2
                    xI = squeeze(I(j1,j2,:));
                    yBm = squeeze(Bm(j1,j2,:));
                    f = isfinite(xI) & isfinite(yBm);
                    Nf = sum(f);
                    if Nf>=2
                        tmpBm(j1,j2,ib) = interp1(xI(f),yBm(f),I0,'linear');
                    elseif (Nf == 1) && any(xI(f)==I0) % singleton, exact match
                        tmpBm(j1,j2,ib(xI(f)==I0)) = yBm(f);
                    end
                end
            end
            traces{i1,i2}.foot_points = Phi_contours(LAT,LON,tmpBm,Bm0,periodic,digits);
            if (verbose>=1) && (now-last_t>1/24/60/60)
                fprintf('Traced foot points for i=(%d,%d): lat = %g, lon = %g)\n',i1,i2,LAT(i1,i2),LON(i1,i2));
                last_t = now;
            end
        end % if ~isempty
    end % for i2
end % for i1

function y = round_to_digits(x,digits)

y = x;
if isfinite(digits) && (digits>0)
    iz = (x ~= 0); % false for nan
    s = sign(y(iz)); % store sign
    y(iz) = abs(y(iz)); % remove sign
    p = floor(log10(y(iz))); % power of 10 represented by y
    p = p-digits+1; % allow for digits
    y(iz) = 10.^p.*round(10.^(-p).*y(iz)); % shift by power of 10
    y(iz) = y(iz).*s; % restore sign
end


function [Lstar,traces] = test

% run a short demo using all 3 grid types on a LAT/LON grid at R=1.0
[LAT,LON] = ndgrid(35:2:80,0:15:359); % finer grid for demo
% [LAT,LON] = ndgrid(40:4:80,0:30:359); % coarser grid for debugging
if true % dipole
    field_model = 'dipole';
    options = {'DIPOLE'};
    kext = '';
    maginput = [];
else % T89, Kp=2
    field_model = 'T89, Kp=2'; %#ok<UNRCH>
    options = {};
    kext = 'T89';
    maginput = onera_desp_lib_maginputs(2); % Kp=2
end
R = ones(size(LAT));
args = {kext,options,'RLL',datenum(2010,1,1),R,LAT,LON,maginput,'R0',0.9,'verbose',inf,'periodic2'};

zvar = {'a0','Bm','K'}; % test all three options for definining mirror point grid
for ivar = 1:length(zvar)
    extras = {}; % extra arguments specific to zvar
    switch(zvar{ivar})
        case 'K' % K grid
            z = [0;0.1;0.3;1.0]; % sqrt(G)*RE
            extras{end+1} = 'G'; % force Gauss in Bunit
            iplot = 1:length(z); % which z to plot
        case 'a0' % a0 grid
            z = 2:2:90; % deg
            iplot = length(z):-2:1; % which z to plot
        case 'Bm' % Bmirror grid
            z = [0.3e3;1e3;3e3;10e3;30e3]; % nT
            iplot = 1:length(z); % which z to plot
    end

    % trace!
    [Lstar,traces] = onera_desp_lib_grid_lstar(args{:},zvar{ivar},z,extras{:});
    for iz = iplot % make contour plots
        figure;
        contourf(LON,LAT,Lstar(:,:,iz),1:0.25:10); % contours in L*
        axis([0 360 0 90]);
        xlabel('Longitude, ^o East');
        ylabel('Latitude, ^o North');
        title(sprintf('%s, %s=%g',field_model,zvar{ivar},z(iz)));
        cb = colorbar('vert');
        ylabel(cb,'L*');
    end
end
