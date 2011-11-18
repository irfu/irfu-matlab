function out = solar_orbiter(varargin)
%SOLAR_ORBITER  common information on Solar Orbiter
%
% SOLAR_ORBITER - return basic information on Solar Orbiter
%
% out=SOLAR_ORBITER(flag) - return specific information defined by flag
%
%  Example:
%   solar_orbiter('antenna')         % display probe information
%   probe=solar_orbiter('antenna');  % return probe information
%
% Some common flags: probe, spacecraft
%
% To see all possibilities execute
%   SOLAR_ORBITER('?')
%
% $Id$

if nargin==0, % return basic information on Solar Orbiter
    fid=fopen(which('solar_orbiter'));
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        b=regexp(tline,'^\s+%-.*','match');
        if numel(b)==1
            disp(b{1});
        end
    end
    %-    Solar Orbiter ESA M1 mission
    %-    Launch 2017.01 (next window 2018.08)
    %-    Cruise phase:   2017.01-2022.12
    %-    Science phase:  2022.12-2024.06
    %-    Extended phase:        -2026.11
    %-    ESA: http://sci.esa.int/solarorbiter
    %-    SWT: http://www.solarorbiter.org/
    %-
    %-    more: help solar_orbiter
    return
end

flag=varargin{1};
irf_units

switch flag
    case '?'      % show all possibel flag options
        fid=fopen(which('solar_orbiter'));
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            b=regexp(tline,'case ''(?<flag>\w*)''\s*[%](?<comment>.*)','names');
            if numel(b)==1
                disp(['solar_orbiter(''' b.flag ''') - ' b.comment]);
            end
        end
    case 'antenna'% antenna properties, can be used for probe current estimates
        probe.type='cylindrical';
        probe.surface='gold';
        probe.radius=0.575*Units.cm;
        probe.length=500*Units.cm;
        probe.cross_section_area=2*probe.radius*probe.length;
        probe.total_area=pi*probe.cross_section_area;
        probe.total_vs_sunlit_area=probe.total_area/probe.cross_section_area;
        probe.capacitance=2*pi*Units.eps0*probe.length/log(probe.length/probe.radius); % assuming length >> radius
        if nargout==0, % display information
            disp('Solar Orbiter has 3 cylindrical antennas.')
            disp(['Antenna surface  = ' probe.surface]);
            disp(['Average radius   = ' num2str(probe.radius/Units.cm,3) ' cm']);
            disp(['Antenna length   = ' num2str(probe.length/Units.cm,3) ' cm']);
            disp(['Antenna cross section area = ' num2str(probe.cross_section_area,3) ' m2']);
            disp(['Antenna total section area = ' num2str(probe.total_area,3) ' m2']);
            disp(['Antenna capacitance = ' num2str(probe.capacitance*1e12,3) ' pF']);
        else % return infomration
            out=probe;
        end
    case 'probe'  % same as 'antenna'
        out=solar_orbiter('antenna');
    case 'spacecraft' % spacecraft properties
        sc.sunlit_area=5.95;
        sc.cross_section_area=sc.sunlit_area;
        sc.total_area=28.11;
        sc.probe_refpot_as_fraction_of_scpot=.2;
        sc.number_of_probes=3;
        sc.antenna_guard_area=0;
        out=sc;
    case 'sc'     % same as 'spacecraft'
        out=solar_orbiter('spacecraft');
    case 'plasma' % typical plasma parameters
        plasma.perihelion=struct('q',[-1 1],'m',[0 1],'n',[100 100],'T',[10 25],'vsc',[4e5 4e5]);
        out=plasma;
    otherwise
        disp(['!!! solar_orbiter: unknown flag ''' lower(flag) '''.'])
        out=[];
end
