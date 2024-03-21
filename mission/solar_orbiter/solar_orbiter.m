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
% Some common flags: antenna, spacecraft, plasma
%
% To see all flags:
%   solar_orbiter('?')

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
Units=irf_units;

switch flag
  case '?'         % Show all possible flag options
    fid=fopen(which('solar_orbiter'));
    while 1
      tline = fgetl(fid);
      if ~ischar(tline), break, end
      b=regexp(tline,'case ''(?<flag>\w*)''\s*[%](?<comment>.*)','names');
      if numel(b)==1
        disp(['solar_orbiter(''' b.flag ''') - ' b.comment]);
      end
    end
  case 'antenna'   % antenna properties, can be used for probe current estimates
    probe.type='cylindrical';
    probe.surface='gold';
    probe.radius=0.575*Units.cm;
    probe.length=500*Units.cm;
    probe.cross_section_area=2*probe.radius*probe.length;
    probe.total_area=pi*probe.cross_section_area;
    probe.total_vs_sunlit_area=probe.total_area/probe.cross_section_area;
    %        probe.capacitance=2*pi*Units.eps0*probe.length/log(probe.length/probe.radius); % assuming length >> radius
    probe.capacitance=30e-12; % 30pF from Chris simulations
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
  case 'probe'     % same as 'antenna'
    out=solar_orbiter('antenna');
  case 'spacecraft'% spacecraft properties
    sc.sunlit_area=5.95;
    sc.cross_section_area=sc.sunlit_area;
    sc.total_area=28.11;
    sc.probe_refpot_as_fraction_of_scpot=.2;
    sc.number_of_probes=3;
    sc.antenna_guard_area=0;
    if nargout==0, % display information
      disp('Solar Orbiter spacecraft properties.')
      disp(['Sunlit area  = ' num2str(sc.sunlit_area,3) ' m^2']);
      disp(['Total area   = ' num2str(sc.total_area,3) ' m^2']);
    else % return infomration
      out=sc;
    end
  case 'sc'        % same as 'spacecraft'
    out=solar_orbiter('spacecraft');
  case 'plasma'    % typical plasma parameters
    plasma.perihelion=struct('q',[-1 1],'m',[0 1],'n',[100 100],'T',[10 25],'vsc',[4e5 4e5]);
    plasma.aphelion  =struct('q',[-1 1],'m',[0 1],'n',[ 10  10],'T',[ 4 10],'vsc',[4e5 4e5]);
    out=plasma;
  case 'orbit'
  case 'plot_sc'   % plot Solar Orbiter coordinate system
    fn=figure;
    set(fn,'color','white'); % white background for figures (default is grey)
    axis(6*[-1 1 -1 1 -1 1]), daspect([1 1 1]);
    hold on, grid on
    ref.origo=[5 0 0];
    ref.scale=3;
    ref.x = [1 0 0];
    ref.y = [0 1 0];
    ref.z = [0 0 1];
    c_eval('arrow3(ref.origo, ref.origo+ref.scale*ref.?,[],1,[],1.2);',{'x','y','z'});
    c_eval('qq=ref.origo+ref.scale*ref.?;text(qq(1),qq(2),qq(3),''\bf?'');',{'x','y','z'});
    qq=ref.origo+ref.scale*1.5*ref.x;
    text(qq(1),qq(2),qq(3),'Sun');
    lAnt = 5;
    lAntBoom = 1; % 1m for antenna boom length
    rAnt1 = [1.642 -0.030 1.120];
    rAnt2 = [1.642 0.760 -1.020];
    rAnt3 = [1.642 -0.760 -1.020];
    rAnt1_vec = [0 0 1];
    angle23 = 110; % angle between Ant2 Ant3 after cinlination change
    angleRad = angle23/2*pi/180; % angle from negative Z axis to V2 antenna
    rAnt2_vec = [0 sin(angleRad) -cos(angleRad)]; % see RPW ANT straylight pdf file
    rAnt3_vec = [0 -sin(angleRad) -cos(angleRad)];
    c_eval('rBoomEnd? =     rAnt? + lAntBoom*rAnt?_vec;',[1 2 3]);
    c_eval(' rAntEnd? = rBoomEnd? +     lAnt*rAnt?_vec;',[1 2 3]);
    c_eval('arrow3(rBoomEnd?,rAnt?,''3d'',0,[],1.5);',[1 2 3]);
    c_eval('arrow3(rBoomEnd?,rAntEnd?,''2p'',0,[],0);',[1 2 3]);
    c_eval('text(rAntEnd?(1),rAntEnd?(2),rAntEnd?(3),''\bfV?'');',[1 2 3]);
  otherwise
    disp(['!!! solar_orbiter: unknown flag ''' lower(flag) '''.'])
    out=[];
end
