function add_position(h, r, varargin)
%ADD_POSITION  add xlabels of s/c position in RE
%
% ADD_POSITION(H,R)  Add position labels to axis with handle H
%		R can be TSeries or matrix with 1st column epochtime, 2-4th column XYZ
%		coordinates and 5th column absolute distance to s/c.
%   'obj'       'Earth', 'Mars', ...
%   See also IRF_TIMEAXIS

% narginchk(2, 2)

Units = irf_units;
obj_radius = Units.RE/1e3;              % Default: Earth radius in [km];
obj_label = '[R_E]';

if ~ishandle(h), error('H is not an axis handle'), end
if isempty(r), irf_log('func','empty position'), return, end

args = varargin;
if isempty(args)
  have_options = 0;
else
  have_options = 1;
end

while have_options
  l = 1;
  switch(lower(args{1}))
    case 'obj'
      l = 2;
      obj = args{2};
  end
  args = args(l+1:end);
  if isempty(args), break, end
end


switch lower(obj)
  case 'mercury'
    obj_radius = Units.Mercury.radius/1e3;
    obj_label = '[R_M]';
  case 'venus'
    obj_radius = Units.Venus.radius/1e3;
    obj_label = '[R_V]';
  case 'earth'
    obj_radius = Units.RE/1e3;
    obj_label = '[R_E]';
  case 'mars'
    obj_radius = Units.Mars.radius/1e3;
    obj_label = '[R_M]';
  case 'jupiter'
    obj_radius = Units.Jupiter.radius/1e3;
    obj_label = '[R_J]';
  case 'saturn'
    obj_radius = Units.Saturn.radius/1e3;
    obj_label = '[R_S]';
  case 'ganymede'
    obj_radius = Units.Ganymede.radius/1e3;
    obj_label = '[R_G]';
  case 'europa'
    obj_radius = Units.Europa.radius/1e3;
    obj_label = '[R_e]';
  case 'callisto'
    obj_radius = Units.Callisto.radius/1e3;
    obj_label = '[R_c]';
  otherwise
    return;
end

if isa(r,'TSeries'), r = irf.ts2mat(r); end
if size(r,2)==4
  r = irf_abs(r);
elseif size(r,2)~=5
  error('R has bad size')
end

TF = issorted(r(:, 1), 'strictascend');             % check time is in ascending order;
if ~TF
  r = unique(r, 'rows');          % remove duplicated data;
  irf.log('warning', 'Time is not in strict ascending order, UNIQUE is used. ');
end

irf_timeaxis(h,'usefig',[r(:,1) r(:,2:end)/obj_radius],...
  {['X ', obj_label], ['Y ', obj_label], ['Z ', obj_label], ['R ', obj_label]});

%%