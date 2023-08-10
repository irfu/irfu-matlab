function [ed,d]=irf_edb(e,b,angle_lim,flag)
%IRF_EDB   Compute Ez under assumption E.B=0 or E.B~=0
%
% [ed,d]=IRF_EDB(e,b,angle_lim) Calculate Ez under assumption E.B=0
%          e - [Ex Ey] or [t Ex Ey] or [t Ex Ey Ez] (Ez whatever)
%          b - [Bx By Bz] or [t Bx By Bz]
%  angle_lim - B angle with respect to the spin plane should be at least
%              angle_lim degrees otherwise Ez is set to 0 or to NaN
%              if flag 'Eperp+NaN' is given.
%         ed - E field output
%          d - B elevation angle above spin plane
%
% [ed,d]=IRF_EDB(e,b,angle_lim,'Epar') Calculate Ez under assumption that
%           the measured electric field along the B projection comes from Eparallel
%  angle_lim - B angle with respect to the spin plane should be less than
%              angle_lim degrees otherwise Ez is set to 0.

if isa(e,'TSeries')
  if ~isa(b,'TSeries')
    errS = 'Both E and B must be of the same class';
    irf.log('critical',errS), error(errS)
  end
  flagTs = true;
elseif isa(b,'TSeries')
  errS = 'Both E and B must be of the same class';
  irf.log('critical',errS), error(errS)
else, flagTs = false;
end

flag_method='E.B=0'; % default method for Ez calculation
defaultValue = 0;
if nargin==0, help irf_edb;return;end
if nargin == 2
  angle_lim=20;
  irf_log('fcal','Using limiting angle of 20 degrees');
end
if nargin==4
  if strcmpi(flag,'epar'), flag_method='Epar';
  elseif strcmpi(flag,'Eperp+NaN'), defaultValue = NaN;
  end
end

if flagTs
  if length(e) ~= length(b)
    irf.log('warn','E and B are not of the same length. Interpolating B.');
    b = b.resample(e.time);
  end
  le = size(e.data,2); lb = size(b.data,2);
  
  if le < 2 || lb < 3, error('E || B has not enough components'), end
  ed = [e.x.data e.y.data e.y.data*defaultValue];
  bd = [b.x.data b.y.data b.z.data];
else
  if size(b,1) ~= size(e,1)
    irf.log('warn','E and B are not of the same length. Interpolating B.');
    b = irf_resamp(b,e);
  end
  le = size(e,2); lb = size(b,2);
  if le < 2
    error('E has not enough components');
  elseif le == 2
    ed=[e(:,1) e(:,2) e(:,1)*defaultValue];
  elseif le == 3
    ed=[e(:,2) e(:,3) e(:,1)*defaultValue];
  else
    ed=[e(:,2) e(:,3) e(:,1)*defaultValue];
  end
  
  if lb == 2
    error('B has not enough components');
  elseif lb == 3
    bd=[b(:,1) b(:,2) b(:,3)];
  elseif lb > 3
    bd=[b(:,2) b(:,3) b(:,4)];
  end
end

switch lower(flag_method)
  case 'ex=0' % Solar Orbiter
    d=atan2d(bd(:,3),bd(:,2));
    ind=find(abs(d)>angle_lim);
    
    ed(:,3) = ed(:,3)*NaN;
    if ~isempty(ind)
      ed(ind,3)=-(ed(ind,2).*bd(ind,2))./bd(ind,3);
    end
    
    irf.log('notice',sprintf(...
      'E.B=0 (Ex=0) was used for %0.0f from %0.0f data points',...
      length(ind),length(d)));
    
  case 'e.b=0'
    % Calculate using assumption E.B=0
    
    d=atan2(bd(:,3),sqrt(bd(:,1).^2+bd(:,2).^2))*180/pi;
    ind=find(abs(d)>angle_lim);
    
    if ~isempty(ind)
      ed(ind,3)=-(ed(ind,1).*bd(ind,1)+ed(ind,2).*bd(ind,2))./bd(ind,3);
    end
    
    irf.log('notice',sprintf(...
      'E.B=0 was used for %0.0f from %0.0f data points',...
      length(ind),length(d)));
    
  case 'epar'
    % Calculate using assumption that E field along the B projection is
    % coming from parallel electric field
    % Ez=(Ex bx + Ey by)*bz/(bx^2+by^2)
    d=atan2(bd(:,3),sqrt(bd(:,1).^2+bd(:,2).^2))*180/pi;
    ind=find(abs(d)<angle_lim);
    
    if ~isempty(ind)
      ed(ind,3)=(ed(ind,1).*bd(ind,1)+ed(ind,2).*bd(ind,2));
      ed(ind,3)=ed(ind,3).*bd(ind,3)./(bd(ind,1).^2+bd(ind,2).^2);
    end
    
    irf.log('notice',...
      'E.B=Emeasured.Bproj_to_spin_plane was used for %0.0f from %0.0f data points',...
      length(ind),length(d));
    
  otherwise
    error('called with unknown input');
end
if flagTs
  ed = irf.ts_vec_xyz(e.time,ed); ed.units = e.units;
  d = irf.ts_scalar(e.time, d); d.units = 'deg';
else
  if le > 2
    ed=[e(:,1) ed];
  end
end