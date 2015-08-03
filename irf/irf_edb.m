function [ed,d]=irf_edb(e,b,angle_lim,flag)
%IRF_EDB   Compute Ez under assumption E.B=0 or E.B~=0
%
% [ed,d]=IRF_EDB(e,b,angle_lim) Calculate Ez under assumption E.B=0
%          e - [Ex Ey] or [t Ex Ey] or [t Ex Ey Ez] (Ez whatever)
%          b - [Bx By Bz] or [t Bx By Bz]
%  angle_lim - B angle with respect to the spin plane should be at least
%              angle_lim degrees otherwise Ez is set to 0 or to NaN
%              if flag 'Eperp+NaN' i sgiven.
%         ed - E field output
%          d - B elevation angle above spin plane
%
% [ed,d]=IRF_EDB(e,b,angle_lim,'Epar') Calculate Ez under assumption that
%           the measured electric field along the B projection comes from Eparallel
%  angle_lim - B angle with respect to the spin plane should be less than
%              angle_lim degrees otherwise Ez is set to 0.
%
% $Id$

%temporary fix, convert TS series to older format, output unchanged
if isa(e,'TSeries'), 
    ttemp = e.time.epochUnix;
    datatemp = double(e.data);
    e = [ttemp, double(datatemp)];
end

if isa(b,'TSeries'), 
    ttemp = b.time.epochUnix;
    datatemp = double(b.data);
    b = [ttemp, datatemp];
end
%end of temporary fix


flag_method='E.B=0'; % default method for Ez calculation
defaultValue = 0;
if nargin==0, help irf_edb;return;end
if nargin == 2,
    angle_lim=20;
    irf_log('fcal','Using limiting angle of 20 degrees');
end
if nargin==4,
    if strcmpi(flag,'epar'),
        flag_method='Epar';
    elseif strcmpi(flag,'Eperp+NaN')
        defaultValue = NaN;
    end
end

le= size(e,2); %
lb= size(b,2); %
if size(b,1) ~= size(e,1),
    irf_log('fcal','E and B are not of the same length. Interpolating B.');
    b=irf_resamp(b,e);
end

if le < 2
    error('E has not enough components');return;
elseif le == 2
    ed=[e(:,1) e(:,2) e(:,1)*defaultValue];
elseif le == 3
    ed=[e(:,2) e(:,3) e(:,1)*defaultValue];
else
    ed=[e(:,2) e(:,3) e(:,1)*defaultValue];
end

if lb == 2
    disp('B has not enough components');return;
elseif lb == 3
    bd=[b(:,1) b(:,2) b(:,3)];
elseif lb > 3
    bd=[b(:,2) b(:,3) b(:,4)];
end


switch lower(flag_method)
    case 'e.b=0'
        % Calculate using assumption E.B=0
        
        d=atan2(bd(:,3),sqrt(bd(:,1).^2+bd(:,2).^2))*180/pi;
        ind=find(abs(d)>angle_lim);
        
        if ~isempty(ind)
            ed(ind,3)=-(ed(ind,1).*bd(ind,1)+ed(ind,2).*bd(ind,2))./bd(ind,3);
        end
        
        tt=sprintf('E.B=0 was used for %0.0f from %0.0f data points',length(ind),length(d));
        irf_log('fcal',tt);
        
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
  
        tt=sprintf('E.B=Emeasured.Bproj_to_spin_plane was used for %0.0f from %0.0f data points',length(ind),length(d));
        irf_log('fcal',tt);
        
    otherwise
        irf_log('fcal','called with unknown input');
        ed=[];d=[];
end
if le > 2
    ed=[e(:,1) ed];
end
