function [y] = R_c_gse2dsc( x, spin_axis, direction, db )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the original function is c_gse2dsc()
% modified by Robert in the summer of -03
% the change is made in the choice of databas: from disco:10 to unix:99
% the opening and closing of the db has been changed by db as input
%
%  Usage:
%  function [out] = c_gse2dsc( inp, spin_axis, [direction],[db])
%  function [out] = c_gse2dsc( inp, [isdat_epoch sc], [direction],[db])
%  function [out] = c_gse2dsc( inp, sc, [direction],[db])
%
%     Convert vector from GSE into DSC reference system.
%     From STAFF manual:
%        inp, out - vectors with 3 components,
%                   inp(:,1) is X,  inp(:,2) is Y ...
%        if more than 3 columns then columns
%                   inp(:,2) is X, inp(:,3) is Y ...
%        spin_axis = vector in GSE or ISDAT epoch.
%        direction = -1 to convert from DSC into GSE.
%        sc        = spacecraft number.
%        db        = isdat database pointer, that is db = Mat_DbOpen(DATABASE)
%
%     Assume the spin orientation does not change significantly during the
%     choosen interval. Only values at start time point is used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug_flag=0;
flag_read_isdat=0;

if nargin <  2, disp('Not enough arguments'); help c_gse2dsc; return; end
if nargin <  3, direction=1;                                          end
if nargin == 4, flag_db=1; else, flag_db=0;                           end

% Spin_axis gives only s/c number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(spin_axis,2) == 1
  ic        = spin_axis;
  t         = x(1,1);
  clear spin_axis;
  spin_axis = [t ic];
end

% Use time to get spin_axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(spin_axis,2) == 2
  t  = spin_axis(1);
  ic = spin_axis(2);
  if exist('./maux.mat','file')
    if debug_flag, disp('Using maux.mat file');end %#ok<UNRCH>
    load maux Epoch__CL_SP_AUX;
    tlat   = (Epoch__CL_SP_AUX-62167219200000)/1000;   % isdat time
    clear Epoch__CL_SP_AUX;
    x1=fromepoch(tlat(1,1));x2=fromepoch(t);
    if x1(1:3) == x2(1:3)  % maux file is from the right day
      tmin = tlat(1);
      tmax = tlat(end);
      if (t > tmin) || (t < tmax)
        eval( av_ssub('load maux sc_at?_lat__CL_SP_AUX sc_at?_long__CL_SP_AUX; lat=sc_at?_lat__CL_SP_AUX; long = sc_at?_long__CL_SP_AUX;', ic) );
        latlong   = av_interp([tlat lat long],t);
      end
    else  % maux file from the wrong day
      disp('c_gse2dsc() OBS!!!  maux.mat from the wrong date');
      disp('                    get the right one or delete the existing one');
      disp('            I am getting attitude data from isdat instead');
      flag_read_isdat=1;
    end
  else
    flag_read_isdat=1;
  end
  if flag_read_isdat==1  % load from isdat satellite ephemeris
    if debug_flag, disp('loading spin axis orientation from isdat database');end %#ok<UNRCH>
    start_time=fromepoch(x(1,1)); % time of the first point
    Dt=600; % 10 min, in file they are saved with 1 min resolution
    if flag_db==0 % open ISDAT database disco:10
      if debug_flag, disp('Starting connection to disco:10');end %#ok<UNRCH>
      db = Mat_DbOpen('disco:20');
    end
    [tlat, lat] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_lat__CL_SP_AUX'], ' ', ' ',' ');
    [tlong, long] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_long__CL_SP_AUX'], ' ', ' ',' ');
    xxx=[double(tlat) double(lat) double(long)];
    if isempty(xxx), y=NaN; return;
    else
      latlong=xxx(1,:);
    end
    if debug_flag, disp(['lat=' num2str(latlong(2)) '  long=' num2str(latlong(3))]); end %#ok<UNRCH>
    if flag_db==0
      Mat_DbClose(db);
    end
  end
  [xspin,yspin,zspin]=sph2cart(latlong(3)*pi/180,latlong(2)*pi/180,1);
  spin_axis=[xspin yspin zspin];
end

spin_axis=spin_axis/norm(spin_axis);
if debug_flag, disp('Spin axis orientation');spin_axis, end %#ok<UNRCH>

lx=size(x,2);
if lx > 3
  inp=x(:,[2 3 4]); % assuming first column is time
elseif lx == 3
  inp=x;
else
  disp('too few components of vector')
  exit
end

Rx=spin_axis(1);
Ry=spin_axis(2);
Rz=spin_axis(3);
a=1/sqrt(Ry^2+Rz^2);
M=[[a*(Ry^2+Rz^2) -a*Rx*Ry -a*Rx*Rz];[0 a*Rz	-a*Ry];[Rx	Ry	Rz]];
Minv=inv(M);

if direction == 1
  out=M*inp';
  out=out';
  if length(out(:,1))==1
    if debug_flag == 1,sprintf('x,y,z = %g, %g, %g [DSC]',out(1), out(2),out(3));end
  end
elseif direction==-1
  out=Minv*inp';
  out=out';
  if length(out(:,1))==1
    if debug_flag == 1, sprintf('x,y,z = %g, %g, %g [GSE]',out(1), out(2),out(3));end
  end
else
  disp('No coordinate transformation done!')
end

y=x;
if lx > 3
  y(:,[2 3 4])=out; % assuming first column is time
else
  y=out;
end


