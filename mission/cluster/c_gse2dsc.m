function [y] = c_gse2dsc( x, spin_axis, direction, db )
% C_GSE2DSC  Convert vector between GSE and DSC reference systems.
% replaced by C_COORD_TRANS
% will be removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%        direction = 1 to convert from GSE into DSC (default).
%                    -1 to convert from DSC into GSE.
%                    2 convert from GSE to DSI
%                    -2 convert from DSI to GSE
%        sc        = spacecraft number.
%        db        = isdat database pointer, that is db = Mat_DbOpen(DATABASE)
%
%     Assume the spin orientation does not change significantly during the
%     choosen interval. Only values at start time point is used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('WARNING!!!!! will be removed. use c_coord_trans instead.')

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
     clear spin_axis;
     if exist('./maux.mat','file')
         irf.log('notice','Loading maux.mat file');
         try c_eval('load maux sc_at?_lat__CL_SP_AUX sc_at?_long__CL_SP_AUX; lat=sc_at?_lat__CL_SP_AUX; long = sc_at?_long__CL_SP_AUX;', ic);
         catch
             irf.log('warning','Loading maux.mat file failed'); flag_read_isdat=1;
         end
         if flag_read_isdat==0 % if reading maux file suceeded
             tmin = lat(1,1);
             tmax = lat(end,1);
             if (t > tmin) || (t < tmax)
                 eval( irf_ssub('load maux sc_at?_lat__CL_SP_AUX sc_at?_long__CL_SP_AUX; lat=sc_at?_lat__CL_SP_AUX; long = sc_at?_long__CL_SP_AUX;', ic) );
                 latlong   = irf_resamp([lat long(:,2)],t);
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
     if flag_read_isdat  % try if there is SP CDF file, otherwise continue to isdat
         cdf_files=dir(['CL_SP_AUX_' irf_time(t,'epoch>yyyymmdd') '*']);
         switch numel(cdf_files)
             case 1
                 cdf_file=cdf_files.name;
                 irf.log('warning',['converting CDF file ' cdf_file ' -> maux.mat']);
                 cdf2mat(cdf_file,'maux.mat');
                 irf.log('warning',['Loading from CDF file:' cdf_file '. Next time will use maux.mat']);
                 c_eval('lat=irf_cdf_read(cdf_file,{''sc_at?_lat__CL_SP_AUX''});',ic);
                 c_eval('long=irf_cdf_read(cdf_file,{''sc_at?_long__CL_SP_AUX''});',ic);
                 if (t > lat(1,1)) && (t < lat(end,1))
                     flag_read_isdat=0;
                     latlong   = irf_resamp([lat long(:,2)],t);
                 end
         end
     end
     if flag_read_isdat  % try if there are auxilaruy files from CAA
       caa_load CL_SP_AUX
       if exist('CL_SP_AUX','var') % succeeded to load CAA data files
         c_eval('lat=getmat(CL_SP_AUX,''sc_at?_lat__CL_SP_AUX'');',ic);
         c_eval('long=getmat(CL_SP_AUX,''sc_at?_long__CL_SP_AUX'');',ic);
         if (t > lat(1,1)-60) && (t < lat(end,1)+60)
           flag_read_isdat=0;
           latlong   = irf_resamp([lat long(:,2)],t);
         end
       end
     end
     if flag_read_isdat  % try if there are SAX variables in mEPH
         if exist('./mEPH.mat','file')
             try
                 c_eval('load mEPH SAX?',ic);
                 c_eval('spin_axis=SAX?;',ic);
                 flag_read_isdat=0;
             catch
                 irf.log('warning','no SAX variable in mEPH, tryinig to read isdat');
             end
         end
     end
     if flag_read_isdat  % load from isdat satellite ephemeris
      irf.log('notice','loading spin axis orientation from isdat database');
       start_time=t; % time of the first point
       Dt=600; % 10 min, in file they are saved with 1 min resolution
        if flag_db==0 % open ISDAT database
					DB_S = c_ctl(0,'isdat_db');
          irf.log('debug',['Starting connection to ' DB_S]);
          db = Mat_DbOpen(DB_S);
        end
        [tlat, lat] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_lat__CL_SP_AUX'], ' ', ' ',' ');
        [~, long] = isGetDataLite( db, start_time, Dt, 'CSDS_SP', 'CL', 'AUX', ['sc_at' num2str(ic) '_long__CL_SP_AUX'], ' ', ' ',' ');
        xxx=[double(tlat) double(lat) double(long)];
        if isempty(xxx), y=NaN; return;
        else
          latlong=xxx(1,:);
        end
        irf.log('notice',['lat=' num2str(latlong(2)) '  long=' num2str(latlong(3))]);
        if flag_db==0
          Mat_DbClose(db);
        end
     end
     if ~exist('spin_axis','var')
         [xspin,yspin,zspin]=sph2cart(latlong(3)*pi/180,latlong(2)*pi/180,1);
         spin_axis=[xspin yspin zspin];
     end
  end

spin_axis=spin_axis/norm(spin_axis);
%if debug_flag, disp('Spin axis orientation');spin_axis, end

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

if direction == 1   % GSE -> DSC
 out=M*inp';
 out=out';
 if length(out(:,1))==1
  irf.log('notice',sprintf('x,y,z = %g, %g, %g [DSC]',out(1), out(2),out(3)));
 end
elseif direction == 2 % GSE -> DSI
 out=M*inp';
 out=out';
 out(:,2)=-out(:,2);
 out(:,3)=-out(:,3);
 if length(out(:,1))==1
  irf.log('notice',sprintf('x,y,z = %g, %g, %g [DSI]',out(1), out(2),out(3)));
 end
elseif direction==-1  % DSC -> GSE
 out=Minv*inp';
 out=out';
 if length(out(:,1))==1
  irf.log('notice',sprintf('x,y,z = %g, %g, %g [GSE]',out(1), out(2),out(3)));
 end
elseif direction==-2   % DSI -> GSE
 inp(:,2)=-inp(:,2);
 inp(:,3)=-inp(:,3);
 out=Minv*inp';
 out=out';
 if length(out(:,1))==1
  irf.log('notice',sprintf('x,y,z = %g, %g, %g [GSE]',out(1), out(2),out(3)));
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


