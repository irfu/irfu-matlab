function magnetic_topology(R1,R2,R3,R4,B1,B2,B3,B4,varargin)
%MAGNETIC_TOPOLOGY - extrapolates and plots magnetic flux tubes topology at a
%given time step
%
%   MAGNETIC_TOPOLOGY extrapolates the magnetic flux tubes by assuming linearly
%   varying magnetic field and plots it for a single time step. The plot is
%   made in the coordinate system of input, so if input is LMN where x is N, y is M and z is L
%   then X axis is N axis, Y axis is M axis, and Z axis is L axis.
%
%   IMPORTANT: Input must either be from just one time
%              step or include the tint variable.
%
%   MAGNETIC_TOPOLOGY(R1,R2,R3,R4,B1,B2,B3,B4,[tint,'BLim', Blim,'Boxsize', boxvalue,'startx', sx,'starty', sy, 'startz', sz,'FluxWidth', Fluxvalue])
%
%
%   INPUT
%     B1..B4    - are B field times series measured by satellites 1..4.
%                 can be in Tseries or vector format (with or without the
%                 time vector).
%     R1..R4    - are position of 1..4 can be in Tseries or vector format (with or without the
%                 time vector)
%     bStart    - B field point from where the interpolation should be
%     calculated in same d units as B1..B4. Give in vector format [x y z].
%     Need to be a value inside the spacecraft tetrahedron.
%     [tint]      - time step to plot given in the format of EpochUnix. If given it must always be given as first variable.
%     ['BLim', Blim]      - removes extrapolated data points with B
%     magnitude larger than Blim. Blim should be given in same units as B1..B4.
%     ['Boxsize', boxvalue]      - give the size of the box where the magnetic
%     flux tubes are plotted use same units as spacecraft position
%     ['startx', sx]      - give the starting points for x can be a
%     single value or a vector ex: 0:2:10
%     ['starty', sy]      - give the starting points for y can be a
%     single value or a vector ex: 0:2:10
%     ['startz', sz]      - give the starting points for z can be a
%     single value or a vector ex: 0:2:10
%     plotting the magnetic flux tubes. Can be given as a single value or
%     as a meshgrid. ex: meshgrid(0:2:10,4,-2:2:10);
%     ['FluxWidth', Fluxvalue]      - give the width of the flux tube (might
%     need to decrease this when you have a lot of flux tubes). Smaller
%     value => thinner tubes.
%
%
%
%
%   OUTPUT
%   Plots the 3D magnetic topology at a specific time.
%
%
% Examples:
%   tint = irf_time('2015-11-30T00:24:22.244870000Z','utc>epoch');
%   irf.magnetic_topology(R1,R2,R3,R4,B1,B2,B3,B4,bStart,tint) - plots the
%   magnetic field topology around the spacecraft tetrahedron center at the
%   specific time given in tint from magnetic field point bStart.
%   irf.magnetic_topology(R1,R2,R3,R4,B1,B2,B3,B4,bStart,tint,'Boxsize',Boxvalue) - plots the
%   magnetic field topology around the spacecraft tetrahedron center at the
%   specific time given in tint and in a box of size +- Boxvalue.
%
%   irf.magnetic_topology(R1,R2,R3,R4,B1,B2,B3,B4,bStart) - plots the
%   magnetic field topology around the spacecraft tetrahedron center at the
%   time of the data. OBS: needs to be data for just one time step.



%% Check TSeries
idR = {R1,R2,R3,R4};
idB = {B1,B2,B3,B4};
timeInd=false; %Is time included in the data
if isa(B1,'TSeries') || isa(R1,'TSeries')
  if size(B1.data,1)>1
    timeInd=true;
  end
  for i=1:4
    if isa(idR{i},'TSeries')
      idR{i} = [idR{i}.time.epochUnix double(idR{i}.data)];

    end
    if isa(idB{i},'TSeries')
      idB{i} =  [idB{i}.time.epochUnix double(idB{i}.data)];
    end
  end

else

  if size(idR{1},2)==4
    timeInd=true;
  end
end
GivenboxWidth=false;
GivenSX=false;
GivenSY=false;
GivenSZ=false;
GiventubeWidth=false;
GivenBLim=false;

if isempty(varargin)
  % Give default values
  SpecificTime=false; %No tint given
end
if isscalar(varargin)
  tint=varargin{1};
  if isa(tint, 'EpochTT')
    tint=tint.epochUnix;
  end
  SpecificTime=true;
elseif ~mod(length(varargin),2)
  % Time to check the combination of elements and change the default values for the given ones.
  for i=1:length(varargin)
    if i+1>length(varargin)
      break;
    else
      %Case 'BLim'
      if strcmp(varargin(i),'BLim')
        if isnumeric(cell2mat(varargin(i+1)))
          Blim=cell2mat(varargin(i+1));
          GivenBLim=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
        %Case 'Boxsize'
      elseif strcmp(varargin(i),'Boxsize')
        if isnumeric(cell2mat(varargin(i+1)))
          boxWidth=cell2mat(varargin(i+1));
          GivenboxWidth=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
        %Case 'startx'
      elseif strcmp(varargin(i),'startx')

        if isnumeric(cell2mat(varargin(i+1)))
          xStart=cell2mat(varargin(i+1));
          GivenSX=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
        %Case 'starty'
      elseif strcmp(varargin(i),'starty')

        if isnumeric(cell2mat(varargin(i+1)))
          yStart=cell2mat(varargin(i+1));
          GivenSY=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
        %Case 'startz'
      elseif strcmp(varargin(i),'startz')

        if isnumeric(cell2mat(varargin(i+1)))
          zStart=cell2mat(varargin(i+1));
          GivenSZ=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
        %Case 'FluxWidth'
      elseif strcmp(varargin(i),'FluxWidth')

        if isnumeric(cell2mat(varargin(i+1)))
          tubeWidth=cell2mat(varargin(i+1));
          GiventubeWidth=true;
        else
          error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
        end
      elseif isnumeric(cell2mat(varargin(i)))
        continue
      else
        error('Unapproved arguments. See usage: help 3D_magnetic_topology')
      end
    end
  end
elseif mod(length(varargin),2)
  % Time to check the combination of elements and change the default values for the given ones.
  if length(varargin)==1 && isa(varargin{1},'Tseries')
    tint=varargin{1};
    if isa(tint, 'EpochTT')
      tint=tint.epochUnix;
    end
    SpecificTime=true;
  elseif length(varargin)==1
    error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
  else
    tint=varargin{1};
    if isa(tint, 'EpochTT')
      tint=tint.epochUnix;
    end
    SpecificTime=true;
    for i=2:length(varargin)
      if i+1>length(varargin)
        break;
      else
        %Case 'BLim'
        if strcmp(varargin(i),'BLim')
          if isnumeric(cell2mat(varargin(i+1)))
            Blim=cell2mat(varargin(i+1));
            GivenBLim=true;

          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
          %Case 'Boxsize'
        elseif strcmp(varargin(i),'Boxsize')
          if isnumeric(cell2mat(varargin(i+1)))
            boxWidth=cell2mat(varargin(i+1));
            GivenboxWidth=true;
          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
          %Case 'startx'
        elseif strcmp(varargin(i),'startx')

          if isnumeric(cell2mat(varargin(i+1)))
            xStart=cell2mat(varargin(i+1));
            GivenSX=true;
          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
          %Case 'starty'
        elseif strcmp(varargin(i),'starty')

          if isnumeric(cell2mat(varargin(i+1)))
            yStart=cell2mat(varargin(i+1));
            GivenSY=true;
          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
          %Case 'startz'
        elseif strcmp(varargin(i),'startz')

          if isnumeric(cell2mat(varargin(i+1)))
            zStart=cell2mat(varargin(i+1));
            GivenSZ=true;
          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
          %Case 'FluxWidth'
        elseif strcmp(varargin(i),'FluxWidth')

          if isnumeric(cell2mat(varargin(i+1)))
            tubeWidth=cell2mat(varargin(i+1));
            GiventubeWidth=true;

          else
            error('Unapproved combination of arguments. See usage: help 3D_magnetic_topology')
          end
        elseif isnumeric(cell2mat(varargin(i)))
          continue
        else
          error('Unapproved arguments. See usage: help 3D_magnetic_topology')
        end
      end
    end
  end
else
  error('Unapproved arguments. See usage: help 3D_magnetic_topology')
end


%% Resample all input vectors to B1 timeline and remove time column
if timeInd
  t = idB{1,1}(:,1);
  idB{1,1}(:,1)=[];
  idB{1,2} = irf_resamp(idB{1,2},t);idB{1,2}(:,1)=[];
  idB{1,3} = irf_resamp(idB{1,3},t);idB{1,3}(:,1)=[];
  idB{1,4} = irf_resamp(idB{1,4},t);idB{1,4}(:,1)=[];
  idR{1,1} = irf_resamp(idR{1,1},t);idR{1,1}(:,1)=[];
  idR{1,2} = irf_resamp(idR{1,2},t);idR{1,2}(:,1)=[];
  idR{1,3} = irf_resamp(idR{1,3},t);idR{1,3}(:,1)=[];
  idR{1,4} = irf_resamp(idR{1,4},t);idR{1,4}(:,1)=[];
else
  t = idB{1,1}(:,1);
  idB{1,1}(:,1)=[];
  idB{1,2}(:,1)=[];
  idB{1,3}(:,1)=[];
  idB{1,4}(:,1)=[];
  idR{1,1}(:,1)=[];
  idR{1,2}(:,1)=[];
  idR{1,3}(:,1)=[];
  idR{1,4}(:,1)=[];
end

Rmean=0.25.*(idR{1,1}+idR{1,2}+idR{1,3}+idR{1,4});
Bmean=0.25.*(idB{1,1}+idB{1,2}+idB{1,3}+idB{1,4});

gradB = c_4_grad(idR{1,1},idR{1,2},idR{1,3},idR{1,4},idB{1,1},idB{1,2},idB{1,3},idB{1,4});

if SpecificTime
  % Picks out the time index closest to the time given in tint
  [~,minpos] = min(abs(t-tint));
  t=t(minpos,1);
  idB{1,1}=idB{1,1}(minpos,:);
  idB{1,2}=idB{1,2}(minpos,:);
  idB{1,3}=idB{1,3}(minpos,:);
  idB{1,4}=idB{1,4}(minpos,:);

  idR{1,1}=idR{1,1}(minpos,:);
  idR{1,2}=idR{1,2}(minpos,:);
  idR{1,3}=idR{1,3}(minpos,:);
  idR{1,4}=idR{1,4}(minpos,:);
  Rmean=Rmean(minpos,:);
  Bmean=Bmean(minpos,:);
  gradB=gradB(minpos,:);
end

gradB = reshape(gradB(1,:),3,3);

if ~GivenboxWidth
  % Give default for box width
  minX = min(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
  maxX = max(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
  boxWidth=maxX-minX;
end
if ~GiventubeWidth
  tubeWidth=100;
end

if ~GivenBLim
  Blim=10000;
end
% Check if boxWidth is given


%Creates the x,y,z values for the box
x=linspace(-boxWidth,+boxWidth,10);
y=linspace(-boxWidth,+boxWidth,10);
z=linspace(-boxWidth,+boxWidth,10);

%Starting points
if ~GivenSX

  xstep=abs(x(1,1)-x(1,end))/10;
  xStart=x(1,2):1.5*xstep:x(1,8);
end
if  ~GivenSY

  ystep=abs(y(1,1)-y(1,end))/10;
  yStart=y(1,2):1.5*ystep:y(1,8);
end
if ~GivenSZ

  zstep=abs(z(1,1)-z(1,end))/10;
  zStart=z(1,2):1.5*zstep:z(1,8);
end


%Check if flag for starting points are given and if they are change to
%those values

% Removes spacecraft center to get seperation distances
idR{1,1}=idR{1,1}-Rmean;
idR{1,2}=idR{1,2}-Rmean;
idR{1,3}=idR{1,3}-Rmean;
idR{1,4}=idR{1,4}-Rmean;

%% Calculate Bfield for each point inside the box volume defined by s/c separation
dR=zeros(1000,3);
B=zeros(1000,3);
for i=1:length(x)
  for j=1:length(y)
    for k=1:length(z)
      ind=sub2ind([10, 10, 10],i,j,k);
      dR(ind,:)=[x(i),y(j),z(k)];
      B(ind,:)=Bmean+(gradB*dR(ind,:)')';
    end
  end
end
Bmag=irf_abs(B,1);
IndexBlim=Bmag>Blim;
B(IndexBlim,:)=NaN;
for l=1:size(dR,1)
  [i,j,k]=ind2sub([10, 10, 10], l);
  Bx(i,j,k)=B(l,1);
  By(i,j,k)=B(l,2);
  Bz(i,j,k)=B(l,3);
end


%% Trace magnetic field lines



[sx,sy,sz] = meshgrid(xStart,yStart,zStart);

% Calculates the field lines
XYZfront=stream3(x,y,z,Bx,By,Bz,sx,sy,sz);
XYZback=stream3(x,y,z,-Bx,-By,-Bz,sx,sy,sz);
% % 2D plane cuts around zero
% [vertices, arrowvertices]=streamslice(x,y,z,Bx,By,Bz,0,0,0,'arrows');

% Go through the field lines going forward and backwards and calculates the magnetic field amplitude at each
% point so that the magnetiude can be used to set the width of the flux
% tubes
k=1;
l=1;
for i=1:size(XYZfront,2)
  if size(XYZfront{i},1)>2 || size(XYZback{i},1)>2
    Bamplfront=[];
    Bamplback=[];
    rlinefront{k}=XYZfront{i};
    rstepSizefront(k,1)=size(XYZfront{i},1);
    rlineback{k}=XYZback{i};
    rstepSizeback(k,1)=size(XYZback{i},1);
    for j=1:rstepSizefront(k,1)
      Blinefront=Bmean+(gradB*XYZfront{i}(j,:)')';
      if irf_abs(Blinefront,1)>Blim
        Bamplfront(j,1)=NaN;
      else
        Bamplfront(j,1)=irf_abs(Blinefront,1);
      end
    end
    for m=1:rstepSizeback(k,1)
      Blineback=Bmean+(gradB*XYZback{i}(m,:)')';
      if irf_abs(Blineback,1)>Blim
        Bamplback(m,1)=NaN;
      else
        Bamplback(m,1)=irf_abs(Blineback,1);
      end
    end

    Bmagfront{k}=Bamplfront;
    Bmagback{l}=Bamplback;
    k=k+1;
    l=l+1;
  end
end



%% Plotting 3D magnetic flux tubes
% S/C position
mms_marker={{'ks','markersize',10},{'rd','markersize',10},...
  {'go','markersize',10,'color',[0 0.6 0]},{'bv','markersize',10}};

SCpos = {idR{1,1},idR{1,2},idR{1,3},idR{1,4}};

irf_plot(1,'newfigure');
hold('on')
for ic=1:4
  % put Cluster markers

  plot3(SCpos{ic}(1),SCpos{ic}(2),SCpos{ic}(3),mms_marker{ic}{:});
end

plot3(0,0,0,'k*') %Center point of tetrahedron

hleg = legend({'MMS1','MMS2','MMS3','MMS4'},'Orientation','horizontal','Location','northoutside','fontsize',11);
hleg.Box = 'off';

daspect([1,1,1])
for i=1:size(rlinefront,2)
  if size(rlinefront{i},1)==1
    continue;
  else
    htubesfront=streamtube({rlinefront{i}},{1./Bmagfront{:,i}},[tubeWidth 20])    ;
    htubesfront.LineStyle='none';
    for j=1:size(htubesfront.XData,1)
      for k=1:size(htubesfront.XData,2)
        dRfront=[htubesfront.XData(j,k), htubesfront.YData(j,k), htubesfront.ZData(j,k)];
        htubesfront.CData(j,k)=irf_abs((Bmean+(gradB*dRfront')'),1);
      end
    end
  end
end

for i=1:size(rlineback,2)
  if size(rlineback{i},1)==1
    continue
  else
    htubesback=streamtube({rlineback{i}},{1./Bmagback{:,i}},[tubeWidth 20])    ;
    htubesback.LineStyle='none';
    for j=1:size(htubesback.XData,1)
      for k=1:size(htubesback.XData,2)
        dRback=[htubesback.XData(j,k), htubesback.YData(j,k), htubesback.ZData(j,k)];
        htubesback.CData(j,k)=irf_abs((Bmean+(gradB*dRback')'),1);
      end
    end
  end
end
view(3)
axis equal
colormap(jet)
hold('off')
h=colorbar;
ylabel(h,'|B| [nT]')
%colormap('jet')
ylabel('Y_{coord} [km]');
zlabel('Z_{coord} [km]');
xlabel('X_{coord} [km]');
if SpecificTime
  tintlab = irf_time(t,'utc');
  title({'3D Magnetic topology',[tintlab(12:23),'UT']});
end

% %% Plotting Streamlines
% fn1=irf_plot(1,'newfigure');
% %set(fn1,'color','white')
% streamline([vertices arrowvertices]);
% view(3)
% axis equal
% ylabel('Y [km]');
% zlabel('Z [km]');
% xlabel('X [km]');
%
% if SpecificTime
% tintlab = irf_time(tint,'utc');
% title({'2D slices of streamlines at tetrahedrons gyrocenter',[tintlab(12:23),'UT']});
% else
%  title('2D slices of streamlines at tetrahedrons gyrocenter')
% end



end

% for k=1:length(rline)
% h=surface([rline{k}(:,1), rline{k}(:,1)],[rline{k}(:,2), rline{k}(:,2)],...
%    [rline{k}(:,3), rline{k}(:,3)],...
%     [Bmagfront{k}, Bmagfront{k}],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',1);
%     %colormap(jet(1:numel(c)))
% end
