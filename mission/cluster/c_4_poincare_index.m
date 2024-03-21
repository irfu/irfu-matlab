function index=c_4_poincare_index(B1,B2,B3,B4)
%C_4_POINCARE_INDEX    Calculates the Poincare index
%
%  C_4_POINCARE_INDEX calculates the Poincare index from B-fields
%  measurements. If index is a nonzero number than there is uneven number
%  of magnetic nulls within tetrahedron and if the index is zero then there
%  is zero or even number of nulls inside the volume.
%
%  index=C_4_POINCARE_INDEX(B1,B2,B3,B4)
%
%  INPUT
%    B1..B4 are magnetic fields by satellites. They are column vectors
%     B1 = [[time] B1x B1y B1z], ...
%
%  OUTPUT
%    index  column vector [[time] index]
%
% Important: This method only works for 3D nulls, i.e. when none of the
% eigenvalues in deltaB vanish.
%
% See Also IRF.SOLIDANGLE
% Reference: Greene 1992 JCP (98) p.194-198
%% Check TSeries
idB = {B1,B2,B3,B4};

if isa(B1,'TSeries')
  for i=1:4
    if isa(idB{i},'TSeries')
      idB{i} =  [idB{i}.time.epochUnix double(idB{i}.data)];
    end
  end
end
%% Check inputs
nCol=size(idB{1,1},2);
nRow=size(idB{1,1},1);
if nCol<3
  error('Time tag must be included in each input vector. Please do so and try again.')
elseif nCol==3
  if size(idB{1,2},1)==nRow && size(idB{1,3},1)==nRow && size(idB{1,4},1)==nRow
    isTimeSpecified = false;
  else
    error('The vectors needs to be the same size or have the time vector in them. See usage: help c_4_poincare_index')
  end
else
  isTimeSpecified = true;
end
if nargin==0
  help c_4_poincare_index;
  return;
elseif nargin < 4
  error('Too few input values. See usage: help c_4_poincare_index')
elseif nargin>4
  error('Too many input values. See usage: help c_4_poincare_index')
end

%% Treat the case when the time column is specified
if isTimeSpecified
  %Each vector contains time tags so all vectors needs to be resampled to
  %establish synchronisation between the S/C's
  idB{1,2}=irf_resamp(idB{1,2},idB{1,1});
  idB{1,3}=irf_resamp(idB{1,3},idB{1,1});
  idB{1,4}=irf_resamp(idB{1,4},idB{1,1});

  %Renaming of vectors for simplicty sake and removing time for calculations
  %using B?
  time=idB{1,1}(:,1); %time= first column of SC1
  idB{1,1}=idB{1,1}(:,2:4); %vec1= 2 to 4th column (Bx,By,Bz) for s/c 1 (if that's the one you placed in position 1)
  idB{1,2}=idB{1,2}(:,2:4); %vec2= 2 to 4th column (Bx,By,Bz)
  idB{1,3}=idB{1,3}(:,2:4); %vec3= 2 to 4th column (Bx,By,Bz)
  idB{1,4}=idB{1,4}(:,2:4); %vec4= 2 to 4th column (Bx,By,Bz)
end

%%
lx=size(idB{1,1},1); %lx becomes the number of rows in B1

%Create the zero matrix that will be used to map each point from xyz to
%By,Bx,Bz space

Map_sc1= zeros(lx,3);
Map_sc2= zeros(lx,3);
Map_sc3= zeros(lx,3);
Map_sc4= zeros(lx,3);
% map the points from x,y,z to magnetic three-dimensional field space (Bx,By,Bz instead of
% x,y,z) (ref. Greene, J.M. 1990)
% A null point in configuration space (xyz satellites) corresponds to the
% origin in M space.

%Mapping of the B fields for each s/c but as a unit vector (length=1)
%Calculate the length of each vector
vl1= sqrt(dot(idB{1,1},idB{1,1},2)); %norm (length of vector)=sqrt(dot(vec1,vec1,2)).
%dot(vec1,vec1,2) treats each row as a vector in the matrix so A=dot(vec1,vec1,2)
% would give A(1,:) dot product of vec1(1,:) and vec1(1,:)
vl2= sqrt(dot(idB{1,2},idB{1,2},2));
vl3= sqrt(dot(idB{1,3},idB{1,3},2));
vl4= sqrt(dot(idB{1,4},idB{1,4},2));

%% Unit vectors is used in solid angle so we need divide the vector with
%its norm (magnitude of the vector/length)
for i=1:3   %Each column on Map_sc1 is given unit vector (first column in 3D divided by the length of vector in the
  % direction of each s/c)
  Map_sc1(:,i)=idB{1,1}(:,i)./vl1; %Needs to use the for loop so that matrix dimensions agrees
  Map_sc2(:,i)=idB{1,2}(:,i)./vl2;
  Map_sc3(:,i)=idB{1,3}(:,i)./vl3;
  Map_sc4(:,i)=idB{1,4}(:,i)./vl4;
end

%% Calculate the solid angles
% Calculate the solid angles of an unit sphere that has taken the sign into
%account to give the number of nulls for each triangles with the direction
%going counterclockwise seen from origin which determines which sc you
%choose as point A,B,C see solidangle for more details
%In a tetrahedron you have four triangles so you need the solid angle for
%each corner
area1=irf.solidangle(Map_sc1,Map_sc2,Map_sc3);
area2=irf.solidangle(Map_sc1,Map_sc4,Map_sc2);
area3=irf.solidangle(Map_sc1,Map_sc3,Map_sc4);
area4=irf.solidangle(Map_sc2,Map_sc4,Map_sc3);
%The Poincar? index is the total area (including each sign) divided by 4pi
index=(area1+area2+area3+area4)/(4*pi);

if isTimeSpecified
  index=[time index];
end
end
