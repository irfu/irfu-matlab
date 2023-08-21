function Nulls=c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,varargin)
%C_4_NULL - finds magnetic nulls within a s/c box
%
%   C_4_NULL finds magnetic nulls inside the square box defined by 4 s/c
%   location. In addition, the properties of magnetic nulls, such as
%   magnetic null type, eigenvalues, eigenvectors .., are calculated.
%		The routine assumes that B-field is varying linearly.
%
%   Nulls=C_4_NULL(R1,R2,R3,R4,B1,B2,B3,B4,arguments)
%
%
%   INPUT
%     B1..B4    - are B field times series measured by satellites 1..4.
%                 That can be in either columns [t Bx By Bz] or Tseries
%     R1..R4    - are position of 1..C, That can be in either columns [t Rx Ry Rz] or Tseries
%     all arguments allowed by C_4_NULL command that can be mixed freely as
%     long as it follows this format for each argument. See Examples below.
%
%     ['strong',currentLim] - returns nulls corresponding with a current
%     larger than currentLim. Default option is all located nulls are returned.
%
%     ['threshold',thresholdValue] - returned are only nulls satisfying
%                 |divB/max(all eigenvalues of gradB)| < thresholdValue.
%                 thresholdValue is given in units of percentage.If
%                 threshold is 100% all points are accepted. Default value
%                 is 40%.
%     ['boxLim',boxLimValue]    - only nulls from box size less and equal to boxLim are returned
%                 (center of box is considered spacecraft tetrahedron center).
%                 Default is to have no limits on box size. Same units as
%                 positional data should be used. If boxLimValue is set to
%                 0 then the maximum and minimum distance in all direction
%                 given by the spacecraft tetrahedron is used.
%
%
%   OUTPUT
%
%   Nulls.t
%   Nulls.R    [xn yn zn] - the position of the null.
%   Nulls.dR   [dxn dyn dzn] - the distance from the
%      spacecraft tetrahedron centrum to the magnetic null.
%   Nulls.j =  [jx jy jz] - the total current for each null.
%         j = curlB/mu0*1E-6 [A/m^2 if B in nT and R in km]
%
%   Nulls.PoincareIndex = [t PI] - the Poincare index for each time tag.
%   Nulls.error = [err] - for each found null returns
%       err = |divB/max(all eigenvalues of gradB)|
%
%   Nulls.gradB  gradB for each time tag in the format given by c_4_grad function
%
%   Nulls.eigenVectors = [eigVect1 eigVect2 eigVect3] contains the
%   eigenvectors for each null, where eigVector1 = [x y z]
%
%   Nulls.eigenValues
%
%   Nulls.eigenvectorCorrespondingToSpine   the number of eigenvector,
%               values can 1..3
%
%   Nulls.Is.A
%   Nulls.Is.As
%   Nulls.Is.B
%   Nulls.Is.Bs
%   Nulls.Is.O
%   Nulls.Is.X
%   Nulls.Is.unknown
%       logical indexes for each type of null. 1(true) is when a null of
%       that type has been found at that time tag.
%
%   Nulls.dB_AB the minimum magnetic field fluctuation that changes the
%   kind of the null from A/As to B/Bs or vice versa.
%   of a null. [nT]
%
% Examples:
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4) - Tries to locate nulls and returns
%   the output for all nulls found. thresholdValue is set to 40%,
%   boxLimValue is no box.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'threshold',thresholdValue) - Tries to
%   locate nulls and returns the output for all nulls found where the
%   threshold is set to thresholdValue. thresholdValue needs to be given in %.
%   boxLimValue is no box.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'boxLim',boxLimValue) - Tries to
%   locate nulls and returns the output for all nulls found where the
%   threshold is set to default. boxLim is set to boxLimValue. boxLim needs
%   to be set with same units as the s/c position.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'strong',currentLim) - Tries to
%   locate nulls with default thresholdValue, no box limit and returns the
%   output for all null corresponding with currents larger than currentLim.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'strong',currentLim,'threshold',thresholdValue) - Tries to
%   locate nulls with a set thresholdValue, no box limit and returns the
%   output for all null corresponding with currents larger than currentLim.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'threshold',thresholdValue,'boxLim',boxLimValue) - Tries to
%   locate nulls with a set thresholdValue and a set box limit and returns the
%   output for all nulls found.
%
%   c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'boxLim',boxLimValue,'strong',currentLim) - Tries to
%   locate nulls with a set box limit and returns the
%   output for all nulls found that has a current larger than currentLim.
%% Check TSeries
idR = {R1,R2,R3,R4};
idB = {B1,B2,B3,B4};

if isa(B1,'TSeries') || isa(R1,'TSeries')
  for i=1:4
    if isa(idR{i},'TSeries')
      idR{i} = [idR{i}.time.epochUnix double(idR{i}.data)];

    end
    if isa(idB{i},'TSeries')
      idB{i} =  [idB{i}.time.epochUnix double(idB{i}.data)];
    end
  end
end

%% Check variables

nCol=size(idB{1,1},2);
if nCol<4
  error('You need the time vector for each input [t x y z] in either Tseries or in vector format')
end
if nargin == 0
  help c_4_null;
  return;
elseif nargin < 8
  error('Too few input values. See usage: help c_4_null')
elseif nargin > 14
  error('Too many input values. See usage: help c_4_null')
end
if isempty(varargin)
  % Gives Default values
  toSaveSpecial=false;
  threshold = 40;
  noBox=true;
  %Check if first variable is one of the argument strings otherwise give
  %error.
elseif ~strcmp(varargin{1},'strong') && ~strcmp(varargin{1},'threshold') && ~strcmp(varargin{1},'boxLim')
  error('Unapproved arguments given. See usage: help c_4_null')
else
  % Gives Default values that can be changed
  toSaveSpecial=false;
  threshold = 40;
  noBox=true;
  % Check values given in varargin
  % All arguments need to include a string+value so if the amount of vargin
  % isn't even then give error
  if length(varargin)==1
    error('Unapproved combination of arguments. See usage: help c_4_null')

  elseif length(varargin)==3
    error('Unapproved combination of arguments. See usage: help c_4_null')
  elseif length(varargin)==5
    error('Unapproved combination of arguments. See usage: help c_4_null')
  else
    % Time to check the combination of elements and change the default values for the given ones.
    for i=1:length(varargin)
      if i+1>length(varargin)
        break;
      else
        %Case 'strong'
        if strcmp(varargin(i),'strong')
          if isnumeric(cell2mat(varargin(i+1)))
            currentLim=cell2mat(varargin(i+1));
            toSaveSpecial=true;
          else
            error('Unapproved combination of arguments. See usage: help c_4_null')
          end
          %Case 'threshold'
        elseif strcmp(varargin(i),'threshold')

          if isnumeric(cell2mat(varargin(i+1)))
            threshold=cell2mat(varargin(i+1));
          else
            error('Unapproved combination of arguments. See usage: help c_4_null')
          end
          %Case 'boxLim'
        elseif strcmp(varargin(i),'boxLim')

          if isnumeric(cell2mat(varargin(i+1)))
            noBox=false;
            boxLim=cell2mat(varargin(i+1));
          else
            error('Unapproved combination of arguments. See usage: help c_4_null')
          end
          % If two strings are after each other give error
        elseif iscellstr(varargin(i)) && iscellstr(varargin(i+1))
          error('Unapproved combination of arguments. See usage: help c_4_null')
          % If to numeric values are after each other give error
        elseif isnumeric(cell2mat(varargin(i)))
          if isnumeric(cell2mat(varargin(i+1)))
            error('Unapproved combination of arguments. See usage: help c_4_null')
          else
            continue;
          end
          % If strings are other than 'threshold', boxLim', and
          % 'strong' give error.
        else
          error('Unapproved arguments. See usage: help c_4_null')
        end
      end
    end
  end
end


%% Resample all input vectors to B1 timeline and remove time column
t = idB{1,1}(:,1);
idB{1,1}(:,1)=[];
nPoints = numel(t);
idB{1,2} = irf_resamp(idB{1,2},t);idB{1,2}(:,1)=[];
idB{1,3} = irf_resamp(idB{1,3},t);idB{1,3}(:,1)=[];
idB{1,4} = irf_resamp(idB{1,4},t);idB{1,4}(:,1)=[];
idR{1,1} = irf_resamp(idR{1,1},t);idR{1,1}(:,1)=[];
idR{1,2} = irf_resamp(idR{1,2},t);idR{1,2}(:,1)=[];
idR{1,3} = irf_resamp(idR{1,3},t);idR{1,3}(:,1)=[];
idR{1,4} = irf_resamp(idR{1,4},t);idR{1,4}(:,1)=[];

%% Calculate 4-s/c measures
% Calculates gradB by assuming linarity between the spacecraft
gradB = c_4_grad(idR{1,1},idR{1,2},idR{1,3},idR{1,4},idB{1,1},idB{1,2},idB{1,3},idB{1,4});
% remove points with NaN magnetic field data
badPoints = any(isnan(gradB(:,2:end)),2);
if any(badPoints)
  t(badPoints)       = [];
  nPoints            = numel(t);
  gradB(badPoints,:) = [];
  idR{1,1}(badPoints,:)=[];
  idR{1,2}(badPoints,:)=[];
  idR{1,3}(badPoints,:)=[];
  idR{1,4}(badPoints,:)=[];
  idB{1,1}(badPoints,:)=[];
  idB{1,2}(badPoints,:)=[];
  idB{1,3}(badPoints,:)=[];
  idB{1,4}(badPoints,:)=[];
end
divB  = c_4_grad(idR{1,1},idR{1,2},idR{1,3},idR{1,4},idB{1,1},idB{1,2},idB{1,3},idB{1,4},'div');
curlB = c_4_grad(idR{1,1},idR{1,2},idR{1,3},idR{1,4},idB{1,1},idB{1,2},idB{1,3},idB{1,4},'curl');
j     = curlB./1.0e3.*1e-9./(4*pi*1e-7); %A/m^2 if B in nT and R in km

%% Calculate the null position
dR  = t*[0 0 0];
dB_AB = t*0;
B41v = idB{1,1}-idB{1,4};
B42v = idB{1,2}-idB{1,4};
B43v = idB{1,3}-idB{1,4};
B12v = idB{1,2}-idB{1,1};
B13v = idB{1,3}-idB{1,1};

%Calculates mean position and magnetic field for the spacecraft tetrahedron
Rmean=0.25.*(idR{1,1}+idR{1,2}+idR{1,3}+idR{1,4});
Bmean=0.25.*(idB{1,1}+idB{1,2}+idB{1,3}+idB{1,4});
% 1) calculate the distance on null from the center of the four satellites

for i=1:nPoints
  deltaBnull = reshape(gradB(i,:),3,3)';
  % distance from tetrahedron center to null
  dR(i,:) = Bmean(i,:)/deltaBnull;
end

R = Rmean - dR; %position of null


% 2) Determinant of gradB
detB = abs(B41v(:,1).*B42v(:,2).*B43v(:,3)...
  + B41v(:,2).*B42v(:,3).*B43v(:,1)...
  + B41v(:,3).*B42v(:,1).*B43v(:,2)...
  - B41v(:,3).*B42v(:,2).*B43v(:,1)...
  - B41v(:,2).*B42v(:,1).*B43v(:,3)...
  - B41v(:,1).*B42v(:,3).*B43v(:,2));

% 3) estimate deltaB_AB
l1vec=cross(B42v,B43v,2);
l2vec=cross(B43v,B41v,2);
l3vec=cross(B41v,B42v,2);
l4vec=cross(B12v,B13v,2);
l1 = sqrt(sum(l1vec.^2,2))./detB;
l2 = sqrt(sum(l2vec.^2,2))./detB;
l3 = sqrt(sum(l3vec.^2,2))./detB;
l4 = sqrt(sum(l4vec.^2,2))./detB;
%Minimum fluctuation to change the sign of volume not satellite
%reference dependent
dB_AB(:,1)=1./(max(([l1 l2 l3 l4]),[],2));

%% Check if box satisfies boxLim conditions and null is within box


%For each eigenvalue corresponding to the tolerance level (the two errors
%less or equal to 40%) break out their corresponding time and dR value
%(the minimum distance from all s/c to the null)
if noBox
  disp('Sorting based on the null located anywhere');
  okNulls    = true(length(R(:,1)));
elseif boxLim==0
  % Calculate the minimum and maximum values for all s/c's in each direction
  % to see the distance among s/c's
  disp('Sorting based on the null inside the size given by s/c tetrahedron');
  minX = min(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
  maxX = max(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
  minY = min(([idR{1,1}(:,2) idR{1,2}(:,2) idR{1,3}(:,2) idR{1,4}(:,2)]),[],2);
  maxY = max(([idR{1,1}(:,2) idR{1,2}(:,2) idR{1,3}(:,2) idR{1,4}(:,2)]),[],2);
  minZ = min(([idR{1,1}(:,3) idR{1,2}(:,3) idR{1,3}(:,3) idR{1,4}(:,3)]),[],2);
  maxZ = max(([idR{1,1}(:,3) idR{1,2}(:,3) idR{1,3}(:,3) idR{1,4}(:,3)]),[],2);
  sortNullDx = R(:,1) >= minX & R(:,1) <= maxX;
  sortNullDy = R(:,2) >= minY & R(:,2) <= maxY;
  sortNullDz = R(:,3) >= minZ & R(:,3) <= maxZ;
  okNulls    = sortNullDx & sortNullDy & sortNullDz;
else
  disp(['Sorting based on the null located inside the box ', num2str(-1*boxLim,3), ' and ',num2str(boxLim,3), ' in each direction']);
  sortNullDx = dR(:,1) >= -1*boxLim & dR(:,1) <= boxLim;
  sortNullDy = dR(:,2) >= -1*boxLim & dR(:,2) <= boxLim;
  sortNullDz = dR(:,3) >= -1*boxLim & dR(:,3) <= boxLim;
  okNulls    = sortNullDx & sortNullDy & sortNullDz;
end

%Removes nulls not within spacecraft box
t(~okNulls)       = [];
gradB(~okNulls,:) = [];
idB{1,1}(~okNulls,:)=[];
idB{1,2}(~okNulls,:)=[];
idB{1,3}(~okNulls,:)=[];
idB{1,4}(~okNulls,:)=[];
dB_AB(~okNulls,:)=[];
R(~okNulls,:)=[];
dR(~okNulls,:)=[];
divB(~okNulls,:)=[];
j(~okNulls,:)=[];
nPoints= numel(t);

%% Calculate eigenvalues and eigenvectors
eigValues   = zeros(nPoints,3);
eigVectors  = zeros(nPoints,9);
eigValErr   = t*NaN;

for i=1:nPoints
  deltaB_null      = reshape(gradB(i,:),3,3);
  [V,D]            = eig(deltaB_null,'vector');
  eigVectors(i,:)  = V(:)';
  eigValues(i,:)   = D';
end

%% Identify points where error larger than threshold
eigValErr(:,2)          = abs(divB./max(abs(eigValues),[],2)) * 100;
badPoints=false(length(eigValErr(:,1)),1);
if threshold == 100
else
  badPoints               = (eigValErr(:,2) > threshold);
end
t(badPoints)       = [];
gradB(badPoints,:) = [];
j(badPoints,:) = [];
R(badPoints,:)=[];
dR(badPoints,:)=[];
dB_AB(badPoints,:)=[];
idB{1,1}(badPoints,:)=[];
idB{1,2}(badPoints,:)=[];
idB{1,3}(badPoints,:)=[];
idB{1,4}(badPoints,:)=[];
eigValues(badPoints,:)  = [];
eigVectors(badPoints,:) = [];
eigValErr(badPoints,:)  = [];


%% Determine spine
% TODO: does this work if sign is equal to zero correctly??
isLambda1Positive     = sign(real(eigValues(:,1))) >= 0;
isLambda2Positive     = sign(real(eigValues(:,2))) >= 0;
isLambda3Positive     = sign(real(eigValues(:,3))) >= 0;

areSignsDifferentL1L2 = (isLambda1Positive ~= isLambda2Positive);
areSignsDifferentL1L3 = (isLambda1Positive ~= isLambda3Positive);
areSignsDifferentL2L3 = (isLambda2Positive ~= isLambda3Positive);

eigenspine1           = areSignsDifferentL1L2 & areSignsDifferentL1L3;
eigenspine2           = areSignsDifferentL1L2 & areSignsDifferentL2L3;
eigenspine3           = areSignsDifferentL1L3 & areSignsDifferentL2L3;


%% Determine the null type by using eigenvalues
%Ideal case
isAllEigenvaluesReal     = abs(max(imag(eigValues),[],2)) == 0;
signOfRealpart           = sign(real(eigValues));
signOfImaginarypart      = sign(imag(eigValues));
nImaginaryNegativeEigenvalues = sum(signOfImaginarypart==-1,2);
nImaginaryPositiveEigenvalues = sum(signOfImaginarypart==+1,2);
nRealNegativeEigenvalues = sum(signOfRealpart==-1,2);
unknowntype              = (nImaginaryNegativeEigenvalues == 3 |sum(signOfImaginarypart,2) ~=0 | nImaginaryPositiveEigenvalues == 3);
aType                    = nRealNegativeEigenvalues == 2 ;
bType                    = nRealNegativeEigenvalues == 1 ;
unknownrealtype          = (nRealNegativeEigenvalues == 3) | (nRealNegativeEigenvalues == 0);
minAbsEigenvalue         = min(abs(eigValues),[],2);
twoDType                 = (minAbsEigenvalue == 0);

Type.A  = ( isAllEigenvaluesReal & aType & ~unknowntype & ~unknownrealtype);
Type.As = (~isAllEigenvaluesReal & aType & ~unknowntype & ~unknownrealtype);
Type.B  = ( isAllEigenvaluesReal & bType & ~unknowntype & ~unknownrealtype);
Type.Bs = (~isAllEigenvaluesReal & bType & ~unknowntype & ~unknownrealtype);
Type.X  = ( isAllEigenvaluesReal & twoDType & ~unknowntype);
Type.O  = (~isAllEigenvaluesReal & twoDType & ~unknowntype);
Type.unknown = (~Type.A & ~Type.B & ~Type.As & ~Type.Bs & ~Type.X & ~Type.O | unknowntype | unknownrealtype);

Nulls.eigenValues=eigValues;
Nulls.eigenVectors=eigVectors;
Nulls.eigenvectorCorrespondingToSpine=[eigenspine1 eigenspine2 eigenspine3];
Nulls.error=eigValErr;

%Calculates Poincare index
index=c_4_poincare_index(idB{1,1},idB{1,2},idB{1,3},idB{1,4});
Nulls.poincareindex       = index;

%Minimum fluctuations
Nulls.dB_AB        = dB_AB;

%Nulltype
Nulls.Is = Type;

%gradB matrix
Nulls.gradB=gradB;
%Current
Nulls.j=j;

%The null positions and time axis
Nulls.R  = R;
Nulls.dR = dR;
Nulls.t  = t;
if sum(Nulls.t)==0
  disp('No Nulls found')
else
  disp('Nulls found')
  if toSaveSpecial
    strongcurrent=(irf_abs(Nulls.j,1)>currentLim);
    if sum(strongcurrent)==0
      disp('No Nulls with strong currents found')
      Nulls.eigenValues=[];
      Nulls.eigenVectors=[];
      Nulls.eigenvectorCorrespondingToSpine=[];
      Nulls.error=[];
      Nulls.poincareindex=[];
      Nulls.dB_AB=[];
      Nulls.Is.A=[];
      Nulls.Is.As=[];
      Nulls.Is.B=[];
      Nulls.Is.Bs=[];
      Nulls.Is.X=[];
      Nulls.Is.O=[];
      Nulls.Is.unknown=[];
      Nulls.gradB=[];
      Nulls.j=[];
      Nulls.R=[];
      Nulls.dR=[];
      Nulls.t=[];
    else
      Nulls.eigenValues=Nulls.eigenValues(strongcurrent,:);
      Nulls.eigenVectors=Nulls.eigenVectors(strongcurrent,:);
      Nulls.eigenvectorCorrespondingToSpine=Nulls.eigenvectorCorrespondingToSpine(strongcurrent,:);
      Nulls.error=Nulls.error(strongcurrent,:);
      Nulls.poincareindex=Nulls.poincareindex(strongcurrent,:);
      Nulls.dB_AB=Nulls.dB_AB(strongcurrent,:);
      Nulls.Is.A=Nulls.Is.A(strongcurrent,:);
      Nulls.Is.As=Nulls.Is.As(strongcurrent,:);
      Nulls.Is.B=Nulls.Is.B(strongcurrent,:);
      Nulls.Is.Bs=Nulls.Is.Bs(strongcurrent,:);
      Nulls.Is.X=Nulls.Is.X(strongcurrent,:);
      Nulls.Is.O=Nulls.Is.O(strongcurrent,:);
      Nulls.Is.unknown=Nulls.Is.unknown(strongcurrent,:);
      Nulls.gradB=Nulls.gradB(strongcurrent,:);
      Nulls.j=Nulls.j(strongcurrent,:);
      Nulls.R=Nulls.R(strongcurrent,:);
      Nulls.dR=Nulls.dR(strongcurrent,:);
      Nulls.t=Nulls.t(strongcurrent,:);
    end
  end
end
end















