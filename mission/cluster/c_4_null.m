function Nulls=c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,varargin)
%C_4_NULL - finds magnetic nulls within a s/c box
%
%   C_4_NULL finds magnetic nulls inside the square box defined by 4 s/c
%   location. In addition, the properties of magnetic nulls, such as
%   magnetic null type, eigenvalues, eigenvectors .., are calculated.
%		The routine assumes that B-field is varying linearly.
%
%   Nulls=C_4_NULL(R1,R2,R3,R4,B1,B2,B3,B4,[threshold],[boxLim])
%
%
%   INPUT
%     B1..B4    - are B field times series measured by satellites C1..C4.
%                 columns are [t Bx By Bz]
%     R1..R4    - are position of C1..C4, columns are [t Rx Ry Rz]
%     threshold - returned are only nulls satisfying
%                 |divB/max(all eigenvalues of gradB)| < threshold
%                 Default threshold value is 40%, threshold is given in percent units.
%                 If threshold is 100% all points are accepted.
%     boxLim    - only nulls from box size less than boxLim are returned.
%                 Default is to have no limits on box size.
%
%   OUTPUT
%   Nulls.(fields) - Structure having different fields characterizing the
%                      found magnetic nulls.
%
%   Nulls.t
%   Nulls.R    [xn yn zn] - the position of the null.
%   Nulls.dR   [dxn dyn dzn] - the distance from the
%      spacecraft tetrahedron centrum to the magnetic null.
%   Nulls.j =  [jx jy jz] - the total current for each null
%   Nulls.PoincareIndex = [t PI] - the Poincar? index for each time tag.
%   Nulls.error = [err] - for each found null returns
%       err = |divB/max(all eigenvalues of gradB)|
%
%   Nulls.gradB  gradB for each time tag in the format given by c_4_grad function.
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
%   of a null.
%

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
	error('You need the time vector for each input [t x y z]')
end
if nargin == 0
	help c_4_null;
	return;
elseif nargin < 8
	error('Too few input values. See usage: help c_4_null')
elseif nargin > 10
	error('Too many input values. See usage: help c_4_null')
end

if isempty(varargin) == true
	threshold = 40;
    noBox=true;
elseif length(varargin)<2
	threshold = varargin{1};
    noBox=true;
else
    threshold = varargin{1};
    boxLim=varargin{2};
    noBox=false;
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
j     = curlB./1.0e3.*1e-9./(4*pi*1e-7);

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
	dR(i,:) = Bmean(i,:)/deltaBnull;
end

R = Rmean - dR;


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
% Calculate the minimum and maximum values for all s/c's in each direction
% to see the distance among s/c's
minX = min(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
maxX = max(([idR{1,1}(:,1) idR{1,2}(:,1) idR{1,3}(:,1) idR{1,4}(:,1)]),[],2);
minY = min(([idR{1,1}(:,2) idR{1,2}(:,2) idR{1,3}(:,2) idR{1,4}(:,2)]),[],2);
maxY = max(([idR{1,1}(:,2) idR{1,2}(:,2) idR{1,3}(:,2) idR{1,4}(:,2)]),[],2);
minZ = min(([idR{1,1}(:,3) idR{1,2}(:,3) idR{1,3}(:,3) idR{1,4}(:,3)]),[],2);
maxZ = max(([idR{1,1}(:,3) idR{1,2}(:,3) idR{1,3}(:,3) idR{1,4}(:,3)]),[],2);

%For each eigenvalue corresponding to the tolerance level (the two errors
%less or equal to 40%) break out their corresponding time and dR value 
%(the minimum distance from all s/c to the null)
disp('Sorting based on the null located within the s/c box made up of the maximum and minimum values for each direction of all satellites');
if noBox
sortNullDx = R(:,1) >= minX & R(:,1) <= maxX;
sortNullDy = R(:,2) >= minY & R(:,2) <= maxY;
sortNullDz = R(:,3) >= minZ & R(:,3) <= maxZ;
okNulls    = sortNullDx & sortNullDy & sortNullDz;   
else
sortsizedX = (maxX-minX) <= boxLim;
sortsizedY = (maxY-minY) <= boxLim;
sortsizedZ = (maxZ-minZ) <= boxLim;
sortNullDx = R(:,1) >= minX & R(:,1) <= maxX;
sortNullDy = R(:,2) >= minY & R(:,2) <= maxY;
sortNullDz = R(:,3) >= minZ & R(:,3) <= maxZ;
okNulls    = sortsizedX & sortsizedY & sortsizedZ & sortNullDx & sortNullDy & sortNullDz;
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
