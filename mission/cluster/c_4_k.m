function [k1,k2,k3,k4]=c_4_k(r1,r2,r3,r4)
%C_4_K calculate reciprocal vectors in barycentric coordinates
%
%  [k1,k2,k3,k4]=C_4_K(r1,r2,r3,r4)
%  r1..r4 are row vectors with the satellite positions,
%  if r1..r4 have more than 3 columns then only columns 2-4 are used,
%  column 1 is assumed to be time and it is added as column 1 to k1..k4
%  r1...r4 should be of the same size and sampled at the same time instants
%
%  [k]=c_4_k(r1,r2,r3,r4) does the same as above but k(1,:,:)=k1;k(2,:,:)=k2;...
%
%  K = C_4_K(R) K and R can be also structures, in this case
%  K.C1=k1,..R.C1=r1,...
%
%  The units of reciprocal vectors are the same as [1/r]
%  If r is in [km] and you want k in [1/m] then you have to divide
%  the obtained values of k by 10^3

%  Reference: ISSI book 14.7
%  k4=cross(r12,r13)/(dot(r14,cross(r12,r13))   r12=r2-r1;
%  k1=cross(r23,r24)/(dot(r21,cross(r23,r24))

%% Prepare input
isOutputStructure = false;
if nargin==0
	help c_4_k;return;
elseif nargin == 1
	R=r1;
	if isstruct(R) && all(isfield(R,{'C1','C2','C3','C4'}))
		isOutputStructure = true;
	else
		errStr = 'wrong syntax';
		irf.log('critica',errStr);
		error(errStr);
	end
elseif nargin == 4
	R.C1 = r1;R.C2 = r2;R.C3 = r3;R.C4 = r4;
else
	errStr = 'wrong syntax';
	irf.log('critica',errStr);
	error(errStr);
end

if size(R.C1,2)>3
	isTimeSpecified=true;
	tVec = R.C1(:,1);
	R.C1(:,1) = [];
	R.C2(:,1) = [];
	R.C3(:,1) = [];
	R.C4(:,1) = [];
else
	isTimeSpecified=false;
end  % check if first column of r1..r4 is time

%% Computation
k.C1 = zeros(size(R.C1));
k.C2 = k.C1; k.C2 = k.C1; k.C4 = k.C1;

id = {'C1','C2','C3','C4','C1','C2','C3'};
for j=1:4
	cc        = cross(R.(id{2+j})-R.(id{1+j}),...
		              R.(id{3+j})-R.(id{1+j}),2);
	dr12      = R.(id{j})-R.(id{1+j});
	denom     = dot(dr12,cc,2);
	k.(id{j}) = [cc(:,1)./denom cc(:,2)./denom cc(:,3)./denom];
end

%% Output
if isTimeSpecified
	for j=1:4
		k.(id{j}) = [tVec k.(id{j})];
	end
end

if isOutputStructure
	k1 = k;
elseif nargout == 1 % output is one big matrix
	K=zeros([4 size(k.C1)]);
	for j=1:4
		K(j,:,:)=k.(id{j});
	end
	k1=K;
elseif nargout==4
	k1=k.C1;k2=k.C2;k3=k.C3;k4=k.C4;
elseif nargout==0
	if size(r1,1)>1; disp('Reciprocal vectors for the first data point');end
	for j=1:4
		strk=['k' num2str(j) '=' num2str(norm(k.(id{j})(1,:)),3)...
			' [ ' num2str(k.(id{j})(1,:)/norm(k.(id{j}(1,:))),...
			' %5.2f') '] '];
		disp(strk);
	end
else
	disp('Check number of output arguments. See usage: help c_4_k');
end
