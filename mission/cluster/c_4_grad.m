function [result,b]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,option)
%C_4_GRAD calculate gradient of physical field using 4 spacecraft technique
%  C_4_GRAD(R1,R2,R3,R4,B1,B2,B3,B4) input positions and magnetic fields
%  values of 4 spacecraft
%  C_4_GRAD(R,B) R and B can be structures R.C1,R.C2,..,B.C1,..
%  C_4_GRAD('R?','B?') implies that R1,R2,..,B1,.. are read in from the calling workspace
%
%  [grad_b,[b]]=C_4_GRAD
%  [grad_b,[b]]=C_4_GRAD(..,'grad')
%  [curl_b,[b]]=C_4_GRAD(..,'curl')
%  [div_b,[b]] =C_4_GRAD(..,'div')
%  [curv_b,[b]]=C_4_GRAD(..,'curvature')
%  [div_T,[b]] =C_4_GRAD(..,'divT')
%  [b_div_b,[b]]=C_4_GRAD(..,'bdivb')
%  [drift_grad,[b]]=C_4_GRAD(..,'drift_grad')
%  [drift_curl,[b]]=C_4_GRAD(..,'drift_curl')
%
%  r1..r4 are row vectors
%         column 1     is time ,if given output will have the 1st column time
%         column 2-4   are satellite position in km
%         if 3 columns, assume x,y,z positions (b should have the same number of rows)
%  b1..b4 are row vectors of physical field (usually magnetc field)
%         column 1     is time, b1 time is used for all interpolation of r1..r4 and b2..b4
%         column 2-4   vector components (if scalar then only one component)
%  curl_b is row vector,
%  div_b  is scalar
%  curv_b is row vector, curvature is defined as (bhat div)bhat, where bhat=norm(b)
%  b_div_b is row vector, defined as (b div) b
%  grad_b structure (answer is tensor if b1..b4 are vectors)
%         first element - row vector with time
%         second element - row tensor
%         grad_b{2}(1,1:3,1:3) - first tensor ..
%  grad_b is row vector if b1..b4 are scalar fields,
%  div_T, divergence of stress tensor, T=BiBj - delta_ij B^2/2
%         T=(div B)B+(curl B)xB
%         div_T is row vector and defined for b1..b4 being vectors
%  drift_grad gradient drift of 1eV positively charged particle in [m/s],
%  assumes that position is given in [km] and magnetic field in [nT]
%  drift_curv curvature drift of 1eV positively charged particle in [m/s],
%  assuming that position is given in [km] and magnetic field in [nT]
%
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17 p. 353

%% Defaults
idC = {'C1','C2','C3','C4'};
doOutputTSeries = false;   % input is not TSeries
%% Check input parameters
if nargin~=9 && nargin~=8 && nargin~=3 && nargin~=2
	disp('Wrong number of input parameters. See usage:');
	help c_4_grad;
	return
end

% first identify what should be calculated
toCalculate='grad'; % default to calculate grad
if nargin==9,
	if ischar(option),
		toCalculate=option;
	end
elseif nargin==8 || nargin==2,
	toCalculate='grad';
elseif nargin==3, % if 3 input arguments, the 3rd is option
	toCalculate=r3;
else
	irf.log('warning','Too many input parameters.');
	return;
end

if nargin==2 || nargin==3, % input is in form 'R?' and 'B?'
	if ischar(r1) && ischar(r2),
		for iC=1:4
			id = idC{iC};
			R.(id) = evalin('caller',irf_ssub(r1,iC));
			B.(id) = evalin('caller',irf_ssub(r2,iC));
            if isa(R.(id),'TSeries')
				doOutputTSeries = true;
				R.(id) = [R.(id).time.epochUnix double(R.(id).data)];
            end
			if isa(B.(id),'TSeries')
				doOutputTSeries = true;
				Bunits = B.(id).units;
				B.(id) = [B.(id).time.epochUnix double(B.(id).data)];
            end
        end
	elseif isstruct(r1) && isstruct(r2)
		R=r1; B=r2;
		for iC=1:4
			id = idC{iC};
			if ~isfield(R,id) || ~isfield(B,id)
				errStr = 'ERROR: input not structures of type R.C1,R.C2,..';
				irf.log('critical',errStr);
				error(errStr);
			else
				if isa(R.(id),'TSeries')
					doOutputTSeries = true;
					R.(id) = [R.(id).time.epochUnix R.(id).data];
				end
				if isa(B.(id),'TSeries')
					doOutputTSeries = true;
					Bunits = B.(id).units;
					B.(id) = [B.(id).time.epochUnix B.(id).data];
				end
			end
		end		
	else
		irf.log('warning','For two input parameters, both should be strings');
		return;
	end
elseif nargin >= 8, 
	R.C1 = r1; R.C2 = r2; R.C3 = r3; R.C4 = r4;
	B.C1 = b1; B.C2 = b2; B.C3 = b3; B.C4 = b4;
end

%%  Check input - vector or scalar
% if time axis supplied, interpolate the vectors
isFieldScalar = false; % default everything false
isFieldVector = false;
isTimeSpecified = false;
if (size(B.C1,2)>=4 || size(B.C1,2)==2) && size(R.C1,2)>3 % input is vector using only columns 2,3,4
	isTimeSpecified = true; % default assume first column is time
	tB = B.C1(:,1); % time vector
	tR = R.C1(:,1); % time vector
	if  size(B.C1,2)>=4,
		isFieldVector = true;
	else
		isFieldScalar = true;
	end
	for iC=1:4
		id=idC{iC};
		ttt      = irf_resamp(B.(id),tB);
		B.(id)   = ttt(:,2:min(4,size(ttt,2))); % remove time column, keep the rest
		ttt      = irf_resamp(R.(id),tR,'spline');
		R.(id)   = ttt(:,2:4);  % remove time column, keep only X,Y,Z coordinates
	end
elseif (size(B.C1,2) == 3 || size(B.C1,2) == 1) ...
		&& size(R.C1,2) == 3 % assume  time not specified
	if ~isequal(numel(B.C1),numel(B.C2),numel(B.C3),numel(B.C4))
		irf.log('critical','ERROR: B input vectors not equal');
		return
	end
	if ~isequal(numel(R.C1),numel(R.C2),numel(R.C3),numel(R.C4))
		irf.log('critical','ERROR: input vectors not equal');
		return
	end
	if ~isequal(size(B.C1,1),size(R.C1,1))
		irf.log('critical','ERROR: R and B not of the same length');
		return
	end
	if size(B.C1,2) == 3 
		isFieldVector   = true;
	else
		isFieldScalar   = true;
	end
end

%% Estimate first reciprical coordinates
%
% because usually r1..r4 is of less time resolution, it is more
% computer friendly first calculate k1..k4 and only after interpolate
% and not the other way around
[K]=c_4_k(R);
%% Do interpolation to b1 time series

b=0.25*B.C1+0.25*B.C2+0.25*B.C3+0.25*B.C4; % estimate mean value of vector or scalar
for iC=1:4
	id=idC{iC};
	if isTimeSpecified
		K.(id) = irf_resamp([tR K.(id)],tB);
		K.(id)(:,1) = []; % remove time column
	end
end

%% Calculate gradient if necessary (grad,curvature) 
if strcmp(toCalculate,'grad')||strcmp(toCalculate,'curvature')||strcmp(toCalculate,'bdivb'),
	if isFieldScalar,  % scalar field, gradient is vector
		gradB=B.C1(:,1)*[0 0 0];
		for iC=1:4
			id=idC{iC};
			gradB = gradB + K.(id).*(B.(id)(:,1)*[1 1 1]);
		end
		result=gradB;
	elseif isFieldVector, % vector field, gradient is matrix 1->(1,1),2->(1,2),3>(1,3),...
		gradB      = B.C1(:,1)*zeros(1,9);
		gradB(:,1) = K.C1(:,1).*B.C1(:,1)+K.C2(:,1).*B.C2(:,1)+K.C3(:,1).*B.C3(:,1)+K.C4(:,1).*B.C4(:,1);%dxBx
		gradB(:,2) = K.C1(:,1).*B.C1(:,2)+K.C2(:,1).*B.C2(:,2)+K.C3(:,1).*B.C3(:,2)+K.C4(:,1).*B.C4(:,2);%dxBy
		gradB(:,3) = K.C1(:,1).*B.C1(:,3)+K.C2(:,1).*B.C2(:,3)+K.C3(:,1).*B.C3(:,3)+K.C4(:,1).*B.C4(:,3);%dxBz
		gradB(:,4) = K.C1(:,2).*B.C1(:,1)+K.C2(:,2).*B.C2(:,1)+K.C3(:,2).*B.C3(:,1)+K.C4(:,2).*B.C4(:,1);%dyBx
		gradB(:,5) = K.C1(:,2).*B.C1(:,2)+K.C2(:,2).*B.C2(:,2)+K.C3(:,2).*B.C3(:,2)+K.C4(:,2).*B.C4(:,2);%dyBy
		gradB(:,6) = K.C1(:,2).*B.C1(:,3)+K.C2(:,2).*B.C2(:,3)+K.C3(:,2).*B.C3(:,3)+K.C4(:,2).*B.C4(:,3);%dyBz
		gradB(:,7) = K.C1(:,3).*B.C1(:,1)+K.C2(:,3).*B.C2(:,1)+K.C3(:,3).*B.C3(:,1)+K.C4(:,3).*B.C4(:,1);%dzBx
		gradB(:,8) = K.C1(:,3).*B.C1(:,2)+K.C2(:,3).*B.C2(:,2)+K.C3(:,3).*B.C3(:,2)+K.C4(:,3).*B.C4(:,2);%dzBy
		gradB(:,9) = K.C1(:,3).*B.C1(:,3)+K.C2(:,3).*B.C2(:,3)+K.C3(:,3).*B.C3(:,3)+K.C4(:,3).*B.C4(:,3);%dzBz
	else
		irf.log('critical','error: input vector is neither scalar or vector');
		return
	end
end

%% Estimate divergence and curl in all cases
if isFieldVector
	divB = B.C1(:,1)*0;
	curlB= B.C1(:,1)*[0 0 0];
	for iC=1:4
		id=idC{iC};
		divB  = divB  + dot(K.(id),B.(id),2);
		curlB = curlB + cross(K.(id),B.(id),2);
	end
end

%% Estimate special cases
switch toCalculate
	case 'grad'
		% gradients are already calculated so do not do anything
		result=gradB;
	case 'curl'
		if isFieldVector
			result=curlB;
		end
	case 'div'
		if isFieldVector
			result=divB;
		end
	case 'bdivb'
		if isFieldVector
			result=[
				dot(b,gradB(:,[1 4 7]),2) ...%2nd col
				dot(b,gradB(:,[2 5 8]),2) ...%3rd col
				dot(b,gradB(:,[3 6 9]),2)]; %4th col
		end
	case 'curvature' %
		for iC=1:4
			id=idC{iC};
			bhat.(id) = irf_norm(B.(id));
		end
		result=c_4_grad(R,bhat,'bdivb');
	case 'drift_grad' % TODO: test that gives correct values
		for iC=1:4
			id=idC{iC};
			babs.(id) = irf_abs(B.(id));
		end
		gradB=c_4_grad(R,babs,'grad');
		result=irf_tappl(irf_vec_x_scal(irf_cross(b,gradB),irf_abs(b,1),-3),'*1e9/1e3');
	case 'drift_curv' % TODO: test that gives correct values
		% curvature drift v=W[eV]/(B^2)* (Bvec x div_par bhat) assuming
		% B[nT],r[km] (*1e9 because of nT and /1e3 because grad in km)
		for iC=1:4
			id=idC{iC};
			bhat.(id) = irf_norm(B.(id));
		end
		curvature=c_4_grad(R,bhat,'bdivb');
		result=irf_tappl(irf_vec_x_scal(cross(b,curvature,2),irf_abs(b,1),-2),'*2*1e9/1e3');
	case 'div_T' % TODO: test that gives correct values
		if isFieldVector
			divT=irf_add(1,irf_vec_x_scal(b,divB),1,cross(curlB,b,2));
			result=divT;
		end
	otherwise
		irf.log('warning','warning: unknown input option');
end

%% Prepare output
if isTimeSpecified % add time if given
	if doOutputTSeries
        torder = 0;
        if length(result(1,:))==3; torder = 1; end
		result = TSeries(EpochUnix(tB),result,'TensorOrder',torder);
		if nargout == 2,
            torder = 0;
            if length(b(1,:))==3; torder = 1; end
			b = TSeries(EpochUnix(tB),b,'TensorOrder',torder);
			result.units = Bunits;
		end
	else
		result = [tB result];
	end
end