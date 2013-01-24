function [result,b]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,option)
%C_4_GRAD calculate gradient of physical field  using 4 spacecraft technique
%  [grad_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4)
%  [grad_b,[b]]=c_4_grad('r?','b?','grad')
%  [curl_b,[b]]=c_4_grad('r?','b?','curl')
%  [div_b,[b]]=c_4_grad('r?','b?','div')
%  [curv_b,[b]]=c_4_grad('r?','b?','curvature')
%  [div_T,[b]]=c_4_grad('r?','b?','divT')
%  [b_div_b,[b]]=c_4_grad('r?','b?','bdivb')
%  [drift_grad,[b]]=c_4_grad('r?','b?','drift_grad')
%  [drift_curl,[b]]=c_4_grad('r?','b?','drift_curl')
%
%  r1..r4 are row vectors
%         column 1     is time
%         column 2-4   are satellite position in km
%         if 3 columns, assume x,y,z positions (b should have the same number of rows)
%  b1..b4 are row vectors of physical field (usually magnetc field)
%         column 1     is time, b1 time is used for all interpolation of r1..r4 and b2..b4
%         column 2-4   vector components (if scalar then only one component)
%  curl_b is row vector,
%         column 1     time
%         column 2-4   curl b
%  div_b  column 1     time
%         column 2     div(B)
%  curv_b is row vector, curvature is defined as (bhat div)bhat, where bhat=norm(b) 
%         column 1     time
%         column 2-4   curv b
%  b_div_b is row vector, defined as (b div) b
%         column 1     time
%         column 2-4   curv b
%  grad_b structure (answer is tensor if b? is vector)
%         first element - row vector with time
%         second element - row tensor
%         grad_b{2}(1,1:3,1:3) - first tensor ..
%  grad_b is row vector if b is scalar,
%         column 1     time
%         column 2-4   grad_b
%  div_T, divergence of stress tensor, T=BiBj - delta_ij B^2/2
%         T=(div B)B+(curl B)xB
%         div_T is row vector and defined for b being vector
%         column 1     time
%         column 2-4   div_T
%  drift_grad gradient drift of 1eV positively charged particle in [m/s]
%  drift_curv curvature drift of 1eV positively charged particle in [m/s]
%
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17 p. 353
%
% $Id$

%%%%%%%%%%%%%%%%%%  Check input parameters %%%%%%%%%%%%%%

if nargin~=9 && nargin~=8 && nargin~=3 && nargin~=2
	disp('Wrong number of input parameters. See usage:');
	help c_4_grad;     
	return
end

% first identify what should be calculated 
flag_option='grad'; % default to calculate grad
if nargin==9,
  if ischar(option),
    flag_option=option;
  end
elseif nargin==8 || nargin == 2,
  flag_option='grad';
elseif nargin==3, % if 3 input arguments, the 3rd is option
  flag_option=r3;
else
  warning('c_4_grad','Too many input parameters.'); return;
end

if nargin==2 || nargin==3, % input is in form 'R?' and 'B?'
	if ischar(r1) && ischar(r2),
		rs = r1;
		bs = r2;
		for cl_id=1:4
			ttt = evalin('caller',irf_ssub(rs,cl_id));
			c_eval('r? =ttt;',cl_id); clear ttt
			ttt = evalin('caller',irf_ssub(bs,cl_id));
			c_eval('b? =ttt;',cl_id); clear ttt
		end
		clear bs rs
	else
		irf_log('fcal','For two input parameters, both should be strings');
		return;
	end
end

%%%%%%%%%%%%%%%%%%  Check input - vector or scalar %%%%%%%%%%%%%%
isTimeSpecified = true; % default assume first column is time
isField='scalar'; % default assume scalar input
if size(b1,2)>=4, % input is vector using only columns 2,3,4
  isField='vector';
elseif size(b1,2) == 3, % assume vector components, time not specified
	isField='vector';
	isTimeSpecified = false;
	if isequal(size(b1,1),size(b2,1),size(b3,1),size(b4,1))
		if isequal(size(r1,1),size(r2,1),size(r3,1),size(r4,1))
			if size(r1,2) == 3 % no time column
				t=1:size(r1,1);
				t=t(:);
				c_eval('r?=[t r?];');
			end
			c_eval('b?=[t b?];');
		else
			irf_log('fcal','ERROR: input vectors not equal');
			return
		end
	else
		irf_log('fcal','ERROR: input vectors not equal');
		return
	end
end

%%%%%%%%%%%%%%%% Estimate first reciprical coordinates %%%%%%%%%%%%%%
%
% because usually r1..r4 is of less time resolution, it is more
% computer friendly first calculate k1..k4 and only after interpolate
% and not the other way around
c_eval('R?=irf_resamp(r?,r1,''spline'');')
[k1,k2,k3,k4]=c_4_k(R1,R2,R3,R4);

%%%%%%%%%%%%%%%% Do interpolation to b1 time series %%%%%%%%%%%%%%%%%%%%%%
c_eval('B?=irf_resamp(b?,b1);');
b=0.25*B1+0.25*B2+0.25*B3+0.25*B4; % estimate mean value of vector or scalar
c_eval('K?=irf_resamp(k?,b1);');

%%%%%%%%%%%%%%%% Calculate gradient if necessary (grad,curvature) %%%%%%%%%%%%%%%%%%%%%%
if strcmp(flag_option,'grad')||strcmp(flag_option,'curvature')||strcmp(flag_option,'bdivb'),
  if strcmp(isField,'scalar'),  % scalar field, gradient is vector
    grad_b=zeros(size(B1,1),4);grad_b(:,1)=b1(:,1);
    c_eval('grad_b(:,2:4)=grad_b(:,2:4)+K?(:,2:4).*repmat(B?(:,2),1,3);');
    result=grad_b;
  elseif strcmp(isField,'vector'), % vector field, gradient is matrix 1->(1,1),2->(1,2),3>(1,3),...
    grad_b=zeros(size(B1,1),1,10);grad_b(:,1)=b1(:,1);
    grad_b(:,2) = K1(:,2).*B1(:,2)+K2(:,2).*B2(:,2)+K3(:,2).*B3(:,2)+K4(:,2).*B4(:,2);
    grad_b(:,3) = K1(:,2).*B1(:,3)+K2(:,2).*B2(:,3)+K3(:,2).*B3(:,3)+K4(:,2).*B4(:,3);
    grad_b(:,4) = K1(:,2).*B1(:,4)+K2(:,2).*B2(:,4)+K3(:,2).*B3(:,4)+K4(:,2).*B4(:,4);
    grad_b(:,5) = K1(:,3).*B1(:,2)+K2(:,3).*B2(:,2)+K3(:,3).*B3(:,2)+K4(:,3).*B4(:,2);
    grad_b(:,6) = K1(:,3).*B1(:,3)+K2(:,3).*B2(:,3)+K3(:,3).*B3(:,3)+K4(:,3).*B4(:,3);
    grad_b(:,7) = K1(:,3).*B1(:,4)+K2(:,3).*B2(:,4)+K3(:,3).*B3(:,4)+K4(:,3).*B4(:,4);
    grad_b(:,8) = K1(:,4).*B1(:,2)+K2(:,4).*B2(:,2)+K3(:,4).*B3(:,2)+K4(:,4).*B4(:,2);
    grad_b(:,9) = K1(:,4).*B1(:,3)+K2(:,4).*B2(:,3)+K3(:,4).*B3(:,3)+K4(:,4).*B4(:,3);
    grad_b(:,10)= K1(:,4).*B1(:,4)+K2(:,4).*B2(:,4)+K3(:,4).*B3(:,4)+K4(:,4).*B4(:,4);
  else
    irf_log('fcal','error: input vector is neither scalar or vector');
    return
  end
end

% estimate divergence and curl in all cases (they are used in most cases)
if strcmp(isField,'vector'), 

  % estimate divergence
  div_b=zeros(size(B1,1),2);div_b(:,1)=b1(:,1);
  c_eval('div_b(:,2)=div_b(:,2)+dot(K?(:,2:4),B?(:,2:4),2);');
  
  % estimate curl  
  curl_b=zeros(size(B1,1),4);curl_b(:,1)=b1(:,1);
  c_eval('curl_b(:,2:4)=curl_b(:,2:4)+cross(K?(:,2:4),B?(:,2:4),2);');
  
end


switch flag_option
  case 'grad'
    % gradients are already calculated so do not do anything
    result=grad_b;
  case 'curl'
    if strcmp(isField,'vector'),
      result=curl_b;
    end
  case 'div'
    if strcmp(isField,'vector'),
      result=div_b;
    end
  case 'bdivb' % TODO: test that gives correct values
    if strcmp(isField,'vector'),
      curv=zeros(size(B1,1),4);curv(:,1)=b1(:,1);
      curv(:,2)=dot(b(:,2:4),grad_b(:,[2 5 8]),2);
      curv(:,3)=dot(b(:,2:4),grad_b(:,[3 6 9]),2);
      curv(:,4)=dot(b(:,2:4),grad_b(:,[4 7 10]),2);
      clear grad_b;
      result=curv; 
    end
  case 'curvature' % 
    c_eval('bhat?=irf_norm(b?);');
    curvature=c_4_grad('r?','bhat?','bdivb');
    result=curvature;
  case 'drift_grad' % TODO: test that gives correct values
    c_eval('babs?=[b?(:,1) irf_abs(b?,1)];');
    gradB=c_4_grad('r?','babs?','grad');
    result=irf_tappl(irf_vec_x_scal(irf_cross(b,gradB),[b(:,1) irf_abs(b,1)],-3),'*1e9/1e3');
  case 'drift_curv' % TODO: test that gives correct values
      % curvature drift v=W[eV]/(B^2)* (Bvec x div_par bhat) assuming
      % B[nT],r[km] (*1e9 because of nT and /1e3 because grad in km)
    c_eval('bhat?=irf_norm(b?);');
    curvature=c_4_grad('r?','bhat?','bdivb');
    result=irf_tappl(irf_vec_x_scal(irf_cross(b,curvature),[b(:,1) irf_abs(b,1)],-2),'*2*1e9/1e3');
  case 'div_T' % TODO: test that gives correct values
    if strcmp(isField,'vector'),
      div_T=irf_add(1,irf_vec_x_scal(b,div_b),1,irf_cross(curl_b,b));
      result=div_T;  
    end
  otherwise
    irf_log('fcal','warning: unknown input option');
end

if ~isTimeSpecified % remove time column if time was not given
	result(:,1)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
%  if nargout==0&size(B1,1)==1,
%         strj=['j= ' num2str(norm(j(1,2:4)),3) ' [ ' num2str(j(1,2:4)/norm(j(1,2:4)),' %5.2f') '] A '];
%         strdivB=['divB= ' num2str(divB(1,2),3) '] A '];
%         disp(strj);disp(strdivB);
%  end

