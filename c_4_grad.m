function [grad_b,b]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,option)
%c_4_j  Calculate gradient using 4 spacecraft technique
%  [grad_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4)  
%  [grad_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,'grad')  
%  [curl_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,'curl')
%  [div_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,'div')  
%  [curv_b,[b]]=c_4_grad(r1,r2,r3,r4,b1,b2,b3,b4,'curvature')  
%
%  r1..r4 are row vectors
%         column 1     is time
%         column 2-4   are satellite position in km
%  b1..b4 are row vectors
%         column 1     is time, b1 time is used for all interpolation of r1..r4 and b2..b4
%         column 2-4   vector components (if scalar then only one component)
%  curl_b is row vector,
%         column 1     time
%         column 2-4   curl b
%  div_b  column 1     time
%         column 2     div(B)
%  grad_b structure (answer is tensor if b? is vector)
%         first element - row vector with time 
%         second element - row tensor 
%         grad_b{2}(1,1:3,1:3) - first tensor ..
%  grad_b is row vector if b is scalar,
%         column 1     time
%         column 2-4   grad_b
%
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17
%
if nargin<8;    warning('c_4_grad','Too few input parameters. See usage:');help c_4_grad;     return;end
if nargin=9,
  if ischar(option),
    flag_option=option;
  end
elseif nargin=8,
  flag_option='grad';
else
  warning('c_4_grad','Too many input parameters.'); return;
end

%%%%%%%%%%%%%%%% Estimate first reciprical coordinates %%%%%%%%%%%%%%
%
% because usually r1..r4 is of less time resolution, it is more
% computer friendly first calculate k1..k4 and only after interpolate
% and not the other way around
for ic=1:4,eval(av_ssub('R?=av_interp(r?,r1,''spline'');',ic)),end
[k1,k2,k3,k4]=c_4_k(R1,R2,R3,R4);

%%%%%%%%%%%%%%%% Do interpolation to b1 time series %%%%%%%%%%%%%%%%%%%%%%
for ic=1:4,eval(av_ssub('B?=av_interp(b?,b1);',ic)),end
for ic=1:4,eval(av_ssub('K?=av_interp(k?,b1);',ic)),end

%%%%%%%%%%%%%%%% Calculate gradient %%%%%%%%%%%%%%%%%%%%%%
if size(b1,2)=2 % scalar field, gradient is vector
  grad_b=zeros(size(B1,1),4);
  flag_vector_input=0;
  c_eval('grad_b(:,2:4)=grad_b(:,2:4)+K?(:,2:4).*repmat(B?(:,2),1,3);',ic)); 
elseif size(b1,2)=4, % vector field, gradient is matrix
  grad_b_temp=zeros(size(B1,1),3,3);
  flag_vector_input=1;
  for j=1:size(B1,1)
    c_eval('grad_b_temp(j,:,:)=grad_b_temp(j,:,:)+K?(j,2:4)''*B?(j,2:4);');
  end
  grad_b={B1(:,1) grad_b_temp};
else
  warning('c_4_grad','input vector is neither scalar or vector');
  return
end

switch flag_option
  case 'grad'
    %
    grad_b=zeros(size(B1,1),4);
    c_eval('grad_b(:,2:4)=grad_b(:,2:4)+K?(:,2:4).*repmat(B?(:,2),3);'); 
  case 'curl'
    if flag_vector_input=1,
      grad_b=zeros(size(B1,1),4);
      c_eval('grad_b(:,2:4)=grad_b(:,2:4)+cross(K?(:,2:4),B?(:,2:4),2);'); 
    end
  case 'div'
    if flag_vector_input=1
      grad_b=zeros(size(B1,1),4);
      c_eval('grad_b(:,2:4)=grad_b(:,2:4)+dot(K?(:,2:4),B?(:,2:4),2);'); 
    end
  case 'curvature'
    if flag_vector_input=1
      grad_b=zeros(size(B1,1),4);
      for j=1:size(B1,1),
        c_eval('grad_b(j,2:4)=grad_b(j,2:4)+B?(j,2:4)*grad_b_temp(j,:,:);'); 
      end
    end
  otherwise
  warning('c_4_grad','unknown input option');
end

% initialize matrix j and divB with right time column
j=B1(:,1:4);j(:,2:4)=0;
divB=B1(:,1:2);divB(:,2)=0;

% Calculate j and divB
for ic=1:4, eval(av_ssub(  'divB(:,2)=divB(:,2)+dot(K?(:,2:4),B?(:,2:4),2);'  ,ic));   end
divB(:,2)=divB(:,2)/1.0e3*1e-9/(4*pi*1e-7); % to get right units

for ic=1:4, eval(av_ssub('j(:,2:4)=j(:,2:4)+cross(K?(:,2:4),B?(:,2:4),2);',ic));   end
j(:,2:4)=j(:,2:4)/1.0e3*1e-9/(4*pi*1e-7);   % to get right units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0&size(B1,1)==1,
       strj=['j= ' num2str(norm(j(1,2:4)),3) ' [ ' num2str(j(1,2:4)/norm(j(1,2:4)),' %5.2f') '] A '];
       strdivB=['divB= ' num2str(divB(1,2),3) '] A '];
       disp(strj);disp(strdivB);
end

