function [k1,k2,k3,k4]=c_4_k(r1,r2,r3,r4)
%C_4_K calculate reciprocal vectors in barycentric coordinates
%
%  [k1,k2,k3,k4]=c_4_k(r1,r2,r3,r4)
%  r1..r4 are row vectors with the satellite positions,
%  if r1..r4 have more than 3 columns then only columns 2-4 are used,
%  column 1 is assumed to be time and it is added as column 1 to k1..k4
%  r1...r4 should be of the same size and sampled at the same time instants
%
%  [k]=c_4_k(r1,r2,r3,r4) does the same as above but k(1,:,:)=k1;k(2,:,:)=k2;...
%
%  The units of reciprocal vectors are the same as [1/r]
%  If r is in [km] and you want k in [1/m] then you have to divide
%  the obtained values of k by 10^3
%
% $Id$

%  Reference: ISSI book 14.7
%  k4=cross(r12,r13)/(dot(r14,cross(r12,r13))   r12=r2-r1;
%  k1=cross(r23,r24)/(dot(r21,cross(r23,r24))

if nargin<4;disp('Not enough arguments. See usage:');help c_4_k;return;end

if size(r1,2)>3, flag_timecolumn=1;col=4; else flag_timecolumn=0;col=3; end  % check if first column of r1..r4 is time

R=zeros(4,size(r1,1),3);k=zeros(size(R));
R(1,:,:)=r1(:,col-2:col);R(2,:,:)=r2(:,col-2:col);R(3,:,:)=r3(:,col-2:col);R(4,:,:)=r4(:,col-2:col);R(5,:,:)=R(1,:,:);R(6,:,:)=R(2,:,:);R(7,:,:)=R(3,:,:);

for j=0:3,
    cc       = squeeze(cross((R(3+j,:,:)-R(2+j,:,:)),(R(4+j,:,:)-R(2+j,:,:)),3));
    dr12     = squeeze(R(1+j,:,:)-R(2+j,:,:));
    if size(dr12,2)==1, dr12=dr12';end  % input is only one time point 
    if size(cc,2)==1, cc=cc';end
    denom    = dot(dr12,cc,2); 
    k(1+j,:,:) = [cc(:,1)./denom cc(:,2)./denom cc(:,3)./denom];
end

%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%
if flag_timecolumn==1;    % add time column from r1
  K=zeros(size(k)+[0 0 1]);
  K(:,:,2:4)=k;
  for j=1:4,  K(j,:,1)=r1(:,1);    end
end

if nargout==1,
   k1=K;
elseif nargout==4,
       k1=squeeze(K(1,:,:));k2=squeeze(K(2,:,:));k3=squeeze(K(3,:,:));k4=squeeze(K(4,:,:));
       if size(k1,2)==1,   k1=k1';k2=k2';k3=k3';k4=k4';    end   % in case only one data point is requested squeeze makes column vector therefore result should be transposed
elseif nargout==0,
       if size(r1,1)>1; disp('Reciprocal vectors for the first data point');end
       for j=1:4,
              strk=['k' num2str(j) '=' num2str(norm(squeeze(k(j,1,:))),3) ' [ ' num2str(k(j,1,:)/norm(squeeze(k(j,1,:))),' %5.2f') '] '];
              disp(strk);
       end
else
       disp('Check number of output arguments. See usage:');help c_4_k;
end


