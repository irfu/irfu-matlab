function [out,l,v]=av_minvar(inp)
% function [out,l,v]=av_minvar(inp)
% minimum variance analysis
% inp, out - column vectors
% if inp has more than 3 columns (for example first column is time) 
% then assume inp(:,2) is X, inp(:,3) is Y, inp(:,3) is Z 
% l - [l1 l2 l3], eigenvalues where l1 is smallest
% v - column vector with corresponding eigenvectors
% v(1,:) - first vector, ...

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_minvar')

ooo=inp;
qx = size(inp(1,:)); lx=qx(2); % the number of vector components
if lx > 3, inp=inp(:,[2  3 4]);
elseif lx < 3,
 disp('not enough components for x vector');
end


inp_m=mean(inp);
Mm2=inp_m([1 2 3 1 1 2]).*inp_m([1 2 3 2 3 3]);
Mm1=mean(inp(:,[1 2 3 1 1 2]).*inp(:,[1 2 3 2 3 3])); % all 6 elements of triagonal matrix
Mm=Mm1-Mm2;
M=[[Mm(1) Mm(4) Mm(5)]' [Mm(4) Mm(2) Mm(6)]' [Mm(5) Mm(6) Mm(3)]'];

%    [V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
%    full matrix V whose columns are the corresponding eigenvectors so
%    that X*V = V*D.
 
[V,L]=eig(M);
ll=[L(1,1) L(2,2) L(3,3)];
v=V';l=ll;
for i=1:2       % third eigenvector obtain from right hand system
 [y,ii]=max(ll);
 l(i)=ll(ii);
 v(i,:)=V(:,ii)';
 ll(ii)=-1e20;
end
v(3,:)=cross(v(1,:),v(2,:));
l(3)=max(ll);
out=(v*inp')';

if nargout < 1,
disp(strcat('l1=',num2str(l(1),3),' v1=[',num2str(v(1,1),3),', ',num2str(v(1,2),3),', ',num2str(v(1,3),3),']'))
disp(strcat('l2=',num2str(l(2),3),' v2=[',num2str(v(2,1),3),', ',num2str(v(2,2),3),', ',num2str(v(2,3),3),']'))
disp(strcat('l3=',num2str(l(3),3),' v3=[',num2str(v(3,1),3),', ',num2str(v(3,2),3),', ',num2str(v(3,3),3),']'))
end

if lx>3, ooo(:,[2 3 4])=out;out=ooo; end
