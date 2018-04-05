function [out,l,v]=irf_minvar(inp,flag)
%IRF_MINVAR minimum variance analysis
%
% [out,l,v]=irf_minvar(inp)
% inp, out - column vectors
% if inp has more than 3 columns (for example first column is time)
% then assume inp(:,2) is X, inp(:,3) is Y, inp(:,3) is Z
% l - [l1 l2 l3], eigenvalues where l1 is smallest
% v - column vector with corresponding eigenvectors
% v(1,:) - first vector, ...
%
% [out,l,v]=irf_minvar(inp,'td') minimum variance under tangential discontinuity
% constraint (minimizing sum of squares of normal component), Eq. 8.17 ISSI book
%
% [out,l,v]=irf_minvar(inp,'<Bn>=0') minimum variance under constraint <Bn>=0
%                                      see Eq. 8.17 ISSI book
%
% See also IRF_MINVAR_GUI, IRF_MINVAR_NEST, IRF_MINVAR_NEST_GUI
% Works with TSeries as input

if nargin==1
    flag='mvar'; % default is to do unconstrained minimum variance
end

rtrnTS = 0;
isaTSeries = isa(inp,'TSeries');
if isaTSeries
    inptemp = inp;
    inp = inptemp.data;
    rtrnTS = 1;
end

ooo=inp;
qx = size(inp(1,:)); lx=qx(2); % the number of vector components
if lx > 3, inp=inp(:,[2  3 4]);
elseif lx < 3
    disp('not enough components for x vector');
end


inp_m=irf.nanmean(inp);
Mm2=inp_m([1 2 3 1 1 2]).*inp_m([1 2 3 2 3 3]);
Mm1=irf.nanmean(inp(:,[1 2 3 1 1 2]).*inp(:,[1 2 3 2 3 3])); % all 6 elements of triagonal matrix

switch lower(flag) % define matrix M for which calculate eigenvalues and vectors
    case {'mvar','<bn>=0'}
        Mm=Mm1-Mm2;
        M=[[Mm(1) Mm(4) Mm(5)]' [Mm(4) Mm(2) Mm(6)]' [Mm(5) Mm(6) Mm(3)]'];
    case 'td'
        Mm=Mm1;
        M=[[Mm(1) Mm(4) Mm(5)]' [Mm(4) Mm(2) Mm(6)]' [Mm(5) Mm(6) Mm(3)]'];
end

%    [V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
%    full matrix V whose columns are the corresponding eigenvectors so
%    that X*V = V*D.

[V,L]=eig(M);
ll=[L(1,1) L(2,2) L(3,3)];
v=V';l=ll;
for i=1:2       % third eigenvector obtain from right hand system
    [~,ii]=max(ll);
    l(i)=ll(ii);
    v(i,:)=V(:,ii)';
    ll(ii)=-1e20;
end
v(3,:)=cross(v(1,:),v(2,:));
l(3)=max(ll);

if strcmpi(flag,'<bn>=0') % <Bn>=0 requires further calculations
    c_eval('inp_mvar_?_mean=mean(dot(inp,repmat(v(?,:),size(inp,1),1)));',1:3); % TODO get rid of c_eval for speed
    % polynom roots
    a=inp_mvar_1_mean.^2+inp_mvar_2_mean.^2+inp_mvar_3_mean.^2;
    b=-inp_mvar_1_mean.^2*(l(2)+l(3))- inp_mvar_2_mean.^2*(l(1)+l(3))- inp_mvar_3_mean.^2*(l(1)+l(2));
    c=inp_mvar_1_mean.^2*l(2)*l(3)+inp_mvar_2_mean.^2*l(1)*l(3)+inp_mvar_3_mean.^2*l(1)*l(2);
    r=roots([a b c]);
    lmin=min(r);
    c_eval('n?=inp_mvar_?_mean/(l(?)-lmin);',1:3); % normal in MVAR reference frame
    nnorm=norm([n1 n2 n3]); %#ok<NASGU>
    c_eval('n?=n?/nnorm;',1:3);
    n=[0 0 0];
    c_eval('n=n+n?*v(?,:);',1:3); % calculate normal in input reference frame
    bn=sum(inp.*repmat(n,size(inp,1),1),2);
    inp_2=inp-repmat(bn,1,3).*repmat(n,size(inp,1),1);
    [~,l,v]=irf_minvar(inp_2,'mvar');
    l(3)=lmin;
    out=(v*inp')';
elseif strcmpi(flag,'td') % 'td' method
    ln=l(3);
    bn=sum(inp.*repmat(v(3,:),size(inp,1),1),2);
    inp_2=inp-repmat(bn,1,3).*repmat(v(3,:),size(inp,1),1);
    [~,l,v]=irf_minvar(inp_2,'mvar');
    l(3)=ln;
    out=(v*inp')';
else % 'mvar' method (default)
    out=(v*inp')';
end

if nargout < 1
    disp(strcat('   max l1=',num2str(l(1),3),'v1=[',num2str(v(1,1),3),', ',num2str(v(1,2),3),', ',num2str(v(1,3),3),']'))
    disp(strcat('interm l2=',num2str(l(2),3),'v2=[',num2str(v(2,1),3),', ',num2str(v(2,2),3),', ',num2str(v(2,3),3),']'))
    disp(strcat('   min l3=',num2str(l(3),3),'v3=[',num2str(v(3,1),3),', ',num2str(v(3,2),3),', ',num2str(v(3,3),3),']'))
end

if lx>3, ooo(:,[2 3 4])=out;out=ooo; end
if rtrnTS, ooo = out; out=inptemp; out.data = ooo; end
