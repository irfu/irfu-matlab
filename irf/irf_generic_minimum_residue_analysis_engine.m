function [L,V,U]=irf_generic_minimum_residue_analysis_engine(varargin)
% IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINGE implements general GMRA solution
%
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS is based on:
%			Sonnerup et al. 2006 JGR
%			Sonnerup et al. 2007 JGR (Correction to Sonnerup et al., 2006 JGR)
%
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINE('Q',Q) calculates
%	 eigenvalues L, eigenvectors V and velocity vector U.
%
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINE('eta',eta,'q',q) calculates
%	 L, V and U given the density of conserved quantity eta and transport
%	 tensor q.
%
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINE(...,'constraint',vector)
%	 calculates eigenvectors using constraint that discontinuity normal is
%	 perpendicular to the constraint vector.
%
%  Inputs:
%       Q - [3,3] matrix
%     eta - [M,3] density of conserved quantity, M is number of points
%       q - [M,3,3] transport tensor
%  vector - [1,3] constraint vector
%
%  Output
%       L - [lmin,lmean,lmax] eigenvalues
%       V - V(:,1) eigenvector corresponding to lmin, V(:,2) - lmean, V(:,3) - lmax
%       U - transport velocity, boundary velocity is Un=U*V(:,1)

%% Defaults
doCalculateVelocity = true;
doCalculateQ        = true;
haveConstraint      = false;
%% Input check
if nargin == 0 && nargout == 0
  help irf_generic_minimum_residue_analysis_engine;
  return;
elseif nargin == 1 && isstruct(varargin{1})
  InputParameters = varargin{1};
  inputParameterFields = fieldnames(InputParameters);
  for j=1:numel(inputParameterFields)
    fieldname = inputParameterFields{j};
    eval([fieldname ' = InputParameters.' fieldname ';']);
  end
elseif nargin > 1
  args = varargin;
  while numel(args)>=2
    switch args{1}
      case 'Q'
        Q = args{2};
        doCalculateVelocity = false;
        doCalculateQ = false;
      case 'eta'
        eta = args{2};
      case 'q'
        q = args{2};
      case 'constraint'
        constraintVector = args{2};
        haveConstraint = true;
      otherwise
        irf.log('critical','unrecognized input');
        return;
    end
    args(1:2)=[];
  end
end
%% Check inputs
if numel(eta) == 1 % scalar
  eta = repmat(eta,size(q,1),1);
end
%% Calculate Q from eta and q
if doCalculateQ
  % U estimate
  % U = <deta dq>/<|deta|^2>
  deta = bsxfun(@minus,eta,irf.nanmean(eta,1)); % Eq. 10
  dq   = bsxfun(@minus,q,irf.nanmean(q,1));
  detadqAver = irf.nanmean(matrix_dot(deta,1,dq,1),1);
  deta2Aver  = irf.nanmean(dot(deta,deta,2),1);
  U = detadqAver/deta2Aver; % Eq. 12

  % Q estimate
  dqdqAver = shiftdim(irf.nanmean(matrix_dot(dq,1,dq,1),1),1);
  detadqAver2Matrix = detadqAver' *detadqAver;
  if eta==0
    Q = dqdqAver; % Eq. 19
  else
    Q = dqdqAver - detadqAver2Matrix / deta2Aver; % Eq. 15b (see correction in Sonnerup 2007)
  end
  % Correct Q by the number of dimensions so that eigenvalues coorespond to
  % the variance of dq
  Q = Q/size(dq,2);
end

%% Check for constraints
if haveConstraint
  P = eye(numel(constraintVector)) - constraintVector(:)*constraintVector(:)'; % Eq 41
  Q = P*Q*P;
end

%% Calculate eigenvalues and eigenvectors from Q

[V,lArray]=eig(Q); % L is diagonal matrix of eigenvalues and V matrix with columns eigenvectors
[L,I] = sort(diag(lArray));
V = V(:,I);

if ~doCalculateVelocity
  U = NaN;
  return;
end

%% Calculate normal velocity

Un = dot(U,V(:,1));

%% Print output
if nargout == 0
  disp(['Eigenvalues: ' sprintf('%7.5f ',L)]);
  disp(vector_disp('N',V(:,1)));
  disp(vector_disp('M',V(:,2)));
  disp(vector_disp('L',V(:,3)));
  disp(vector_disp('U',U,'km/s'));
  disp(['Un = ' num2str(Un,3),' km/s']);
end

%% Define output
if nargout == 0
  clear L V U;
end

%% Functions
function out = vector_disp(vectSymbol,vect,vectUnit)
if nargin==2, vectUnit ='';end
out = sprintf(['|' vectSymbol '| = %8.4f ' vectUnit ...
  ', ' vectSymbol ' = [ %8.4f %8.4f %8.4f ] ' vectUnit],...
  norm(vect),vect(1),vect(2),vect(3));

function out = matrix_dot(inp1,ind1,inp2,ind2)
% MATRIX_DOT summation over one index multiplication
%
% MATRIX_DOT(inp1,ind1,inp2,ind2)
% inp1,inp2 are the matrixes and summation is over dimensions (ind1+1) and
% (ind2+1). +1 because first dimension is always time.
szinp1 = size(inp1);ndimsInp1 = ndims(inp1)-1;
szinp2 = size(inp2);ndimsInp2 = ndims(inp2)-1;
szout1 = szinp1; szout1(ind1+1)=[];
szout2 = szinp2; szout2([1 ind2+1])=[];
szout = [szout1 szout2];
out = zeros(szout);
if ndimsInp1 == 1
  if ndimsInp2 == 1
    out = sum(inp1.*inp2,2);
  elseif ndimsInp2 == 2 && ind2 == 1
    for jj = 1:szinp2(3)
      for kk = 1:szinp1(2)
        out(:,jj) = out(:,jj) + inp1(:,kk).*inp2(:,kk,jj);
      end
    end
  elseif ndimsInp2 == 2 && ind2 == 2
    for jj = 1:szinp2(2)
      for kk = 1:szinp1(2)
        out(:,jj) = out(:,jj) + inp1(:,kk).*inp2(:,jj,kk);
      end
    end
  else
    error('Not yet implemented'); % not implemented
  end
elseif ndimsInp1 == 2 && ndimsInp2 == 1 && ind1 == 1
  for jj = 1:szinp1(2)
    for kk = 1:szinp2(2)
      out(:,jj) = out(:,jj) + inp1(:,kk,jj).*inp2(:,kk);
    end
  end
elseif ndimsInp1 == 2 && ndimsInp2 == 1 && ind1 == 2
  for jj = 1:szinp1(2)
    for kk = 1:szinp2(2)
      out(:,jj) = out(:,jj) + inp1(:,jj,kk).*inp2(:,kk);
    end
  end
elseif ndimsInp1 == 2 && ndimsInp2 == 2 && ind1 == 1 && ind2 == 1
  for jj = 1:szinp1(3)
    for kk = 1:szinp2(3)
      for ss = 1:szinp1(ind1+1)
        out(:,jj,kk) = out(:,jj,kk) + inp1(:,ss,jj).*inp2(:,ss,kk);
      end
    end
  end
else
  error('Not yet implemented'); % not implemented
end
