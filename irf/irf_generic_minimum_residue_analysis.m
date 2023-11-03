function [L,V,U]=irf_generic_minimum_residue_analysis(varargin)
% IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS implements different GMRA methods: MVAB, MFR,...
%
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS is based on:
%			Sonnerup et al. 2006 JGR
%			Sonnerup et al. 2007 JGR
%
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS(...)
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('MVAB',B)   Minimum Variance Analysis
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('MFR' ,E,B) Minimum Faraday Residue
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('MFRV',V,B) Minimum Faraday Residue, E=-vxB
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('MMR' ,n,V) Conservation of mass
%
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS(..,'TD',V) Tangential Discontinuity
%  to vector V, most often V is magnetic field B.
%
%  Inputs:
%       B - [M,3] magnetic field [nT]   matrix where M is number of points
%       E - [M,3] electric field [mV/m]
%       V - [M,3] velocity field [km/s]
%       n - [M,1] density        [cc]
%  Output
%       L - [lmin,linterm,lmax] eigenvalues
%       V - V(:,1) eigenvector corresponding to lmin, V(:,2) - linterm, V(:,3) - lmax
%       U - transport velocity, boundary velocity is Un=U*V(:,1)
%
%  See also IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINE

%% Defaults
% Define Levi Civita for 3D
leviCivita3D = zeros(3,3,3);
leviCivita3D([8 12 22]) = 1;
leviCivita3D([6 16 20]) = -1;
doConstraint = false;

%% Input check
if nargin == 0 && nargout == 0
  help irf_generic_minimum_residue_analysis;
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
      case {'MVAB'} % Minimum Variance Analysis of B
        B = args{2};
        args(1:2)=[];
        eta = 0;
        q=reshape(repmat(B,3,1),[size(B) 3]);

      case 'MFR' % Minimum Faraday Residue
        B = args{3};
        E = args{2};
        args(1:3)=[];
        eta = B;
        q = zeros([size(E,1) size(E,2) size(E,2)]);
        for mm = 1:size(q,1)
          for ii = 1:size(q,2)
            for jj = 1:size(q,2)
              for kk = 1:size(q,2)
                q(mm,ii,jj) = q(mm,ii,jj) + leviCivita3D(ii,jj,kk)*E(mm,kk);
              end
            end
          end
        end

      case 'MFRV' % Minimum Faraday Residue
        B = args{3};
        V = args{2};
        args(1:3)=[];
        E=-irf_cross(V*1e3,B*1e-9)*1e-3; % mV/m
        [L,V,U] = irf_generic_minimum_residue_analysis('MFR',E,B,args{:});

      case 'MMR'
        n = args{2};
        V = args{3};
        args(1:3)=[];
        eta = n;
        q = bsxfun(@times,n,V);
        q = reshape(q,size(q,1),1,size(q,2));

      case 'TD'
        v = args{2};
        args(1:2)=[];
        vAverage = irf.nanmean(v,1);
        nConstraint = vAverage/norm(vAverage);
        doConstraint = true;

      otherwise
        irf.log('critical','unrecognized input');
        return;
    end

  end
end

%% Run the analysis
if ~exist('L','var')
  if doConstraint
    [L,V,U] = irf_generic_minimum_residue_analysis_engine('eta',eta,'q',q,'constraint',nConstraint);
  else
    [L,V,U] = irf_generic_minimum_residue_analysis_engine('eta',eta,'q',q);
  end
end
%% Calculate normal velocity

X3 = V(:,1)';
Un = dot(U,X3);

%% Print output
if nargout == 0
  disp(['Eigenvalues: ' sprintf('%7.3d ',L)]);
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

