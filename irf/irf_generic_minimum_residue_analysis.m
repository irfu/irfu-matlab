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
%  Inputs:
%       B - [M,3] magnetic field [nT]   matrix where M is number of points 
%       E - [M,3] electric field [mV/m] 
%       V - [M,3] velocity field [km/s] 
%       n - [M,1] density        [cc]
%  Output
%       L - [lmin,linterm,lmax] eigenvalues
%       V - V(:,1) eigenvector corresponding to lmin, V(:,2) - linterm, V(:,3) - lmax
%
%  See also IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS_ENGINE

%% Defaults
% Define Levi Civita for 3D
leviCivita3D = zeros(3,3,3);
leviCivita3D([8 12 22]) = 1;
leviCivita3D([6 16 20]) = -1;

%% Input check
if nargin == 0 && nargout == 0,
	help irf_generic_minimum_residue_analysis;
	return;
elseif nargin == 1 && isstruct(varargin{1}),
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
				eta = 0;
				q=zeros([size(B) 3]);
				q=reshape(repmat(B,3,1),[size(B) 3]);
				[L,V,U] = irf_generic_minimum_residue_analysis_engine('eta',eta,'q',q);
				
			case 'MFR' % Minimum Faraday Residue 
				B = args{3};
				eta = B; 
				E = args{2};
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
				[L,V,U] = irf_generic_minimum_residue_analysis_engine('eta',eta,'q',q);
		
			case 'MFRV' % Minimum Faraday Residue
				B = args{3};
				V = args{2};
				E=-irf_cross(V*1e3,B*1e-9)*1e-3; % mV/m
				[L,V,U] = irf_generic_minimum_residue_analysis('MFR',E,B);
				
			case 'MMR'
				n = args{2};
				V = args{3};
				eta = n;
				q = bsxfun(@times,n,V);
				q = reshape(q,size(q,1),1,size(q,2));
				[L,V,U] = irf_generic_minimum_residue_analysis_engine('eta',eta,'q',q);
				
			otherwise
				irf.log('critical','unrecognized input');
				return;
		end
		args(1:2)=[];
	end
end

%% Print output
if nargout == 0
	disp(['Eigenvalues: ' sprintf('%7.3f ',L)]);
	disp(vector_disp('N',V(:,1)));
	disp(vector_disp('M',V(:,2)));
	disp(vector_disp('L',V(:,3)));
	disp(vector_disp('U',U,'km/s'));
end

%% Define output
if nargout == 0,
	clear L V U;
end

%% Functions
function out = vector_disp(vectSymbol,vect,vectUnit)
if nargin==2, vectUnit ='';end
	out = sprintf(['|' vectSymbol '| = %7.3f ' vectUnit ...
		', ' vectSymbol ' = [ %7.3f %7.3f %7.3f ] ' vectUnit],...
		norm(vect),vect(1),vect(2),vect(3));
	
