function [L,V,U]=irf_generic_minimum_residue_analysis(varargin)
% IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS implements different GMRA methods: MVAB, MFR,...
% 
%  IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS is based on Sonnerup et al. 2006 JGR
%
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('Q',Q) calculates
%	 eigenvalues L, eigenvectors V and velocity vector U.
% 
%	 [L,V,U]=IRF_GENERIC_MINIMUM_RESIDUE_ANALYSIS('eta',eta,'q','q') calculates
%	 L, V and U given the density of conserved quantity eta and transport
%	 tensor q. 
% 
%  Inputs:
%       Q - [3,3] matrix
%     eta - [M,3] density of conserved quantity, M is number of points
%       q - [M,3,3] transport tensor
%  Output 
%       L - [lmin,lmean,lmax] eigenvalues
%       V - V(:,1) eigenvector corresponding to lmin, V(:,2) - lmean, V(:,3) - lmax

%% Defaults
doCalculateVelocity = true;

%% Input check
if nargin == 0 && nargout == 0,
	help irf_generic_minimum_variance_analysis;
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
			case {'Q'}
				Q = args{2};
				doCalculateVelocity = false;
			case {'eta'}
				eta = args{2};
			case {'q'}
				q = args{2};
			otherwise
				irf.log('critical','unrecognized input');
				return;
		end
		args(1:2)=[];
	end
end
%% Calculate eigenvalues and eigenvectors from Q


[V,lArray]=eig(Q); % L is diagonal matrix of eigenvalues and V matrix with columns eigenvectors
L=[lArray(1,1) lArray(2,2) lArray(3,3)];

if ~doCalculateVelocity
	U = NaN;
	return;
end

%% Calculate velocity
