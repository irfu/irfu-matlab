classdef plasma
	% LP.PLASMA Define plasma models
	%
	% Plasma = LP.PLASMA describes plasma model consisting of several plasma
	% components where each component is characterized by charge size and
	% sign, density, mass of particles, temperature and drift velocity.
	% Plasma properties can be a single number applicable to all plasma
	% components or a vector of the length equal to the number of plasma
	% components.
	%
	%   Plasma.name   % string
	%   Plasma.qe     % number of elementary charge
	%   Plasma.n      % [m^-3]
	%   Plasma.mp     % in proton mass, if mp=0, it gives electron mass.
	%   Plasma.T      % eV
	%   Plasma.v      % [m/s]
	%
	% Example: (e- proton plasma, 1cc, Te=1eV,Tp=10eV)
	%   Plasma    = lp.plasma
	%   Plasma.n  = 1e6;
	%   Plasma.qe = [-1 1];
	%   Plasma.T  = [1 10];
	
% TODO include also magnetic field and anisotropies

	properties
		name   % string
		qe     % number of elementary charge
		n      % [m^-3]
		mp     % in proton mass, if mp=0, it gives electron mass.
		T      % eV
		v      % [m/s]
	end
	properties (Dependent)
		m
		q
	end
	
	methods
		function m = get.m(Plasma)
			m =  Plasma.mp*1.6726e-27;
			m(m==0) = 9.1094e-31;
		end
		function q = get.q(Plasma)
			q =  Plasma.qe*1.6022e-19;
		end
	end
	
end

