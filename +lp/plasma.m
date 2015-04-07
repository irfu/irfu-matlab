classdef plasma
	%LP.PLASMA Define plasma models
	
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

