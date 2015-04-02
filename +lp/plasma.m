classdef plasma
	%LP.PLASMA Define plasma models
	
	properties
		q
		n
		mp
		T
		v
	end
	properties (Dependent)
		m
	end
	
	methods
		function m = get.m(Plasma)
			m =  Plasma.mp*1.6726e-27;
			m(m==0) = 9.1094e-31;
		end
	end
	
end

