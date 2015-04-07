classdef lprobe
	%LP.LPROBE Class of Langmuir probes
	
	properties
		name
		surface
		radiusSphere
		radiusWire     % if radius wire is two numbers than use stazer formula
		lengthWire
	end
	properties (Dependent)
		type
		Area           % structure with fields total, sunlit, totalVsSunlit, sunlitVsTotal
		capacitance
	end

	methods
		
		function type = get.type(Lp)
			if ~isempty(Lp.radiusSphere) && ~isempty(Lp.radiusWire) && ~isempty(Lp.lengthWire),
				type = 'sphere+wire';
			elseif ~isempty(Lp.radiusWire) && ~isempty(Lp.lengthWire),
				type = 'wire';
			elseif ~isempty(Lp.radiusWire),
				type = 'sphere';
			else
				type =[];
			end
		end
		
		function Area = get.Area(Lp)
			Area.total = 0;
			Area.sunlit = 0;
			if isnumeric(Lp.radiusSphere) && ~isempty(Lp.radiusSphere),
				Area.total  = Area.total  + 4*pi*Lp.radiusSphere^2;
				Area.sunlit = Area.sunlit + pi*Lp.radiusSphere^2;
			end
			if isnumeric(Lp.radiusWire) && ~isempty(Lp.radiusWire) ...
					&& isnumeric(Lp.lengthWire) && ~isempty(Lp.lengthWire),
				Area.total  = Area.total  + pi*Lp.radiusWire*Lp.lengthWire;
				Area.sunlit = Area.sunlit + Lp.radiusWire*Lp.lengthWire;
			end
			Area.totalVsSunlit = Area.total / Area.sunlit;
			Area.sunlitVsTotal = Area.sunlit / Area.total;
		end
		
		function capacitance = get.capacitance(Lp)
			cSphere  = irf_estimate('capacitance_sphere',Lp.radiusSphere);
			if isnumeric(Lp.radiusWire) && isnumeric(Lp.lengthWire),
				if Lp.lengthWire > 10*Lp.radiusWire
					cWire    = irf_estimate('capacitance_wire',  Lp.radiusWire,Lp.lengthWire);
				elseif Lp.lengthWire > Lp.radiusWire
					cWire    = irf_estimate('capacitance_cylinder',  Lp.radiusWire,Lp.lengthWire);
				else
					irf.log('error','estiamte of capacitance for cylinder requires length > radius');
					cWire = [];
				end
			end
			capacitance = sum([ cSphere cWire]);
		end
		
	end
	
end

