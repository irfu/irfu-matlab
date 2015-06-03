classdef lprobe
	%LP.LPROBE Class of Langmuir probes
	
	properties
		name
		surface
		radiusSphere
		radiusWire     % if radius wire is two numbers than use stazer formula
		lengthWire
		surfacePhotoemission = []; % if not specified, obtain from lp.photocurrent
	end
	properties (Dependent)
		type
		Area           % structure with fields total, sunlit, totalVsSunlit, sunlitVsTotal,sphere,cylinder
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
			areaSphereSunlit = 0;
			areaSphereTotal  = 0;
			areaWireSunlit   = 0;
			areaWireTotal    = 0;
			if isnumeric(Lp.radiusSphere) && ~isempty(Lp.radiusSphere),
				areaSphereSunlit = pi*Lp.radiusSphere^2;
				areaSphereTotal  = 4*areaSphereSunlit;
			end
			if isnumeric(Lp.radiusWire) && ~isempty(Lp.radiusWire) ...
					&& isnumeric(Lp.lengthWire) && ~isempty(Lp.lengthWire),
				areaWireSunlit = Lp.radiusWire*Lp.lengthWire;
				areaWireTotal  = pi*areaWireSunlit;
			end
			
			Area.sphere = areaSphereTotal;
			Area.wire   = areaWireTotal;
			Area.total  = areaWireTotal + areaSphereTotal;
			Area.sunlit = areaWireSunlit + areaSphereSunlit;
			
			Area.totalVsSunlit = Area.total / Area.sunlit;
			Area.sunlitVsTotal = Area.sunlit / Area.total;
			
		end
		
		function capacitance = get.capacitance(Lp)
			cWire = 0;
			cSphere  = irf_estimate('capacitance_sphere',Lp.radiusSphere);
			if isnumeric(Lp.radiusWire) && any(Lp.radiusWire) ...
					&& isnumeric(Lp.lengthWire) && any(Lp.lengthWire),
				if Lp.lengthWire > 10*Lp.radiusWire
					cWire    = irf_estimate('capacitance_wire',  Lp.radiusWire,Lp.lengthWire);
				elseif Lp.lengthWire > Lp.radiusWire
					cWire    = irf_estimate('capacitance_cylinder',  Lp.radiusWire,Lp.lengthWire);
				else
					irf.log('critical','estiamte of capacitance for cylinder requires length > radius');
					cWire = [];
				end
			end
			capacitance = sum([ cSphere cWire]);
		end

		function surfacePhotoemission = get.surfacePhotoemission(Lp)
			if any(strcmp(Lp.surface,'user defined')) || ~isempty(Lp.surfacePhotoemission)
				surfacePhotoemission = Lp.surfacePhotoemission;
			else
				surfacePhotoemission = lp.photocurrent(1,0,1,Lp.surface);
			end
		end

	end
	
end

