classdef lprobe
	%LP.LPROBE Class of Langmuir probes
	
	properties
		name
		type
		surface
		radiusSphere
		radiusWire     % if radius wire is two numbers than use stazer formula
		lengthWire
	end
	properties (Dependent)
		Area           % structure with fields total, sunlit
		capacitance
	end
	%     data.probe.surface='themis';
%     set(data.inp.probe.surface,'Value',find(strcmp('cluster',lp.photocurrent))+1);
%     set(data.inp.probe.length_value,'style','text','string','');
%     set(data.inp.probe.radius_value,'string','4');
%     set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells

	methods
		function Area = get.Area(Lp)
			Area.total = 0;
			Area.sunlit = 0;
			if isnumeric(Lp.radiusSphere),
				Area.total  = Area.total  + 4*pi*Lp.radiusSphere^2;
				Area.sunlit = Area.sunlit + pi*Lp.radiusSphere^2;
			end
			if isnumeric(Lp.radiusWire) && isnumeric(Lp.lengthWire),
				Area.total  = Area.total  + pi*Lp.radiusWire*Lp.lengthWire;
				Area.sunlit = Area.sunlit + Lp.radiusWire*Lp.lengthWire;
			end
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

