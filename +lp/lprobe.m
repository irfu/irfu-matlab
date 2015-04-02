classdef lprobe
	%LP.LPROBE Class of Langmuir probes
	
	properties
		name
		type
		surface
		radiusSphere
		radiusWire
		lengthWire
	end
	properties (Dependent)
		capacitance
	end
	%     data.probe.surface='themis';
%     set(data.inp.probe.surface,'Value',find(strcmp('cluster',lp.photocurrent))+1);
%     set(data.inp.probe.length_value,'style','text','string','');
%     set(data.inp.probe.radius_value,'string','4');
%     set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells

	methods
		function capacitance = get.capacitance(Lp)
			cSphere  = irf_estimate('capacitance_sphere',Lp.radiusSphere);
			cWire    = irf_estimate('capacitance_wire',  Lp.radiusWire,Lp.lengthWire);
			capacitance = sum([ cSphere cWire]);
		end
	end
	
end

