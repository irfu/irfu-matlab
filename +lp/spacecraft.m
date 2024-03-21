classdef spacecraft
  %LP.SPACECRAFT spaceraft class

  properties
    name
    areaTotal
    areaSunlit
    areaSunlitGuard
    surface
    surfacePhotoemission = []; % if not specified, obtain from lp.photocurrent
    probeRefPotVsSatPot
    probeDistanceToSpacecraft
    nProbes
    probe
  end
  properties (Dependent)
    areaTotalVsSunlit
  end
  % 	    data.probe.type='spherical';
  %     set(data.inp.probe.type,'Value',1);
  %     data.probe.surface='themis';
  %     set(data.inp.probe.surface,'Value',find(strcmp('cluster',lp.photocurrent))+1);
  %     set(data.inp.probe.length_value,'style','text','string','');
  %     set(data.inp.probe.radius_value,'string','4');
  %     set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells
  %     set(data.inp.sc.total_area_value,'string','25.66');
  %     set(data.inp.sc.sunlit_area_value,'string','3.87');
  %     set(data.inp.sc.antenna_guard_area_value,'string','0.039');
  %     set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
  %     set(data.inp.sc.number_of_probes_value,'string','4');
  %     set(data.inp.sc.probe_distance_to_spacecraft_value,'string','44');
  %     set(data.inp.Rsun_value,'string','1');
  %     set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
  %     data.probe.total_vs_sunlit_area=4;
  %     set(data.inp.n_value,'string','1');
  %     set(data.inp.T_value,'string','100 500');

  methods

    function surfacePhotoemission = get.surfacePhotoemission(Lp)
      if any(strcmp(Lp.surface,'user defined')) || ~isempty(Lp.surfacePhotoemission)
        surfacePhotoemission = Lp.surfacePhotoemission;
      else
        surfacePhotoemission = lp.photocurrent(1,0,1,Lp.surface);
      end
    end

  end

end

