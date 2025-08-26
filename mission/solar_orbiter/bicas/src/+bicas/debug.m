%
% Code to help with debugging, in particular for quickly creating standard
% plots for inspecting common variables/objects. This code is such that it
% should be easy to add and remove calls to it from arbitrary places in the
% BICAS proper source code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef debug
  % TODO-DEC: How should one handle NaN/FV?
  %   NOTE: Plotting lines is bad for NaN/FV.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Method for creating a very simple plot of the content of the object in a
    % separate figure.
    %
    function Fig = plot_QRCBM(Qrcsm, tt2000Ar, figName)
      assert(isa(Qrcsm, "bicas.proc.QrcbMap"))
      bicas.utils.assert_ZV_Epoch(tt2000Ar)
      assert(isstring(figName))
      assert(numel(tt2000Ar) == Qrcsm.nRecords)

      Fig = figure('WindowState', 'maximized', "Name", figName);
      tiledlayout(numel(Qrcsm.qrcidAr'), 1, "TileSpacing", "compact", "Padding", "none");
      for i = 1:numel(Qrcsm.qrcidAr')
        qrcid = Qrcsm.qrcidAr(i);
        h = nexttile;
        plot(tt2000Ar, Qrcsm.get(qrcid), '.-');
        h.YLim = [-0.05, 1.05];
        grid on

        legend(irf.graph.escape_str(qrcid))
      end
    end



    function Fig = plot_VDC_EDC_FPA(VDC_Fpa, EDC_Fpa, tt2000Ar, figName)
      % PROPOSAL: Separate functions for VDC, EDC.
      %   NOTE: Different legends.
      FV = single(-1);

      assert(isa(VDC_Fpa, "bicas.utils.FPArray"))
      assert(isa(EDC_Fpa, "bicas.utils.FPArray"))
      bicas.utils.assert_ZV_Epoch(tt2000Ar)
      assert(isstring(figName))

      Fig = figure("WindowState", "maximized", "Name", figName)
      tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "none");

      nexttile
      plot(tt2000Ar, VDC_Fpa.array(FV), ".-")
      ylabel("VDC")
      legend(["V1", "V2", "V3"])

      nexttile
      plot(tt2000Ar, EDC_Fpa.array(FV), ".-")
      ylabel("EDC")
      legend(["V12", "V13", "V23"])
    end



  end    % methods(Static)



end
