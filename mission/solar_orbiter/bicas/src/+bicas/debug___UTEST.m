%
% matlab.unittest automatic test code for bicas.debug.
%
% NOTE: Tests plotting, i.e. (1) can only check that code does not crash, (2)
% will create and close figure windows.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef debug___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_plot_QRCBM(T)
      Qrcbm = bicas.proc.QrcbMap(5);
      Qrcbm.add("QRCID_1", logical([0; 1; 1; 0; 0]))
      Qrcbm.add("QRCID_2", logical([1; 1; 1; 1; 1]))
      Qrcbm.add("QRCID_3", logical([0; 0; 0; 0; 0]))

      Fig = bicas.debug.plot_QRCBM(int64([1:5]' * 1e9), Qrcbm, "TEST FIGURE NAME");
      close(Fig)
    end



    function test_plot_VDC_EDC_FPA(T)
      DATA_AR = single([...
        0 0   0;
        1 2   2;
        2 NaN 4;
        3 4   5;
        2 NaN 2]);

      VDC_Fpa = bicas.utils.FPArray(DATA_AR,   'FILL_VALUE', single(NaN));
      EDC_Fpa = bicas.utils.FPArray(DATA_AR+3, 'FILL_VALUE', single(NaN));
      tt2000Ar = int64([1:5]' * 1e9);

      Fig = bicas.debug.plot_VDC_EDC_FPA(VDC_Fpa, EDC_Fpa, tt2000Ar, "TEST FIGURE NAME");
      close(Fig)
    end



  end    % methods(Test)



end
