classdef test_irf_plot < matlab.unittest.TestCase
  properties
    T
    t
    x
    y
    TS1
    TS2
    TS3
  end
  methods (TestClassSetup)
    function createMap(testCase)
      testCase.T = EpochTT('2002-03-04T09:30:00Z'):.2:EpochTT('2002-03-04T10:30:00Z'); % define time line as EpochTT object
      testCase.t = testCase.T.tts - testCase.T.tts(1);                     % define relative time in s from start
      testCase.x = exp(0.001*(testCase.t)).*sin(2*pi*testCase.t/180);      % define function x(t)=exp(0.001(t-to))*sin(t-to)
      testCase.TS1 = irf.ts_scalar(testCase.T,testCase.x);                 % define scalar TSeries object
      testCase.y = exp(0.001*(testCase.t)).*cos(2*pi*testCase.t/180);	     % y(t)=exp(0.001(t-to))*cos(t-to)
      testCase.TS2 = irf.ts_vec_xy(testCase.T,[testCase.x testCase.y]);
      testCase.TS3 = irf.ts_vec_xy(testCase.T,[.8*testCase.x 1.1*testCase.y]);
    end
  end

  methods (Test)

    function testPlotTwoPanels(testCase)
      import matlab.unittest.diagnostics.Diagnostic;
      import matlab.unittest.diagnostics.FigureDiagnostic;

      h = irf_plot(2);
      testCase.addTeardown(@close, h.Parent);
      hca = irf_panel('panel A');
      irf_plot(hca,testCase.TS1);
      hca = irf_panel('panel B');
      irf_plot(hca,testCase.TS2);
      % Now we log it for fun and for profit.
      testCase.log(3, ...
        Diagnostic.join('Please confirm two panels with data, both with shape of y=exp(0.001x)*sin(x).', ...
        FigureDiagnostic(gcf)));
    end

    function testPlotTwoTSeries(testCase)
      import matlab.unittest.diagnostics.Diagnostic;
      import matlab.unittest.diagnostics.FigureDiagnostic;

      irf_plot({testCase.TS1,testCase.TS2});
      % Now we log it for fun and for profit.
      testCase.addTeardown(@close, gcf);
      testCase.log(3, ...
        Diagnostic.join('Please confirm two panels of data, first panel one line and second panel lines of data.', ...
        FigureDiagnostic(gcf)));
    end

    function testPlotTwoTSeriesCompared(testCase)
      import matlab.unittest.diagnostics.Diagnostic;
      import matlab.unittest.diagnostics.FigureDiagnostic;

      irf_plot({testCase.TS1,testCase.TS2},'comp');
      testCase.addTeardown(@close, gcf);
      % Now we log it for fun and for profit.
      testCase.log(3, ...
        Diagnostic.join('Please confirm there is one panel of data with one overplotted line.', ...
        FigureDiagnostic(gcf)));
    end

    function testPlotTwoTSeriesComparedSubpanels(testCase)
      import matlab.unittest.diagnostics.Diagnostic;
      import matlab.unittest.diagnostics.FigureDiagnostic;

      h = irf_plot(2);
      hca = irf_panel('panel A');
      irf_plot(hca,{testCase.TS2.x,testCase.TS3.x},'comp');
      hca = irf_panel('panel B');
      irf_plot(hca,{testCase.TS2.y,testCase.TS3.y},'comp');
      testCase.addTeardown(@close, h.Parent);
      % Now we log it for fun and for profit.
      testCase.log(3, ...
        Diagnostic.join('Please confirm there are two panels with two lines of data in each.', ...
        FigureDiagnostic(gcf)));
    end

  end
end

