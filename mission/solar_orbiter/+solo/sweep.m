classdef sweep
  %SOLO.SWEEP  Class representing BIAS sweep
  %   Class created from a sweep CDF file

  properties
    do
    I
    V
    time
    fName
  end

  methods
    function obj = sweep(fName)
      %UNTITLED2 Construct an instance of this class
      %   Detailed explanation goes here
      [~,obj.fName,~] = fileparts(fName);
      obj.do = dataobj(fName);
      obj.I = struct('p1',[],'p2',[],'p3',[]);
      obj.V = struct('p1',[],'p2',[],'p3',[]);

      V=get_ts(obj.do,'V');
      obj.time = irf.tint(V.time);
      I=get_ts(obj.do,'BIAS_SWEEP_CURRENT');
      AntFlag = get_ts(obj.do,'ANT_FLAG');

      for iProbe=1:3
        ii = AntFlag.data==iProbe;
        dataI = I.data(ii,iProbe);
        dataV = V.data(ii,iProbe);
        %Remove initial jumps
        dataI(1:2) = []; dataV(1:2) = [];

        obj.I.(sprintf('p%d',iProbe)) = dataI;
        obj.V.(sprintf('p%d',iProbe)) = dataV;
      end
    end

    function [iPhoto,iPhotoStd] = i_photo(obj)
      %I_PHOTO  Compute photo-saturation current
      %
      % [iPhoto,iPhotoStd] = i_photo(obj)
      iPhoto = NaN(1,3); iPhotoStd = NaN(1,3);
      for iProbe=1:3
        PHOTO_LIM = [-55 -5];
        dataI = obj.I.(sprintf('p%d',iProbe));
        dataV = obj.V.(sprintf('p%d',iProbe));

        %Determine photo saturation current
        iiPhoto = dataV>PHOTO_LIM(1) & dataV<PHOTO_LIM(2);
        if any(iiPhoto)
          iPhoto(iProbe) = median(dataI(iiPhoto));
          iPhotoStd(iProbe) = std(dataI(iiPhoto));
          if (iPhotoStd(iProbe) > 5*1000) || (iPhoto(iProbe) > 0), iPhoto(iProbe) = NaN; end
        else, iPhotoStd(iProbe) = NaN;
        end
      end
    end

    function [iMinRes,RMinRes] = min_resistance(obj)
      %MIN_RESISTANCE  Find point of min resistance (where we must bias)
      %
      % [iMinRes,RMinRes] = min_resistance(obj)

      %%%%% XXX placeholder %%%

      iMinRes = []; RMinRes = [];
    end

    function [iPhotoTS,iPhotoStdTS] = i_photo_ts(obj)
      %I_PHOTO  Get photo-saturation current as TS obj
      %
      % [iPhotoTS,iPhotoStdTS] = i_photo_ts(obj)

      [iPhoto,iPhotoStd] = i_photo(obj);

      iPhotoTS = irf.ts_scalar(obj.time.start,iPhoto);
      iPhotoTS.units = 'nA';
      iPhotoTS.siConversion = '1e-9>A';
      iPhotoTS.name = 'I saturation';

      iPhotoStdTS = irf.ts_scalar(obj.time.start,iPhotoStd);
      iPhotoStdTS.units = 'nA';
      iPhotoStdTS.siConversion = '1e-9>A';
      iPhotoStdTS.name = 'I saturation (STD)';
    end

    function hOut = plot_raw(obj,hardcopyFlag)
      %PLOT_RAW  Plot raw sweep
      %
      %  h = plot_raw(obj,hardcopyFlag)

      if nargin<2, hardcopyFlag = false; end

      Vvar=get_ts(obj.do,'V');
      Ivar=get_ts(obj.do,'BIAS_SWEEP_CURRENT');

      irf_figure(93485,2,'reset');
      h = irf_plot({Ivar*1e-3,Vvar});
      ylabel(h(1),'I [\mu A]')
      title(h(1),obj.fName,'Interpreter','none')

      if hardcopyFlag

        set(gcf,'InvertHardcopy','off','PaperPositionMode','auto')
        print(gcf,'-dpng','-painters','-r150',[obj.fName '.png'])
      end
      if nargout, hOut = h; end
    end

    function hOut = plot_iv(obj,hardcopyFlag)
      %PLOT_IV  Plot IV sweep
      %
      %  h = plot_iv(obj,hardcopyFlag)

      if nargin<2, hardcopyFlag = false; end

      [iPhoto,iPhotoStd] = i_photo(obj);

      h = irf_figure(93486,1,'reset');

      for iProbe=1:3
        plot(h,obj.V.(sprintf('p%d',iProbe)), obj.I.(sprintf('p%d',iProbe))/1000,'.-')
        legS = sprintf('I_{sat}%d = %.1f uA (std = %.1f)',iProbe,...
          iPhoto(iProbe)/1000, iPhotoStd(iProbe)/1000);
        irf_legend(legS ,[0.8, 0.4-0.1*iProbe]);
        switch iProbe
          case 1, hold(h,'on')
          case 3, hold(h,'off')
        end
      end
      ylabel(h,'I [\mu A]')
      xlabel(h,'U [V]')
      grid(h,'on')
      set(h,'XLim',[-70 70])
      legend(h,'V1','V2','V3','Location','northwest')
      title(h(1),obj.fName,'Interpreter','none')

      if hardcopyFlag
        set(gcf,'InvertHardcopy','off','PaperPositionMode','auto')
        print(gcf,'-dpng','-painters','-r150',[obj.fName '_IV.png'])
      end

      if nargout, hOut = h; end
    end

  end
end

