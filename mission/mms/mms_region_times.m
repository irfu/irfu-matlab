function out = mms_region_times(varargin)
%MMS_REGION_TIMES Return times when MMS is in SW, MSH, or other
%
% Input:
%     Tint - Time interval over which regions are identified, TSeries or
%     string.
%
% Options:
%     SWVx - set solar wind speed threshold
%     SWB - set maximum B for solar wind threshold
%     SWn - set density maximum density threshold or range (ne)
%     SWTi - set maximum ion temperature threshold
%     SWTe - set maximum electron temperature threshold
%     MSPn - set maximum density of magnetosphere or whatever
%     plot - set to 1 to plot overview figure (otherwise no figure)
%     table - set to 1 to save times and region flags as a table
%
% Output:
%     out - TSeries with times flags corresponding to the start of a new
%     region (flag values: 0 - Magnetosheath, 1 - Solar wind, 2 - other, usually magnetosphere)
%
% Notes: function uses 4 spacecraft average
%
% Example:
%     roiList = mms_sci_roi();
%     outTS = mms_region_times(roiList(14,:),'plot',1,'table',1,'SWTe',30)

% Check the input
if nargin == 0
  help irf_region_intervals;
  return;
end

% Flags
plotfig = false;
savetable = false;

% Default Solar wind conditions
SWVx = -300;
SWn = [2 10];
SWB = 10;
SWTi = 75;
SWTe = 20;
MSPn = 2;

Tint=varargin{1};
if isa(Tint,'char')
  Tint = irf.tint(Tint);
end

args=varargin(2:end);
if numel(args)>0
  haveoptions=1;
else
  haveoptions=0;
end

while haveoptions
  l = 2;
  switch(lower(args{1}))
    case lower('SWVx')
      if numel(args)>1 && isnumeric(args{2})
        SWVx = args{2};
      end
    case lower('SWB')
      if numel(args)>1 && isnumeric(args{2})
        SWVx = args{2};
      end
    case lower('SWn')
      if numel(args)>1 && isnumeric(args{2})
        if numel(args{2})== 2
          SWn = sort(args{2});
        elseif numel(args{2})== 1
          SWn = args{2};
        else
          irf.log('warning','SWn format not recognized. Using default');
        end
      end
    case lower('SWTi')
      if numel(args)>1 && isnumeric(args{2})
        SWTi = args{2};
      end
    case lower('SWTe')
      if numel(args)>1 && isnumeric(args{2})
        SWTe = args{2};
      end
    case lower('MSPn')
      if numel(args)>1 && isnumeric(args{2})
        MSPn = args{2};
      end
    case 'plot'
      if numel(args)>1 && isnumeric(args{2})
        if args{2} > 0
          plotfig = true;
        end
      end
    case 'table'
      if numel(args)>1 && isnumeric(args{2})
        if args{2} > 0
          savetable = true;
        end
      end
    otherwise
      irf.log('warning',['Unknown flag:' args{1}])
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), haveoptions=0; end
end

ic = 1:4;
c_eval('B? = mms.get_data(''B_dmpa_srvy'',Tint,?);',ic);
if isempty(B1)
  c_eval('B? = mms.get_data(''B_dmpa_dfg_srvy_l2pre'',Tint,?);',ic)
end
if isempty(B1)
  c_eval('B? = mms.get_data(''B_dmpa_dfg_srvy_ql'',Tint,?);',ic)
end
c_eval('Vifpi? = mms.get_data(''Vi_dbcs_fpi_fast_l2'',Tint,?);',ic)
if isempty(Vifpi1)
  c_eval('Nifpi? = mms.get_data(''Ni_fpi_ql'',Tint,?);',ic)
  c_eval('Vifpi? = mms.get_data(''Vi_dbcs_fpi_ql'',Tint,?);',ic)
  c_eval('Tifpi? = mms.get_data(''Ti_dbcs_fpi_ql'',Tint,?);',ic)
  c_eval('Nefpi? = mms.get_data(''Ne_fpi_ql'',Tint,?);',ic)
  c_eval('Vefpi? = mms.get_data(''Ve_dbcs_fpi_ql'',Tint,?);',ic)
  c_eval('Tefpi? = mms.get_data(''Te_dbcs_fpi_ql'',Tint,?);',ic)
else
  c_eval('Nifpi? = mms.get_data(''Ni_fpi_fast_l2'',Tint,?);',ic)
  c_eval('Tifpi? = mms.get_data(''Ti_dbcs_fpi_fast_l2'',Tint,?);',ic)
  c_eval('Vefpi? = mms.get_data(''Ve_dbcs_fpi_fast_l2'',Tint,?);',ic)
  c_eval('Nefpi? = mms.get_data(''Ne_fpi_fast_l2'',Tint,?);',ic)
  c_eval('Tefpi? = mms.get_data(''Te_dbcs_fpi_fast_l2'',Tint,?);',ic)
end
c_eval('Tifpi? = Tifpi?.trace/3;',ic)
c_eval('Tefpi? = Tefpi?.trace/3;',ic)

% Load E
Exy1 = mms.get_data('E2d_dsl_edp_fast_ql',Tint,1);
Ex = irf.ts_scalar(Exy1.time,Exy1.data(:,1));

% Resample to spin resolution
epoch20 = fix(Tint.start.epochUnix/60)*60:5:ceil(Tint.stop.epochUnix/60)*60;
EpochS = EpochUnix(epoch20);
c_eval('B? = B?.resample(EpochS,''median'');',ic)
%c_eval('Nifpi? = Nifpi?.resample(EpochS,''median'');',ic)
%c_eval('Vifpi? = Vifpi?.resample(EpochS,''median'');',ic)
%c_eval('Tifpi? = Tifpi?.resample(EpochS,''median'');',ic)
%c_eval('Tefpi? = Tefpi?.resample(EpochS,''median'');',ic)
%c_eval('Nefpi? = Nefpi?.resample(EpochS,''median'');',ic)
%c_eval('Vefpi? = Vefpi?.resample(EpochS,''median'');',ic)
c_eval('Nifpi? = Nifpi?.resample(EpochS);',ic)
c_eval('Vifpi? = Vifpi?.resample(EpochS);',ic)
c_eval('Tifpi? = Tifpi?.resample(EpochS);',ic)
c_eval('Tefpi? = Tefpi?.resample(EpochS);',ic)
c_eval('Nefpi? = Nefpi?.resample(EpochS);',ic)
c_eval('Vefpi? = Vefpi?.resample(EpochS);',ic)

% Take averages
B = (B1+B2+B3+B4)/4;
c_eval('idxnotnan? = ~isnan(Nifpi?.data);',ic)
idxnotnan = idxnotnan1+idxnotnan2+idxnotnan3+idxnotnan4;
c_eval('Nifpi?.data(~idxnotnan?) = 0.0;',ic);
c_eval('Vifpi?.data(~idxnotnan?,:) = 0.0;',ic);
c_eval('Tifpi?.data(~idxnotnan?) = 0.0;',ic);
c_eval('Nefpi?.data(~idxnotnan?) = 0.0;',ic);
c_eval('Vefpi?.data(~idxnotnan?,:) = 0.0;',ic);
c_eval('Tefpi?.data(~idxnotnan?) = 0.0;',ic);
Nifpi = irf.ts_scalar(EpochS,(Nifpi1.data+Nifpi2.data+Nifpi3.data+Nifpi4.data)./idxnotnan);
Vifpi = irf.ts_vec_xyz(EpochS,(Vifpi1.data+Vifpi2.data+Vifpi3.data+Vifpi4.data)./[idxnotnan idxnotnan idxnotnan]);
Tifpi = irf.ts_scalar(EpochS,(Tifpi1.data+Tifpi2.data+Tifpi3.data+Tifpi4.data)./idxnotnan);
Nefpi = irf.ts_scalar(EpochS,(Nefpi1.data+Nefpi2.data+Nefpi3.data+Nefpi4.data)./idxnotnan);
Vefpi = irf.ts_vec_xyz(EpochS,(Vefpi1.data+Vefpi2.data+Vefpi3.data+Vefpi4.data)./[idxnotnan idxnotnan idxnotnan]);
Tefpi = irf.ts_scalar(EpochS,(Tefpi1.data+Tefpi2.data+Tefpi3.data+Tefpi4.data)./idxnotnan);

% Find solar wind and magnetosheath change times
if length(SWn) < 2
  idxSW1 = Nifpi.data < SWn(1);
else
  idxSW1 = Nifpi.data>SWn(1) & Nifpi.data < SWn(2);
end
idxSW2 = Vifpi.x.data < SWVx;
idxSW3 = B.abs.data < SWB;
idxSW4 = Tifpi.data < SWTi;
idxSW5 = Tefpi.data < SWTe;
idxSW = idxSW1+idxSW2+idxSW3+idxSW4+idxSW5;
idxSWint = idxSW;
idxSWint = irf.ts_scalar(EpochS,idxSWint);
idxSW = idxSW > 3.5;
idxSWf = idxSW;
idxchange = diff(idxSW);
idxcpos = find(abs(idxchange) > 0.5);
numpntthres = 20;
while min(diff(idxcpos)) < numpntthres
  [numpoints,rmpoint] = min(diff(idxcpos));
  startidx = idxcpos(rmpoint);
  idxSWf(startidx:startidx+numpoints) = idxSWf(startidx-1);
  idxcpos(rmpoint)= [];
  idxcpos(rmpoint)= [];
end
if ~isempty(idxcpos)
  if idxcpos(1) < numpntthres
    idxSWf(1:idxcpos(1)) = idxSWf(idxcpos(1)+1);
    idxcpos(1) = [];
  end
  %if ~isempty(idxcpos)
  if length(idxSWf)-idxcpos(end) < numpntthres
    idxSWf(idxcpos(end):end) = idxSWf(idxcpos(end));
    idxcpos(end) = [];
  end
  %end
end

% Find magnetosheath and magnetosphere or whatever change times
idxMS = ~idxSWf & Nifpi.data < MSPn;
idxchangeMSP = diff(idxMS);
idxcposMSP = find(abs(idxchangeMSP) > 0.5);
while min(diff(idxcposMSP)) < numpntthres
  [numpoints,rmpoint] = min(diff(idxcposMSP));
  startidx = idxcposMSP(rmpoint);
  idxMS(startidx:startidx+numpoints) = idxMS(startidx-1);
  idxcposMSP(rmpoint)= [];
  idxcposMSP(rmpoint)= [];
end
if ~isempty(idxcposMSP)
  if idxcposMSP(1) < numpntthres
    idxMS(1:idxcposMSP(1)) = idxMS(idxcposMSP(1)+1);
    idxcposMSP(1) = [];
  end
  if ~isempty(idxcposMSP)
    if length(idxMS)-idxcposMSP(end) < numpntthres
      idxMS(idxcposMSP(end):end) = idxMS(idxcposMSP(end));
      idxcposMSP(end) = [];
    end
  end
end

idxMS = irf.ts_scalar(EpochS,2*single(idxMS));
idxSWf = irf.ts_scalar(EpochS,single(idxSWf));
idxfinal = irf.ts_scalar(B.time,abs(idxSWf.data+idxMS.data-2));
idxt = abs(diff(idxfinal.data)) > 0.5;
if sum(idxt) > 0.5
  times = [(Tint(1)+-900); EpochS(idxt)];
  idxx = idxfinal.data([true; idxt]);
  out = irf.ts_scalar(times,idxx);
else
  times = (Tint(1)+-900);
  idxx = idxfinal.data(true);
  out = irf.ts_scalar(times,idxx);
end

if plotfig
  h=irf_plot(8,'newfigure');
  %h=irf_figure(540+ic,8);
  xSize=750; ySize=750;
  set(gcf,'Position',[10 10 xSize ySize]);

  xwidth = 0.86;
  ywidth = 0.115;
  set(h(1),'position',[0.10 0.99-ywidth xwidth ywidth]);
  set(h(2),'position',[0.10 0.99-2*ywidth xwidth ywidth]);
  set(h(3),'position',[0.10 0.99-3*ywidth xwidth ywidth]);
  set(h(4),'position',[0.10 0.99-4*ywidth xwidth ywidth]);
  set(h(5),'position',[0.10 0.99-5*ywidth xwidth ywidth]);
  set(h(6),'position',[0.10 0.99-6*ywidth xwidth ywidth]);
  set(h(7),'position',[0.10 0.99-7*ywidth xwidth ywidth]);
  set(h(8),'position',[0.10 0.99-8*ywidth xwidth ywidth]);

  Ball = irf.ts_scalar(B.time,[B.data B.abs.data]);
  h(1)=irf_panel('B');
  irf_plot(h(1),Ball);
  ylabel(h(1),'B (nT)','Interpreter','tex');
  irf_legend(h(1),{'B_{x}','B_{y}','B_{z}','|B|'},[0.1 0.12])
  irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)

  h(2)=irf_panel('V');
  irf_plot(h(2),Vifpi);
  ylabel(h(2),'V_{i} (km s^{-1})','Interpreter','tex');
  irf_legend(h(2),{'V_{x}','V_{y}','V_{z}'},[0.1 0.12])
  irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)

  h(3)=irf_panel('N');
  irf_plot(h(3),Nifpi);
  hold(h(3),'on')
  irf_plot(h(3),Nefpi,'b');
  hold(h(3),'off')
  ylabel(h(3),'n (cm^{-3})','Interpreter','tex');
  irf_legend(h(3),{'n_i','n_e'},[0.1 0.12])
  irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)

  h(4)=irf_panel('Ti');
  irf_plot(h(4),Tifpi);
  ylabel(h(4),'T_{i} (eV)','Interpreter','tex');
  irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)

  h(5)=irf_panel('Te');
  irf_plot(h(5),Tefpi);
  ylabel(h(5),'T_{e} (eV)','Interpreter','tex');
  irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)

  h(6)=irf_panel('Ex');
  irf_plot(h(6),Ex);
  ylabel(h(6),'E_x (mV m^{-1})','Interpreter','tex');
  irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)

  h(7)=irf_panel('SWF1');
  irf_plot(h(7),idxSWint);
  irf_zoom(h(7),'y',[0 5.5]);
  ylabel(h(7),'SW flag1','Interpreter','tex');
  irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)

  %h(7)=irf_panel('SWF2');
  %irf_plot(h(7),idxSWf);
  %irf_zoom(h(7),'y',[0 1.2]);
  %ylabel(h(7),'SW flag2','Interpreter','tex');
  %irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)

  h(8)=irf_panel('idxfinal');
  irf_plot(h(8),idxfinal);
  hold(h(8),'on')
  irf_plot(h(8),out,'ro');
  hold(h(8),'off')
  irf_zoom(h(8),'y',[0 2.2]);
  ylabel(h(8),'Flag','Interpreter','tex');
  irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12)
  irf_legend(h(8),'0 - other, 1 - SW, 2 - MSH',[0.5 0.1],'color','k','fontsize',12)

  irf_plot_axis_align(h(1:8));
  irf_zoom(h(1:8),'x',Tint);
  set(h(1:8),'fontsize',12);
end

if savetable
  fid = fopen(['mms_region_times_' irf_fname(Tint.start.epochUnix,3) '_v0.0.0.txt'],'w');
  fprintf(fid,'%s\n','% UTC time [yyyy-mm-ddThh:mm:ss.mmmuuunnnZ]     Region Flag [0 - other, 1 - SW, 2 - MSH], TAB Separated');
  for ii = 1:length(out.data)
    fprintf(fid,'%s\t%.0f\n',out.time(ii).utc(),out.data(ii));
  end
  fclose(fid);
end

end