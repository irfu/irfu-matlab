function [newTS,colors]=irf_plot_linecolor_2TS(inputTS,colorTS,values)
% Takes the input TS, splits it into different TSeries
% based on how colorTS relates to data value. Note that inputTS and colorTS
% must share time line. (Resampled)
% Only use on 1-component TSeries.
% To plot, use: (Or see mms.Example_MMS_colorcode_plots)
% vals=linspace(min(colorTS.data),max(colorTS.data),5);
%[newTS,colors]=my.irf_plot_linecolor_2TS(Es34.x,colorTS,vals)
% charvec=[];
% charvec{1,1}='';%%
% %
% h=irf_plot(1,'newfigure');
% hold(h(1),'on');
% for ii=1:length(newTS)
% irf_plot(h(1),newTS{ii},'color',colors(ii,:),'linewidth',2)
% if ii<length(newTS)
% charvec{ii+1,1}=num2str(vals(ii));
% else
% charvec{ii+1}=[];
% end
% end
% colormap(colors)
% cbar=colorbar;
% cbar.Ticks = [0:1/size(colors,1):1];
% cbar.TickLabels = charvec;
%
% Questions? Contact konrad.steinvall@irfu.se

data = inputTS.data;
colorTSdata=colorTS.data;
for ivals = 1:length(values)+1

  if ivals==1
    interval = [0,values(ivals)];

  elseif ivals==length(values)+1
    interval = [values(end),inf];

  else
    interval = [values(ivals-1),values(ivals)];

  end

  temp=NaN(size(data));
  crit = logical(((colorTSdata>=interval(1)).*(colorTSdata<interval(2))));
  temp(crit)=data(crit);
  newTS{ivals}=irf.ts_scalar(inputTS.time,temp);
end

% Pick any colormap
colors = parula(length(newTS));


end