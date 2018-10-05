function r_lap_plot(data,desc)
%R_LAP_PLOT  plot RPC-LAP data
%
% R_LAP_PLOT(data,[sweep_desc])
%
% See also: R_LAP_READ_PDS

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


fName = data.name;
switch fName(21:23)
  case 'PSD'
    mesType = 'P';
  otherwise
    flagIV = fName(25);
    flagPs = fName(26);
    mesType = fName(27);
    
    switch flagIV
      case {'I','V','B','A','U','P'}
      otherwise, error('unrecoglized IV')
    end
    switch flagPs
      case {1,2,3}
        flagP = str2double(fName(26));
      case {'S','H'}
      otherwise, error('unrecoglized Probe')
    end
    switch mesType
      case {'L','H','S','D','W','C','O'}
      otherwise, error('unrecoglized Measurement Type')
    end
end

switch mesType 
  case {'L','H','D','W','C','O','P'}
    f = fieldnames(data);
    fields = setdiff(f,{'f','name','t','tstop','obt','obtstop','qne','qiphs','qui','qte1','qte2','qvphk','q'});
    fields = [fields; {'q'}];
    nPanels = length(fields);
    
    h = irf_plot(nPanels);
    for i=1:nPanels
      hca = irf_panel(fields{i});
      l = fields{i}; coef = 1;
      switch l(1)
        case 'v', l = [l ' [V]'];
        case 'i', l = [l ' [nA]']; coef = 1e9;
        case 'n', l = [l ' [cm^{-3}]'];
        case 't', l = [l ' [eV]'];
        case 'u', l = [l ' [m/s]'];  
        case 'q', l = 'quality';
        case 'p', l = 'frequency [Hz]';  
        otherwise, error('eee')
      end
      if l(1)=='f' % PSD
        irf_spectrogram(hca,data)
      else
        irf_plot(hca,[data.t coef*double(data.(fields{i}))]);
      end
      ylabel(hca,l)
      if i==1, title(hca,data.name,'interpreter','none'), end
    end
    irf_plot_ylabels_align(h)
    if data.t(end)>data.t(1)
      irf_zoom(h,'x',[data.t(1) data.t(end)])
    end
  case 'S'
    switch flagIV
      case 'B'
        h = irf_plot(1);
        hca = irf_panel('sweep desc');
        plot(hca,data.dt,data.di)
        ylabel(hca,'v [V]')
        xlabel(hca,'time [s]')
        title(hca,data.name,'interpreter','none')
      case 'I'
        h = irf_plot(2);
        hca = irf_panel('sweep');
        data.p = data.sweep*1e9;
        data.f = desc.di;
        data.f_label = 'v [V]';
        data.p_label = 'i [nA]';
        [hca,hcb] = irf_spectrogram(hca,data);
        title(hca,data.name,'interpreter','none')
        p = get(hca,'Position');
        
        hca = irf_panel('quality');
        pp = get(hca,'Position');
        set(hca,'Position',[pp(1),pp(2),p(3),pp(4)])
        irf_plot(hca,[data.t data.q]);
        ylabel(hca,'quality')
        irf_zoom(h,'x',[data.t(1) data.t(end)])
    end 
end

end