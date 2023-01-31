function r_lap_plot(data,desc)
%R_LAP_PLOT -- plot a Rosetta RPC-LAP PSA/PDS-formated data file
%
% R_LAP_PLOT(data,[sweep_desc])
%
% If data is I1S or I2S, desc can be used for giving B1S or B2S.
% If data is ASW, then passing any second argument (e.g. 0 or 'semla')
%   results in plotting all quality values; otherwise no qvals are plotted
%   for ASW data (as yielding too many panels).
% For all 1D data, passing a second argument means data points will be
%   plotted as disconnected dots (default is connected lines).
%
% See also: R_LAP_READ_PDS

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
% Additions in March 2019 by Anders E.

if(nargin < 2)
  pstyle = 'b';  % Defines plot style
else
  pstyle = 'b.';
end

fName = data.name;

if length(fName) == 28
  
  error('Geometry files not handled')
  
else % go on
  
  flagIV = fName(end-6);
  
  switch fName((end-10):(end-8))
    case 'PSD'
      mesType = 'P';
    otherwise
      flagPs = fName(end-5);
      mesType = fName(end-4);
      
      switch flagIV
        case {'I','V','B','A','U','P','E','N'}
        otherwise, error('unrecognized IV')
      end
      switch flagPs
        case {'1','2','3'}
          flagP = str2double(flagPs);
        case 'F'  % EFL data
          flagPs = '3';
          flagP = 3;
        case {'S','H','P'}
        otherwise, error('unrecognized Probe')
      end
      switch mesType
        case {'L','H','S','D','W','C','O','G'}
        otherwise, error('unrecognized Measurement Type')
      end
  end
  
  switch mesType
    case {'L','H','D','W','C','O','P'}
      f = fieldnames(data);
      % Here we list what not to plot:
      if((mesType == 'W' || flagIV =='N') && nargin > 1)
        % Then plot all qvals with ASW data
        fields = setdiff(f,{'f','name','t','tstop','obt','obtstop','qval','g'},'stable');
      else
        fields = setdiff(f,{'f','name','t','tstop','obt','obtstop','qne','qiphs','qui','qte1','qte2','qvphk','qval','g'},'stable');
      end
      fields = [fields; {'g'}];  % Adds qflag as last panel
      nPanels = length(fields);
      
      h = irf_plot(nPanels);
      for i=1:nPanels
        hca = irf_panel(fields{i});
        l = fields{i}; coef = 1; lb = l;
        switch l(1)
          case 'v', l = [l ' [V]'];
          case 'i', l = [l ' [nA]']; coef = 1e9;
          case 'n', l = [l ' [cm^{-3}]'];
          case 't', l = [l ' [eV]'];
          case 'u', l = [l ' [m/s]'];
          case 'q', l = 'qval';
          case 'g', l = 'qflag';
          case 'p', l = 'frequency [Hz]';
          case 'e', l = [l ' [mV/m]'];
          case 's'
          case 'c'
          otherwise, error('eee')
        end
        if l(1)=='f' % PSD
          irf_spectrogram(hca,data)
        else
          pdat = coef*double(data.(fields{i}));
          irf_plot(hca,[data.t pdat],pstyle);
          if data.t(end)>data.t(1)
            irf_zoom(h,'x',[data.t(1) data.t(end)]);
          end
          if(lb(1) == 'q')
            % Fixed range [0,1] for qval
            irf_zoom(hca,'y',[-0.05 1.05]);
          elseif any(~isnan(pdat))
            % For other put 5% margins on y axis
            range = max(pdat)-min(pdat);
            mid = (max(pdat)+min(pdat))/2;
            if(range > 1e-13)
              pspan = mid+1.05*range*[-1 1]/2;
            else
              pspan = [mid-1 mid+1];
            end
            irf_zoom(hca,'y',pspan);
          end
        end
        ylabel(hca,l)
        if i==1, title(hca,data.name,'interpreter','none'), end
      end
      irf_plot_ylabels_align(h);
    case 'S'
      switch flagIV
        case 'B'
          h = irf_plot(1);
          hca = irf_panel('sweep desc');
          plot(hca,data.dt,data.di,pstyle)
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
          [hca,~] = irf_spectrogram(hca,data);
          title(hca,data.name,'interpreter','none')
          p = get(hca,'Position');
          
          hca = irf_panel('quality');
          pp = get(hca,'Position');
          set(hca,'Position',[pp(1),pp(2),p(3),pp(4)])
          irf_plot(hca,[data.t data.q],pstyle);
          ylabel(hca,'quality')
          irf_zoom(h,'x',[data.t(1) data.t(end)])
      end
  end
  
end