function irf_pl_ebsp_new(ebsp,params)



if nargin==1 || isempty(params)
  params = {{'bb_xxyyzzss',4},...
    {'ee_ss'},...
    {'dop'},{'planarity'},{'ellipticity'},{'k_tp'},...
    {'pf_rtp',[1 2]}};
end

fieldsEBSP = fields(ebsp);
IGNORE_FIELDS = {'t','f','fac'};
fieldsPlottable = setxor(fieldsEBSP,IGNORE_FIELDS);

plotFields = ''; plotComps = ''; nPanels = 0;
for idx = 1:length(params)
  p = params{idx};
  if ~ischar(p{1})
    error('invalid FIELD_NAME for parameter #%d, expecting string',idx)
  end
  
  if isempty(intersect(p(1),fieldsPlottable)), continue, end
  plotFields = [plotFields p(1)]; %#ok<AGROW>
  if length(p)>1
    if ~isnumeric(p{2})
      error('invalid COMP for parameter %s, expecting numerical array',p{1})
    end
    plotComps = [plotComps p(2)]; %#ok<AGROW>
    nComps = length(p{2});
  else
    nComps = size(ebsp.(p{1}),3);
    plotComps = [plotComps {1:nComps}]; %#ok<AGROW>
  end
  nPanels = nPanels + nComps;
end
    
h = irf_plot(nPanels);

sr = struct('t',ebsp.t,'f',ebsp.f);


for idxField = 1:length(plotFields)
  for idxComp = 1:length(plotComps{idxField})
    field = plotFields{idxField}; comp = plotComps{idxField}(idxComp);
    [paramStr,compStr] = GetCompStrings();
    if isempty(compStr), panelStr = paramStr;
    else panelStr = [paramStr '_' compStr]; 
    end
    hca = irf_panel(panelStr);
    sr.p = UpdateUnits(ebsp.(field)(:,:,comp));
    [sr.plot_type,sr.p_label] = GetPlotTypeLabel();
    irf_spectrogram(hca,sr)
    set(hca,'YScale','log')
    SetCaxis()
  end
end

irf_zoom(h,'x',ebsp.t([1 end])')

  function [f,c] = GetCompStrings
    a = tokenize(field,'_');
    f = a{1}; c = '';
    if length(a) ==1, return, end 
    r = a{2};
    if length(unique(r)) == length(r), c = r(comp);
    else c = r((comp-1)*2+(1:2));
    end
    if strcmpi(c,'ss'), c = 'sum'; end
  end
  function [t,s] = GetPlotTypeLabel  
    t = 'lin';
    switch compStr
      case {'r','x','y','z','xx','yy','zz','sum'}
        s = ['log(' paramStr '_{' upper(compStr) '}) \newline ' GetUnits()];
        t = 'log';
      case 't'
        s = ['\Theta_{' paramStr '} [deg]'];
      case 'p'
        s = ['\Phi_{' paramStr '} [deg]'];
      otherwise
        s = paramStr;
    end
    function s = GetUnits
      switch paramStr
        case 'bb'
          s = '[nT^2/Hz';
        case 'ee'
          s = '[(mV/m)^2/Hz]';
        otherwise
          s = '[\mu W/m^2 Hz]^{1/2}';
      end
    end
  end
  function SetCaxis
    switch paramStr
      case {'dop','dop2d','planarity'}
        caxis(hca,[0 1])
      case 'ellipticity'
        caxis(hca,[-1 1])
      otherwise
        % do nothing
    end
    switch compStr
      case 't'
        caxis(hca,[0 180])
      case 'p'
        caxis(hca,[-180 180])
      otherwise
        % do nothing
    end
  end

  function a=UpdateUnits(a)
    switch compStr
      case {'t','p'}
        a = a*180/pi; % to degrees
    end
  end
end
    


