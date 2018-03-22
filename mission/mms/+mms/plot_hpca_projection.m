function [hax, hcb] = plot_hpca_projection(varargin)
%MMS.PLOT_HPCA_PROJECTION: plot HPCA projection on a specified plane.
% ! The purpose of this function is to make quick look of HPCA data. 
% ! If you decide to use HPCA projection/PAD results, please contac HPCA
% ! team to check the results, or ask for official figures.
%
%   [hax, hcb] = mms.plot_hpca_projection(hca, dist, tint, startaz, aze, elevation, ienergy, 'Opt1', OptVal1,...);
%
% Input:
%       dist - ion PSD or flux; [nt, npo16, ner63], looking direction;
%       tint - time interval for projection, can be smaller than 10 sec;
%       startaz - start index of azimuthal angle; [nt, 1], (0 - 15);
%       aze - azimuthal angle per energy; [nT, naz16, npo16, ner63];
%       elevation - polar angle; [npo16, 1];
%       ienergy - ion energy level; [ner63, 1];
%
% Options:
%       'vectors' - Nx2 cell array with 1x3 vector in first column and textlabel in second column
%       'xyz' - 3x3 matrix with [x;y;z]. z is normal to the plotted plane and
%               x and y are made orthogonal to z and each other if they are 
%               not already. If you want to plot different planes you have to
%               rotate this matrix -> [y;z;x] -> [z;x;y]
%       'clim' - [cmin cmax], colorbar limits in logscale
%       'vlim' - vmax in km/s, zoom in to XLim = YLim = vlim*[-1 1]
%       'elevationlim' - in degrees [0 90], limiting elevation angle  
%                        above/below projection plane to include in projection
%       'elim' - [emin, emax] energy range for projection;
%
% Example:
%    plot_hpca_projection(hca, HepPSD, tint, startaz, aze, elevation, ...
%        ienergy, 'xyz', xyz, 'elevationlim', elevlim, 'vlim', vlim, ...
%        'vlabel', vlabels);
% History:
%   1. v1 on 2016-11-20;  

    % 1. get data
    [ax,args,nargs] = axescheck(varargin{:});
    dist = args{1};     tint = args{2};         startaz = args{3};
    aze = args{4};      elevation = args{5};    ienergy = args{6};
    args = args(7: end);
    doFlipX = 0;
    doFlipY = 0;
    have_vlabels = 0;
    have_vectors = 0;
    have_clim = 0;
    have_vlim = 0;          % for plot
    have_elim = 0;          % for data
    limElevation = 20;
    correctForBinSize = 0;
    includescpot = 0;
    x = [1 0 0]; y = [0 1 0]; z = [0 0 1]; % default vectors
    naz = 16;       % azimuthal angle #
    nAZ = 32;       % projection azimuthal angle #;
    
    % 2. check species
    species = '0';
    if ~isempty(strfind(dist.name, 'hplus')),  species = 'H+'; end
    if ~isempty(strfind(dist.name, 'heplus')),  species = 'He+'; end
    if ~isempty(strfind(dist.name, 'heplusplus')),  species = 'He++'; end
    if ~isempty(strfind(dist.name, 'oplus')),  species = 'O+'; end    
    if strcmp(species, '0');     error('Can''t recognize species');  end
    distunits = dist.units;        % These used are assumed if data is not PDist format
    
    % 3. handle options 'xyz', xyz, 'elevationlim', elevlim, 'vlim', vlim, 'vlabel',vlabels
    if nargs > 6, have_options = 1; end
    while have_options
        l = 1;
        switch(lower(args{1}))   
            case 'vectors'
                l = 2;
                vectors = args{2};
                have_vectors = 1;
            case 'xyz'
                l = 2;
                coord_sys = args{2};
                x = coord_sys(1,:)/norm(coord_sys(1,:));
                y = coord_sys(2,:)/norm(coord_sys(2,:));
                z = coord_sys(3,:)/norm(coord_sys(3,:));
                z = cross(x,y); z = z/norm(z);
                y = cross(z,x); y = y/norm(y);    
                if abs(acosd(y*(coord_sys(2,:)/norm(coord_sys(2,:)))'))>1 
                    irf.log('warning',['y (perp1) changed from [' num2str(coord_sys(2,:)/norm(coord_sys(2,:)),'% .2f') '] to [' num2str(y,'% .2f') '].']);
                end
                if abs(acosd(x*(coord_sys(1,:)/norm(coord_sys(1,:)))'))>1 
                    irf.log('warning',['x (perp2) changed from [' num2str(coord_sys(1,:)/norm(coord_sys(1,:)),'% .2f') '] to [' num2str(x,'% .2f') '].']);
                end
            case 'clim'
                l = 2;
                clim = args{2};
                have_clim = 1;
            case 'vlim'
                l = 2;
                vlim = args{2};
                have_vlim = 1;
            case 'elim'
                l = 2;
                elim = args{2};
                have_elim = 1;                 
            case 'elevationlim'
                l = 2;
                limElevation = args{2}; 
            case 'vlabel'
                l = 2;
                have_vlabels = 1;
                vlabels = args{2};   
            case 'usebincorrection'      
                l = 2;
                correctForBinSize = args{2};
            case 'flipx'
                doFlipX = 1;
            case 'flipy'
                doFlipY = 1;  
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
    end
    if have_vlabels
        vlabelx = vlabels{1};
        vlabely = vlabels{2};
        vlabelz = vlabels{3};
    else
        vlabelx = ['v_{x=[' num2str(x,'% .2f') ']}'];
        vlabely = ['v_{y=[' num2str(y,'% .2f') ']}'];
        vlabelz = ['v_{z=[' num2str(z,'% .2f') ']}'];
    end    
       
    % 4. set Tint and TLIM(dist & startaz)
    nstartaz = find(startaz.data == 0); nstartaz = nstartaz(1);
    start_time_match = aze.time(1).epochUnix + (-dist.time(nstartaz).epochUnix);
    if start_time_match < 0.001
        irf.log('warning','start times of aze and dist match.');
    else
        error('Error: start times of aze and dist don''t match!');
    end
    nstopaz = find(startaz.data == 15); nstopaz = nstopaz(end);
    stop_time_match = (nstopaz - nstartaz + 1) - length(aze.time) * naz;
    if stop_time_match == 0
        irf.log('warning','stop times of aze and dist match.');
    else
        error('Error: stop times of aze and dist don''t match!');        
    end
    Tint = irf.tint(startaz.time(nstartaz), startaz.time(nstopaz));
    dist = dist.tlim(Tint);
    startaz = startaz.tlim(Tint);
    tt = dist.time;
    
    % 5. data dimension
    % 5.1. data index
    npo = length(elevation.data);
    ner = length(ienergy.data);
    nT = length(aze.time);
    nt = length(dist.time); 
    if nt ~= nT * naz, error('Error: aze and dist don''t match!'); end
    % 5.2. reshape aze to Az & Po matrix
    Az = permute(aze.data, [2, 1, 3, 4]);           % [nt, npo, ner]
    Az = reshape(Az, nt, npo, ner);    
    Po = repmat(elevation.data, 1, nt, ner);        % [nt, npo, ner]
    Po = permute(Po, [2, 1, 3]);
    xx = sind(Po) .* cosd(Az);
    yy = sind(Po) .* sind(Az);
    zz = cosd(Po);
    % 5.3. data within tint;
    dist = dist.tlim(tint);
    itstart = find(tt.epochUnix >= tint.start.epochUnix);
    itstop = find(tt.epochUnix <= tint.stop.epochUnix);
    itstart = itstart(1);       itstop = itstop(end);
    if or(abs(tt(itstart).epochUnix - dist.time(1).epochUnix) > 0.0001, ...
            abs(tt(itstop).epochUnix - dist.time(end).epochUnix) > 0.0001)
        error('Error: time within tint don''t match!'); 
    end                     % check dist.tlim(tint) is consistent with up codes.
    Az = Az(itstart: itstop, :, :);         
    Po = Po(itstart: itstop, :, :);
    xx = xx(itstart: itstop, :, :);
    yy = yy(itstart: itstop, :, :);
    zz = zz(itstart: itstop, :, :);
    % 5.4. i energy edges
    ienergy = ienergy.data;
    dE = median(diff(log10(ienergy)))/2;
    ienergyEdges = 10.^(log10(ienergy)-dE);
    ienergyEdges = [ienergyEdges; 10.^(log10(ienergy(end))+dE)];    

    % 6. compute projection Matrix
    % 6.0. define variables    
    FF = zeros(length(dist.time), nAZ, ner); % azimuthal, energy
    edgesAZ = linspace(0, 2*pi, nAZ+1);
    % 6.1. Transform into different coordinate system
    xX = reshape(xx, size(xx, 1) * size(xx, 2) * size(xx, 3) , 1);
    yY = reshape(yy, size(yy, 1) * size(yy, 2) * size(yy, 3), 1);
    zZ = reshape(zz, size(zz, 1) * size(zz, 2) * size(zz, 3), 1);
    newTmpX = [xX yY zZ]*x';            % [nt, npo16, ner63]
    newTmpY = [xX yY zZ]*y';
    newTmpZ = [xX yY zZ]*z';
    newX = reshape(newTmpX, size(xx, 1), size(xx, 2), size(xx, 3));
    newY = reshape(newTmpY, size(yy, 1), size(yy, 2), size(yy, 3));
    newZ = reshape(newTmpZ, size(zz, 1), size(zz, 2), size(zz, 3));
    elevationAngle = atan(newZ./sqrt(newX.^2+newY.^2));
    planeAz = atan2(newY,newX);         
    planeAz(planeAz < 0) = planeAz(planeAz < 0) + 2 * pi;
    geoFactorElev = cos(elevationAngle);
    geoFactorBinSize = 1;
    % if correctForBinSize, geoFactorBinSize = sin(POL);    
    % else, geoFactorBinSize = 1; end 
    
    % 6.2. collect data in newdistribution function FF
    C = dist.data .* geoFactorElev.*geoFactorBinSize;
    %C = dist.data;    
    C(abs(elevationAngle)>limElevation*pi/180) = NaN;
    if have_elim
        [~, ielim] = min(abs(ienergy - elim(1))); 
        [~, jelim] = min(abs(ienergy - elim(2)));
        irf.log('warning', ['PSD/pflux plot projection from ' ...
            num2str(ienergy(ielim), '%.2f') ' [eV] to ' ...
            num2str(ienergy(jelim), '%.2f') ' [eV].']); 
        if ielim > 1; C(:, :, 1: ielim-1) = NaN; end
        if jelim < ner; C(:, :, jelim + 1: ner) = NaN; end   
    end
    FF = zeros(nAZ, ner);
    for iAZ = 1: nAZ
        tmp = C;
        tmp(planeAz < edgesAZ(iAZ)) = NaN;
        tmp(planeAz>edgesAZ(iAZ+1)) = NaN;
        tmpp = squeeze(irf.nansum(tmp, 1));
        tmpp = squeeze(irf.nansum(tmpp, 1));
        FF(iAZ, :) = tmpp';
    end
    
    % 7. plot
    % 7.0 units & energy table
    Units = irf_units; % Use IAU and CODATA values for fundamental constants.
    switch(species)
        case 'H+'
            m = Units.mp;
        case 'He+'
            m = Units.mp * 4;
        case 'He++'
            m = Units.mp * 4;
        case 'O+'
            m = Units.mp * 16;
    end
    speedTable = real(sqrt((ienergyEdges)*Units.e*2/m)*1e-3); % km/s
    rE = speedTable;
    plX = rE*cos(edgesAZ);
    plY = rE*sin(edgesAZ);
    FF(FF == 0) = NaN; % set to white the zero points
    % 7.1 plot
    irf.log('warning','Please verify that you think the projection is done properly!');    
    hs = surf(ax,plX,plY,plY*0,log10(FF'));  
    vUnitStr= '[km/s]'; 
    if have_vlim, ax.YLim = vlim*[-1 1]; ax.XLim = ax.YLim; end
    hcb = colorbar('peer',ax);
    hcb.YLabel.String = ['log_{10} f_i [' distunits ']'];     
    % 7.2. axis label
    ax.XLabel.String = [vlabelx ' ' vUnitStr];
    ax.YLabel.String = [vlabely ' ' vUnitStr];
    ax.ZLabel.String = [vlabelz ' ' vUnitStr];
    axis(ax,'square')
    axes(ax)
    view([0 0 1])
    % 7.2 doFlipX
    if doFlipX
        cbPos = hcb.Position;
        axPos = ax.Position;
        view(ax,180,90)
        hcb.Position = cbPos;
        ax.Position = axPos;
    end
    % 7.3. title
    ax.Box = 'on';
    titleString = {[species ' : ' irf_time(tt(itstart).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' ~ ' ...
        irf_time(tt(itstop).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm')]};
    ax.Title.String = titleString;   
    shading(ax,'flat');
    if have_clim, ax.CLim = clim; end
    % Plot vectors
    hold(ax,'on');
    while have_vectors  
        if isempty(vectors), break, end  
        vecHat = vectors{1,1}/norm(vectors{1,1});
        vecHat = vecHat;
        vecTxt = vectors{1,2};
        scale = ax.YLim(2)/3;
        scaleQuiver = 1.5*scale;   
        quiver3(ax,0,0,0,(vecHat*x'),(vecHat*y'),abs(vecHat*z'),scaleQuiver,'linewidth',1,'linestyle','-','color','k')
        scaleTxt = 1.9*scale;
        text(double(scaleTxt*vecHat*x'),double(scaleTxt*vecHat*y'),1,vecTxt,'fontsize',10)
        vectors = vectors(2:end,:);  
        axis(ax,'square')
    end
    hold(ax,'off');    
    
end
%%
