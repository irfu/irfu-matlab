function [PAdistspec, emin, emax] = get_hpca_pad(varargin)
%MMS.GET_HPCA_PAD: compute HPCA pitch angle distribution.
%
%   [PDdistspec] = get_hpca_pad(dist, startaz, aze, elevation, ienergy, dmpaB, 'Opt1', OptVal1);
% Output:
%       PAdistspec - PAD spectrum
%       emin - minimum energy
%       emax - maximum energy
% Input:
%       dist - ion PSD or flux; [nt, npo16, ner63], looking direction;
%       startaz - start index of azimuthal angle; [nt, 1], (0 - 15);
%       aze - azimuthal angle per energy; [nT, naz16, npo16, ner63];
%       elevation - polar angle; [npo16, 1];
%       ienergy - ion energy level; [ner63, 1];
%       dmpaB - B in dmpa coordinate
%
% Options:
%       'elim' - [emin, emax], energy range for PAD;
%
% Example:
%   PAdistspec = mms.get_hpca_pad(dist, startaz, aze, elevation, ienergy, dmpaB, 'elim', [500, 3000]);
%
% History:
%   1. v1 on 2016-11-18;

    % 1. get data
    [~, args, nargs] = axescheck(varargin{:});
    dist = args{1};             % [nt, npo16, ner63], looking direction of particles
    startaz = args{2};          % [nt, 1]
    aze = args{3};              % [nT, naz16, npo16, ner63]
    elevation = args{4};        % 
    ienergy = args{5};    
    dmpaB = args{6};
    elim = [ienergy.data(1), ienergy.data(end)];
    if nargs == 8
        if strcmp(args{7}, 'elim')
        elim = args{8};
        end
    end
    naz = 16;       % azimuthal angle #
       
    % 2. set tint and TLIM(dist & startaz)
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
    tint = irf.tint(startaz.time(nstartaz), startaz.time(nstopaz));
    dist = dist.tlim(tint);
    startaz = startaz.tlim(tint);
    
    % 3. compute PAD
    % 3.1. pitchangle
    anglevec = 15:15:180;      % Default pitch angles. 15 degree angle widths
    dangle = median(diff(anglevec))*ones(1,length(anglevec));
    pitchangle = anglevec-dangle/2;
    
    % 3.2. data dimension
    npo = length(elevation.data);
    ner = length(ienergy.data);
    nT = length(aze.time);
    nt = length(dist.time);
    if nt ~= nT * naz, error('Error: aze and dist don''t match!'); end
    tt = dist.time;

    % 3.3. reshape aze to Az
    Az = permute(aze.data, [2, 1, 3, 4]);           % [nt, npo, ner]
    Az = reshape(Az, nt, npo, ner);    
    Po = repmat(elevation.data, 1, nt, ner);        % [nt, npo, ner]
    Po = permute(Po, [2, 1, 3]);
    xx = sind(Po) .* cosd(Az);
    yy = sind(Po) .* sind(Az);
    zz = cosd(Po);   

    % 3.4. make B vector: [nt, npo, ner]    
    t0 = dist.time.epochUnix;     t0start = t0(1);
    t0 = t0 - t0start;  t1 = interp(t0, ner) + t0start;     t1TT = EpochUnix(t1);
    tmp1 = irf.ts_scalar(t1TT, zeros(length(t1TT), 1));
    dmpaBr = dmpaB.resample(tmp1);      dmpaBr = dmpaBr / dmpaBr.abs;
    Bx = dmpaBr.x.data;     Bx = reshape(Bx, ner, nt);      Bx = Bx';
    By = dmpaBr.y.data;     By = reshape(By, ner, nt);      By = By';
    Bz = dmpaBr.z.data;     Bz = reshape(Bz, ner, nt);      Bz = Bz';
    Bx = repmat(Bx, 1, 1, npo);     Bx = permute(Bx,[1 3 2]);       % [nt, npo, ner]    
    By = repmat(By, 1, 1, npo);     By = permute(By,[1 3 2]);    
    Bz = repmat(Bz, 1, 1, npo);     Bz = permute(Bz,[1 3 2]);    
    thetab = acosd(xx .* Bx + yy .* By + zz .* Bz);     
    
    % 3.5. select dist for PAD
    c_eval('dist? = dist.data;', 1: length(anglevec));
    dist1(thetab > anglevec(1)) = NaN;
    for jj = 2: (length(anglevec) - 1)
        c_eval('dist?(thetab < (anglevec(?)-dangle(?))) = NaN;',jj);
        c_eval('dist?(thetab > anglevec(?)) = NaN;',jj); 
    end
    c_eval('dist?(thetab < (anglevec(end)-dangle(?))) = NaN;', length(anglevec)); 
    c_eval('dist? =  squeeze(irf.nanmean(dist?,2));', 1: length(anglevec));     % [nt, npo, ner] --> [nt, ner]
    % average among energy dimension
    [~, ielim] = min(abs(ienergy.data - elim(1)));      emin = ienergy.data(ielim);
    [~, jelim] = min(abs(ienergy.data - elim(2)));      emax = ienergy.data(jelim);
    irf.log('warning', ['PSD/pflux pitch angle distribution from ' ...
        num2str(ienergy.data(ielim), '%.2f') ' [eV] to ' ...
        num2str(ienergy.data(jelim), '%.2f') ' [eV].']);
    c_eval('dist? =  squeeze(irf.nanmean(dist?(:, ielim: jelim),2));', 1: length(anglevec));     % [nt, ner] --> [nt]
    PAdist = dist1;
    for ii = 2: length(anglevec)
        c_eval('PAdist = cat(2, PAdist, dist?);',ii);      % iPADHep: [nt, ner, npitcha12]
    end
    
    % 3.6. make spectrum
    if ~isempty(strfind(dist.name, 'flux')),   distname = 'log(dF)'; end
    if ~isempty(strfind(dist.name, 'phase_space_density')),   distname = 'log(PSD)'; end
    if ~isempty(strfind(dist.name, 'hplus')),  species = 'H+'; end
    if ~isempty(strfind(dist.name, 'heplus')),  species = 'He+'; end
    if ~isempty(strfind(dist.name, 'heplusplus')),  species = 'He++'; end
    if ~isempty(strfind(dist.name, 'oplus')),  species = 'O+'; end
    PAdistspec = struct('t', tt.epochUnix);
    PAdistspec.p = double(PAdist);
    PAdistspec.p_label = {distname,dist.units};
    PAdistspec.f_label = {[species ' PAD']};
    PAdistspec.f = pitchangle;
end
%%