function res = RBSP_get_data(varName,dataDir,tint,varargin)
%RBSP_get_data load RBSP data from specified data directory
%
% res = RBSP_get_data(varName,dataDir,tint)
%   Read variable varName from RBSP CDF files located in dataDir directory.
%   Variable names has to be in the following structure:
%   RBSPX_instr_var_lvl
%   X is either A or B for which probe
%   instr is the name of the instrument. Implementation started: REPT, MAGEIS, HOPE, EMFISIS. To be
%   implemented: ICE, EFW, RPS
%   var is the name of the data product to be loaded (below is a list of
%   implemented data products)
%   lvl is the level of the data product. At the moment only l3 is
%   supported.
%   dataDir is a string for the directory containing the data
%   tint is the time interval in EpochTT format
%   flags:
%   fetch_online (default 0), if set to 1 downloads to dataDir the cdf
%   files from CDAWEB if file was not found.
%   add_ancillary (default 1), if set to 0 the code doesn't load additional
%   data suchs as L_str, L, I, B_Eq, Position,...
%   Variable names that are implemented are:
%   REPT and MAGEIS:
%   FEDU, FPDU
%   HOPE:
%   Counts_E_Omni, Counts_E, Counts_He_Omni, Counts_He, Counts_O_Omni,
%   Counts_O, Counts_P_Omni, Counts_P, FEDO, FEDU, FPDO, FPDU, FODO, FODU,
%   FHEDO, FHEDU, Dens_e_200, Dens_he_30, Dens_o_30, Dens_p_30,
%   Ion_density,Tpar_e_200,Tpar_he_30,Tpar_o_30,Tpar_p_30, Tperp_e_200,
%   Tperp_he_30, Tperp_o_30, Tperp_p_30, e_Dens_flag, e_T_flag,
%   e_Tperp_flag, e_nonzero, he_Dens_flag, he_T_flag, he_Tperp_flag,
%   he_nonzero, he_o_ratio, he_p_ratio, o_Dens_flag, o_T_flag,
%   o_Tperp_flag, o_nonzero, o_p_ratio, p_Dens_flag, p_T_flag,
%   p_Tperp_flag, p_nonzero
%   EMFISIS:
%   Mag_coord_res: Bfield vector, coord is the coordinate
%   system (GEI,GEO,GSE,GSM,SM), and res is the resolution of the
%   measurement (1sec,4sec,hires)

% Written by Ahmad Lalti

Var = varargin;
fetch_online = 0;%default don't attempt to download files from the web
add_anc = 1;%default load ancillary
while ~isempty(Var)
  flag = Var{1};
  switch flag
    case 'fetch_online'
      fetch_online = Var{2};
      Var(1:2) = [];
    case 'add_ancillary'
      add_anc = Var{2};
      Var(1:2) = [];
    otherwise
      error(['undefined input flag: ' flag])

  end
end


splt = split(varName,'_');

RBSPID = lower(splt{1});
instr = lower(splt{2});
levelID = lower(splt{end});
k = 3;
while k<length(splt)
  if k == 3
    VAR = splt{k};
  else
    VAR = [ VAR '_' (splt{k})];
  end
  k = k+1;
end
clear k

switch lower(instr)
  case 'rept'
    rel = 'rel03';
    instrumentID = ['ect', filesep, instr, filesep, 'sectors', filesep, rel];
    cdfname = [RBSPID '_' rel '_ect-rept-sci-' levelID '_'];
    return_type = 'rept-PA';
  case 'mageis'
    rel = 'rel04';
    instrumentID = ['ect', filesep, instr ,filesep, 'sectors', filesep, rel];
    cdfname = [RBSPID '_' rel '_ect-mageis-' levelID '_'];
    return_type = 'mageis-PA';
  case 'hope'
    rel = 'rel04';
    %include a condition to check if the file is in moments or pitchangle
    %folder
    options_dist = {'Counts_E_Omni','Counts_E','Counts_He_Omni',...
      'Counts_He','Counts_O_Omni','Counts_O','Counts_P_Omni','Counts_P',...
      'FEDO','FEDU','FPDO','FPDU','FODO','FODU','FHEDO','FHEDU'};
    options_moms = {'Dens_e_200','Dens_he_30','Dens_o_30','Dens_p_30',...
      'Ion_density','Tpar_e_200','Tpar_he_30','Tpar_o_30','Tpar_p_30',...
      'Tperp_e_200','Tperp_he_30','Tperp_o_30','Tperp_p_30',...
      'e_Dens_flag','e_T_flag','e_Tperp_flag','e_nonzero','he_Dens_flag',...
      'he_T_flag','he_Tperp_flag','he_nonzero','he_o_ratio',...
      'he_p_ratio','o_Dens_flag','o_T_flag','o_Tperp_flag','o_nonzero',...
      'o_p_ratio','p_Dens_flag','p_T_flag','p_Tperp_flag','p_nonzero'};
    condist = strcmpi(options_dist,VAR);condist = condist(condist);
    conmoms = contains(options_moms,VAR);
    if max(conmoms) == 0
      conmoms = 0;
    else
      VAR = options_moms{conmoms};
      conmoms = conmoms(conmoms);
    end
    if condist
      instrumentID = ['ect', filesep, instr, filesep, 'pitchangle', filesep, rel];
      cdfname = [RBSPID '_' rel '_ect-hope-pa-' levelID '_'];
      return_type = 'hope-PA';
    elseif conmoms
      instrumentID = ['ect', filesep, instr, filesep, 'moments', filesep, rel];
      cdfname = [RBSPID '_' rel '_ect-hope-mom-' levelID '_'];
      return_type = 'hope-TS';
    else
      error(['Unkown variab: ' VAR])
    end

  case 'emfisis'
    vat = split(VAR,'_');
    VAR = lower(vat{1});
    resolution = lower(vat{3});
    coord = lower(vat{2});
    instrumentID = ['emfisis', filesep, 'magnetometer', filesep, resolution, filesep, coord];
    cdfname = [RBSPID(1:end-1) '-' RBSPID(end) '_magnetometer_' resolution '-' coord '_emfisis-' levelID '_'];
    return_type = 'emfisis_bfield';
  otherwise
    error(['The instrument ' instr ' is not implemented'])
end


tint = epochUnix(tint);
timeVecStart = fromepoch(tint(1));
epochFileStart = toepoch([timeVecStart(1:3) 0 0 0]);
timeVecEnd = fromepoch(tint(end));
epochFileEnd = toepoch([timeVecEnd(1:3) 0 0 0]);
res = [];
while true
  %%
  fileToRead = find_file_to_read();
  %%
  if ~isempty(fileToRead)
    irf.log('notice',['reading ' fileToRead])
    varTmp = get_ts();
    if ~isempty(varTmp)
      if isempty(res)
        res =  varTmp.tlim(EpochUnix(tint));
      else
        % join time intervals

        res = combine_ts(res,varTmp);

      end
    else
      return
    end
  end
  epochFileStart = epochFileStart + 3600*24;
  if epochFileStart>=epochFileEnd, break, end
  timeVecStart = fromepoch(epochFileStart);
end


  function res = find_file_to_read
    res = '';
    fullPath = [dataDir, filesep, RBSPID, filesep, levelID, filesep, instrumentID, filesep...
      num2str(timeVecStart(1))];
    files = dir([fullPath, filesep, cdfname epoch2yyyymmdd(epochFileStart) '*']);

    if isempty(files) && fetch_online
      files = download_from_web();
    end
    if ~isempty(files)
      maxVer = 0; fileIDx = 1;
      for iFile = 1:length(files)
        ver = str2double(files(iFile).name(end-5:end-4));
        if ver>maxVer
          maxVer = ver;
          fileIDx = iFile;
        end
      end
      res = [fullPath, filesep, files(fileIDx).name];
    end
    function res = download_from_web()
      res = [];
      cdaweb = 'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/';
      full_url = [cdaweb RBSPID '/' levelID '/' instrumentID '/' num2str(timeVecStart(1)) '/'];
      html = webread(full_url);
      filePattern = 'href="([^"]+)"';
      fileList = regexp(html, filePattern, 'tokens');
      fileList = [fileList{:}];
      fileList = convertCharsToStrings(fileList);
      substr = [cdfname epoch2yyyymmdd(epochFileStart)];
      ix = (contains(fileList,substr));

      if isempty(ix)
        disp('No file downloaded')
        return
      else
        full_url = [full_url char(fileList(ix))];
        dest = [dataDir, filesep, RBSPID, filesep, levelID, filesep, instrumentID, filesep, ...
          num2str(timeVecStart(1))];
        if ~exist(dest,'dir')
          mkdir(dest)
        end
        websave([dest, filesep, char(fileList(ix))],full_url);
        res = dir([dest, filesep, char(fileList(ix))]);
      end

    end
  end

  function res = get_ts
    [load_cdf, cdf_info] = spdfcdfread(fileToRead);
    var_names = cdf_info.Variables(:,1);
    switch return_type
      case 'rept-PA'
        %protons or electrons?
        species = lower(VAR(2));
        if strcmpi(species,'e')
          s = 'electrons';
        elseif strcmpi(species,'p')
          s = 'ions';
        end
        %time
        if strcmpi(species,'e')
          time_str = 'epoch';
        elseif strcmpi(species,'p')
          time_str = 'epoch_prot';
        end
        ixt = strcmpi(var_names,time_str);
        time = load_cdf{ixt};
        time = EpochUnix(toepoch(datevec(time)));%convert from matlab datenum to EpochUnix
        ix_tlim = epochUnix(time)>=(tint(1)) & epochUnix(time)<=(tint(2));
        %distribution
        load_var = spdfcdfread(fileToRead,'Structure',true,'variable',VAR);
        dist = load_var.Data;
        dist = permute(dist,[3,2,1]);%dimensions represent time, energy, and PA
        dist_units = load_var.Attributes.UNITS;
        formatted_units = fix_units(dist_units,'!n','!e');
        %energy
        ixe = strcmpi(var_names,[VAR '_Energy']);
        E = load_cdf{ixe};
        %pitch angle
        ixpa = strcmpi(var_names,[VAR '_Alpha']);
        Alpha = load_cdf{ixpa};
        %build the pdist object
        res = PDist(EpochTT(time(ix_tlim)),dist(ix_tlim,:,:),'pitchangle',E,Alpha);
        %additional info
        res.species = s;
        res.name = varName;
        res.units = formatted_units;
        res.userData.GlobalAttributes = cdf_info.GlobalAttributes;
        res.userData.VariableAttribute = load_var.Attributes;
        %add anciellary
        ixtt = strcmpi(var_names,[VAR '_Energy_DELTA_minus']);
        res.ancillary.delta_energy_minus = load_cdf{ixtt};
        ixtt = strcmpi(var_names,[VAR '_Energy_DELTA_plus']);
        res.ancillary.delta_energy_plus = load_cdf{ixtt};
        ixtt = strcmpi(var_names,[VAR '_Alpha_DELTA']);
        res.ancillary.Alpha_DELTA = load_cdf{ixtt};
        if add_anc
          ixtt = strcmpi(var_names,'L_star');
          res.ancillary.L_star = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'L');
          res.ancillary.L = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'I');
          res.ancillary.I = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'B_Calc');
          res.ancillary.B_Calc = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'B_Eq');
          res.ancillary.B_Eq = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'MLT');
          res.ancillary.MLT = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'MLAT');
          res.ancillary.MLAT = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'Position');
          res.ancillary.Position = load_cdf{ixtt}(ix_tlim,:);
        end
      case 'mageis-PA'
        %protons or electrons?
        species = lower(VAR(2));
        if strcmpi(species,'e')
          s = 'electrons';
        elseif strcmpi(species,'p')
          s = 'ions';
        end
        %time
        if strcmpi(species,'e')
          time_str = 'epoch';
        elseif strcmpi(species,'p')
          time_str = 'fpdu_epoch';
        end
        ixt = strcmpi(var_names,time_str);
        time = load_cdf{ixt};
        time = EpochUnix(toepoch(datevec(time)));%convert from matlab datenum to EpochUnix
        ix_tlim = epochUnix(time)>=(tint(1)) & epochUnix(time)<=(tint(2));
        %distribution
        load_var = spdfcdfread(fileToRead,'Structure',true,'variable',VAR);
        dist = load_var.Data;
        dist = permute(dist,[3,2,1]);%dimensions represent time, energy, and PA
        dist_units = load_var.Attributes.UNITS;
        formatted_units = fix_units(dist_units,'!n','!e');
        %energy
        ixe = strcmpi(var_names,[VAR '_Energy']);
        E = load_cdf{ixe};
        %pitch angle
        ixpa = strcmpi(var_names,[VAR '_Alpha']);
        Alpha = load_cdf{ixpa};
        %build the pdist object
        res = PDist(EpochTT(time(ix_tlim)),dist(ix_tlim,:,:),'pitchangle',E,Alpha);
        %additional info
        res.species = s;
        res.name = varName;
        res.units = formatted_units;
        res.userData.GlobalAttributes = cdf_info.GlobalAttributes;
        res.userData.VariableAttribute = load_var.Attributes;
        %add anciellary
        ixtt = strcmpi(var_names,[VAR '_Energy_DELTA_minus']);
        res.ancillary.delta_energy_minus = load_cdf{ixtt};
        ixtt = strcmpi(var_names,[VAR '_Energy_DELTA_plus']);
        res.ancillary.delta_energy_plus = load_cdf{ixtt};
        if add_anc
          ixtt = strcmpi(var_names,'L_star');
          res.ancillary.L_star = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'L');
          res.ancillary.L = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'I');
          res.ancillary.I = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'B_Calc');
          res.ancillary.B_Calc = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'B_Eq');
          res.ancillary.B_Eq = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'MLT');
          res.ancillary.MLT = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'MLAT');
          res.ancillary.MLAT = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'Position');
          res.ancillary.Position = load_cdf{ixtt}(ix_tlim,:);
        end
      case 'hope-PA'
        %which species
        species = split(VAR,'_');

        if length(species)>1
          if length(species) == 3
            o_or_u = 'o';
          else
            o_or_u = 'u';
          end
          species = species{2};
          if strcmpi(species,'e')
            s = 'electrons';
          elseif strcmpi(species,'p')
            s = 'protons';
          elseif strcmpi(species,'he')
            s = 'helium';
          elseif strcmpi(species,'o')
            s = 'oxygen';
          end
        else
          species = species{1};
          o_or_u = species(end);
          species(1) = [];
          species(end-1:end) = [];
          if strcmpi(species,'e')
            s = 'electrons';
          elseif strcmpi(species,'p')
            s = 'protons';
          elseif strcmpi(species,'he')
            s = 'helium';
          elseif strcmpi(species,'o')
            s = 'oxygen';
          end
        end
        %time
        if strcmpi(s,'electrons')
          time_str = 'Epoch_Ele';
        else
          time_str = 'Epoch_Ion';
        end
        ixt = strcmpi(var_names,time_str);
        time = load_cdf{ixt};
        time = EpochUnix(toepoch(datevec(time)));%convert from matlab datenum to EpochUnix
        ix_tlim = epochUnix(time)>=(tint(1)) & epochUnix(time)<=(tint(2));
        %distribution
        load_var = spdfcdfread(fileToRead,'Structure',true,'variable',VAR);
        dist = load_var.Data;
        if strcmpi(o_or_u,'u')
          dist = permute(dist,[3,2,1]);%dimensions represent time, energy, and PA
        end
        dist_units = load_var.Attributes.UNITS;
        formatted_units = fix_units(dist_units,'!N','!E');
        %energy
        if strcmpi(s,'electrons')
          ixe = strcmpi(var_names,'hope_energy_ele');
        else
          ixe = strcmpi(var_names,'hope_energy_ion');
        end
        E = load_cdf{ixe};
        %pitch angle
        if strcmpi(o_or_u,'u')
          ixpa = strcmpi(var_names,'pitch_angle');
          Alpha = load_cdf{ixpa};
        end
        %build the pdist object
        if strcmpi(o_or_u,'u')
          res = PDist(EpochTT(time(ix_tlim)),dist(ix_tlim,:,:),'pitchangle',E,Alpha);
        elseif strcmpi(o_or_u,'o')
          res = PDist(EpochTT(time(ix_tlim)),dist(ix_tlim,:,:),'omni',E);
        end
        %additional info
        res.species = s;
        res.name = varName;
        res.units = formatted_units;
        res.userData.GlobalAttributes = cdf_info.GlobalAttributes;
        res.userData.VariableAttribute = load_var.Attributes;
        %add anciellary
        if add_anc
          if strcmpi(species,'electrons')
            ixtt = strcmpi(var_names,'Energy_Ele_DELTA');
            res.ancillary.delta_energy = load_cdf{ixtt};
          else
            ixtt = strcmpi(var_names,'Energy_Ion_DELTA');
            res.ancillary.delta_energy = load_cdf{ixtt};
          end
          if strcmpi(s,'electrons')
            sp = 'Ele';
          else
            sp = 'Ion';
          end
          ixtt = strcmpi(var_names,['L_star' '_' sp]);
          res.ancillary.L_star = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['L' '_' sp]);
          res.ancillary.L = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['I' '_' sp]);
          res.ancillary.I = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['B_Calc' '_' sp]);
          res.ancillary.B_Calc = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['B_Eq' '_' sp]);
          res.ancillary.B_Eq = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['MLT' '_' sp]);
          res.ancillary.MLT = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,['Position' '_' sp]);
          res.ancillary.Position = load_cdf{ixtt}(ix_tlim,:);
        end
      case 'hope-TS'
        %which species
        species = split(VAR,'_');
        str = {'e','he','o','p','Ion'};
        con = 1;
        q = 1;
        while con
          con = ~max(strcmpi(str{q},species));
          if con == 0
            species = str{q};
          end
          q = q+1;
          if q>5
            error('species not defined')
          end
        end

        if strcmpi(species,'e')
          s = 'electrons';
        elseif strcmpi(species,'p')
          s = 'protons';
        elseif strcmpi(species,'he')
          s = 'helium';
        elseif strcmpi(species,'o')
          s = 'oxygen';
        elseif strcmpi(species,'ion')
          s = 'ion';
        end

        %time
        if strcmpi(s,'electrons')
          time_str = 'Epoch_Ele';
        else
          time_str = 'Epoch_Ion';
        end
        ixt = strcmpi(var_names,time_str);
        time = load_cdf{ixt};
        time = EpochUnix(toepoch(datevec(time)));%convert from matlab datenum to EpochUnix
        ix_tlim = epochUnix(time)>=(tint(1)) & epochUnix(time)<=(tint(2));
        %data
        load_var = spdfcdfread(fileToRead,'Structure',true,'variable',VAR);
        dat = load_var.Data;
        dist_units = load_var.Attributes.UNITS;
        formatted_units = fix_units(dist_units,'!n','!u');

        %build the TSeries
        res = TSeries(time(ix_tlim),dat(ix_tlim,:));
        %additional info
        res.name = varName;
        res.units = formatted_units;
        res.userData.GlobalAttributes = cdf_info.GlobalAttributes;
        res.userData.VariableAttribute = load_var.Attributes;
        %add anciellary
        if add_anc
          ixt_anc = strcmpi(var_names,'Epoch_Ion');
          time_anc = load_cdf{ixt_anc};
          time_anc = EpochUnix(toepoch(datevec(time_anc)));%convert from matlab datenum to EpochUnix
          ix_tlim_anc = epochUnix(time_anc)>=(tint(1)) & epochUnix(time_anc)<=(tint(2));
          ancillary.Time_anc = time_anc(ix_tlim_anc);
          ixtt = strcmpi(var_names,'L_star');
          ancillary.L_star = load_cdf{ixtt}(ix_tlim_anc);
          ixtt = strcmpi(var_names,'L');
          ancillary.L = load_cdf{ixtt}(ix_tlim_anc);
          ixtt = strcmpi(var_names,'I');
          ancillary.I = load_cdf{ixtt}(ix_tlim_anc);
          ixtt = strcmpi(var_names,'B_Eq');
          ancillary.B_Eq = load_cdf{ixtt}(ix_tlim_anc);
          ixtt = strcmpi(var_names,'MLT');
          ancillary.MLT = load_cdf{ixtt}(ix_tlim_anc);
          ixtt = strcmpi(var_names,'Position');
          ancillary.Position = load_cdf{ixtt}(ix_tlim_anc,:);
          res.userData.ancillary = ancillary;
        end
      case 'emfisis_bfield'


        %time

        ixt = strcmpi(var_names,'Epoch_centered');
        if isempty(load_cdf{ixt})
          ixt = strcmpi(var_names,'Epoch');
        end
        time = load_cdf{ixt};
        time = EpochUnix(toepoch(datevec(time)));%convert from matlab datenum to EpochUnix
        ix_tlim = epochUnix(time)>=(tint(1)) & epochUnix(time)<=(tint(2));
        %data
        if strcmpi(VAR,'mag')
          VAR = 'Mag';
        else
          error(['Variable ' VAR ' does not exist.'])
        end
        load_var = spdfcdfread(fileToRead,'Structure',true,'variable',VAR);
        dat = load_var.Data;
        B_units = load_var.Attributes.UNITS;
        formatted_units = fix_units(B_units,'!n','!u');

        %build the TSeries
        res = irf.ts_vec_xyz(time(ix_tlim),dat(ix_tlim,:));
        %additional info
        res.name = varName;
        res.units = formatted_units;
        res.userData.GlobalAttributes = cdf_info.GlobalAttributes;
        res.userData.VariableAttribute = load_var.Attributes;
        %add anciellary
        if add_anc

          ixtt = strcmpi(var_names,'range_flag');
          ancillary.range_flag = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'calState');
          ancillary.calState = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'magInvalid');
          ancillary.magInvalid = load_cdf{ixtt}(ix_tlim);
          ixtt = strcmpi(var_names,'magFill');
          ancillary.magFill = load_cdf{ixtt}(ix_tlim);
          if ~strcmpi(resolution,'hires')
            ixtt = strcmpi(var_names,'rms');
            ancillary.rms = load_cdf{ixtt}(ix_tlim);
          end
          ixtt = strcmpi(var_names,'coordinates');
          ancillary.coordinates = load_cdf{ixtt}(ix_tlim,:);
          res.userData.ancillary = ancillary;
        end
    end

  end
  function formatted_units = fix_units(cdf_units,del1,del2)

    ut = split(cdf_units,del1);

    for i = 1:length(ut)
      if isempty(ut{i}); continue;end
      utt = split(ut{i},del2);
      bse{i} = utt{1};
      if length(utt)~=1
        exp(i) = str2double(utt{2});
      else
        exp(i) = 1;
      end
      clear utt
    end

    %set numerator
    num = '';
    if isempty(exp(exp>=1))
      num = '1';
    elseif length(exp(exp>=1)) == 1
      num = [bse{exp >= 1} '^' num2str(exp(exp>=1))];
    else
      expt = (find(exp>=1));
      for i = 1:length(expt)
        if i == 1
          num = ['(' bse{expt(i)} '^' num2str(exp(expt(i))) ')'];
        else
          num = [num(1:end-1) ' ' bse{expt(i)} '^' num2str(exp(expt(i))) num(end)];
        end
      end
    end

    %set denominator
    den = '';
    if isempty(exp(exp<0))
      den = '';
      formatted_units = num;
      return
    elseif length(exp(exp<1)) == 1
      den = bse{exp<1};
      formatted_units = [num '/' den '^' num2str(exp(exp<1))];
      return
    else
      expt = (find(exp<1));
      for i = 1:length(expt)
        if i == 1
          den = ['(' bse{expt(i)} '^' num2str(exp(expt(i))) ')'];
        else
          den = [den(1:end-1) ' ' bse{expt(i)} '^' num2str(exp(expt(i))) den(end)];
        end
      end
      formatted_units = [num '/' den];
    end

  end
  function res = combine_ts(inp1,inp2)
    if strcmpi(return_type,'rept-PA') || strcmpi(return_type,'mageis-PA')

      %%combine inp1 and inp2 to 1 PDist object also combine data in
      %%ancillary

      t1 = epochUnix(inp1.time);
      t2 = epochUnix(inp2.time);
      tn = [t1; t2];
      [Tn, itn, ~] = unique(tn);

      dat1 = inp1.data;
      dat2 = inp2.data;
      datn = [dat1; dat2];
      Datn = datn(itn,:,:);

      anc1 = struct2cell(inp1.ancillary);
      anc2 = struct2cell(inp2.ancillary);
      ancn = struct;

      flds = fields(inp1.ancillary);
      for i = 1:length(flds)
        if i<4
          ancn.(flds{i}) = anc1{i};
        else
          tmp1 = anc1{i};
          tmp2 = anc2{i};
          tmpn = [tmp1; tmp2];
          ancn.(flds{i}) = tmpn(itn,:);

        end
      end

      res = PDist(EpochUnix(Tn),Datn,'pitchangle',inp1.depend{1}',inp1.depend{2}');
      res.species = inp1.species;
      res.ancillary = ancn;
      res.name = inp1.name;
      res.units = inp1.units;
      res.userData = inp1.userData;


    end

  end
end
