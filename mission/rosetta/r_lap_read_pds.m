function data = r_lap_read_pds(fPath)
%R_LAP_READ_PDS  load RPC-LAP pds file
%
% data = R_LAP_READ_PDS(fName)
%
% fName - name of a .TAB file. LBL files are ignored.
%
% Example:
% data = R_LAP_READ_PDS('LAP_20160401_000000_617_I1H.TAB');
%
% See also: R_LAP_PLOT

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


%%

dummy = -1000;  % Put the current fill value in LAP PDS archive here!

[FILEPATH,NAME,EXT] = fileparts(fPath);

fName = [NAME EXT];
data.name = fName;

switch NAME((end-6):(end-4))
  case 'PSD'
    fSuf = 'PSD';
    frFile = NAME; frFile((end-6):(end-4)) = 'FRQ';
    freq = load([ FILEPATH filesep frFile EXT]);
    fmt = ['%s %s %f %f %d %f %f' repmat(' %f',1,length(freq))];
  otherwise
    fSuf = NAME((end-2):end);
    switch fSuf
      case 'ASW'
        fmt = '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %d';
        % fmt = '%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d';
        fields = {'obt','ne','qne','iphs','qiphs','ui','qui','te1','qte1','te2','qte2','vphk','qvphk','g'};
        % fields = {'obt','obtstop','ne','qne','iphs','qiphs','ui','qui','te1','qte1','te2','qte2','vphk','qvphk','g'};
      case 'USC'
        fmt = '%s %f %f %f %d %d';
        fields = {'obt','vsc','qvsc','src','g'};
      case 'PHO'
        fmt = '%s %f %f %f %d';
        fields = {'obt','iph','qiph','g'};
      otherwise
        flagIV = NAME(end-2);
        flagPs = NAME(end-1);
        switch flagPs
          case 'F'  % Means EFL data
            flagP = 3;
          case 'P'  % Means NPL data
            flagP = 3; % Given in source variable
          case {'1','2','3'}
            flagP = str2double(flagPs);
          otherwise, error('unrecognized Probe')
        end
        mesType = NAME(end);
        switch flagIV
          case {'I','V','B','A','E','N'}
          otherwise, error('unrecognized IV')
        end
        switch mesType
          case {'L','H','S','D'}
          otherwise, error('unrecognized Measurement Type')
        end

        switch flagIV
          case {'I','V','E'}
            switch mesType
              case {'L','H'}
                if flagP<3
                  fmt = '%s %f %f %f %f';
                  fields = {'obt','i','v','g'};
                elseif(flagP==3)
                  if(flagIV == 'E')
                    fmt = '%s %f %f %f %f';
                    fields = {'obt','e12','config','g'};
                  else
                    fmt = '%s %f %f %f %f %f';
                    if flagIV == 'V', fields = {'obt','i1','i2','v12','g'};
                    elseif flagIV == 'I', fields = {'obt','i12','v1','v2','g'};
                    else, error('should not be here')
                    end
                  end
                else, error('should not be here')
                end
              case {'D'}
                fmt = '%s %f %f %f %f %f %f';
                fields = {'obt','i','istd','v','vstd','g'};
              case {'S'}
                fileID = fopen(fPath);
                l = fgetl ( fileID );
                fclose(fileID);
                n = length(strfind(l,',')) - 4; % N of steps in a sweep
                fmt = [ '%s %s %f %f %f' repmat(' %f',1,n)];
              otherwise, error('should not be here')
            end
          case {'B'}
            fmt = '%f %f';
            fields = {'dt','di'};
          case {'N'}
            fmt = '%s %f %f %f %f %f';
            fields = {'obt','ne','qne','src','g'};
        end
    end
end
%%
fileID = fopen(fPath);
C = textscan(fileID,fmt,'Delimiter',',');
fclose(fileID);

%%
switch fSuf
  case {'ASW','NPL'}
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    %         times = cell2mat(C{2}); times(:,end+1)='Z'; t = iso2epoch(times);
    %         data.tstop = t;
    for i=1:length(fields)
      data.(fields{i}) = C{i+1};
      % data.(fields{i}) = C{i+2};
      data.(fields{i})(data.(fields{i})==dummy) = NaN;
    end
  case 'PHO'
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    for i=1:length(fields)
      data.(fields{i}) = C{i+1};
    end
    data.iph(data.iph==dummy) = NaN;
  case 'USC'
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    for i=1:length(fields)
      data.(fields{i}) = C{i+1};
    end
    data.vsc(data.vsc==dummy) = NaN;
    %       case {'USC','PHO'}
    %         times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    %         data.t = t;
    %         for i=1:length(fields)
    %             data.(fields{i}) = C{i+1};
    %         end

  case 'PSD'
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    times = cell2mat(C{2}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.tstop = t;
    data.obt = C{3};
    data.obtstop = C{4};
    data.g = C{5};
    data.ip = C{6};
    data.vp = C{7};
    data.p = cell2mat(C(8:end));
    %data.sweep(data.sweep==dummy) = NaN;
    data.f = freq;
  otherwise
    switch flagIV
      case {'I','V','E'}
        times = cell2mat(C{1});
        times(:,end+1)='Z';
        t = iso2epoch(times);
        switch mesType
          case {'L','H','D'}
            data.t = t;
            for i=1:length(fields)
              data.(fields{i}) = C{i+1};
            end
          case 'S'
            data.t = t;
            times = cell2mat(C{2});
            times(:,end+1)='Z';
            t = iso2epoch(times);
            data.tstop = t;
            data.obt = C{3};
            data.obtstop = C{4};
            data.f = C{5};
            data.sweep = cell2mat(C(6:end));
            data.sweep(data.sweep==dummy) = NaN;
        end
      case {'B'}
        for i=1:length(fields)
          data.(fields{i}) = C{i};
        end
    end
end
end