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
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


%%
[FILEPATH,NAME,EXT] = fileparts(fPath);

fName = [NAME EXT];
data.name = fName;

switch fName(21:23)
  case 'PSD'
    fSuf = 'PSD';
    frFile = NAME; frFile(21:23) = 'FRQ';
    freq = load([ FILEPATH filesep frFile EXT]);
    fmt = [ '%s %s %f %f %d %f %f' repmat(' %f',1,length(freq))];
  otherwise
    
    fSuf = fName(25:27);
    switch fSuf
      case 'ASW'
        fmt = '%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d';
        fields = {'obt','obtstop','ne','qne','iphs','qiphs','ui','qui','te1','qte1','te2','qte2','vphk','qvphk','q'};
      case 'USC'
        fmt = '%s %f %f %f %d';
        fields = {'obt','vsc','qvsc','q'};
      case 'PHO'
        fmt = '%s %f %f %f %d';
        fields = {'obt','iph','qiph','q'};
      otherwise
        flagIV = fName(25);
        flagP = str2double(fName(26));
        mesType = fName(27);
        
        switch flagIV
          case {'I','V','B','A'}
          otherwise, error('unrecoglized IV')
        end
        switch flagP
          case {1,2,3}
          otherwise, error('unrecoglized Probe')
        end
        switch mesType
          case {'L','H','S','D'}
          otherwise, error('unrecoglized Measurement Type')
        end
        
        switch flagIV
          case {'I','V'}
            switch mesType
              case {'L','H'}
                if flagP<3
                  fmt = '%s %f %f %f %f';
                  fields = {'obt','i','v','q'};
                elseif flagP==3
                  fmt = '%s %f %f %f %f %f';
                  if flagIV == 'V', fields = {'obt','i1','i2','v12','q'};
                  elseif flagIV == 'i', fields = {'obt','i12','v1','v2','q'};
                  else, error('should not be here')
                  end
                else, error('should not be here')
                end
              case {'D'}
                fmt = '%s %f %f %f %f %f %f';
                fields = {'obt','i','istd','v','vstd','q'};
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
        end
    end
end
%%
fileID = fopen(fPath);
C = textscan(fileID,fmt,'Delimiter',',');
fclose(fileID);

%%
switch fSuf
  case 'ASW'
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    times = cell2mat(C{2}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.tstop = t;
    for i=1:length(fields)
      data.(fields{i}) = C{i+2};
    end
  case {'USC','PHO'}
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t;
    for i=1:length(fields)
      data.(fields{i}) = C{i+1};
    end   
  case 'PSD'
    times = cell2mat(C{1}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.t = t; 
    times = cell2mat(C{2}); times(:,end+1)='Z'; t = iso2epoch(times);
    data.tstop = t;
    data.obt = C{3};
    data.obtstop = C{4};
    data.q = C{5};
    data.im = C{6};
    data.vm = C{7};
    data.p = cell2mat(C(8:end));
    %data.sweep(data.sweep==-1000) = NaN;
    data.f = freq;
  otherwise
    
    switch flagIV
      case {'I','V'}
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
            data.q = C{5};
            data.sweep = cell2mat(C(6:end));
            data.sweep(data.sweep==-1000) = NaN;
        end
      case {'B'}
        for i=1:length(fields)
          data.(fields{i}) = C{i};
        end
    end
end