function out = psp(varargin)
%PSP  common information on Parker Solar Probe
%
% PSP - return basic information on Solar Orbiter
%
% PSP(flag) - return specific information defined by flag or keyword
% including "flag" in its text. Keywords and different useful info 
% are at the end of the file.
%
%  Example:
%   psp enc        % display encounter information
%


if nargin==0 % return basic information on Parker Solar Probe
  psp ?
  return
end

flag=varargin{1};
Units=irf_units;

switch flag
    case '?'         % Show all possible flag options
      content = fileread( which('psp') ) ;
      keyWords = regexp(content,'(?<key>##\w+)','names');
      keywords = unique(squeeze(struct2cell(keyWords)));
      for i=1:numel(keywords)
        disp(['<a href="matlab: psp ' keywords{i} '">' keywords{i} '</a> ' ]);
      end
        fid=fopen(which('psp'));
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            b=regexp(tline,'case ''(?<flag>\w*)''\s*[%](?<comment>.*)','names');
            if numel(b)==1
                disp(['<a href="matlab: psp ' b.flag '">' b.flag '</a> : ' b.comment]);
            end
        end
        fclose(fid);
	otherwise
        psp_txt(lower(flag));
end
end

function psp_txt(txt)
  fid=fopen(which('psp'));
  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    b=regexp(tline,'^%%%%.*','match'); % find where start txt paragraphs
    if numel(b)==1, break; end
  end

  while 1 % look for keywords
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    b=regexpi(tline,['^%%\s.*' txt '.*'],'match'); % find where start txt paragraphs
    if numel(b)==1
      disp(b{1}(2:end));
      while 1
        tline = fgetl(fid);
        if ~ischar(tline) || isempty(tline), break, end
        disp(tline(2:end));
      end
    end
  end
end
%%%%
%% ##basic info
% Parker Solar Probe NASA Living With a Star mission
% 

%% ##Encounters ##Venus
% List of encounters with ENLIL plots <a href="https://sppgway.jhuapl.edu/encounters#enc1">www</a>
%   V1 2018-10-03/08:44
% Enc1 2018-11-06/03:27  35.6 Rs  0.166 AU
% Enc2 2019-04-04/22:39  35.6 Rs  0.166 AU
% Enc3 2019-09-01/17:50  35.6 Rs  0.166 AU (no SPC data)
%   V2 2019-12-26/18:14
% Enc4 2020-01-29/09:37  27.8 Rs  0.129 AU
% Enc5 2020-06-07/08:23  27.8 Rs  0.129 AU
%   V3 2020-07-11
% Enc6 2020-09-27/09:16  20.3 Rs  0.094 AU
% Enc7 2021-01-17/17:40  20.3 Rs  0.094 AU
%
% 

%% ##Payload
% FIELDS paper: <a href="https://doi.org/10.1002/2016JA022344">Malaspina 2016</a>





