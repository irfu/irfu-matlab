function [varsb, freq] = c_efw_burst_param(filename)
%C_EFW_BURST_PARAM  EFW internal burst parameters
%
% [PARAMS, FREQ] = C_EFW_BURST_PARAM(FILENAME)
%

fid=fopen(filename,'rb'); % opens the binary file for reading
if fid==-1
  error(['Can not find burst file ' filename]);
end

fseek(fid,128,0); % skip first 128 bytes

burstreadbytes=44;

[~,name,ext] = fileparts(filename);
fns = length(name) + length(ext);

% run on intel little endian read
data(fns+1:fns+burstreadbytes,1) = fread(fid, burstreadbytes, ...
  'uint8=>uint8', 'ieee-be'); % start at 18 due to filname byte length (17)
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting frequency and number of parameters    %%%%%%%%%%%%%%%%%%%%%
%%%%%                                                  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting the information from the data and makes the corresponding number smaller.
fff = bitand(data(19),7);
ss = bitand(data(19),48);
SS = bitshift(ss,-4);

switch SS
  case 0
    switch fff
      case 0
        output=[450  8];
      case 1
        output=[900 8];
      case 2
        output=[2250 8];
      case 3
        output=[4500 8];
      case 4
        output=[9000 4];
      case 7
        output=[25000 2];
      otherwise
        output=[0 0];
    end

  case 1
    if fff==5
      output=[18000 2];
    else
      error('bad fff');
    end

  case 2
    switch fff
      case 4
        output=[9000 8];
      case 5
        output=[18000 2];
      otherwise
        output=[0 0];
    end

  case 3
    switch fff
      case 0
        output=[450 16];
      case 1
        output=[900 16];
      case 2
        output=[2250 16];
      case 3
        output=[4500 16];
      case 4
        output=[9000 8];
      case 5
        output=[18000 4];
      case 6
        output=[36000 2];
      otherwise
        output=[0 0];
    end
end

freq = output(1);

i=1;
iii=2;
ii=54;

varsb = {};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Getting information about the burst data    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
madc0={ 'V1L' 'V1M' 'V1H' 'V1U' 'V3L' 'V3M' 'V3H' 'V3U'...
  'V12M' 'V43H' 'SCX' 'SCZ' 'BAD' };
madc1={ 'V2L' 'V2M' 'V2H' 'V2U' 'V4L' 'V4M' 'V4H' 'V4U'...
  'V43M' 'V12H' 'SCY' 'BP12' 'BAD' };

while length(data)>=ii && data(ii)~=63

  adc0 = bitand(data(ii),15);
  adc1 = bitand(data(ii),240);
  adc1 = bitshift(adc1,-4);
  adc1 = bitor(bitand(adc0,8),bitand(adc1,7));

  if adc0>11 || adc0<0
    adc0=11;
  end
  varsb{i} = madc0{adc0+1}; %#ok<AGROW>

  if adc1>11 || adc1<0
    adc1=11;
  end
  varsb{iii} = madc1{adc1+1}; %#ok<AGROW>

  i=i+2;
  ii=ii+1;
  iii=iii+2;

end
