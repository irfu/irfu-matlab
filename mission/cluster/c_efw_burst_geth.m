function OUT=c_efw_burst_geth(filename)
%C_EFW_BURST_GETH Get the Cluster EFW internal burst header
%
% HEADERTIME=C_EFW_BURST_GETH('FILENAME') reads the EFW burst header
% from the file 'FILENAME', and transforms the relevant times into
% seconds.
%
% HEADERTIME(1) is the time for playback according to the EFW clock.
% HEADERTIME(2) is the time for playback according to the S/C clock.
% HEADERTIME(3) is the time for the data according to the EFW clock.
%
% Note! This does not always work if the computer running this
% program is not the computer generating the index files due to
% byte ordering problems.
%

% By Anders Tjulin, last update 4/9-2003


%%% Get the data in the burstfile

% b-stands for big endian, which is the case if
% burst files are generated on SUN/Sparc.
%fid = fopen(filename,'r','b');
fid = fopen(filename,'r','l'); % l is for DB/amd64
data1 = fread(fid,'uint');
fclose(fid);

fid = fopen(filename,'r','l');
data2 = fread(fid,'uchar');
fclose(fid);

%%% Get the relevant parts of the header

if data2(129)~=178 % Check that it is burst data
  OUT=[NaN;NaN;NaN];
else
  
  % The playback time according to EFW [s]
  
  playbackEFWclock=data1(1)+data1(2)/1000000000;
  
  % The playback time according to the spacecraft clock (epoch)
  
  playbackSATclock=data1(3)+data1(4)/1000000000;
  
  % The starting time for data collecting according to EFW [s]
  
  starttimeEFWclock=((((((data2(143)*256)+data2(142))*256)+data2(141))*256+data2(140))*256+data2(139))/1000+11/450;
  
  OUT = [playbackEFWclock;playbackSATclock;starttimeEFWclock];
  
end

