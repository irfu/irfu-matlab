function B = get_B_from_ISDAT(spacecraft,start_time, Dt)
%
%Input:
% spacecraft-which spacecraft to get the data from: range 1-4.
% start_time-the time when to start the reading [yyyy mm dd hh mm ss].
% Dt-the duration of the reading in seconds.
%
%Output:
% B - Magneticfield data [ time|x-component | y-component | z-component ]
%
%Descrition of the function:
%
%Using:
% Mat_DbOpen()
% database unix:99 !not a function
% isGetDataLite()
% R_c_despin()
% R_c_gse2dsc()
% the differance in these files compared to c_despin() and c_gse2dsc(), is the choice of databas 'unix:99, instead of 'disco:10'
% MatDbClose()
% create_download_blocks
%
%Work method:
% Open the databas, does the three readings and then closes the databas
%
%Error: [-1 0 0 0] is returned if no data was loaded
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
B = 0;
block_size = 15*60;

disp( ' ');
disp(['get b from isdat, cl ' int2str(spacecraft) ' start time: ' R_datestring(start_time) ' duration: ' int2str(Dt) ' s']);
disp( ' ');

download_blocks = create_download_blocks(start_time, Dt, block_size);
[i_end ,c ] = size(download_blocks);
db = Mat_DbOpen('unix:99');
for i = 1:i_end
  start_time = fromepoch(download_blocks(i,1));
  Dt = download_blocks(i,2);

  fgm_t = -1;
  fgm_data = 0;

  [fgm_t, fgm_data] = isGetDataLite(db, start_time, Dt, 'Cluster', num2str(spacecraft),'fgm', 'Bprim', ' ', ' ', ' ');
  fgm_data=double(fgm_data);

  [a, b] = size(fgm_data);
  [c, d] = size(fgm_t);

  % to solve the situation if there are only zeros... error in download
  if a == 0 && b == 0 && c == 0 && d == 0
    fgm_t = -1;
  end

  if fgm_t(1) ~= -1 && a == 3 && b > 1

    %to reduce sensitivity in R_c_despin in the aspect of start and end points
    % the first point and the last points are thrown away
    if i == 1
      fgm_t = fgm_t(6:c,:);
      fgm_data = fgm_data(:,6:b);
    end
    %to reduce sensitivity in R_c_despin in the aspect of start and end points
    if i == i_end
      fgm_t = fgm_t(1:c-6,:);
      fgm_data = fgm_data(:,1:b-6);
    end

    [a, b] = size(fgm_data);
    [c, d] = size(fgm_t);

    srB = [double(fgm_t) double(real(fgm_data))'];
    dB = R_c_despin(srB, spacecraft, 'sat',db);
    B_t = R_c_gse2dsc(dB, spacecraft,-1,db);

  else
    disp(['did not load data from ' R_datestring(start_time) ' dt: ' int2str(Dt) ]);
    B_t = 0;
  end

  B = add_A2M(B, B_t);

  clear fgm_t fgm_data srB dB

end

Mat_DbClose(db);

%if the hole intervall is empty
if B == 0
  B = [-1 0 0 0];
end

