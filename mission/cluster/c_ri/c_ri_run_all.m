function c_ri_run_all(run_steps,st_m, et_m, min_angle, min_ampl, period, d2MP, psw)
%
% c_ri_run_all
% c_ri_run_all(run_steps)
% c_ri_run_all(run_steps,st,et)
% c_ri_run_all(run_steps,st_m, et_m, min_angle, min_ampl, period)
%
%Example
% c_ri_run_all([1 1 1 1 1 1 1], [2002 02 03 0 0 0], [2002 07 08 0 0 0], 150, 5, 3, 3,2);
% does all steps from [2002 02 03 0 0 0] to [2002 07 08 0 0 0], with
% minum angle at 150 degrees, ,minumum amplitud 5 nT, the period to class to events as one
% is 3 seconds and the maximum distance to MP is 3 Re and the solarwind preassure i
% 2 nPa.
%
%Input:
% st_m -[yyyy | mm | dd | hh | mm | ss] - can be in matrix
% et_m -[yyyy | mm | dd | hh | mm | ss] - can be in matrix
% min_angle -[150|165], must have same number of rows as
%            "st_m"
% min_ampl -[5|10], must have same number of rows as
%            "st_m"
% period -time distane, in seconds
% d2MP -distance to MP, in Re
% psw -solarwind preassure, in nPa
% run_steps - [ 0 | 1 | 1 | 1 | 0 | 0 | 0]
%       if one element is zero then the step will be jumped
% 		the steps are:
%		1) calculating the MP-crossings
%		2) Geting the data for each crossing
%		3) Converting the data into matlab format
%		4) Preprocessing the data. ! this step cannot be jumped	!
%		5) Caluculating the angles for the preprocessed data
%		6) Classing the angles as events
%   	7) Filtering the events (reducing the numbers)
% 		   and turn the events into a ascii file and a jpeg figure
%Output:
% saved to file:
% p_MP = '/share/robert/MP2/';  -the MP-crossings
% p_Bt = '/share/robert/B_temp2/'; -the satellite data in binary format
% p_Bd = '/share/robert/Bd2/';   -the satellite data in matlabformat
% p_Bp = '/share/robert/Bp2/';   -the preprocessed data
% p_A = '/share/robert/A/';      -the calculated angles
% p_E = '/share/robert/E/';      -the found events
% p_R ='/share/robert/R/';       -the events in ascii and a figure
%
% File types:
% Files with MP-crossings
% MP_20010201_00:00:00_to_20010708_00:00:00.mat
%
% Files with B-data in ascii format one for each satellite.
% "b" means burst mode and "n" means normal mode
% Ba_010613_F125419_T125434_b.01
% Ba_010613_F125419_T125434_b.02
% Ba_010613_F125419_T125434_b.03
% Ba_010613_F125419_T125434_b.04
%
% Files where the B-data from the cluster satellites are in one file
% and are matlab compatible.
% Bm_020620_F052501_T052506_b.mat
%
% Files where the B-data have been preprocessed into
% matrixes with the same timeline and the same size
% Bp_020302_F030001_T040720_b.ma_020302_F030001_T040720_b.mat
%
% Files where the angles have been calculated
% Ap_020302_F030001_T040720_b.mat
%
% Files that contains the events
% E_F20020302t030000_T20020302t040000.mat
%
% Result file, the processed events written to file
% R_F20020302t030000_T20020302t040000.txt
%
% Figure files of the processed events
% F_20020302t033512.jpg
%
%Descrition of the function:
% Finds the times when the satellites are within +-d2MP, downloads the data
% for this period. Then the angles are calculated and if the angles and the
% amplitude are larger than a certian threshold. This time is classed as
% an event. The events are reduced by classing two events as the same events
% if the timedifferance is less than period.
%
%Using:
% c_ri_auto_event_search
% c_ri_run_get_B
% c_ri_four_B_files_2_one
% c_ri_load_B_for_preprocess
% c_ri_run_calc_angles_w_pre
% c_ri_run_class_angle_as_event
% c_ri_run_events_into_pictures
%
%Work method:
% Calls the functions one at the time for each timeintervall
%
%Error:
% There are serveral calls to unix using ls and the ouput is written to file.
% (This file is not removed) and end up in the output file. This system can cause
% an error if ls produces no result. Then no file will be written and the error will
% say that it could not find the file. Or there might be a problem if the cataloge
% contains to many files.
%
% The filtering alogoithm can cause caining, meaing that several different events will
% be classed as one event if the time period is to large.
%
% This program produces lots of files, the reson this files are keept is to be able to
% only do one step.
%
% To change download source for binare data.
% The function "c_ri_run_get_B" load the binary data from "/data/cluster/DDS/" this can be changed
% in the function "create_file" and "c_ri_get_B".
%
% To work properly following shell variables should be defined
% export FGMPATH=/share/isdat_files/cal/fgm
% export SATTPATH=/data/cluster/DDS
% export ORBITPATH=/data/cluster/DDS
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03
%
% See also c_ri.m
%

%--------------------- the beginning --------------------------
if  nargin == 0
  %This is where you write the matrix with the timeintervalls
  st_m =[2002 02 02 0 0 0; 2001 02 01 0 0 0];
  et_m =[2002 07 09 0 0 0; 2001 07 08 0 0 0];
  min_angle(1:2) = [150 170];
  min_ampl(1:2) = 5;
  period(1:2) = 3;
  d2MP = 3;
  psw =2;
  run_steps(1:7) = [1 1 1 1 1 1 1];
end

if  nargin == 1
  %This is where you write the matrix with the timeintervalls
  st_m =[2002 02 02 0 0 0; 2001 02 01 0 0 0];
  et_m =[2002 07 09 0 0 0; 2001 07 08 0 0 0];
  min_angle(1:2) = [150 170];
  min_ampl(1:2) = 5;
  period(1:2) = 3;
  d2MP = 3;
  psw =2;
end


if nargin == 3
  [r , c] = size(st_m);
  min_angle(1:r) = 150;
  min_ampl(1:r) = 5;
  period(1:r) = 3;
  d2MP = 3;
  psw =2;
end

if nargin == 6
  [r , c] = size(st_m);
  d2MP = 3;
  psw =2;
end

if exist('.c_ri_parameters.mat','file')
  load .c_ri_parameters.mat;
else
  %-----sets the paths ---------------------
  % paths for saving data
  % to use the current directory, path=[pwd '/'];
  %The MP-crossing info
  p_MP = '/share/robert/MP2/';
  %temporary files. one file for each satellite
  p_Bt = '/share/robert/B_temp2/';
  % data file. B1..B4 in each file
  p_Bd = '/share/robert/Bd2/';
  % Pre_processed data. B1...B4 in each file
  p_Bp = '/share/robert/Bp2/';
  % angles are calculated
  p_A = '/share/robert/A/';
  % Events are saved
  p_E = '/share/robert/E/';
  %The result in ascii files and jpeg pictures
  p_R ='/share/robert/R/';
  %-----------------------------------------------
end

path_ok='disp([])';
while ~strcmp(path_ok,'c')
  eval(path_ok);
  disp('=========== Path information =========');
  disp(['MP-crossings          > p_MP = ''' p_MP ''';']);
  disp(['data in binary format > p_Bt = ''' p_Bt ''';']);
  disp(['data in matlab format > p_Bd = ''' p_Bd ''';']);
  disp(['preprocessed data     > p_Bp = ''' p_Bp ''';']);
  disp(['calculated angles     > p_A  = ''' p_A ''';']);
  disp(['found events          > p_E  = ''' p_E ''';']);
  disp(['events ASCII, figures > p_R  = ''' p_R ''';']);
  disp('======================================');
  disp('To change enter new value, e.g. >p_A=[pwd ''/''];  or >p_R=''/share/tmp'';');
  disp('To continue >c');
  path_ok=input('>','s');
  if exist('.c_ri_parameters.mat','file')
    try save -append .c_ri_parameters.mat p_MP p_Bt p_Bd p_Bp p_A p_E p_R;
    catch, disp('Paths changes valid only for this run!');
    end
  else
    try save .c_ri_parameters.mat p_MP p_Bt p_Bd p_Bp p_A p_E p_R;
    catch, disp('Paths changes valid only for this run!');
    end
  end

end
[i_end,c] = size(st_m);

for i = 1:i_end
  st = st_m(i,:);
  et = et_m(i,:);

  %step one
  if run_steps(1) == 1
    c_ri_auto_event_search(st,et,d2MP,psw,p_MP);
  end

  %step two
  if run_steps(2) == 1
    c_ri_run_get_B(st,et,p_MP, p_Bt);
  end

  %step three
  if run_steps(3) == 1
    c_ri_four_B_files_2_one(p_Bt,p_Bd,st,et);
  end

  %step four
  if run_steps(4) == 1
    c_ri_load_B_for_preprocess(st,et,p_Bd, p_Bp);
  end

  %step five
  if run_steps(5) == 1
    c_ri_run_calc_angles_w_pre(st,et, p_Bp, p_A);
  end

  %step  sex
  if run_steps(6) == 1
    m_ang = min_angle(i);
    m_ampl = min_ampl(i);
    c_ri_run_class_angle_as_event(st,et,m_ang,m_ampl,p_A,p_E);
  end

  %step seven
  if run_steps(7) == 1
    per = period(i);
    c_ri_run_events_into_pictures(st,et,p_MP,p_Bp,p_E,p_R, per);
  end

end
