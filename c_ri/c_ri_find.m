function [events_tint]=c_ri_find(run_steps,st_m, et_m, min_angle, min_ampl, period, d2MP, psw)
%
% c_ri_find(run_steps,st_m, et_m, min_angle, min_ampl, period)
%
%Example
% c_ri_find([1 1 1 1 1 1 1], [2002 02 03 0 0 0], [2002 07 08 0 0 0], 150, 5, 3, 3,2);
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
%		2) obtaining angles for potential events
%		3) Classing the angles as events
%   4) Filtering the events (reducing the numbers), get event data (not ready yet)
% 		   and turn the events into a ascii file and a jpeg figure
%Output:
%  save to mMP variables
% saved to file:
% p_E = './E/';      -the found events
% p_R ='./R/';       -the events in ascii and a figure
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
% Adopted from c_ri_run_all

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
  run_steps(1:4) = [1 1 1 1];
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

if exist('.c_ri_parameters.mat'),
  load .c_ri_parameters.mat;
else
  %-----sets the paths ---------------------
  % paths for saving data
  % to use the current directory, path=[pwd '/'];
  % Events are saved
  p_E = './E/';
  %The result in ascii files and jpeg pictures
  p_R ='./R/';
  %-----------------------------------------------
end

path_ok='disp([])';
while ~strcmp(path_ok,'c'),
 eval(path_ok);
 disp('=========== Path information =========');
 disp(['found events          > p_E  = ''' p_E ''';']);
 disp(['events ASCII, figures > p_R  = ''' p_R ''';']);
 disp('======================================');
 disp('To change enter new value, e.g. >p_A=[pwd ''/''];  or >p_R=''/share/tmp'';');
 disp('To continue >c');
 path_ok=input('>','s');
 if exist('.c_ri_parameters.mat'),
  try save -append .c_ri_parameters.mat p_E p_R;
  catch disp('Paths changes valid only for this run!');
  end
 else
  try save .c_ri_parameters.mat st_m et_m min_angle min_ampl period d2MP psw run_steps;
  catch disp('Paths changes valid only for this run!');
  end
 end

end
try save -append .c_ri_parameters.mat p_E p_R;
catch disp('Input parameters not saved');
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
