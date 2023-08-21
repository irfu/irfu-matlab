function c_ri_calc_angle_w_preprocess(path_input,filename,path_output)
%
% c_ri_calc_angle_w_preprocess(path_input,filename,path_output)
%
%Input:
% path_input -
% filename -
% path_output -
%
%Output:
% Saves to file: /path_output/Ap_020302_F030100_T045500_b.mat
%
%Descrition of the function:
% Calculates the angels of the loaded file
%
%Using:
% c_ri_calc_angles_and_ampl
%
%Work method:
%
%Error:
% Matrixes must have the same size. Checks the size. No file written
% if there is size not is correct.
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

p_and_f = sprintf('%s%s',path_input,filename);
disp(['loading: ' p_and_f]);

load(p_and_f)

[b1_max,c] = size(B1);
[b2_max,c] = size(B2);
[b3_max,c] = size(B3);
[b4_max,c] = size(B4);

%only if they all have the same size this will work
if b1_max == b2_max && b1_max == b3_max && b1_max == b4_max

  [angles, ampl] = c_ri_angles_and_ampl(B1,B2,B3,B4);

  p_and_f = sprintf('%sA%s',path_output,filename(2:length(filename)));
  save(p_and_f,'angles','ampl')

  if angles == 0
    disp('angles and ampl = 0 ')
  end

end

