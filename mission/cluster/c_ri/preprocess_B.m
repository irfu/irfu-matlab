function [B1_c,B2_c,B3_c,B4_c] = preprocess_B(B1,B2,B3,B4,mode)
%
%Input:
% B1, B2, B3, B4 -matrices with B-data
% mode -burst or normal (b/n)
%
%Output:
% B1,B2,B3,B4 -preprocessed data. All the B-matrices have the same size,
%             start and end time.
%
%Descrition of the function:
%
%Using:
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
[B1_max,col] = size(B1);
[B2_max,col] = size(B2);
[B3_max,col] = size(B3);
[B4_max,col] = size(B4);

sample_period = c_ri_calc_dt(B1,1000,mode);

if sample_period == -1
  %calculating the sample fq
  if B1_max ~= 1 && B2_max ~= 1 && B3_max ~= 1 && B4_max ~= 1
    per1 = (B1(B1_max,1) - B1(1,1))/(B1_max-1);
    per2 = (B2(B2_max,1) - B2(1,1))/(B2_max-1);
    per3 = (B3(B3_max,1) - B3(1,1))/(B3_max-1);
    per4 = (B4(B4_max,1) - B4(1,1))/(B4_max-1);
  else
    per1 = -1;
    per2 = -1;
    per3 = -1;
    per4 = -1;
  end

  % The one with the highest sample fq should have lovest packet loss.
  sample_period = min([per1 per2 per3 per4]);
end

if B1(1,1) == -1
  B1_1 = NaN;
  B1_m = NaN;
else
  B1_1 = B1(1,1);
  B1_m = B1(B1_max,1);
end

if B2(1,1) == -1
  B2_1 = NaN;
  B2_m = NaN;
else
  B2_1 = B2(1,1);
  B2_m = B2(B2_max,1);
end

if B3(1,1) == -1
  B3_1 = NaN;
  B3_m = NaN;
else
  B3_1 = B3(1,1);
  B3_m = B3(B3_max,1);

end

if B4(1,1) == -1
  B4_1 = NaN;
  B4_m = NaN;
else
  B4_1 = B4(1,1);
  B4_m = B4(B4_max,1);

end

%finding start and end times
start_time = min([B1_1 B2_1 B3_1 B4_1]);
end_time = max([B1_m B2_m B3_m B4_m]);

if ~isnan(start_time) && ~isnan(end_time)

  time_line = (start_time:sample_period:end_time)';

  B1_c = c_ri_comple(time_line,B1, sample_period);
  B2_c = c_ri_comple(time_line,B2, sample_period);
  B3_c = c_ri_comple(time_line,B3, sample_period);
  B4_c = c_ri_comple(time_line,B4, sample_period);

else
  B1_c = [-1 0 0 0];
  B2_c = [-1 0 0 0];
  B3_c = [-1 0 0 0];
  B4_c = [-1 0 0 0];

end
