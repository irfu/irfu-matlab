function B = load_file(f_path,name)
%
% B = load_file(f_path,name)
%
%Input:
% f_path -path to file
% name -file to load. this file should contain B-data
%
%Output:
% B -the B-matrix [time,Bx,By,Bz]
%
%Descrition of the function:
%
%Using:
% timestr2epoch
% hhmmss2epoch
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

B =[-1 0 0 0];
S=0; 

p_and_f = sprintf('%s%s',f_path,name);
disp(['loading: ' name])
%tic
fp = fopen(p_and_f);
[S,count] = fscanf(fp,'%4d%1s%2d%1s%2d%1s%2d%1s%2d%1s%6f%2s%f%f%f',[15,inf]);
fclose(fp);

[r_s, c_s] = size(S);
if r_s == 15
temp_B = S([1 3 5 7 9 11 13 14 15],:)';
B = [toepoch(temp_B(:,1:6)) temp_B(:,7:9)]; 
%[time, Bx, By, Bz] = textread(p_and_f,'%s %f %f %f');
%time_of_load = toc;
%disp(['loading took ' int2str(time_of_load) ' seconds'])
end
