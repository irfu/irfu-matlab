function number = Ut2number(t1,length, t2,sr) 
%
%Input: 
% t1 is an epoch,lengt is the number of samples. 
% t2 is a Ut and sr is the samplerate
%
%Output:
% The number of steps from time t1 to time t2 with the samplerate sr 
% The number is in the range 1-length
%
%Descrition of the function:
% This function determins the number of steps between to times when a 
% certain samplerate i used. 
% The lower time in epoch and the higher time in Ut. 
%
%Using:
% toepoch
% 
%Work method:
%
%Error:
% There is an errorcheck to see if the returning number is smaller than 1.
% Then 1 is returned.
% There is an errorcheck to see if the time i longer that the number of
% samples. Then the last sample i returned.
%
%Discription of variables:
% t1 in epoch
% t2 in Ut
% sr in samples per second
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
t2 = toepoch(t2);
t_delta = t2-t1;
number = round(t_delta*sr);

if number < 1
number = 1;
end

if number > length
number = length;
end
