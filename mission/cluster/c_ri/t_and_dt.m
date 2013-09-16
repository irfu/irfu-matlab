function t_dt = t_and_dt(t,dt) 
%
%Input:
% t -containing the time t in [yyyy mm dd hh mm ss.mss]
% dt - containing the duration dt of each time, in seconds
%
%Output:
% t_dt -one matrix with [t |dt ]
% where t are in epochs and dt in seconds
%
%Descrition of the function:
% Make one single matris of the time t and dt. Converts the time t
% into epochs
%
%Using:
% toepoch
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
[rows,col] = size(t);
for n = 1:rows
t_dt(n,1) = toepoch(t(n,:));
end
t_dt(:,2) = dt;
