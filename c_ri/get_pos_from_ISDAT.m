function P_GSE = get_pos_from_ISDAT(spacecraft,start_time, Dt)
%
%Input:
% spacecraft-which spacecraft to get the data from: range 1-4.
% start_time-the time when to start the reading [yyyy mm dd hh mm ss].
% Dt-the duration of the reading in seconds.
%
%Output:
% P_GSE - The position in GSE [ time|x-component | y-component | z-component ]
%
%Descrition of the function:
%
%Using:
% Mat_DbOpen()
% isGetDataLite()
% Mat_DbClose()
% database unix:99 !not a function! 
%
%Work method:
% Open the databas, does the three readings and then closes the databas 
%
%Error: [-1,0, 0, 0] is returned, if it wasnt possible to load the data
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
clear P_GSE;
db = Mat_DbOpen('disco:10');

[tp, p_gse] = isGetDataLite(db, start_time, Dt, 'Cluster', num2str(spacecraft),'ephemeris', 'position', ' ', ' ', ' ');
Mat_DbClose(db);
p_gse=double(p_gse);

%if there is an error in the reading
if isempty(tp)
tp = -1;
p_gse = [0 0 0]';
end
 
P_GSE(:,1) = tp;
P_GSE(:,2:4) = p_gse';




