function timetable = get_timetable(from,to,p_f) 
%
%Input:
% p_f - the path and filename to the file with time of beginnings and end.
%
%Output:
% timetable - [start_time | end_time]
%
%Descrition of the function:
% Finds a row starting with B and saves the time in start.
% finds a row starting with E and saves the time in end.
%
%Using:
% 
%Work method:
%
%Error:
% timetable will be 0 if no timetable is found
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

fp = fopen(p_f);
timetable = 0;

i = 0;
while feof(fp) == 0
    line = fgetl(fp);
    if(line(1) == 'B')
		s_time = line2time(line);
    end
    if(line(1) == 'E')
        e_time = line2time(line);
		timetable = addTime2Table(timetable,from,to,s_time,e_time);
    end
end
fclose(fp);



%converts the line with the time in a string to time in epoch
function time = line2time(line)
e =length(line)-5;
t_line = line(3:e);

yyyy_mm_dd = t_line(1:10);
mm_dd_yyyy(4:5) = yyyy_mm_dd(9:10);
mm_dd_yyyy(1:2) = yyyy_mm_dd(6:7);
mm_dd_yyyy(7:10) = yyyy_mm_dd(1:4);
mm_dd_yyyy(3) = '-';
mm_dd_yyyy(6) = '-';

yyyy_mm_dd = datevec(mm_dd_yyyy);

hr_mm_ss = t_line (12:length(t_line));
hr_mm_ss = datevec(hr_mm_ss);

time_UT = [yyyy_mm_dd(1:3) hr_mm_ss(4:6)];
time = toepoch(time_UT);

function timetable = addTime2Table(timetable,from,to,s_time, e_time)

if to >= s_time & e_time >= from
s_t = max([from s_time]);
e_t = min([to e_time]);

t = [s_t e_t];
timetable = add_A2M(timetable,t);

end

