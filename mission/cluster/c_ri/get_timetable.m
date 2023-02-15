function timetable = get_timetable(from,to,p_f)
%GET_TIMETABLE  read data from a table file
%
% timetable = get_timetable(from,to,p_f)
%
%	Finds a row starting with B and saves the time in start.
%	finds a row starting with E and saves the time in end.
%
% Input:
%	from, to - time in epoch
%	p_f - path+filename of the file with time of beginnings and end.
%
% Output:
%	timetable - [start_time | end_time]
%
% See also CREATE_FILE, CREATE_TIMETABLE
%

%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
narginchk(3,3)
if ~(isnumeric(from) && from >0), error('FROM must be epoch'), end
if ~(isnumeric(to) && to>from), error('TO must be larger then FROM'), end
if ~ischar(p_f), error('P_F must be a file path'), end

fp = fopen(p_f);
timetable = 0;

while feof(fp) == 0
  line = fgetl(fp);
  if isempty(line)
    disp('empty line')
    continue
  end
  switch line(1)
    case 'B'
      s_time = line2time(line);
    case 'E'
      e_time = line2time(line);
      timetable = addTime2Table(timetable,from,to,s_time,e_time);
    otherwise
      %disp('FGM get_timetable: bogus line')
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

if to >= s_time && e_time >= from
  s_t = max([from s_time]);
  e_t = min([to e_time]);
  t = [s_t e_t];
  timetable = add_A2M(timetable,t);
end

