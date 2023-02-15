function date_string=epoch2yyyymmdd(isdat_epoch)
%EPOCH2YYYYMMDD   Convert ISDAT epoch to YYYYMMDD
%
% date_string=epoch2yyyymmdd(isdat_epoch)
% epoch2yyyymmdd convert isdat epoch to string YYYYMMDD
%

t=fromepoch(isdat_epoch);
date_string=sprintf('%04d%02d%02d',t(1),t(2),t(3));

