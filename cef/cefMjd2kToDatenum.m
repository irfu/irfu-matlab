%
% Conversion between modified julian seconds to matlab internal date format.
% Modified julian seconds since 00:00 UT 2000.
% and Matlab format since: 00:00 Jan 1 0000.

function dnum=cefMjd2kToDatenum(mjd2k)

     dnum=mjd2k + 730486;
