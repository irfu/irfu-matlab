%
% Conversion between modified julian seconds to matlab internal date format.
% Modified julian seconds since 00:00 UT on Jan 1, 1950.
% and Matlab format since: 00:00 Jan 1 0000.

function mjs=cefDatenumToMjs(date)


    mjs= round((date-712224)*86400);
