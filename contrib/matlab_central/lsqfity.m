% lsqfity.m                                      by:  Edward T Peltzer, MBARI
%                                                revised:  2016 Mar 17.
%
% M-file to calculate a "MODEL-1" least squares fit.
%
%     The line is fit by MINIMIZING the residuals in Y only.
%
%     The equation of the line is:     Y = my * X + by.
%
%     Equations are from Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       pp: 104, 108-109, 199.
%
%     Data are input and output as follows:
%
%         [my,by,ry,smy,sby] = lsqfity(X,Y)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%
%             my    =    slope
%             by    =    y-intercept
%             ry    =    correlation coefficient
%             smy   =    standard deviation of the slope
%             sby   =    standard deviation of the y-intercept

% Source: Monterey Bay Aquarium Research Institute
% https://www.mbari.org/index-of-downloadable-files/

function [my,by,ry,smy,sby]=lsqfity(X,Y)

% Determine the size of the vector

n = length(X);

% Calculate the sums

Sx = sum(X);
Sy = sum(Y);
Sx2 = sum(X.^2);
Sxy = sum(X.*Y);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sx * Sy;
den = n * Sx2 - Sx^2;

% Calculate my, by, ry, s2, smy and sby

my = num / den;
by = (Sx2 * Sy - Sx * Sxy) / den;
ry = num / (sqrt(den) * sqrt(n * Sy2 - Sy^2));

diff = Y - by - my .* X;

s2 = sum(diff .* diff) / (n-2);
smy = sqrt(n * s2 / den);
sby = sqrt(Sx2 * s2 / den);
