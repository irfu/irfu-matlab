% lsqfitx.m                                      by:  Edward T Peltzer, MBARI
%                                                revised:  2016 Mar 17.
%
% M-file to calculate a "MODEL-1" least squares fit.
%
%     The line is fit by MINIMIZING the residuals in X only.
%
%     The equation of the line is:     Y = mx * X + bx.
%
%     Equations are modified from those in Bevington & Robinson (1992)
%       Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
%       pp: 104, 108-109, 199.
%
%     Data are input and output as follows:
%
%         [mx,bx,rx,smx,sbx] = lsqfitx(X,Y)
%
%             X     =    x data (vector)
%             Y     =    y data (vector)
%
%             mx     =    slope
%             bx     =    y-intercept
%             rx     =    correlation coefficient
%             smx    =    standard deviation of the slope
%             sbx    =    standard deviation of the y-intercept

function [mx,bx,rx,smx,sbx]=lsqfitx(X,Y)

% Determine the size of the vector

n = length(X);

% Calculate the sums

Sx = sum(X);
Sy = sum(Y);
Sx2 = sum(X.^2);
Sxy = sum(X.*Y);
Sy2 = sum(Y.^2);

% Calculate re-used expressions

num = n * Sxy - Sy * Sx;
den = n * Sy2 - Sy^2;

% Calculate m, a, rx, s2, sm, and sb

mxi = num / den;
a = (Sy2 * Sx - Sy * Sxy) / den;
rx = num / (sqrt(den) * sqrt(n * Sx2 - Sx^2));

diff = X - a - mxi .* Y;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sa = sqrt(Sy2 * s2 / den);

% Transpose coefficients

mx = 1 / mxi;
bx = -a / mxi;

smx = mx * sm / mxi;
sbx = abs(sa / mxi);
