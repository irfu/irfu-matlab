function out = hanning(n)
out = .5*(1 - cos(2*pi*(1:n)'/(n+1)));