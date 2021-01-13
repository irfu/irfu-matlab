function [x_reduced, y_reduced, c_reduced] = reduce_spec_to_width(x, y, c, width, lims)

% x should be n x 1 matrix
% y should be vector or m x n matrix
% c should be m x n matrix
%
% x_reduced is 2*witdh x 1 - matrix if reduction done
%
% [x_reduced, y_reduced, c_reduced] = reduce_spec_to_width(x, y, c, width, lims)
% 
% (This function is primarily used by LinePlotReducer, but has been
%  provided as a stand-alone function outside of that class so that it can
%  be used potentially in other projects. See help LinePlotReducer for 
%  more.)
%
% Reduce the data contained in x and y for display on width pixels,
% stretching from lims(1) to lims(2).
%
% For example, if x and y are 20m points long, but we only need to plot
% them over 500 pixels from x=5 to x=10, this function will return 1000
% points representing the first point of the window (x=5), the last point
% of the window (x=10), and the maxima and minima for each pixel between 
% those two points. The result is that x_reduced and y_reduced can be 
% plotted and will look exactly like x and y, without all of the waste of 
% plotting too many points.
%
% [xr, yr] = reduce_to_width(x, y, 500, [5 10]); % Reduce the data.
%
% plot(xr, yr); % This contains many fewer points than plot(x, y) but looks
%                 the same.
%
% Tucker McClure
% Copyright 2013, The MathWorks, Inc.

    % We'll need the first point to the left of the limits, the first point
    % to the right to the right of the limits, and the min and max at every
    % pixel inbetween. That's 1 + 1 + 2*(width - 2) = 2*width total points.
    n_points = 2*width;
    
    % Find the starting and stopping indices for the current limits.
    % Rename the column. This actually makes stuff substantially
    % faster than referencing x(:, k) all over the place. On my
    % timing trials, this was 20x faster than referencing x(:, k)
    % in the loop below.
    
    % when interpolating, start point moves to start of bin and end point
    % moves to the end of the previous bin. when analyzing each bin one
    % needs to take care of both end points without corresponding start
    % points and startpoints.   
    % each spectra point consists of start and stop time
    % indSpectra = ceil(indX/2);
    
    % Map the lower and upper limits to indices.
    [nSpecLines,nx] = size(c);
    if isinf(lims(1)), lims(1) = x(1);end
    if isinf(lims(2)), lims(2) = x(end);end
    
    % check if y is matrix of the same size as c
    isYMatrix = all(size(y) == size(c));
    if ~isYMatrix, y_reduced = y;end % assume y is vector 
    
    [~,indLowerLimit]      = binary_search(x, lims(1), 1, nx);
    if x(indLowerLimit-mod(indLowerLimit,2)) > lims(2) % no data in the interval
      x_reduced = []; 
      y_reduced = [];
      c_reduced = [];
      return;
    else
      [indUpperLimit,~] = binary_search(x, lims(2), indLowerLimit, nx);
    end
    
    % no reduction needed
    if indUpperLimit-indLowerLimit <= n_points
      indInterval = (indLowerLimit-(mod(indLowerLimit,2)==0)) : ...
                    (indUpperLimit+(mod(indUpperLimit,2)==1));
        x_reduced = x(indInterval);
        c_reduced = c(:,indInterval);
        if isYMatrix
          y_reduced = y(:,indInterval);
        end
        return;
    end
    
    % Reduce the data to the new axis size.
    % x_reduced = nan(n_points, 1);  % TODO swap xx to x_reduced
    c_reduced = nan(nSpecLines,n_points);
    if isYMatrix, y_reduced = c_reduced; end
      
    
    % Make the windows mapping to each pixel.
    x_divisions = linspace(lims(1), ...
      lims(2), ...
      width + 1 );
    xx = [x_divisions(1:end-1) ; x_divisions(2:end)];
    xx = xx(:)';

    % Create a place to store the indices we'll need.
    %indices = [lower_limit, zeros(1, n_points-2), upper_limit];
    indices = nan(1, n_points);
    
    % For each pixel...
    right = indLowerLimit-1;
    
    for z = 1:width % go through bins
      
      left              = right;
      isSpectraCrossingBinBoundary = (mod(left,2) == 1);
      
      if left == indUpperLimit % empty bins at the end
        if isSpectraCrossingBinBoundary % open spectra until end
          indices(2*z-1)=left;
          xx(2*z:end-1) = NaN;
        else  % no open spectra, remove all last points
           xx(2*z-1:end) = NaN;
        end
        break;
      end
      
       [right,~]         = binary_search(x, ...
          x_divisions(z+1), ...
          left, indUpperLimit);
        
        if x(right) > x_divisions(z+1) % data gap
          if isSpectraCrossingBinBoundary % open spectra
            indices(2*z-1)=left;
          else  % no open spectra, remove the bin
            xx(2*z-1:2*z) = NaN;
          end
          right = left;
          continue;
        end
             
        if isSpectraCrossingBinBoundary % open spectra until end
          indices(2*z-1)=left;
        else  % no open spectra, remove all last points
           indices(2*z-1)=left+1;
        end
    end
    
    x_reduced = xx;
    c_reduced(:,isfinite(indices)) = c(:,indices(isfinite(indices)));
    
    indEmpty = isnan(x_reduced);
    x_reduced(indEmpty)=[];
    c_reduced(:,indEmpty)=[];
    if isYMatrix
      y_reduced(:,isfinite(indices)) = ...
        y(:,indices(isfinite(indices)));
      y_reduced(:,indEmpty)=[];
      y_reduced(:,2:2:end)=y_reduced(:,1:2:end);
    end
    
end

% Binary search to find boundaries of the ordered x data.
function [L, U] = binary_search(x, v, L, U)
    while L < U - 1                 % While there's space between them...
        C = floor((L+U)/2);         % Find the midpoint
        if x(C) < v                 % Move the lower or upper bound in.
            L = C;
        else
            U = C;
        end
    end
end

% value index 

