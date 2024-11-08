function [x_reduced, y_reduced] = reduce_to_width(x, y, width, lims)

% x should be n x 1 matrix
% x_reduced is 2*witdh x 1 - matrix if reduction done
%
% [x_reduced, y_reduced] = reduce_to_width(x, y, width, lims)
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
% x can be n-by-1 or n-by-m with n samples of m columns (that is, there can
% be 1 x for all y or 1 x for each y.
%
% y must be n-by-m with n samples of m columns
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
    % If the data is already small, there's no need to reduce.

    % Find the starting and stopping indices for the current limits.
    % Rename the column. This actually makes stuff substantially
    % faster than referencing x(:, k) all over the place. On my
    % timing trials, this was 20x faster than referencing x(:, k)
    % in the loop below.
    
    % Map the lower and upper limits to indices.
    nx = size(x, 1);
    if isinf(lims(1)), lims(1) = x(1);end
    if isinf(lims(2)), lims(2) = x(end);end
    
    [~,lower_limit]      = binary_search(x, lims(1), 1, nx);
    if x(lower_limit) > lims(2)
      x_reduced = []; % no data
      y_reduced = [];
      return;
    else
      [upper_limit,~] = binary_search(x, lims(2), lower_limit, nx);
    end
    
    if upper_limit-lower_limit <= n_points
        x_reduced = x(lower_limit:upper_limit);
        y_reduced = y(lower_limit:upper_limit);
        return;
    end
    
    % Reduce the data to the new axis size.
    x_reduced = nan(n_points, 1);
    y_reduced = nan(n_points, size(y, 2));
      
    
    % Make the windows mapping to each pixel.
    x_divisions = linspace(lims(1), ...
      lims(2), ...
      width + 1);
    xx = linspace(  lims(1), lims(2) ,n_points+2);
    xx = xx(2:end-1)';

    for k = 1:size(y, 2)

  
        % Create a place to store the indices we'll need.
        %indices = [lower_limit, zeros(1, n_points-2), upper_limit];
        indices = nan(1, n_points);
        
        % For each pixel...
        right = lower_limit;
        for z = 1:width
            
            % Find the window bounds.
            left              = right;
            if left == upper_limit,continue; end % empty bins at the end
            if x(left) < x_divisions(z+1) && x(upper_limit) >= x_divisions(z+1)
            [right,~]         = binary_search(x, ...
                                               x_divisions(z+1), ...
                                               left, upper_limit);
            elseif x(upper_limit) < x_divisions(z+1) % end of time series
              right = upper_limit;
            else % empty bin
              if isfinite(x(right)) % remove the point if the last point in the previous bin is finite number
                xx(2*z-1:2*z) = NaN;
              end
              continue 
            end
            if left == right,continue; end % last empty bins
            % Get the indices of the max and min.
            yt = y(left:right, k);
            [~, max_index]     = max(yt);
            [~, min_index]     = min(yt);
            % Record those indices.
            indices(2*z-1:2*z) = sort([min_index max_index]) + left - 1;
        end

        % Sample the original x and y at the indices we found.
  %      x_reduced(:, k) = xt(indices);
  %      y_reduced(:, k) = y(indices, k);
      x_reduced(:,k) = xx;
      y_reduced(isfinite(indices), k) = y(indices(isfinite(indices)), k);
      
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
