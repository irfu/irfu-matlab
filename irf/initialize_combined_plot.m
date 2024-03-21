function [h1, h2] = initialize_combined_plot(nTimePanels,nRows,nCols,totCols,sorting)
% INITIALIZE_COMBINED_PLOT - Combination of irf time panels and subplot(nRows,nCols,i_subplot).
%
%   [h1,h2] = INITIALIZE_COMBINED_PLOT(nTimePanels,nRows,nCols,totCols OR timeseriesFraction,sorting);
%     Input:
%       nTimePanels - number of time panels (placed to the left)
%       nRows - number of rows
%       totCols - total number of columns (used to scale the time panels,
%                 the nRows leftmost are used for the subplots)
%       timeseriesFraction - give the fraction of the figure to be occupied
%                            by the time series panels
%       nCols - rightmost number of columns to make plots in
%       sorting - 'vertical' or ' 'horizontal' (default), sorts axis
%                 handles in given direction
%     Output:
%       h1 - axes handles to time panels: h1 = irf_plot(nTimePanels);
%       h2 - axes handles to subplots
%
%   Examples:
%     [h1,h2] = INITIALIZE_COMBINED_PLOT(8,2,2,0.4,'vertical')
%     [h1,h2] = INITIALIZE_COMBINED_PLOT(8,2,2,3,'vertical')

scrsz = get(groot,'ScreenSize');
%figure('Position',scrsz)
h1 = irf_plot(nTimePanels);
if isempty(sorting)
  sorting = 'horizontal';
end

if totCols<1 % gives the fraction of the plot to be taken up by the timeseries
  for ii = 1:nTimePanels % move timeseries panels
    space = totCols;
    h1(ii).Position(3) = space*0.8;
    h1(ii).Position(1) = h1(ii).Position(1)*0.4;
  end

  % Make "square" plots
  isub = 0;
  h2left = [];
  totCols = 100;

  if strcmp('horizontal',sorting)
    for ii = 1:nRows
      for jj = (totCols-nCols+1):totCols
        isub = isub + 1;
        h2(isub) = subplot(nRows,totCols,jj+(ii-1)*totCols);
        h2left = [h2left; h2(isub).Position(1)];
      end
    end
  else % sorts the handles vertically
    for jj = (totCols-nCols+1):totCols
      for ii = 1:nRows
        isub = isub + 1;
        h2(isub) = subplot(nRows,totCols,jj+(ii-1)*totCols);
        h2left = [h2left; h2(isub).Position(1)];
      end
    end
  end

  % Move panels so they have the appropriate width
  h1right = h1(1).Position(1)+h1(1).Position(3)+0.03;
  x0 = totCols;
  isub = 0;
  for ii = 1:nCols
    for jj = 1:nRows
      isub = isub + 1;
      h2(isub).Position(1) = h1right*1.1 + (1-h1right)/nCols*(ii-1)*0.9; % set the left edge
      h2(isub).Position(3) = (1-h1right)/nCols*0.7; % set the width

    end
  end

else % gives the the "total number of columns" for scaling
  % Adjust the width of the irf_plot panels
  for ii = 1:nTimePanels
    space = 1-(nCols + 1)/(totCols);
    h1(ii).Position(3) = space*0.5;
    h1(ii).Position(1) = h1(ii).Position(1)*0.4;
  end

  % Make the subplot panels
  isub = 0;
  h2left = [];

  if strcmp('horizontal',sorting)
    for ii = 1:nRows
      for jj = (totCols-nCols+1):totCols
        isub = isub + 1;
        h2(isub) = subplot(nRows,totCols,jj+(ii-1)*totCols);
        h2left = [h2left; h2(isub).Position(1)];
      end
    end
  else % sorts the handles vertically
    for jj = (totCols-nCols+1):totCols
      for ii = 1:nRows
        isub = isub + 1;
        h2(isub) = subplot(nRows,totCols,jj+(ii-1)*totCols);
        h2left = [h2left; h2(isub).Position(1)];
      end
    end
  end

  % Adjust irf_plot panels again
  h2left = min(h2left);
  for ii = 1:nTimePanels
    space = 1-(nCols + 1)/(totCols);
    h1(ii).Position(3) = (h2left-h1(ii).Position(1))*0.8;
  end
end




