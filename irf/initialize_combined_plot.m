function [h1, h2] = initialize_combined_plot(type,nTimePanels,nRows,nCols,totCols,sorting)
% INITIALIZE_COMBINED_PLOT - Combination of irf time panels and subplot(nRows,nCols,i_subplot).
%
%   [h1,h2] = INITIALIZE_COMBINED_PLOT(type,nTimePanels,nRows,nCols,totCols OR timeseriesFraction,sorting);
%     Input:
%       type - specifies how the time- and non-time panels should be
%              organized with resect to each other: 'leftright' or 'topbottom'
%       nTimePanels - number of time panels (placed to the left or top)
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


switch type
  case 'leftright'
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
      for jj = (totCols-nCols+1):totCols
        for ii = 1:nRows
          isub = isub + 1;
          h2(isub) = subplot(nRows,totCols,jj+(ii-1)*totCols);
          h2left = [h2left; h2(isub).Position(1)];
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

      if strcmp('horizontal',sorting)
        h2 = reshape(h2,nRows,nCols);
        h2 = h2';
        h2 = h2(:);
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
  case 'topbottom'
    if totCols<1 % gives the fraction of the plot to be taken up by the timeseries
      space = totCols;
      height = space*0.8 - 0.1;
      ybot = 1-height - 0.1;
      ytop = ybot + height;
      panel_height = height/nTimePanels;
      for ii = 1:nTimePanels % move timeseries panels
        h1(ii).Position(2) = ybot + ((nTimePanels-ii))*panel_height;
        h1(ii).Position(4) = panel_height;
      end

      % Make "square" plots
      isub = 0;
      h2left = [];
      %totCols = 100;
      totRows = 100;
      space = 1-totCols;
      height = space*0.9;
      ybot = 0.1;
      ytop = ybot + height;
      panel_height = height/nRows;
      width = h1(1).Position(3);
      space_x_panels = 0.05;
      space_y_panels = 0.05;
      if nCols == 1
        width_panel = (width/nCols);
      else
        width_panel = (width/nCols)*0.8;
      end

      if nRows == 1
        panel_height = panel_height;
      else
        panel_height = panel_height*0.8;
      end

      if nCols == 1
        x_spacing = 0;
      else
        x_spacing = (width/(nCols-1))*(1-0.8);
      end
      if nRows == 1
        y_spacing = 0;
      else
        y_spacing = (panel_height/(nRows-1))*(1-0.8);
      end
      xleft = h1(1).Position(1);
      xright = xleft + width;

      for ii = nRows:-1:1
        for jj = 1:nCols
          isub = isub + 1;
          pos1 = xleft + (jj-1)*(width_panel+x_spacing);
          pos2 = ybot + (ii-1)*(panel_height+y_spacing);
          pos3 = width_panel;
          pos4 = panel_height;
          h2(isub) = axes('position',[pos1 pos2 pos3 pos4]);
          h2left = [h2left; h2(isub).Position(1)];
        end
      end

      if strcmp('vertical',sorting)
        h2 = reshape(h2,nCols,nRows);
        h2 = h2';
        h2 = h2(:);
      end
    end
end
end




