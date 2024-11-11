%% Examples of LinePlotReducer
%
% These cells show several different ways to use a LinePlotReducer to plot
% large amount of data gracefully. In some cases, LinePlotExplorer objects
% are created as well to allow the user to quickly examine a plot. It's
% recommended that this script *not* be executed in its entire, but rather
% cell-by-cell. This can be most easily done by clicking into a cell and
% pressing ctrl+enter to run only that cell. All cells require the
% Initialize cell to have been run, leaving |t| and |x| in the workspace.
%
% This script also acts as a small unit test of the LinePlotReducer.
%
% --- Change Log ---
%
% 2015-05-31: Minor spelling changes.
%
% 2014-06-04: Tests having multiple LinePlotReducers in a single axes.
% Changes copyright 2014 Tucker McClure.
% 
% 2013-03-15: Original. Copyright 2013, The MathWorks, Inc.
%
% ---
%
% Copyright 2015, The MathWorks, Inc. and Tucker McClure

%% Make sure user is running one cell at a time.
a = questdlg(['This script will perform numerous actions with quite ' ...
              'a bit of data. It''s recommended to run this script '...
              'one cell at a time. Are you sure you wish to continue?'],...
             'Continue?', ...
             'Yes', 'No', 'No');
if ~strcmp(a, 'Yes')
    return;
end

%% Initialize
% Create the data and set some options.

clc;

% Create the data.
n = 20e6 + randi(1000);                          % Number of samples
t = sort(100*rand(1, n));                       % Non-uniform sampling
x = [sin(0.10 * t) + 0.05 * randn(1, n); ...
     cos(0.43 * t) + 0.001 * t .* randn(1, n); ...
     round(mod(t/10, 5))];
x(:, t > 40 & t < 50) = 0;                      % Drop a section of data.
x(randi(numel(x), 1, 20)) = randn(1, 20);       % Emulate spikes.

% Ignore certain code analyzer warnings.
%#ok<*NASGU>

%% default
% Show the LinePlotReducer at the simplest level and add a LinePlotExplorer
% to explore the new plot.
o = LinePlotReducer(x');
lpe = LinePlotExplorer(gcf());

%% plot
% Pass |plot| into the LPR and specify x and y values.
lpr = LinePlotReducer(@plot, t', x');
lpe = LinePlotExplorer(gcf(), 0, 100);

%% plot
% Try it with matrix inputs.
lpr = LinePlotReducer(@plot, t', [x' x' x']);

%% plot
% See if it handles linespec values gracefully.
lpr = LinePlotReducer(t, x(1, :), 'b:', ...
                      t, x(2, :), 'g', ...
                      t, x(3, :), 'r--*');

%% plot
% What about a matrix with a linespec?
o = LinePlotReducer(@plot, t, x, ':');

%% plot
% Lots of data with no linespec?
o = LinePlotReducer(@plot, t, x(1, :), t, x(2, :), t, x(3, :));

%% semilogx
% What about other plot options, like semilogx.
o = LinePlotReducer(@semilogx, t, abs(x));

%% semilogy
% ... and semilogy (this might look a little funny).
o = LinePlotReducer(@semilogy, t, abs(x));

%% stairs
% Only a single line can be plotted at a time for |stairs| due to some
% subtleties about how |stairs| accepts its arguments (they can't be
% different sizes, and the line reduction algorithm may make each line a
% different size).
o = LinePlotReducer(@stairs, t, x(1, :), 'b:');

%% stairs
% If we try to use multiple lines, LinePlotReducer will provide some
% feedback and an error. We'll catch the error and carry on. This cell just
% tests that we generate the correct error.
fprintf(['We expect the following to fail with a message that the LPR ' ...
         'can''t manage the stairs function.']);
try
    o = LinePlotReducer(@stairs, t, x, 'b:');
catch er
    fprintf('The LinePlotReducer fails correctly.\n');
end

%% plotyy
% It works with two axes as well, such as are created by plotyy.
o = LinePlotReducer(@plotyy, t, x(1, :), t, x(2, :));
lpe = LinePlotExplorer(gcf(), 0, 100);

%% subplot
% Subplots work just fine.

h1 = subplot(2, 1, 1);
LinePlotReducer(t, x);
h2 = subplot(2, 1, 2);
LinePlotReducer(t, 1-2*x);

%% linkaxes
% We can link the axes and zoom on one subplot and see the LinePlotReducer
% manage the second plot too.
linkaxes([h1 h2]);

%% Can we add a LPE?
% We can explore a subplot.
lpe = LinePlotExplorer(gcf(), 0, 100);

%% Test attached callbacks.
% The LinePlotExplorer uses these callback functions. If we want something
% else to happen for each callback, we can attach a new function. When the
% LinePlotExplorer is done, it will call whatever we attach.
lpe.AttachCallback('WindowButtonDownFcn',   'disp(''down'');', ...
                   'WindowButtonUpFcn',     'disp(''up'');', ...
                   'WindowButtonMotionFcn', 'disp(''moving'');', ...
                   'WindowScrollWheelFcn',  @(~, ~) disp('scrolling'));

%% Ok, but can we stop the LPE?
% Note this doesn't stop the callbacks that we attached. It just stops the
% scrolling-panning behavior.
Stop(lpe);

% Remove the callbacks we set.
set(gcf(), 'WindowButtonDownFcn',   [], ...
           'WindowButtonUpFcn',     [], ...
           'WindowButtonMotionFcn', [], ...
           'WindowScrollWheelFcn',  []);


%% Takeover
% A LPR can take over from an existing plot.
clf();
h = plot(t, x(1, :));
lpr = LinePlotReducer(h);

%% Takeover
% A LPR can take over from multiple existing plots all at once.
h = plot(t, x);
lpr = LinePlotReducer(h);

%% Can we add a LPE?
% This still works fine with everything else.
lpe = LinePlotExplorer(gcf(), 0, 100);

%% Takeover (specifiying x and y)
% We can also specify (potentially bigger) x and y values when we let it 
% take over the existing plot.
h = plot([0 100], [0 0], 'r:', [0 100], [0 0], 'g:', [0 100], [0 0], 'b:');
lpr = LinePlotReducer(h, t, x);
lpe = LinePlotExplorer();

%% Function interface
% Not everyone wants to use an object for the interface, so here's a simple
% function. It takes the same arguments as LinePlotReducer, but returns the
% plot handles, just like a regular plot function would.
h = reduce_plot(t, x);
set(h, 'Color', 'r');

%% Rows vs. columns
% Let's make sure it works if we mix up rows and columns. Here, t is a row,
% but x.' has data in columns.
reduce_plot(t, x.');

%%
% If we have an implied x, can it figure out that this is a single series?
reduce_plot(x(1, :));

%%
% If we have an implied x, can it figure this out too?
reduce_plot(x(1, :).');

%%
% If we have an implied x, can it find the line spec?
reduce_plot(x(1, :).', 'r');

%%
% What if t is a column but x has data in rows?
reduce_plot(t.', x);

%%
% What about multiple LinePlotReducers in a single axis?

clf();
hold on;
LinePlotReducer(t, x(1, :), 'r');
LinePlotReducer(t, x(2, :), 'b');
LinePlotReducer(t, x(3, :), 'g');
LinePlotExplorer();
hold off;
