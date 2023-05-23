%
% Miscellaneous utility functions. Mostly to collect small shared functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils   % < handle

    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Generate text string with information on data source and when the plot
        % was generated.
        %
        function str = generate_data_source_info()
            dateStr = char(datetime("now", "Format", "uuuu-MM-dd"));
            str = sprintf( ...
                [ ...
                    'Swedish Institute of Space Physics, Uppsala (IRFU), %s. ', ...
                    'Data available at http://soar.esac.esa.int/' ...
                ], ...
                dateStr ...
            );
        end



        % ~Utility function that removes duplicated code from plot functions.
        %
        function ensure_axes_data_tick_margins(hAxesArray)
            assert(isa(hAxesArray, 'matlab.graphics.axis.Axes'))

            for i = 1:numel(hAxesArray)
                h = hAxesArray(i);

                h.YLim = solo.qli.ensure_data_tick_margins(...
                    [min(h.YTick), max(h.YTick) ], ...
                    [    h.YLim(1),    h.YLim(2)] ...
                );
            end
        end



        function filename = get_plot_filename(Tint)
            assert(isa(Tint, 'EpochTT'))
            assert(length(Tint) == 2)

            ett1 = Tint(1);
            utcStr1 = ett1.utc;
            utcStr1 = utcStr1(1:13);
            utcStr1([5,8])=[];

            ett2 = Tint(end);
            utcStr2 = ett2.utc;
            utcStr2 = utcStr2(1:13);
            utcStr2([5,8])=[];

            filename = [utcStr1,'_',utcStr2,'.png'];
        end



    end    % methods(Static)



end
