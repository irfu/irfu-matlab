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



    end    % methods(Static)



end
