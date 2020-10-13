%
% Create human-readable multi-line string to represent SETTINGS. Meant for
% logging and printing to stdout.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-22
%
function str = sprint_SETTINGS(SETTINGS)
    
    % PROPOSAL: Make hierarchy visually clearer?!!! Should then have help from data structure itself.
    %
    % PROPOSAL: Print more information.
    %   PROPOSAL: First column of characters to represent how values have been overridden and where.
    %       Ex: C    OUTPUT_CDF.WRITE_POLICY        % C  = Config file (once)
    %       Ex: CC   OUTPUT_CDF.WRITE_POLICY        % CC = Config file twice (or more)
    %       Ex: C AA OUTPUT_CDF.WRITE_POLICY        % C  = Config file (once), overridden by CLI argument (twice or more)
    %   PROPOSAL: Extra indented rows when default value is overridden.
    %       NOTE: Must be clear which values are actually used. Not to be confused with overridden values.
    %       OUTPUT_CDF.WRITE_POLICY = <value actually used>
    %           Default      = ...
    %           Config file  = ...
    %           Config file  = ...   % Second value in config file.
    %           CLI argument = ...
    %       CON: Slightly less clear when not having a separat column for overriding.
    %       PRO: Han handle any situation of overriding.
    %
    
    % IMPLEMENTATION NOTE: Only prints "Settings" as a header (not "constants")
    % to indicate/hint that it is only the content of the "SETTINGS" variables,
    % and not of bicas.constants.
    str = sprintf([...
        '\n', ...
        'SETTINGS\n', ...
        '========\n']);
    
    % Values seem sorted from the method, but sort again just to be sure.
    keyList      = sort(SETTINGS.get_keys());   
    lengthMaxKey = max(cellfun(@length, keyList));
    
    
    
    for iKey = 1:length(keyList)
        key   = keyList{iKey};
        valueStructArray = SETTINGS.get_final_value_array(key);
        %value = valueStructArray(end).value;
        nValues = numel(valueStructArray);
        
        %======================================================================
        % Derive value strings for all historical values: present and previous
        % ones
        %======================================================================
        strValueList = {};   % Must be reset for every key.
        for iVs = 1:nValues
            value = valueStructArray(iVs).value;
            
            if ischar(value)
                
                strValue = ['"', value, '"'];
                
            elseif isnumeric(value)
                
                EJ_library.assert.vector(value)
                if isscalar(value)
                    strValue = sprintf('%d', value);
                else
                    strArray = EJ_library.str.sprintf_many('%d', value);
                    strValue = sprintf('[%s]', strjoin(strArray, ', '));
                end

            elseif iscell(value)

                EJ_library.assert.vector(value)
                strValueCa = {};
                for i = 1:numel(value)
                    cellValue = value{i};
                    if isnumeric(cellValue) && isscalar(cellValue)
                        strValueCa{i} = sprintf('%g', cellValue);
                    elseif ischar(cellValue)
                        strValueCa{i} = sprintf('"%s"', cellValue);
                    else
                        error(...
                            'BICAS:sprintf_settings:IllegalCodeConfiguration', ...
                            'Can not print setting for log since cell array component is neither scalar numeric nor string.')
                    end
                end
                strValue = sprintf('{%s}', strjoin(strValueCa, ', '));

            else

                error(...
                    'BICAS:sprintf_settings:Assertion', ...
                    ['SETTINGS value (overriden or not) for key="%s" has illegal MATLAB class.', ...
                    ' It is neither char, numeric, nor 1D cell array.'], key)

            end
            strValueList{iVs} = strValue;
            clear strValue value
        end
        
        valueStatusStr   = EJ_library.utils.translate({...
            {'default'},            '  --';
            {'configuration file'}, '(conf)';
            {'CLI arguments'},      '(CLI)'}, ...
            valueStructArray(end).valueSource, ...
            'BICAS:sprintf_settings:Assertion', ...
            'Illegal setting value source');
        
        str = [str, sprintf(...
            ['%-6s  %-', int2str(lengthMaxKey),'s = %s\n'], ...
            valueStatusStr, key, strValueList{end})];
        
        % Print previous values.
        %     if nValues > 1
        %         for iVs = 1:nValues
        %             str = [str, sprintf(['            %s = %s\n'], valueStructArray(iVs).valueSource, strValueList{iVs})];
        %         end
        %     end
        
    end
    
    str = [str, newline];
    str = [str, sprintf('Explanations for leftmost column above:\n')];
    str = [str, sprintf('---------------------------------------\n')];
    str = [str, sprintf('  --   = Default value\n')];
    str = [str, sprintf('(conf) = Value comes from configuration file\n')];
    str = [str, sprintf('(CLI)  = Value comes from CLI argument\n')];
    
    str = [str, newline];
    
end
