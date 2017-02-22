% Constants - Singleton class for global constants used by BICAS.
% 
% Class for storing 
% 1) settings (variable values) which could reasonably be set via some user interface (CLI, configuration file), and
% 2) settings and constants which could reasonably be printed for the user.
%
% The class completely encapsulates the settings. Settings are identfied by arbitrary strings. Settings can only be
% obtained and set via methods.
% 
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-22
%
classdef settings < handle
% BOGIQ: 
%
% PROPOSAL: Rename SIMPLE_DEMUXER. SIMPLE_CALIBRATION?
% PROPOSAL: Validation code.
%   PROPOSAL: Check that master CDFs exist, that paths exist.
%       CON: Makes sense to make that kind of check here?!
    
    properties(Access=private)
        map = containers.Map('KeyType', 'char', 'ValueType', 'any')    % Bad name. "storage"? "SETTINGS.g"?
    end
    
    %###################################################################################################################
    
    methods(Access=public)
        
        function obj = settings
            
            %-------------------------------------------------------------------------------------
            % Values common to multiple constants
            % -----------------------------------
            % Only used INDIRECTLY and only INTERNALLY to set the values of the "real" constants.
            %-------------------------------------------------------------------------------------
            D = [];
            D.AUTHOR_NAME  = 'Erik P G Johansson';
            D.AUTHOR_EMAIL = 'erik.johansson@irfu.se';
            D.INSTITUTE    = 'IRF-U';
            %D.SWD_OUTPUT_RELEASE_VERSION = '01';  % For the S/W descriptor output CDFs' release version. Unknown what a sensible value is.



            obj.set('AUTHOR_NAME',  D.AUTHOR_NAME);
            obj.set('AUTHOR_EMAIL', D.AUTHOR_EMAIL);
            obj.set('INSTITUTE',    D.INSTITUTE);
            obj.set('MASTER_CDFS_RELATIVE_DIR', 'data');    % Location of master CDF files. Relative to the software directory structure root.
            
            % Value that shows up in EOut dataset GlobalAttributes.Calibration_version.
            % String value.
            obj.set('CALIBRATION_VERSION', '0.1; Only proportionality constants i.e. no voltage offset tables, no transfer functions; No bias currents');
        
            
            
            %===========================================================================================================
            % Various S/W descriptor release data for the entire software (not specific outputs)
            % ----------------------------------------------------------------------------------
            % EXCEPTION TO VARIABLE NAMING CONVENTION: Field names are used for constructing the JSON object struct and
            % can therefore NOT follow variable naming conventions without modifying other code.
            %===========================================================================================================
            obj.set('SWD_IDENTIFICATION.project',     'ROC-SGSE');
            obj.set('SWD_IDENTIFICATION.name',        'BICAS');
            obj.set('SWD_IDENTIFICATION.identifier',  'ROC-SGSE-BICAS');
            obj.set('SWD_IDENTIFICATION.description', 'BIAS Calibration Software (BICAS) which derives the BIAS L2S input signals (plus some) from the BIAS L2R output signals.');
            %
            obj.set('SWD_RELEASE.version',      '0.1.0');
            obj.set('SWD_RELEASE.date',         '2017-02-22');
            obj.set('SWD_RELEASE.author',       D.AUTHOR_NAME);
            obj.set('SWD_RELEASE.contact',      D.AUTHOR_EMAIL);
            obj.set('SWD_RELEASE.institute',    D.INSTITUTE);
            obj.set('SWD_RELEASE.modification', 'No modification (initial release)');
            %
            obj.set('SWD_ENVIRONMENT.executable', 'roc/bicas');     % Relative path to BICAS executable. See RCS ICD.
            
            
            
            % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
            obj.set('STDOUT_PREFIX', 'STDOUT: ');
        
            % Parameters influencing how JSON objects are printed with function JSON_object_str.
            obj.set('JSON_OBJECT_STR.INDENT_SIZE',     4);
            obj.set('JSON_OBJECT_STR.VALUE_POSITION', 15);
        
            % The epoch for ACQUISITION_TIME.
            % The time in UTC at which ACQUISITION_TIME is [0,0].
            % Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
            % PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
            obj.set('ACQUISITION_TIME_EPOCH_UTC', [2000,01,01, 12,00,00, 000,000,000]);

            obj.set('INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID',       0);   % Require input CDF Global Attribute "DATASET_ID"       to match the expected value.
            obj.set('INPUT_CDF_ASSERTIONS.STRICT_SKELETON_VERSION', 1);   % Require input CDF Global Attribute "Skeleton_version" to match the expected value.
            obj.set('INPUT_CDF_ASSERTIONS.MATCHING_TEST_ID',        0);   % Require Test_id to be identical for all input CDF datasets.
            obj.set('OUTPUT_CDF.SET_TEST_ID',  1);           % Set CDF GlobalAttribute "Test_id". ROC DFMD says that it should really be set by ROC.            
            obj.set('OUTPUT_CDF.DATA_VERSION', '01');        % Set CDF GlobalAttribute "Data_version". ROC DFMD says it should be updated in a way which can not be automatized?!!! Set here for now.

            obj.set('PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION', 1);
            
            % zVariables which are still empty after copying data into the master CDF assigned a correctly sized array
            % with fill values. This should only be necessary for S/W modes with incomplete processing.
            obj.set('OUTPUT_CDF.EMPTY_ZVARIABLES_SET_TO_FILL', 0);
            
            obj.set('LOGGING.MAX_UNIQUES_PRINTED', 5);    % When logging contents of matrix/vector, maximum number of unique values printed before switching to shorter representation (min-max range)
            obj.set('LOGGING.IRF_LOG_LEVEL', 'notice');   % Log level for "irf.log".
            
            
            
            %=====================================================================
            % Define constants relating to interpreting LFR datasets
            % ------------------------------------------------------
            % F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz. Names are LFR's naming.
            %=====================================================================
            obj.set('LFR.F0', 24576);  % = 6 * 4096
            obj.set('LFR.F1',  4096);
            obj.set('LFR.F2',   256);
            obj.set('LFR.F3',    16);
        
            %========================================================
            % Constants for how the "simple demuxer" calibrates data
            %========================================================
            obj.set('SIMPLE_DEMUXER.ALPHA',           1/17);
            obj.set('SIMPLE_DEMUXER.BETA',               1);
            obj.set('SIMPLE_DEMUXER.GAMMA_HIGH_GAIN',  100);
            obj.set('SIMPLE_DEMUXER.GAMMA_LOW_GAIN',     5);   % NOTE/POSSIBLE BUG: Uncertain which value is high-gain, and low-gain.
                        
        end
        
        
        
        function value = get(obj, key)
            if ~obj.map.isKey(key)
                error('BICAS:settings:Assertion:IllegalArgument', 'There is no setting "%s".', key)
            end
            
            value = obj.map(key);
        end
        
        
        
        function keyList = get_keys(obj)
            keyList = obj.map.keys;
        end
        
        
        
        %function value = get_struct(obj, keyRoot)
        %    keyList = obj.map.keys;
        %    
        %    for iKey = 1:length(keyList)
        %        for 
        %    end
        %end

        
        
        % Modify settings
        %
        % ARGUMENTS
        % =========
        % ModifiedSettingsAsStrings : containers.Map with
        %   keys   = Recursive struct names / settings names
        %   values = Settings values as strings.
        %
        %
        % NOTE: This function is only supposed to be called only once, and as soon as possible after the constants
        % object has been initialized.
        %
        % IMPLEMENTATION NOTE: This function can NOT be trivially merged with the constructor since
        % (1) "constants" have to be initialized before parsing CLI arguments (for S/W modes).
        % (2) "constants"/settings can be modified by the CLI arguments.
        %
        function modify_settings(obj, ModifiedSettingsAsStrings)
            
            keysList = ModifiedSettingsAsStrings.keys;
            for iModifSetting = 1:length(keysList)
                key = keysList{iModifSetting};
                newValueAsString = ModifiedSettingsAsStrings(key);
                
                % Use old value to convert string value to appropriate MATLAB class.
                oldValue = obj.get(key);
                if isnumeric(oldValue)
                    newValue = str2double(newValueAsString);
                elseif ischar(oldValue)
                    newValue = newValueAsString;
                else
                    error('BICAS:constants:Assertion:ConfigurationBug', 'Can not handle the MATLAB class=%s of internal setting "%s".', class(oldValue), key)
                end
            
                % Overwrite old setting.
                obj.set(key, newValue);
            end
            
        end
        
    end
    
    %###################################################################################################################
    
    methods(Access=private)
        
        function set(obj, key, value)
            obj.map(key) = value;
        end
        
    end
    
end

