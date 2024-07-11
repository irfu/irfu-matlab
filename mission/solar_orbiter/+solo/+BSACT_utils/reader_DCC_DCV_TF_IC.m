%
% Class that collects calibration data and the corresponding metadata from
% BSACT. An object can and should contain either DCC, DCV, _or_ TF data.
% Metadata are collected by parsing the testlogbook*.txt files.
%
%
% NOTES
% =====
% IMPLEMENTATION NOTE: The current implementation really "only" collects a list
% of metadata and paths to the corresponding calibration table files. When this
% metadata has been retrieved, the object is no longer needed. The class could
% therefor in principle be replaced by a function. However, the thought behind
% the current design is to make it possible to easily modify it in the future so
% that data can be loaded/retrieved through the class. This design may be
% changed in the future though.
% NOTE: The implementation does not force one to use either DCC, DCV, TF, or IC
% test in a given object. One can mix them, but that form uf usage makes little
% sense and is not intended.
% NOTE: There is a slight difference between BSACT TFs (without inverted diffs).
%   2016 June: Has TF phase shift ~180 degrees at low frequencies.
%   2016 July: Has TF phase shift   ~0 degrees at low frequencies.
% NOTE: The class name is chosen to reflect the types of calibration data that
% it may contain in anticipation of eventually creating another analogous class
% for other calibration data.
% NOTE: For DCC log files: Does not read the BIAS output channel (LFR_1/2/3)
% since requires reading "Header1" rows. Should be the same as the antenna
% channels (which are read).
%
%
% RATIONALE
% =========
% Loading all calibration files rather than just their fits, or the relevant
% subsets could be useful:
% (1) Can be used for automatically deriving smaller sets of calibration data to
%     actually be used by BICAS. This includes deriving more "refined" (more
%     granular), smaller sets of calibration data to actually be used by BICAS:
%       Example 1: Select AC (spectrum) transfer functions (TF) based on
%       combination of "channel" and high gain/low gain, rather than just one
%       for AC low gain and AC high gain as expected.
%       Example 2: Interpolate table instead of linear function (offset & slope;
%       y=k*x+m).
%       Example 3: Derive fits (DC) from calibration tables directly rather than
%       using separately tabulated offsets and slopes.
% (2) Can be used for automatically converting calibration data to to other file
%     formats (CDF, XML etc).
% (3) Manually analyzing the calibration data itself (plotting; making one's own
%     fits using other table columns; comparing channels, checking temperature
%     dependence).
%
%
% VARIABLE NAMING CONVENTIONS / DEFINITIONS OF TERMS
% ==================================================
% BSACT             BIAS StandAlone Calibration Tables. The raw files generated
%                   at the BIAS standalone calibration, in
%                   particular those taken (1) June 2016 and (2) July 2016.
%
% CTable            Calibration table. Equal to the contents of one BIAS
%                   standalone calibration file.
%
% metadata          Refers to the metadata for a specific calibration
%                   table/CTable, e.g. MUX mode, input channels,
%                   AC low/high gain, diff/single, stimuli.
%
% Input channel     One or two integers (row vector; values 1-3), signifying a
% (inputChNbrs)     single-probe "channel" or a diff-probe "channel", i.e. which
%                   signal is coming from the ANTENNA(S).
%                       Example: [2] refers to V2_LF.
%                       Example: [2,3] refers to (V2_LF-V3_LF).
%                   NOTE: Does NOT refer to the signal itself (voltage), only
%                   which channel out of several.
%                   NOTE: The probe order is always increasing.
%
% Output channel    Integer (1-5) representing one of the "channels" at
% (outputChNbr)     the BIAS-to-LFR/TDS boundary., Also known as LFR_1...LFR_5,
%                   TDS_1...TDS_3, and BIAS_1...BIAS_5.
%                   NOTE: Does NOT refer to the signal itself (voltage), only
%                   which channel out of several.
%
% Test ID:          Integer defined in the standalone calibration data (raw
% (testIdNbr)       format). For a specific temperature and type of calibration
%                   data, it identifies a specific calibration file table with
%                   specific input/output channels (or TC-bias and bias
%                   currents), stimuli, MUX mode (+latching relay setting). Test
%                   IDs are identical for the same type of calibration data at
%                   different temperatures.
%
% RPS               Radians per second
% --
% EOO               Excluding output offset (offset removed)
% IOO               Including output offset (offset not removed)
% BST               Before stimuli
% AST               After stimuli (after voltage drop due to stimuli)
%
%
% Types of calibration data
% -------------------------
% DCC (DC_CURRENT)            For determining offset and slope for bias currents.
% DCV (DC_VOLTAGE)            For determining offset and slope for DC voltages.
% TF  (TRANSFER_FUNCTION)     Transfer functions (amplification and phase shift
%                             as a function of frequency).
% IC  (BIAS_DC_INTERNAL_CAL)  For determining the internal calibration resistance.
%
%
% Direction, what is input/output
% -------------------------------
% The naming convention here is to use the flow of information in the physical
% instrument to determine what is "input" and "output". Therefor, for voltages,
% BICAS calculates the "input" signals from the "output" signals.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-12-11
%
classdef reader_DCC_DCV_TF_IC < handle
  % BOGIQ
  % =====
  % PROPOSAL: Force user to only use one data type: DCC, DCV, or TF tests.
  %   PROPOSAL: Internal flag that selects which.
  %   PROPOSAL: Subclasses for each case.
  % PROPOSAL: Implement pre-loading of data.
  %   CON: Code does not know whether to load DCV or TF.
  %
  % PROPOSAL: Split into multiple classes, one per calibration data type.
  %   NOTE: Want to reuse "add_test_directory".
  %
  % PROPOSAL: Rename to exclude "DCC_DCV_TF_IC".
  %   PRO: Already includes all data types so mentioning them is superfluous.
  %
  % PROPOSAL: Replace class with functions for retrieving metadata. If one also
  %           wants functionality for loading all calibration data at once (and
  %           saving it memory), or caching, then that should be done by other
  %           function/class which calls those functions.
  %
  % PROPOSAL: Change shortening: EOO-->OOR (Output offset removed), IOO-->OOK (Output offset kept)
  %
  % PROPOSAL: Read temperature which seems contained in all (?) calibration table files.
  %   Ex: "header.reg2 23.96   :Ambient temperature in C"
  %
  % PROPOSAL: Functions for reading files in separate file (or static functions in separate class).
  % PROPOSAL: Function for recognizing and parsing data/CTable filenames.
  % Eliminate cTableFilesPattern, testLogbookFile.
  %   NOTE: Compare solo.adm.dsfn.parse_dataset_filename().
  %   PROBLEM: May not be able to eliminate testLogbookFile since there are
  %            multiple logbooks and no known algorithm for selecting the correct
  %            one.
  %
  % PROPOSAL: Replace term "Output channel" with "BLTS", iBlts as in BICAS.
  %
  % PROPOSAL: read_TF_calib_file() should return struct(s) without redundancy:
  %   Either "human units", "mathematical units", or two separate structs with
  %   each.



  properties(Access=private)
    %doPreload    = [];
    metadataList = [];
  end



  %###########################################################################



  methods(Access=public)

    % CONSTRUCTOR
    function obj = reader_DCC_DCV_TF_IC(varargin)
      % NOTE 2020-11-16: Old code that uses BSACT calls this constructor
      % with a path argument.
      %             if nargin > 1
      %                 error('BICAS:reader_DCC_DCV_TF_IC:IllegalArgument', 'Illegal number of arguments')
      %             end
      %
      %             obj.doPreload = ismember('preload', varargin);
    end



    function metadataList = add_test_directory(obj, ...
        dataType, cTableFilesPattern, testLogbookFile, mebTemperatureCelsius)

      % Add directory with tests (calibration tables). Automatically parse the
      % test logbook and associate each calibration table file with the
      % corresponding metadata, mostly from the logbook.
      %
      % Example: TEMP25C/4_4_DC_VOLTAGE_TEST/ or 4-5_TRANSFER_FUNCTION/.
      %
      % ARGUMENTS
      % =========
      % cTableFilesPattern
      %       Path as a sprintf pattern describing all test files with %02i or
      %       %03i representing the test ID nbr. NOTE: Short test ID numbers
      %       have to be precedeed by zeroes (due to filenaming convention).
      %       Example: '4-5_TRANSFER_FUNCTION/SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FS0_PAFM.txt'
      % testLogbookFile
      %       Path to testlogbook* file that describes and enumerates the
      %       files referred to by cTableFilesPattern.
      % mebTemperatureCelsius
      %       The MEB temperature at which the tests are made.

      testLogbookRowList = irf.fs.read_text_file(...
        testLogbookFile, '\r?\n');
      metadataList = solo.BSACT_utils.parse_testlogbook_DCC_DCV_TF_IC(...
        testLogbookRowList, dataType);

      % TODO-NI: Necessary to use special function here? Can replace call with one-liner?
      metadataList = irf.ds.merge_structs(...
        metadataList, struct('mebTempCelsius', mebTemperatureCelsius));



      % IMPLEMENTATION NOTE: Does not read entire files since does not
      % know whether they are DCC, DCV, TF, or IC files.
      [metadataList.filePath] = deal([]);
      for i = 1:numel(metadataList)
        filePath = sprintf(cTableFilesPattern, metadataList(i).testIdNbr);

        % ASSERTION: Check that files exist.
        if ~exist(filePath, 'file')
          error('BICAS:read_text_file:Assertion', ...
            'Can not find file "%s".', filePath)
        end

        % NOTE: Adding field to metadataList.
        metadataList(i).filePath = filePath;
      end



      obj.metadataList = [obj.metadataList, metadataList];
    end



    function metadataList = get_metadataList(obj)
      % Return the internal list of metadata for every calibration table file
      % that has been registered by this instance of this class.
      metadataList = obj.metadataList;
    end

  end    % methods(Access=public)



  %###########################################################################



  methods(Static, Access=public)

    function Data = read_DCC_calib_file(filePath)
      % Read one BSACT DCC table (text file) into a struct.

      assert(ischar(filePath), ...
        'BICAS:reader_DCC_DCV_TF_IC:IllegalArgument', ...
        'Argument is not a string.');

      Data = solo.BSACT_utils.read_BSACT_file(filePath, {...
        'currentMicroSAmpere', ...     % Set current/design current. Exactly proportional to digital current.
        'currentDigital', ...          % Signed integer value representing current. If not identical to TM,
        ...                            % then at least very-very similar (signed/unsigned?).
        'currentAAmpere', ...          % Measured bias current. NOTE: Opposite sign compared to set current.
        'inputVoltageVolt', ...        % Measured bias voltage.
        'controlVoltageVolt'});        % Measured voltage multiplied by factor, to approximate what is is sent to
      % LFR/TDS. Not meant to be used for actual calibration. Only a reality check.
    end



    function Data = read_DCV_calib_file(filePath)
      % Read one BSACT DCV table (text file) into a struct.
      %
      %
      % Column naming convention to distinguish columns with very similar
      % meanings (and values):
      % See comments at top of file.

      assert(ischar(filePath), ...
        'BICAS:reader_DCC_DCV_TF_IC:IllegalArgument', ...
        'Argument is not a string.')

      Data = solo.BSACT_utils.read_BSACT_file(filePath, {...
        'inputBstVolt', ...
        'outputIooVolt', ...
        'currentAAmpere', ...
        'inputAstVolt', ...
        'outputEooVolt', ...
        'hkInput1AstVolt', ...
        'hkInput2AstVolt', ...
        'hkInput3AstVolt'});
    end



    function [Data] = read_TF_calib_file(filePath)
      % Read one BSACT TF table (text file) into a struct.

      assert(ischar(filePath), ...
        'BICAS:reader_DCC_DCV_TF_IC:IllegalArgument', ...
        'Argument is not a string.')

      Data = solo.BSACT_utils.read_BSACT_file(filePath, {...
        'freqHz', ...
        'gainEnergyDb', ...
        'phaseShiftDeg'});

      % IMPLEMENTATION NOTE: Not returning information in same struct
      % since information redundancy can be dangerous.
      %             % Add fields to Data, fields which are likely to be used for
      %             % plotting.
      %             [Data.freqRps, Data.z] = irf.utils.convert_TF_human2math(...
      %                 Data.freqHz, ...
      %                 Data.gainEnergyDb, ...
      %                 Data.phaseShiftDeg);
    end



    function Data = read_IC_calib_file(filePath)
      % Read one BSACT IC table (text file) into a struct.

      assert(ischar(filePath), ...
        'BICAS:reader_DCC_DCV_TF_IC:IllegalArgument', ...
        'Argument is not a string.')

      Data = solo.BSACT_utils.read_BSACT_file(filePath, {...
        'currentMicroSAmpere', ...   % Set current/design current. Exactly proportional to digital current.
        'currentDigital', ...        % Signed integer value representing current. If not identical to TM,
        ...                          % then at least very-very similar (signed/unsigned?).
        'outputIooVolt', ...
        'outputEooVolt', ...
        'hkInput1Volt', ...          % Slightly unsure, but should be measured voltage in HK data, in Volt. Compare DCV columns.
        'hkInput2Volt', ...
        'hkInput3Volt'});

    end

  end    % methods(Static, Access=public)



  %###########################################################################
end
