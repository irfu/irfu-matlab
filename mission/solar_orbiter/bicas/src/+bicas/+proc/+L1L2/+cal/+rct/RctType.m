%
% Abstract class of which instances of subclasses represent one RCT TYPE (*not*
% the data stored in an RCT). Subclasses should therefore only need to be
% instantiated once, in principle, and should NOT contain any actual RCT data.
%
% NOTE: BICAS may load multiple RCTs for the same RCT type.
%
% IMPLEMENTATION NOTES
% ====================
% * Subclasses effectively collect code associated with each RCT type so that
%   RCTs can processed without knowing the type (to some extent).
% * Class can not contain map to singleton objects of subclasses since MATLAB
%   prevents that due recursive definitions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) RctType
% PROPOSAL: Singleton classes for each type of RCT (not storing the data).  -- ALMOST IMPLEMENTED
%           Share a common superclass. Can contain:
%   (1) methods read_*_RCT()
%   (2) methods modify_*_data()
%   (3) methods log_*_RCTs()
%   (4) filenameRegexpSettingKey (constant)
%   (5) RCTID (constant)
%   PRO: Can replace bulk of
%        bicas.proc.L1L2.cal.rct.fs,
%        bicas.proc.L1L2.cal.rct.typeproc.
%   PRO: Can replace structs created by
%        bicas.proc.L1L2.cal.rct.typeproc.init_RCT_TYPES_MAP.entry().
%   NEED: There is a need to be able to load the "plain" content of RCTs,
%         without modification.
%       Ex: Best practice w.r.t. modularization. Reusability.
%       Ex: Loading and analyzing RCTs outside of BICAS.
%       PROPOSAL: Can still have separate methods for (1) reading raw RCT, and
%                 (2) modifying the raw RCT data.
%
% PROPOSAL: Field for RCTID.
%
% PROPOSAL: Create singleton classes of subclasses and store a map with RCTID
%           keys to them.
%   PROPOSAL: Assert that filenameRegexpSettingKey values are unique.
%
% PROPOSAL: Classes for RCT data (not RCT type).
%   PRO: BIAS data has many fields.
%   PRO: More well-defined data structs.
%   PRO: Automatic assertions.
%   CON: Structs are modified when cal.m uses them, i.e. one could just as well
%        have classes for the format cal.m uses. ==> Too many classes.
%   PROPOSAL: Convert subclasses to stores of RCT data too.
%       CON: Can have multiple non-BIAS RCTs loaded. Multiple instances of same
%            RCT type has no meaning.
%
% PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Assert CDF skeleton/master version number.
% PROPOSAL: Assert skeleton/master.
%   PRO: Can give better error message when reading the wrong RCT.
%
% PROPOSAL: Assert/warn (depending on setting?) when CDF metadata imply that the RCT zVariables have the wrong units.
% PROPOSAL: Use utility function for reading every zVariable.
%   PROPOSAL: Assert units from zVar attributes.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(Abstract, Constant, GetAccess=public)
        % Settings key for value that defines the regular expression that is
        % used for finding the corresponding RCT(s).
        filenameRegexpSettingKey
    end
    properties(GetAccess=public, Constant)

        % Minimum number of expected entries in tabulated transfer functions in
        % RCTs.
        TF_TABLE_MIN_LENGTH = 10;

        % LL = Log Level
        RCT_DATA_LL = 'debug';
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Abstract)



        % Read RCT file.
        %
        %
        % DESIGN INTENT
        % =============
        % Should be implemented so that no calibration data is modified/added
        % to/removed from. The returned data structures reflect the content of
        % the RCTs, not necessarily the data used by BICAS. Modification of data
        % should be done elsewhere, in particular modifications of transfer
        % functions, e.g. extrapolation, cut-offs, inversions.
        % --
        % NOTE: BIAS & LFR RCTs: contain FTFs which are not inverted in this code.
        %       TDS RCTs:        contain ITFs.
        % NOTE: Code still converts RCT TFs slightly:
        %   frequency      : Hz    --> rad/s
        %   phase+amplitude: degrees,dimensionless real value --> Z (complex number)
        %
        [RctData] = read_RCT(filePath);

        % Modify the data structure read by bicas.proc.L1L2.cal.rct.read_RCT()
        % to a data structure that BICAS can use.
        %
        % IMPLEMENTATION NOTE: There is a need to distinguish between (1) the
        % data structures in RCT files, which one may want to inspect manually,
        % or log, and should be quite analogous to the RCT file content, and (2)
        % the calibration data data structures which are convenient for BICAS to
        % use.
        [RctData] = modify_RCT_data(RctData);

        % Custom logging of modified RCT data.
        log_RCT(RctData, L);



    end    % methods(Static)



end
