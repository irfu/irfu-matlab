
%
% EXPERIMENTAL
%
% Find relevant CDFs and search for specific conditions that indicate that
% something is wrong with the processing inside the RPWI (not the GS pipeline).
% Print/log warning for specific conditions.
% 
%
% NOTES
% =====
% * Code is partly a proof-of-concept. It tests whether it makes sense to use
%   the GS pipeline CDFs to automatically look for common well-defined states
%   that indicate warnings or errors in the instrument itself (not problems
%   related to the TM as such), e.g. for end2end tests. It is not obvious how
%   many conditions it should search for or how sophisticated the code should be,
%   e.g. w.r.t. logging.
% * Code can potentially be amended with more tests.
% * Should maybe not be implemented in MATLAB, but in Python.
% * Since the code relies on
%
%
% ARGUMENTS
% =========
% dirPath
%       Path to root directory. Code will search all subdirectories for relevant
%       JUICE/RPWI GS TM-to-L1a CDF files assuming filenaming conventions.
% reportErrorCounterNonzero
%       Scalar logical. Whether to report to when error counters are not zero.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function validate_instrument_processing(dirPath, reportErrorCounterNonzero)
    % PROPOSAL: Better name.
    %   TM, L1a, CDFs, datasets, data
    %   validate, check
    %   instrument
    %   processing
    %       CON: Does not have to do with TM-to-L1a processing or GS.
    %   validate_instrument_processing
    %   
    % PROPOSAL: Implement in Python.
    %   PRO: Could load MIB with pre-existing code.
    %   CON: Harder to plot.
    %       CON-PROPOSAL: Once a problem has been identified, can then just as
    %                     well plot outside the script, i.e. in MATLAB.
    %       CON-PROPOSAL: (Should) use RPWIQuicklooks.
    %
    % PROPOSAL: Print timestamps for when error counters increment.
    %   CON: Makes code too complicated. Detracts from core purpose.
    % PROPOSAL: Add plots somehow for determining when error counters increment.
    %   How plot time axis?!
    %   CON: Makes code too complicated. Detracts from core purpose.
    %
    % PROPOSAL: Log to file.
    %
    % NOTE: Has no access to MIB data here, as opposed to in the GS pipeline.
    %       Can therefore not look up e.g. packet type descriptions.

    assert(ischar(dirPath))
    
    validate_error_packets(dirPath)
    validate_error_counters(dirPath, reportErrorCounterNonzero)
end



% Search for 5.2, 5.3, and 5.4 packets.
function validate_error_packets(dirPath)
    PPD_FILENAME_RE       = 'JUICE_LU_RPWI-PPD_.*.bin.cdf';
    SERVICE_TYPE          = 5;
    SERVICE_SUBTYPE_ARRAY = [2,3,4];
    
    DirInfoArray = find_files(dirPath, PPD_FILENAME_RE);
    
    for i = 1:numel(DirInfoArray)
        cdfPath = fullfile(DirInfoArray(i).folder, DirInfoArray(i).name);
        do      = dataobj(cdfPath);
        
        spidCa              = cellstr(read_DO_field(do, {'SPID'}));
        serviceTypeArray    =         read_DO_field(do, {'DFH_SERVICE_TYPE'});
        serviceSubtypeArray =         read_DO_field(do, {'DFH_SERVICE_SUBTYPE'});
        
        fprintf('cdfFile = "%s"\n', cdfPath)
                
        b = (serviceTypeArray == SERVICE_TYPE) & ismember(serviceSubtypeArray, SERVICE_SUBTYPE_ARRAY);
        iEPackets = find(b);
        
        if ~isempty(iEPackets)
            nErrorPackets = numel(iEPackets);
            
            [~, jUniqueEPacketsArray] = unique(spidCa(iEPackets));
            iUniqueEPacketsArray = iEPackets(jUniqueEPacketsArray);
            
            log_warning(0, 'Found %i error packet(s) (5.2, 5.3, 5.4). Unique packet types:', ...
                nErrorPackets)
            for iUniqueEPacket = iUniqueEPacketsArray(:)'
                st = serviceTypeArray(   iUniqueEPacket);
                ss = serviceSubtypeArray(iUniqueEPacket);
                spid = spidCa{           iUniqueEPacket};
                log_warning(1, '%i.%i: SPID=%s', st, ss, spid)
            end
        end
    end
    
end



function validate_error_counters(dirPath, reportErrorCounterNonzero)
    PPTD_HK64_FILENAME_RE = 'JUICE_LU_RPWI-PPTD-LWYHK[01]0064_.*\.cdf';
    ERROR_CORE0_MPN_CA    = {'LWT0343D', 'LWT0456B'};
    ERROR_CORE1_MPN_CA    = {'LWT0343E', 'LWT0456C'};
    
    DirInfoArray = find_files(dirPath, PPTD_HK64_FILENAME_RE);

    for i = 1:numel(DirInfoArray)
        cdfPath = fullfile(DirInfoArray(i).folder, DirInfoArray(i).name);
        do = dataobj(cdfPath);
        
        fprintf('cdfFile = "%s"\n', cdfPath)
        
        errorCounter0Array = read_DO_field(do, ERROR_CORE0_MPN_CA);
        errorCounter1Array = read_DO_field(do, ERROR_CORE1_MPN_CA);
        
        validate_error_counter_array(0, errorCounter0Array, reportErrorCounterNonzero)
        validate_error_counter_array(1, errorCounter1Array, reportErrorCounterNonzero)
    end
end



% Validate one array of error counter values.
function validate_error_counter_array(iCounter, counterArray, reportErrorCounterNonzero)
    assert(isscalar(reportErrorCounterNonzero) & islogical(reportErrorCounterNonzero))
    
    strUniqueValues = sprintf('Unique values: %s', num2str(unique(counterArray(:)')));

    n = numel(unique(counterArray));
    if ~ismember(n, [0, 1])
        log_warning(0, 'Core%i error counter values are NOT CONSTANT. %s', ...
            iCounter, strUniqueValues)
    end
    
    if reportErrorCounterNonzero
        if ~all(counterArray == 0)
            log_warning(0, 'Core%i error counter values are not equal to zero. %s', ...
                iCounter, strUniqueValues)
        end
    end
end



function DirInfoArray = find_files(rootDirPath, rePattern)
    DirInfoArray = dir(fullfile(rootDirPath, '**/*.cdf'));
    
    bIsFile = ~[DirInfoArray.isdir];
    
    ca = regexp({DirInfoArray.name}, rePattern);
    bFilenameMatch = ~cellfun(@isempty, ca);
    
    DirInfoArray = DirInfoArray(bIsFile & bFilenameMatch);
end



% Extract zVariable from dataobj without knowing the exact field name.
%
% zvNameCandidatesCa
%       1D CA of ZV names. Exactly one should be valid.
function value = read_DO_field(do, zvNameCandidatesCa)
    assert(isa(do, 'dataobj'))
    
    bArray = ismember(zvNameCandidatesCa, fieldnames(do.data));
    assert(sum(bArray) == 1, 'Can not find exactly one matching field.')
    
    fieldName = zvNameCandidatesCa{bArray};
    value = do.data.(fieldName).data;
end



% Print a one-row warning to stdout.
function log_warning(iLevel, s, varargin)
    assert(iLevel >= 0)
    
    if iLevel == 0
        s = ['WARNING: ', s];
    end
    
    s = [repmat('    ', 1, iLevel+1), s, '\n'];
    fprintf(s, varargin{:})
end