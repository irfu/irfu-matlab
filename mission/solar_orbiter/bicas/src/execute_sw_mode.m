% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
% input_files : containers.Map with
%    keys   = dataset IDs
%    values = paths to input files.
%
function execute_sw_mode(i_sw_mode, input_files, output_dir, sw_root_path)
%
% QUESTION: How verify dataset ID and dataset version against constants?
%    NOTE: Need to read cdf first.
%    NOTE: Need S/W mode.
%
% PROPOSAL: In principle, if using the MMS data manager concept, then the code does not need
%       to explicitly check for input files, only ask for a set of datasets.
%    PROBLEM: How check input dataset version? Need to look up S/W mode and constants.
%       Could in principle be different for different modes(?!)
%       QUESTION: Check where? In data_manager.set_input? In execute_sw_mode.
%          data_manager does logically need some knowledge of dataset versions since its algorithms
%          are tailored for different dataset versions.
%       --
%       PROPOSAL: Set input data not one cdf at a time, but one S/W mode at a time.
%
%
%
% QUESTION: What should be the relationship between data manager and S/W modes really?
% Should data manager check anything?
%

global ERROR_CODES
C = bicas_constants;

% QUESTION: Move this check to execute_sw_mode/bicas?!
irf.log('n', sprintf('Output directory = "%s"', output_dir));       
if ~exist(output_dir, 'dir')
    errorp(ERROR_CODES.PATH_NOT_FOUND, 'Output directory "%s" does not exist.', output_dir)
end



dm = data_manager();
for dataset_ID = input_files.keys
    input_file = input_files(dataset_ID{1});
    dm.set_input_cdf(dataset_ID{1}, input_file);
end

end
