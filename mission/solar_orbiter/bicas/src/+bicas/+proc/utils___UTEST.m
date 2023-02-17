%
% matlab.unittest automatic test code for bicas.proc.utils().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08 from older test code.
%
classdef utils___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_set_struct_field_rows(testCase)

            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = bicas.proc.utils.set_struct_field_rows(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            %===================================================================
            % VERY INCOMPLETE TEST SUITE.

            test({...
                struct('asd', reshape([1:24], [4,3,2])), ...
                struct('asd', reshape([1:12], [2,3,2])), 2:3}, {...
                struct('asd', reshape([1,1,2,4, 5,3,4,8, 9,5,6,12,   13,7,8,16, 17,9,10,20, 21,11,12,24], 4,3,2)) ...
                });

            test({...
                struct('asd', [1;2;3;4;5]), ...
                struct('asd', [8;9]), [4,3]}, {...
                struct('asd', [1;2;9;8;5])});
        end



        function test_convert_matrix_to_cell_array_of_vectors(testCase)

            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            %===================================================================
            test({zeros(0,1), zeros(0,1)}, {cell(0,1)});
            test({[1,2,3,4,5], [3]}, {{[1,2,3]}});
            test({[1,2,3,4,5; 6,7,8,9,0], [3; 2]}, {{[1,2,3]; [6,7]}});
        end



        function test_convert_cell_array_of_vectors_to_matrix(testCase)

            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            %===================================================================
            test({zeros(0,1), zeros(0,1)}, {cell(0,1)});
            test({[1,2,3,4,5], [3]}, {{[1,2,3]}});
            test({[1,2,3,4,5; 6,7,8,9,0], [3; 2]}, {{[1,2,3]; [6,7]}});

        end



        function test_set_NaN_after_snapshots_end(testCase)

            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = bicas.proc.utils.set_NaN_after_snapshots_end(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            %===================================================================
            test({ones(0,4),              ones(0,1)},   {ones(0,4)});
            test({[0,1,2],                [3]},   {[0,1,2]});
            test({[0,1,2,3,4; 5,6,7,8,9], [2;4]}, {[0,1,NaN,NaN,NaN; 5,6,7,8,NaN]});
        end
        
        
        
    end    % methods(Test)



end
