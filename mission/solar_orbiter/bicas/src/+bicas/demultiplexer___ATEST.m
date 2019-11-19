% Automatic test code for bicas.demultiplexer.
%
% Very basic tests at this stage. Could be improved but unsure how much is meaningful.
function demultiplexer___ATEST

new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.demultiplexer.main, inputs, outputs));
tl = {};

V1   = 10;
V2   = 11;
V3   = 12;
V12  = V1-V2;
V13  = V1-V3;
V23  = V2-V3;
V12a = 45-56;
V13a = 45-69;
V23a = 56-69;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function AsrSamplesVolt = ASR_samples(varargin)
        assert(nargin == 9)
        AsrSamplesVolt = struct(...
            'V1_LF',     as(varargin{1}, V1), ...
            'V2_LF',     as(varargin{2}, V2), ...
            'V3_LF',     as(varargin{3}, V3), ...
            'V12_LF',    as(varargin{4}, V12), ...
            'V13_LF',    as(varargin{5}, V13), ...
            'V23_LF',    as(varargin{6}, V23), ...
            'V12_LF_AC', as(varargin{7}, V12a), ...
            'V13_LF_AC', as(varargin{8}, V13a), ...
            'V23_LF_AC', as(varargin{9}, V23a));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function BltsSrcArray = BLTS_src_array(categoryArray, antennasArray)
        assert( numel(categoryArray) == numel(antennasArray) )
        
        for i =1:numel(categoryArray)
            BltsSrcArray(i) = bicas.BLTS_src_dest(...
                categoryArray{i}, ...
                antennasArray{i});
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if 1
    tl{end+1} = new_test(...
        {0, true, {V1, V12, V23, V12a, V23a}}, ...
        {BLTS_src_array(...
            {'DC single', 'DC diff', 'DC diff', 'AC diff', 'AC diff'}, ...
            {[1], [1 2], [2 3], [1 2], [2 3]}), ...
        ASR_samples(1,1,1, 1,1,1, 1,1,1)});
end

if 1
    tl{end+1} = new_test(...
        {1, false, {V2, V3, V23, V13a, V23a}}, ...
        {BLTS_src_array(...
            {'DC single', 'DC single', 'DC diff', 'AC diff', 'AC diff'}, ...
            {[2], [3], [2 3], [1 3], [2 3]}), ...
        ASR_samples(0,1,1, 0,0,1, 1,1,1)});
end

EJ_library.atest.run_tests(tl)
end



 % as = assign. Effectively implements ~ternary operator + constant (NaN).
function V = as(v,V)
if v; V = V;
else  V = NaN;
end
end
