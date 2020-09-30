%
% Automatic test code for bicas.demultiplexer.
%
% Very basic tests at this stage. Could be improved but unsure how much is
% meaningful.
%
function demultiplexer___ATEST
    main___ATEST
    complement_ASR___ATEST
end



function complement_ASR___ATEST
    newTest    = @(inputFieldsCa, outputFieldsCa) (erikpgjohansson.atest.CompareFuncResult(...
        @new_test, ...
        {inputFieldsCa}, {}));
    
    tl = {};
    
    % TODO: dlrUsing12
    tl{end+1} = newTest({'dcV1', 19, 'dcV12', 27, 'dcV23', 33,    'acV12', 54, 'acV23', 75});    % mux=0, dlrUsing12=1
    tl{end+1} = newTest({'dcV1', 19, 'dcV13', 27, 'dcV23', 33,    'acV13', 54, 'acV23', 75});    % mux=0, dlrUsing12=0
    tl{end+1} = newTest({'dcV2', 19, 'dcV3',  27, 'dcV23', 19-27, 'acV12', 54, 'acV23', 75});    % mux=1
    tl{end+1} = newTest({'dcV1', 2   'dcV2',  7,  'dcV3',  32,    'acV12', 74, 'acV23', 85});    % mux=4
    
    erikpgjohansson.atest.run_tests(tl)
    
    
    
    function new_test(inputFieldsCa)
        A = bicas.demultiplexer.complement_ASR( struct(inputFieldsCa{:}) );
        
        % Test all possible relationsships.
        %
        % NOTE: Implicitly asserts that all fields are present.
        % NOTE: Must account for that some fields may be NaN, and therefore can
        % not be checked against relations.
        assert_relation(A.dcV1,  A.dcV12, A.dcV2 )
        assert_relation(A.dcV1,  A.dcV13, A.dcV3 )
        assert_relation(A.dcV2,  A.dcV23, A.dcV3 )
        assert_relation(A.dcV13, A.dcV12, A.dcV23)    % DC. All diffs
        %
        assert_relation(A.acV13, A.acV12, A.acV23)    % AC. All diffs
    end
end



% Local utility function.
function assert_relation(A, B, C)
    b = ~isnan(A) & ~isnan(B) & ~isnan(C);
    assert(all(A(b) == B(b) + C(b)))
end



function main___ATEST()
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
    %###################################################################
    function AsrSamplesVolt = ASR_samples(varargin)
        assert(nargin == 9)
        AsrSamplesVolt = struct(...
            'dcV1',  as(varargin{1}, V1), ...
            'dcV2',  as(varargin{2}, V2), ...
            'dcV3',  as(varargin{3}, V3), ...
            'dcV12', as(varargin{4}, V12), ...
            'dcV13', as(varargin{5}, V13), ...
            'dcV23', as(varargin{6}, V23), ...
            'acV12', as(varargin{7}, V12a), ...
            'acV13', as(varargin{8}, V13a), ...
            'acV23', as(varargin{9}, V23a));
    end
    %###################################################################
    function BltsSrcArray = BLTS_src_array(categoryArray, antennasArray)
        assert( numel(categoryArray) == numel(antennasArray) )
        
        for i =1:numel(categoryArray)
            BltsSrcArray(i) = bicas.BLTS_src_dest(...
                categoryArray{i}, ...
                antennasArray{i});
        end
    end
    %###################################################################
    
    
    
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


% Local utility function.
% as = assign. Effectively implements ~ternary operator + constant (NaN).
function V = as(b,V)
    if b; V = V;
    else  V = NaN;
    end
end
