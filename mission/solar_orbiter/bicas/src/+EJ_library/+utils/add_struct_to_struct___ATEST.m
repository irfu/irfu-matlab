function add_struct_to_struct___ATEST
    
    AB  = struct('a', 1, 'b', 2);
    CD  = struct('c', 3, 'd', 4);
    BC  = struct('b', 8, 'c', 9);
    ES  = struct();
    
    % E = Error
    % O = Overwrite
    % N = Do nothing
    SET_EEEE = struct('noStructs', 'Error',      'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Error');
    SET_OEEE = struct('noStructs', 'Overwrite',  'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Error');
    SET_NEEE = struct('noStructs', 'Do nothing', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Error');
    SET_OEEO = struct('noStructs', 'Overwrite',  'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Overwrite');
    SET_OEER = struct('noStructs', 'Overwrite',  'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Recurse');
    SET_EEER = struct('noStructs', 'Error',      'aIsStruct', 'Error', 'bIsStruct', 'Error', 'abAreStructs', 'Recurse');
    
    NewTest = @(inputs, expOutputsException)  (EJ_library.atest.CompareFuncResult(@EJ_library.utils.add_struct_to_struct, ...
        inputs, expOutputsException));
    
    tl = {};



    % No overlap
    tl{end+1} = NewTest({ES, CD}, {CD});
    tl{end+1} = NewTest({AB, CD, SET_OEEE}, {struct('a', 1, 'b', 2, 'c', 3, 'd', 4)});
    tl{end+1} = NewTest({AB, CD},           {struct('a', 1, 'b', 2, 'c', 3, 'd', 4)});    % Test default DFP.
    
    % Test struct arrays
    A   = EJ_library.utils.empty_struct([9,0], 'a', 'b');
    B   = EJ_library.utils.empty_struct([9,0],      'b', 'c');
    tl{end+1} = NewTest({A, B, SET_OEEE}, 'MException');
    
    % Overlap
    tl{end+1} = NewTest({AB, BC, SET_OEEE}, {struct('a', 1, 'b', 8, 'c', 9)});
    tl{end+1} = NewTest({AB, BC, SET_NEEE}, {struct('a', 1, 'b', 2, 'c', 9)});
    tl{end+1} = NewTest({AB, BC, SET_EEEE}, 'MException');
    
    % Overwrite overlap (note recursive).
    A = struct('a', 1.1, 'b', 2,         'S1', struct('d', 4.1, 'e', 5        ));
    B = struct('a', 1.2,         'c', 3, 'S1', struct('d', 4.2,         'f', 6));
    C = struct('a', 1.2, 'b', 2, 'c', 3, 'S1', struct('d', 4.2,         'f', 6));
    tl{end+1} = NewTest({A, B, SET_OEEO}, {C});    
    
    % Overwrite overlaps (recursive).
    A = struct('a', 1.1, 'b', 2,         'S1', struct('d', 4.1, 'e', 5        ));
    B = struct('a', 1.2,         'c', 3, 'S1', struct('d', 4.2,         'f', 6));
    C = struct('a', 1.2, 'b', 2, 'c', 3, 'S1', struct('d', 4.2, 'e', 5, 'f', 6));
    tl{end+1} = NewTest({A, B, SET_OEER}, {C});
    
    % No overlap
    A = struct('b', 2,         'S1', struct('e', 5        ));
    B = struct(        'c', 3, 'S1', struct(        'f', 6));
    C = struct('b', 2, 'c', 3, 'S1', struct('e', 5, 'f', 6));
    tl{end+1} = NewTest({A, B, SET_EEER}, {C});
    
    %tl{end+1} = NewTest({}, {});
    %tl{end+1} = NewTest({}, {});
    
    
    
    EJ_library.atest.run_tests(tl);
end
