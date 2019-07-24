%
% Takes a structure "a" and adds to it the fields of a
% structure "b". If a field exists in both structures,
% then DuplicateFieldBehaviours determines what happens.
% Can be recursive.
%
% NOTE: The operation is NOT symmetric in "a" and "b".
%
%
% ARGUMENTS
% =========
% MANDATORY:
%   a                          : A struct. The fields of "b" will be added to this struct.
%   b                          : A struct.
% OPTIONAL:
%   DuplicateFieldBehaviours : Struct array with four components, each one representing how to react for different cases of
%                              duplicate fields, depending on which of the two fields are structs themselves.
%                              Each field in this argument is a string with one of the following values:
%                                   "Overwrite"  : The field in a is overwritten by the field in b.
%                                   "Do nothing"
%                                   "Recurse"    : Only valid for two structs.
%                                   "Error"      : Default value.
%       .noStructs
%       .aIsStruct
%       .bIsStruct
%       .bothAreStructs
%
%
% Initially created 201x-xx-xx (=<2017) by Erik P G Johansson.
%
function a = add_struct_to_struct(a, b, DuplicateFieldBehaviours)
% PROPOSAL: Function pointer for what to do for duplicate fields?!
% TODO-DECISION: Over-engineered?
%   PRO: Too many duplicate field behaviours.
%   PRO?: Distinction aIsStruct/bIsStruct is unnecessary.
% PROPOSAL: Option for what to do if two fields are present but have different (non-struct) values.

% TODO-DECISION: How determine recursion depth?
%   PROPOSAL: Function pointer argument?!!
%   PROPOSAL: (Flag) always/never recurse structs (if present in both a and b).
%
% PROPOSAL: Use EJ_library.utils.interpret_settings_args.
%   CON: Can not use it since it calls this function! Would get infinite recursion!!!


% NOTE: In reality designed for merging "compatible" data structures (types are compatible; field names may or may not overlap).
% ==> Does not need so many cases.
% ==> Ambiguous if a struct within a struct of settings is to be regarded as a value or as another container for settings (not clear where to end recursion).
% ==> Ambiguous whether to copy/overwrite or not. ==> Operation can not be automatically selected by algorithm.
%
% PROPOSAL: Want to be able to think of structs (recursive) as sets, with or without overlap (intersection).
%   TODO-DECISION: How handle "unequal" intersections (same-named fields).
%       Ex: Different field value types (struct, double, string, cell, object, ...).
%       Ex: Different field values (of same type).
%   TODO-DECISION: How determine recursion depth?
%   PROPOSAL: Comparison with diffdn somehow?
%       CON: Recursion depth is well defined for diffdn.
%   


    % DFB = Duplicate field behaviour.
    DEFAULT_DFB = struct('noStructs', 'Error', 'aIsStruct', 'Error', 'bIsStruct', 'Error', 'twoStructs', 'Error');   % DFB = Duplicate Field Behaviours
    
    
    
    if nargin == 2
        DuplicateFieldBehaviours = DEFAULT_DFB;
    end
    
    % ASSERTIONS
    if ~ismember(nargin, [2,3])
        error('Wrong number of arguments.')
    elseif ~isstruct(a)
        a
        error('a is not a structure.')
    elseif ~isstruct(b)
        b
        error('b is not a structure.')
    elseif ~isstruct(DuplicateFieldBehaviours)
        DuplicateFieldBehaviours
        error('DuplicateFieldBehaviours is not a struct.')
    end

        
    %===========================
    % Iterate over fields in b.
    %===========================
    bFieldNamesList = fieldnames(b);
    for i = 1:length(bFieldNamesList)
        fieldName = bFieldNamesList{i};

        if isfield(a, fieldName)
            %===========================================================
            % CASE: Duplicate fields (field exists in both structures).
            %===========================================================
            afv = a.(fieldName);
            bfv = b.(fieldName);

            %nFieldStructs = isstruct(a.(fieldName)) + isstruct(b.(fieldName));   % Of the duplicate fields, determine how many are structures.
            % 0 = two non-structs
            % 1 = struct + non-struct
            % 2 = struct + struct
            
%             switch nFieldStructs
%                 case 0
%                     behaviour = DuplicateFieldBehaviours.noStructs;
%                 case 1
%                     behaviour = DuplicateFieldBehaviours.oneStruct;
%                 case 2
%                     behaviour = DuplicateFieldBehaviours.twoStructs;
%                 otherwise
%                     error('Should never reach this code.')
%             end

            bothAreStructs = false;
            if ~isstruct(afv) && ~isstruct(bfv)
                behaviour = DuplicateFieldBehaviours.noStructs;
            elseif isstruct(afv) && ~isstruct(bfv)
                behaviour = DuplicateFieldBehaviours.aIsStruct;
            elseif ~isstruct(afv) && isstruct(bfv)
                behaviour = DuplicateFieldBehaviours.bIsStruct;
            else
                behaviour = DuplicateFieldBehaviours.bothAreStructs;
                bothAreStructs = true;
            end

            if strcmp(behaviour, 'Error')
                error('Structures share identically named fields "%s".', fieldName)
            elseif strcmp(behaviour, 'Overwrite')
                a.(fieldName) = b.(fieldName);
            elseif strcmp(behaviour, 'Do nothing')
                % Do nothing.
            elseif strcmp(behaviour, 'Recurse') && (bothAreStructs)
                % NOTE: Recursive call. Needs the original DuplicateFieldBehaviours.
                a.(fieldName) = EJ_library.utils.add_struct_to_struct(a.(fieldName), b.(fieldName), DuplicateFieldBehaviours);
            else
                error('Can not interpret string value behaviour="%s" for this combination of field values.', behaviour)
            end
        else
            %===========================================
            % CASE: Field exists in "b", but not in "a"
            %===========================================
            a.(fieldName) = b.(fieldName);     % NOTE: Not overwrite field, but create new field in a.
        end
    end
end
