%
% Takes a structure "A" and adds to it the fields of a structure "B". If a field
% exists in both structures, then Settings determines what happens. Can be
% recursive.
%
%
% NOTES
% =====
% NOTE: The operation is NOT symmetric in "A" and "B".
% NOTE: The algorithm can ALMOST be generalized to struct arrays. The conceptual
%   cruxes are:
%   (1) a non-empty struct array field can be struct or non-struct for different
%       components,
%   (2) an empty struct array field does not have a class (neither struct
%       or non-struct)
%   A user could iterate over struct array components and call the function for
%   each component separately.
% --
% NOTE: irf.ds.merge_structs is similar.
%
%
% ARGUMENTS
% =========
% A         : Struct array. The fields of "B" will be added to this struct.
% B         : Struct array. Same size as A.
%             NOTE: Struct arrays only work if algorithm does not need to read
%             field values, ie. no collisions.
% varargin  : As interpreted by irf.utils.interpret_settings_args.
%       noStructs,
%       aIsStruct,
%       bIsStruct,
%       abAreStructs :
%           Above settings are string constants. They represent how to react for
%           different cases of duplicate fields, depending on which of the two
%           fields are structs themselves. Each field in this argument is a
%           string with one of the following values:
%              "Overwrite"  : The field in "A" is overwritten by the field in "B".
%              "Do nothing" :
%              "Recurse"    : Only valid if both fields are structs. Apply the
%                             function on the two fields.
%              "Error"      : Default value.
%       .onlyBField : String constant. 'Copy', 'Error'. What to do if field is
%                     present in B, but not A.
%
%
% RETURN VALUE
% ============
% A : Struct. Corresponds to argument "A", with the fields of "B" added to it.
%
%
% Initially created 201x-xx-xx (=<2017) by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function A = add_struct_to_struct(A, B, varargin)
%
% PROPOSAL: Function pointer for what to do for duplicate fields?!
% TODO-DEC: Over-engineered?
%   PRO: Too many duplicate field behaviours.
%   PRO?: Distinction aIsStruct/bIsStruct is unnecessary.
% PROPOSAL: Option for what to do if two fields are present but have different (non-struct) values.

% TODO-DEC: How determine recursion depth?
%   PROPOSAL: Function pointer argument?!!
%   PROPOSAL: (Flag) always/never recurse structs (if present in both "A" and "B").
%
% PROPOSAL: Use irf.utils.interpret_settings_args.
%   CON: Can not use it since it calls this function! Would get infinite recursion!!!
%
% PROPOSAL: Policy alternative for requiring a field in "A" to overwrite.
%   Ex: Useful for irf.utils.interpret_settings_args.
%   Ex: BICAS: write_dataset_CDF: init_modif_dataobj
%   PROPOSAL: Settings.onlyBfield = 'Error'
%   CON: Caller can assert that A fieldnames are superset of B fieldsnames.
%       CON: Not trivial expression.
%       CON: Not recursive.
%
% NOTE: In reality designed for merging "compatible" data structures (types are compatible; field names may or may not overlap).
% ==> Does not need so many cases.
% ==> Ambiguous if a struct within a struct of settings is to be regarded as a
%     value or as another container for settings (not clear where to end recursion).
% ==> Ambiguous whether to copy/overwrite or not. ==> Operation can not be automatically selected by algorithm.
%
% PROPOSAL: Want to be able to think of structs (recursive) as sets, with or without overlap (intersection).
%   TODO-DEC: How handle "unequal" intersections (same-named fields).
%       Ex: Different field value types (struct, double, string, cell, object, ...).
%       Ex: Different field values (of same type).
%   TODO-DEC: How determine recursion depth?
%   PROPOSAL: Comparison with diffdn somehow?
%       CON: Recursion depth is well defined for diffdn.
%
% PROPOSAL: Permit struct arrays, but restrict policies (for struct arrays only) to 'Error' and 'Overwrite'.
% PROPOSAL: Should not support struct arrays. Should use irf.ds.merge_structs instead.



% DFB = Duplicate Field Policy.
DEFAULT_SETTINGS = struct(...
  'noStructs',    'Error', ...
  'aIsStruct',    'Error', ...
  'bIsStruct',    'Error', ...
  'abAreStructs', 'Error', ...
  'onlyBField',   'Copy');
Settings = irf.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})

% ASSERTIONS
assert(isstruct(A),               'A is not a structure.')
assert(isstruct(B),               'B is not a structure.')
assert(isstruct(Settings),        'Settings is not a struct.')
%     assert(isscalar(A),               'A is not scalar.')
%     assert(isscalar(B),               'B is not scalar.')
assert(all(size(A) == size(B)))
%assert(isequal(size(A), size(B)), 'A and B have different array sizes.')

%===========================
% Iterate over fields in B.
%===========================
bFieldNamesList = fieldnames(B);
for i = 1:length(bFieldNamesList)
  fieldName = bFieldNamesList{i};

  if isfield(A, fieldName)
    %===========================================================
    % CASE: Duplicate fields (field exists in both structures).
    %===========================================================
    afv = A.(fieldName);
    bfv = B.(fieldName);

    abAreStructs = false;
    if ~isstruct(afv) && ~isstruct(bfv)
      behaviour = Settings.noStructs;

    elseif isstruct(afv) && ~isstruct(bfv)
      behaviour = Settings.aIsStruct;

    elseif ~isstruct(afv) && isstruct(bfv)
      behaviour = Settings.bIsStruct;

    else
      behaviour = Settings.abAreStructs;
      abAreStructs = true;
    end

    if     strcmp(behaviour, 'Error')
      error('Structures share identically named fields "%s".', ...
        fieldName)
    elseif strcmp(behaviour, 'Overwrite')
      A.(fieldName) = B.(fieldName);
    elseif strcmp(behaviour, 'Do nothing')
      % Do nothing.
    elseif strcmp(behaviour, 'Recurse') && (abAreStructs)
      %====================================================
      % NOTE: RECURSIVE CALL. Needs the original Settings.
      %====================================================
      A.(fieldName) = irf.ds.add_struct_to_struct(...
        A.(fieldName), B.(fieldName), Settings);
    else
      error(['Can not interpret string value behaviour="%s" for', ...
        ' this combination of field values.'], behaviour)
    end
  else
    %===========================================
    % CASE: Field exists in "B", but not in "A"
    %===========================================
    switch(Settings.onlyBField)
      case 'Copy'
        % NOTE: Not overwrite field, but create new field in A.
        [A.(fieldName)] = deal(B.(fieldName));
      case 'Error'
        error('Field B.%s exists but not A.%s.', fieldName, fieldName)
      otherwise
        error('Illegal setting onlyBField=%s', Settings.onlyBField)
    end

  end
end
end
