% IRF.DATATYPES show information on common data types
%
% Common data types
% -----------------
% DataObject     - the same as dataset in CAA, equal to cdf file representation in matlab
% VariableStruct - Structure including all variable data in original format and metadata,
%	usually extracted from DataObject.
% VariableIrf    - Variable as a simple structure in a default format for irfu-matlab
%	Variable.t        - TimeArray
%	Variable.[data]   - data matrix of dimension [t x dim1 x dim2]... (typical example energy spectrograms)
%	Variable.[unit,label,dimunit,dimlabel,dimvec]
%	Variable.[vec]    - data matrix of dimension [t x ndim], (for example ndim=3 for field vector)
%	Variable.[abs]    - absolute value of vector in case exists Variable.vec
% variableMat    - matrix where first column is time and the other columns are data [DEPRECATED]

help irf.datatypes;
