function b = subsref(a,index)
%SUBSREF Define field name indexing for dataobj objects
b = builtin('subsref',a,index);
