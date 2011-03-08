function b = subsref(a,index)
%SUBSREF Define field name indexing for dataobj objects

switch index.type
  case '()'
    switch index.subs{:}
      case 1,
        b=a.GlobalAttributes;
      case 2,
        b=a.VariableAttributes;
      otherwise
        error('Index out of range')
    end
  case '.'
    b=eval(['a.' index.subs]);
  case '{}'
    error('Cell array indexing not supported by dataobj objects')
end