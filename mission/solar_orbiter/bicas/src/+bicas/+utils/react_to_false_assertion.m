function react_to_false_assertion(giveError, msg)
% Function for either giving a warning, or an error depending on a setting (presumably a global setting).

if giveError
    error('BICAS:Assertion', msg)
else
    irf.log('w', ['FALSE ASSERTION: ', msg])
end

end