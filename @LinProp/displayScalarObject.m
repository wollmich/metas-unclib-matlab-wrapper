function displayScalarObject(obj)
    name = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
    fprintf('    %s:\n    %s\n\n', name, LinProp.toUncCharColumn(get_value(obj), get_stdunc(obj)));
    
    stack = dbstack();
    
    if numel(stack) == 1
        if obj.IsComplex
            linkStrReal = LinProp.getMethodLink('displayBudget', 'Budget of real part', inputname(1), class(obj), 'real');
            linkStrImag = LinProp.getMethodLink('displayBudget', 'Budget or imag part', inputname(1), class(obj), 'imag');
            linkStr = [linkStrReal, ', ', linkStrImag];
        else
            linkStr = LinProp.getMethodLink('displayBudget', 'Budget', inputname(1), class(obj), '');
        end
        
        methodsStr = sprintf('<a href="matlab:methods(''%s'')">Methods</a>',class(obj));

        if startsWith(get(0, 'Format'), 'long')
            methodsStr = [LinProp.getMethodLink('displayInFormat', 'in Format short', inputname(1), class(obj), 'short'), ', ', methodsStr];
        else
            methodsStr = [LinProp.getMethodLink('displayInFormat', 'in Format long', inputname(1), class(obj), 'long'), ', ', methodsStr];
        end
        
        fprintf('Show %s, %s.\n',linkStr, methodsStr);
    end
end
