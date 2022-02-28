function displayScalarObject(obj)
    name = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
    fprintf('    %s:\n    %s\n\n', name, LinProp.toUncCharColumn(get_value(obj), get_stdunc(obj)));
    
    stack = dbstack();
    
    if numel(stack) == 1
        if obj.IsComplex
            linkStrReal = LinProp.getMethodLink('displayContributions(%s, ''real'')', 'Contributions of real part', inputname(1), class(obj));
            linkStrImag = LinProp.getMethodLink('displayContributions(%s, ''imag'')', '... imag part', inputname(1), class(obj));
            linkStr = [linkStrReal, ', ', linkStrImag];
        else
            linkStr = LinProp.getMethodLink('displayContributions(%s)', 'Contributions', inputname(1), class(obj));
        end
        
        methodsStr = sprintf('<a href="matlab:methods(''%s'')">Methods</a>',class(obj));

        if startsWith(get(0, 'Format'), 'long')
            methodsStr = [LinProp.getMethodLink('displayInFormat(%s, ''short'')', 'in Format short', inputname(1), class(obj)), ', ', methodsStr];
        else
            methodsStr = [LinProp.getMethodLink('displayInFormat(%s, ''long'')', 'in Format long', inputname(1), class(obj)), ', ', methodsStr];
        end
        
        fprintf('Show %s, %s.\n',linkStr, methodsStr);
    end
end
