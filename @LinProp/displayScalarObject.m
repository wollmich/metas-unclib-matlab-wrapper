function displayScalarObject(obj)
    name = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
    fprintf('    %s:\n    %s\n\n', name, LinProp.toUncCharColumn(get_value(obj), get_stdunc(obj)));
    
    stack = dbstack();
    
    if numel(stack) == 1
        if obj.IsComplex
            linkStrReal = sprintf('<a href="%s">Budget of real part</a>',getBudgetLink(inputname(1), class(obj), 'real'));
            linkStrImag = sprintf('<a href="%s">Budget or imag part</a>',getBudgetLink(inputname(1), class(obj), 'imag'));
            linkStr = [linkStrReal, ', ', linkStrImag];
        else
            linkStr = sprintf('<a href="%s">Budget</a>',getBudgetLink(inputname(1), class(obj), ''));
        end

        methodsStr = sprintf('<a href="matlab:methods(''%s'')">Methods</a>',class(obj));

        fprintf('Show %s, %s.\n',linkStr, methodsStr);
    end
end

function link=getBudgetLink(varName, class, part)

    if nargin < 3
        partStr = '';
    else
        partStr = sprintf(', ''%s''', part);
    end

    link=sprintf(['matlab:if exist(''%s'', ''var'')&&isa(%s, ''%s''), ', ...
        'displayBudget(%s%s); ', ...
        'else, ', ...
        'fprintf(''Unable to display budget. %s refers to a deleted object.\\n''), ', ...
        'end'], ...
        varName,varName,class,varName,partStr,varName);
end