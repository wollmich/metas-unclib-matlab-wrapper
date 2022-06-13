function displayFooter(obj, inputname, matrixDisplay)

    if isempty(matlab.mixin.CustomDisplay.getDetailedFooter(obj))
        return;
    end
    if nargin < 3
        matrixDisplay = false;
    end

    links = '';

    if startsWith(get(0, 'Format'), 'long')
        links = [methodLink('displayInFormat(%s, ''short'')', 'in Format short', inputname, class(obj)), links];
    else
        links = [methodLink('displayInFormat(%s, ''long'')', 'in Format long', inputname, class(obj)), links];
    end

    global UncLibMatrixDisplay
    if matrixDisplay
        if strcmp(UncLibMatrixDisplay, 'separate')
            links = [methodLink('LinProp.setMatrixDisplay(''combined'');display(%s)', 'with combined uncertainties', inputname, class(obj)), links];
        else
            links = [methodLink('LinProp.setMatrixDisplay(''separate'');display(%s)', 'with separate uncertainties', inputname, class(obj)), links];
        end
    end
    
    fprintf('Show %s<a href="matlab:methods(''%s'')">Methods</a>.\n', links, class(obj));

end

function link = methodLink(method, text, varName, class)

    methodCall = sprintf(method, varName);

    link = sprintf('matlab:if exist(''%s'', ''var'')&&isa(%s, ''%s''), %s; else, fprintf(''Unable to display variable. %s refers to a deleted object.\\n''), end', ...
        varName, varName, class, methodCall, varName);

    link = sprintf('<a href="%s">%s</a>, ', link, text);
end
