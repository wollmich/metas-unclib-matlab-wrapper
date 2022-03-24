function displayFooter(obj, inputname)

    if isempty(matlab.mixin.CustomDisplay.getDetailedFooter(obj))
        return;
    end

    links = '';

    if isscalar(obj)
        if obj.IsComplex
            linkStrReal = methodLink('displayContributions(%s, ''real'')', 'Contributions of real part', inputname, class(obj));
            linkStrImag = methodLink('displayContributions(%s, ''imag'')', '... imag part', inputname, class(obj));
            links = [linkStrReal, linkStrImag];
        else
            links = methodLink('displayContributions(%s)', 'Contributions', inputname, class(obj));
        end
    end

    if startsWith(get(0, 'Format'), 'long')
        links = [methodLink('displayInFormat(%s, ''short'')', 'in Format short', inputname, class(obj)), links];
    else
        links = [methodLink('displayInFormat(%s, ''long'')', 'in Format long', inputname, class(obj)), links];
    end

    fprintf('Show %s, <a href="matlab:methods(''%s'')">Methods</a>.\n', links, class(obj));

end

function link = methodLink(method, text, varName, class)

    methodCall = sprintf(method, varName);

    link = sprintf('matlab:if exist(''%s'', ''var'')&&isa(%s, ''%s''), %s; else, fprintf(''Unable to display variable. %s refers to a deleted object.\\n''), end', ...
        varName, varName, class, methodCall, varName);

    link = sprintf('<a href="%s">%s</a>, ', link, text);
end