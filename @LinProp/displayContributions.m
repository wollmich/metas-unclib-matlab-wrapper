function displayContributions(obj, part)

    varName = inputname(1);

    if numel(obj) ~= 1
        error('Can only show budget for scalars');
    end
    
    if obj.IsComplex
        if nargin < 2
            error('Input argument ''part'' not specified.');
        end
        
        switch (part)
            case 'real'
                obj = real(obj);
                varCall = 'real(%s)';
            case 'imag'
                obj = imag(obj);
                varCall = 'imag(%s)';
            otherwise
                error('Input argument ''part'' must be ''real'' or ''imag''.');
        end
    else
        varCall = '%s';
    end

    switch (obj.NetObject.Dependencies.Length)
        case 0
            fprintf('\n  Variable %s has no uncertainty contributions.\n', sprintf(varCall, varName));
            return;
        case 1
            fprintf('\n  Variable %s has 1 uncertainty contribution.\n', sprintf(varCall, varName));
        otherwise
            fprintf('\n  Variable %s has %i uncertainty contributions.\n', sprintf(varCall, varName), obj.NetObject.Dependencies.Length);
    end
    
    if exist('LoadVNATools', 'file')
        LoadVNATools();
        tree = Metas.UncLib.LinProp.UncBudget.ComputeTreeUncBudget(obj.NetObject, @Metas.Vna.Data.UncIdDefs.GetInfluenceInfo2);
    else
        tree = Metas.UncLib.LinProp.UncBudget.ComputeTreeUncBudget(obj.NetObject);
    end
    
    commonPath = strings(0, 0);
    while tree.Length == 1 && tree(1).SubItems.Length > 0
        % Here we have to grow the array in a loop. We cannot know the
        % depth and it is only a few interations.
        commonPath(end+1) = string(tree(1).ShortDescription); %#ok<AGROW>
        tree = tree(1).SubItems;
    end
    commonPath = char(strjoin(commonPath, ' → '));
    
    fprintf('  Top level contributions are:\n\n');
    
    
    description = strings(1, tree.Length);
    componenent = zeros(1, tree.Length);
    percentage  = zeros(1, tree.Length);
    
    for ii = 1:tree.Length
        description{ii} = char(tree(ii).ShortDescription);
        componenent(ii) = tree(ii).UncComponent;
        percentage(ii)  = tree(ii).UncPercentage;
    end

    nTableRows = numel(description);
    maxDescriptionLength = max(cellfun(@numel, description));

    description = sprintf(['%-' num2str(maxDescriptionLength) 's'], description);
    description = reshape(description, [], nTableRows);

    componenent = LinProp.toCharColumn(componenent);
    componenent = reshape(componenent, [], nTableRows);

    percentage = sprintf('%6.2f%%\n', percentage);
    percentage = reshape(percentage, [], nTableRows);
    
    % Idea: Make display depend on number of lines and number of
    % dependencies
    
    spacing1 = 8;
    spacing2 = 4;

    if ~isempty(commonPath)
         description = [repmat(' ', 6, size(description, 2)); description];
    else
         description = [repmat(' ', 4, size(description, 2)); description];
    end
    
    % The description must have a minimum width, so the headers do not
    % overlap too much.
    descriptionMinWidth = max(25, length(commonPath));
    if size(description, 1) < descriptionMinWidth
        description = [description; repmat(' ', descriptionMinWidth - size(description, 1), size(description, 2))];
    end
    
    text = [...
        description; ...
        repmat(' ', spacing1, nTableRows); ...
        componenent; ...
        repmat(' ', spacing2, nTableRows); ...
        percentage ...
    ];
    fprintf('<strong>    Description</strong>%*s<strong>Component</strong>    <strong>Percentage</strong>\n', ...
        size(description, 1) - numel('    Description') - numel('Component') + size(componenent, 1) + spacing1, '');
    if ~isempty(commonPath)
        fprintf('    %s\n', commonPath);
    end
    fprintf('%s\n', text);
    
    methodCall = sprintf('unc_budget(%s)', varCall);
    fprintf('%s.\n', methodLink(methodCall, 'Open full budget', varName, class(obj)));
    
end

function link = methodLink(method, text, varName, class)

    methodCall = sprintf(method, varName);

    link = sprintf('matlab:if exist(''%s'', ''var'')&&isa(%s, ''%s''), %s; else, fprintf(''Unable to display variable. %s refers to a deleted object.\\n''), end', ...
        varName, varName, class, methodCall, varName);

    link = sprintf('<a href="%s">%s</a>', link, text);
end