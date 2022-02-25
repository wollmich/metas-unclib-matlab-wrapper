function displayBudget(obj, part)

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
                varName = sprintf('real(%s)', varName);
            case 'imag'
                obj = imag(obj);
                varName = sprintf('imag(%s)', varName);
            otherwise
                error('Input argument ''part'' must be ''real'' or ''imag''.');
        end
    end

    tree = Metas.UncLib.LinProp.UncBudget.ComputeTreeUncBudget(obj.NetObject);
    if tree.Length == 0
        fprintf('\n  Variable %s has no uncertainty contributions.\n', varName);
        return;
    end
    
    fprintf('\n  Budget of %s:\n', varName);
    
    [description, componenent, percentage] = parseSubTree(tree);

    nTableRows = numel(description);
    maxDescriptionLength = max(cellfun(@numel, description));

    description = sprintf(['%-' num2str(maxDescriptionLength) 's'], description);
    description = reshape(description, [], nTableRows);

    componenent = LinProp.toCharColumn(componenent);
    componenent = reshape(componenent, [], nTableRows);

    percentage = sprintf('%6.2f%%\n', percentage);
    percentage = reshape(percentage, [], nTableRows);

    spacing = 4;

    text = [...
        description; ...
        repmat(' ', spacing, nTableRows); ...
        componenent; ...
        repmat(' ', spacing, nTableRows); ...
        percentage ...
    ];
    fprintf('<strong>    Description</strong>%*s<strong>Component</strong>    <strong>Percentage</strong>\n', ...
        size(description, 1) - numel('    Description') - numel('Component') + size(componenent, 1) + spacing, '');
    fprintf('%s\n', text);
end
    
function [description, componenent, percentage] = parseSubTree(tree, level)
    description = strings(1, 0);
    componenent = zeros(1, 0);
    percentage = zeros(1, 0);
    
    if nargin == 1
        level = 4;
    end

    for ii = 1:tree.Length
        
        [subTest, subComponent, subPercentage] = parseSubTree(tree(ii).SubItems, level+2);
        
        description = [...
            description, ...
            sprintf('%s%s', repmat(' ', 1, level), tree(ii).ShortDescription), ...
            subTest...
        ];
        
        componenent = [...
            componenent, ...
            tree(ii).UncComponent, ...
            subComponent...
        ];
        
        percentage = [...
            percentage, ...
            tree(ii).UncPercentage, ...
            subPercentage...
        ];
    end
end