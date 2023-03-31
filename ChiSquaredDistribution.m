function d = ChiSquaredDistribution(k)
% Chi-squared Distribution
% Michael Wollensack METAS - 16.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.ChiSquaredDistribution(k);
end
