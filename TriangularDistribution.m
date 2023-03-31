function d = TriangularDistribution(a, b)
% Triangular Distribution
% Michael Wollensack METAS - 03.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.TriangularDistribution(a, b);
end
