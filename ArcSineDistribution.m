function d = ArcSineDistribution(a, b)
% Arc Sine Distribution
% Michael Wollensack METAS - 03.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.ArcSineDistribution(a, b);
end
