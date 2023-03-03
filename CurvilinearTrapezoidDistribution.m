function d = CurvilinearTrapezoidDistribution(a, b, d)
% Curvilinear Trapezoid Distribution
% Michael Wollensack METAS - 03.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.CurvilinearTrapezoidDistribution(a, b, d);
end
