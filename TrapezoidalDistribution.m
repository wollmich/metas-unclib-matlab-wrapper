function d = TrapezoidalDistribution(a, b, beta)
% Trapezoidal Distribution
% Michael Wollensack METAS - 03.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.TrapezoidalDistribution(a, b, beta);
end
