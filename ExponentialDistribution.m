function d = ExponentialDistribution(mu)
% Exponential Distribution
% Michael Wollensack METAS - 03.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.ExponentialDistribution(mu);
end
