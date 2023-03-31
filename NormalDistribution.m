function d = NormalDistribution(mu, sigma)
% Normal Distribution
% Michael Wollensack METAS - 01.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.NormalDistribution(mu, sigma);
end
