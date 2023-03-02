function d = DistributionFromSamples(seed, samples)
% Distribution From Samples
% Michael Wollensack METAS - 01.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.DistributionFromSamples(seed, samples);
end
