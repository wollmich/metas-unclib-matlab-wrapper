function d = RandomChoicesFromSamples(seed, samples)
% Random Choices From Samples
% Michael Wollensack METAS - 06.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.RandomChoicesFromSamples(seed, samples);
end
