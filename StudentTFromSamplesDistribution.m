function d = StudentTFromSamplesDistribution(samples)
% Student T From Samples Distribution
% Michael Wollensack METAS - 06.03.2023

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('MCProp');
end

d = Metas.UncLib.Core.Unc.StudentTFromSamplesDistribution(samples);
end
