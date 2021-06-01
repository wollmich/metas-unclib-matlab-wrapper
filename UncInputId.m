function id = UncInputId(varargin)
% Uncertainty Input Id
% Michael Wollensack METAS - 28.02.2013

global UncLibModul;
if ~ischar(UncLibModul)
    UncPropLoadNETAssemblies('LinProp');
end

if (nargin == 1)
    id = Metas.UncLib.Core.Unc.InputId(varargin{1});
else
    id = Metas.UncLib.Core.Unc.InputId.NewInputId();
end

end
