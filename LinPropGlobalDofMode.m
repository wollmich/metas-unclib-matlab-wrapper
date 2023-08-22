function x = LinPropGlobalDofMode(varargin)

UncPropLoadNETAssemblies('LinProp');
if nargin > 0
    if strcmp(varargin{1}, 'WelchSatterthwaite')
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.DofMode', Metas.UncLib.LinProp.Misc.DofModeType.WelchSatterthwaite);
    elseif strcmp(varargin{1}, 'StudentT')
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.DofMode', Metas.UncLib.LinProp.Misc.DofModeType.StudentT);
    else
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.DofMode', varargin{1});
    end
end
x = Metas.UncLib.LinProp.Misc.Global.DofMode;
