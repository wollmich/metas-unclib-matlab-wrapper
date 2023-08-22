function x = LinPropGlobalFromSamplesMode(varargin)

UncPropLoadNETAssemblies('LinProp');
if nargin > 0
    if strcmp(varargin{1}, 'Dof')
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.FromSamplesMode', Metas.UncLib.LinProp.Misc.FromSamplesModeType.Dof);
    elseif strcmp(varargin{1}, 'ExpandInputCovariance')
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.FromSamplesMode', Metas.UncLib.LinProp.Misc.FromSamplesModeType.ExpandInputCovariance);
    else
        NET.setStaticProperty('Metas.UncLib.LinProp.Misc.Global.FromSamplesMode', varargin{1});
    end
end
x = Metas.UncLib.LinProp.Misc.Global.FromSamplesMode;
