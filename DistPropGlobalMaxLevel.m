function x = DistPropGlobalMaxLevel(varargin)

UncPropLoadNETAssemblies('DistProp');
if nargin > 0
    NET.setStaticProperty('Metas.UncLib.DistProp.Misc.Global.MaxLevel', varargin{1});
end
x = Metas.UncLib.DistProp.Misc.Global.MaxLevel;
