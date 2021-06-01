function x = MCPropGlobalN(varargin)

UncPropLoadNETAssemblies('MCProp');
if nargin > 0
    NET.setStaticProperty('Metas.UncLib.MCProp.Misc.Global.n', varargin{1});
end
x = Metas.UncLib.MCProp.Misc.Global.n;
