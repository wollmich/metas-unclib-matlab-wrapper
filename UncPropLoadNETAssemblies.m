function UncPropLoadNETAssemblies(modul)
% UncPropLoadNETAssemblies
% Michael Wollensack METAS - 28.02.2013

global UncLibModul;
if strcmp(UncLibModul, modul) == 0
    UncLibModul = modul;
    % 64 bit or 32 bit MATLAB?
    v = computer;
    if strcmp(v, 'PCWIN')
      d = winqueryreg('HKEY_LOCAL_MACHINE', 'Software\Microsoft\.NETFramework\AssemblyFolders\Metas.UncLib');
    elseif strcmp(v, 'PCWIN64')
      d = winqueryreg('HKEY_LOCAL_MACHINE', 'Software\Wow6432Node\Microsoft\.NETFramework\AssemblyFolders\Metas.UncLib');
    else
      error('Metas.UncLib is not installed.')
    end

    disp(['Loading .NET Assembly: ' d 'Metas.UncLib.Core.dll'])
    NET.addAssembly([d 'Metas.UncLib.Core.dll']);
    disp(['Loading .NET Assembly: ' d 'Metas.UncLib.' modul '.dll'])
    NET.addAssembly([d 'Metas.UncLib.' modul '.dll']);
    disp(['Loading .NET Assembly: ' d 'Metas.UncLib.Optimization.dll'])
    NET.addAssembly([d 'Metas.UncLib.Optimization.dll']);
end
