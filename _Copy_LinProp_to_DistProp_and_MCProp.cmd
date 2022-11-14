powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'DistProp' -replace 'DistProp only', 'LinProp only' | Out-File -encoding default @DistProp\DistProp.m"
powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'MCProp' -replace 'MCProp only', 'LinProp only' | Out-File -encoding default @MCProp\MCProp.m"

powershell -Command "Get-ChildItem -Path +LinProp *.m -recurse | Foreach-Object { (gc +LinProp\$_) -replace 'LinProp', 'DistProp' | Out-File -encoding default +DistProp\$_"}"
powershell -Command "Get-ChildItem -Path +LinProp *.m -recurse | Foreach-Object { (gc +LinProp\$_) -replace 'LinProp', 'MCProp' | Out-File -encoding default +MCProp\$_"}"

pause
