powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'DistProp' | Out-File -encoding default @DistProp\DistProp.m"
powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'MCProp' | Out-File -encoding default @MCProp\MCProp.m"

pause
