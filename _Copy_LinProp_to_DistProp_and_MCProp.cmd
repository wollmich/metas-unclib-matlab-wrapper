powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'DistProp' -replace 'DistProp only', 'LinProp only' | Out-File -encoding default @DistProp\DistProp.m"
powershell -Command "(gc @LinProp\LinProp.m) -replace 'LinProp', 'MCProp' -replace 'MCProp only', 'LinProp only' | Out-File -encoding default @MCProp\MCProp.m"

pause
