@setx AnimRenduDep "%CD%"
@if errorlevel 1 @c:\windows\system32\setx AnimRenduDep "%CD%"
@if errorlevel 1 (echo "Critical error: unable to configure environment...") else echo "Configuration done"
@pause
