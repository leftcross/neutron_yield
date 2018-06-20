@echo off
set JASROOT=$~dp0\..
set JASJAR="jas-core-3.2.7.jar"
java -Xmx1024m -Dapplication.home=%JASROOT% %JASJVM_ARGS% -jar %JASROOT%/lib/%JASJAR% %*