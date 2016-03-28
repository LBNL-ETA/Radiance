@echo off
setlocal
set PATH="c:\program files\python 2.7";%PATH%
set PYTHONPATH="c:\program files\python 2.7\lib"
python run_all.py %*
if errorlevel 9009 echo you do not have python in your PATH
REM color 00 causes this script to exit with non-zero exit status
if errorlevel 1 color 00
endlocal
