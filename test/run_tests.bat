@echo off
setlocal
set PATH=%PATH%;c:\programme\python
set PYTHONPATH=%PYTHONPATH%;h:\src\lang\Python-2.2.1\Lib
python -c "import run_all_exe; run_all_exe.main()" %*
if errorlevel 9009 echo you do not have python in your PATH
REM color 00 causes this script to exit with non-zero exit status
if errorlevel 1 color 00
endlocal
