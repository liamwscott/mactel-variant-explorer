@echo off
REM ==========================================================================
REM Run MacTel Explorer (Windows).bat - Windows double-click launcher.
REM
REM Double-click this file to start the MacTel Variant Explorer. The app opens
REM in your default web browser. Close this window to stop the app.
REM ==========================================================================

title MacTel Variant Explorer
color 1F

echo  ============================================================
echo.
echo                  MacTel Variant Explorer
echo.
echo      Starting up - the app will open in your browser.
echo      Keep this window open while you use the app;
echo      close it when you are finished to stop the app.
echo.
echo  ============================================================
echo.

REM Move into the folder this script lives in, so R finds app.R and data\.
cd /d "%~dp0"

REM Locate Rscript: PATH first, then the usual install locations.
set "RSCRIPT="
where Rscript >nul 2>nul
if %errorlevel%==0 (
  set "RSCRIPT=Rscript"
) else (
  for /d %%V in ("%ProgramFiles%\R\R-*") do if exist "%%V\bin\Rscript.exe" set "RSCRIPT=%%V\bin\Rscript.exe"
)

if not defined RSCRIPT (
  echo R is not installed, or could not be found.
  echo Please install R from https://cran.r-project.org first, then run this launcher again.
  echo.
  pause
  exit /b 1
)

echo Using R at: %RSCRIPT%

REM If a previous instance is still holding the port, stop it so this launch
REM can start cleanly (otherwise R fails with "address already in use").
for /f "tokens=5" %%P in ('netstat -ano ^| findstr ":7766 " ^| findstr LISTENING') do (
  echo Stopping a previous instance still using port 7766...
  taskkill /F /PID %%P >nul 2>nul
)

"%RSCRIPT%" launch.R

if %errorlevel% neq 0 (
  echo.
  echo The app exited with an error ^(code %errorlevel%^).
  pause
)
