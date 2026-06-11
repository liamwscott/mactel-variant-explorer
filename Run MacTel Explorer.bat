@echo off
REM ==========================================================================
REM Run MacTel Explorer.bat - Windows double-click launcher.
REM
REM Double-click this file to start the MacTel Variant Explorer. The app opens
REM in your default web browser. Close this black window to stop the app.
REM ==========================================================================

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
"%RSCRIPT%" launch.R

if %errorlevel% neq 0 (
  echo.
  echo The app exited with an error ^(code %errorlevel%^).
  pause
)
