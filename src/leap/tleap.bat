@echo off

set AMBERHOME=%~dp0
set AMBERHOME=%AMBERHOME:~0,-1%\..

%AMBERHOME%\bin\teLeap -I%AMBERHOME%/dat/leap/prep -I%AMBERHOME%/dat/leap/lib -I%AMBERHOME%/dat/leap/parm -I%AMBERHOME%/dat/leap/cmd %*