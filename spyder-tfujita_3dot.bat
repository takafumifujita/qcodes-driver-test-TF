@echo off

set USERNAME=tfujita
set QTTUSERDIR=D:\users\%USERNAME%
set QCODESFRONTEND=spyder

cd D:/Users/%USERNAME%
start cmd.exe /k "activate tfujita_3dot"

echo "Starting spyder" 
SET PYTHONPATH=%PYTHONPATH%;%QTTUSERDIR%
call activate tfujita_3dot

spyder --show-console --new-instance -w d:\users\%USERNAME%
call deactivate tfujita_3dot