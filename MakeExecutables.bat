:: This file is part of EMBL-HLA-Submission.
::
:: EMBL-HLA-Submission is free software: you can redistribute it and/or modify
:: it under the terms of the GNU Lesser General Public License as published by
:: the Free Software Foundation, either version 3 of the License, or
:: (at your option) any later version.
::
:: EMBL-HLA-Submission is distributed in the hope that it will be useful,
:: but WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
:: GNU Lesser General Public License for more details.
::
:: You should have received a copy of the GNU Lesser General Public License
:: along with EMBL-HLA-Submission. If not, see <http://www.gnu.org/licenses/>.

:: Version 1.0

:: This bat file is intended to create an executable for the windows environment.  
:: It uses Anaconda for python 2.7 to keep track of packages. 

:: See the file README.MD for how to set up your anaconda environment.

SET CodePath=src
SET BinPath=..\bin
SET SpecFile=AlleleSubInstallerOptions_Windows.spec
SET CondaEnvironment=AlleleSubEnvironment

:: Run Pyinstaller to create executables
cd %CodePath%
activate %CondaEnvironment% && pyinstaller %SpecFile% --distpath %BinPath% && deactivate
