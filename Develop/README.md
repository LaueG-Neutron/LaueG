A snapshot of the files in my C:\DevApps\LaueG folder which I used to develop the software, minus many temporary files

See the "Documentation" on this Git has more information on Developer and Release versions of LaueG

The file _CONTENTS.txt has information on individual files and folders

In particular, the _Release\ folder is used to create new installations and contains:
```	Source\				Fortran & Scilab source files, XML and some files copied from ..\ and .\
	WindowsBuild\		Various misc. files needed to make an installation
	BuildFiles\			All files needed to make installations for Windows/Mac/Linux
	Installs\			Installation executable for Windows made using Inno Setup

	_Notes.txt			Earlier notes on building a new *.exe installation file
	SaveSource.bat		Copies source files from ..\ to Source\
	MakeBuilds.bat		Makes a build in BuildFiles\ from ..\, Source\ & WindowsBuild\
	WindowsInstall.iss	Uses Inno Setup to make a Windows installation file from BuildFiles\

The *.bat files are meant to be run manually, not by LaueG
