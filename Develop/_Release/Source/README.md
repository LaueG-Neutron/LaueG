I zipped the files/folders and removed the originals

The BuildFiles\ subfolder is used in the creation of installation executables

Folder Contents of Source\BuildFiles\*.7z files:

	Source\				Fortran & Scilab source files, XML and some files copied from ..\ and .\
	WindowsBuild\		Various misc. files needed to make an installation
	BuildFiles\			All files needed to make installations for Windows/Mac/Linux
	Installs\			Installation executable for Windows made using Inno Setup

	_Notes.txt			This file
	SaveSource.bat		Copies source files from ..\ to Source\
	MakeBuilds.bat		Makes a build in BuildFiles\ from ..\, Source\ & WindowsBuild\
	WindowsInstall.iss	Uses Inno Setup to make a Windows installation file from BuildFiles\

The *.bat files are meant to be run within LaueG using the "My Version" commands
