===================================================================
Setup for building new jar files
===================================================================

I assume the jar files are being built in D:\Work\Scilab\Rebuilds\
for the new Scilab version and then the final results are saved
in D:\Work\Scilab\Rebuilds\_Changes\*.*.*\


Delete files/folders in D:\Work\Scilab\Rebuilds\ except for:
	__ReadMe__.txt
	_ant_setup.bat
	_apach-ant-*\
	_Changes\	

Copy the modules\ and thirdparty\ folders from the Scilab executable folder to Rebuilds\

Extract from the top folder of the source-code Zip file, and put in Rebuilds\
	build.incl.xml
	build.qa.incl.xml
	scilab.properties.in
	scilab-lib.properties.in
	scilab-lib-doc.properties.in
Remove the .in suffix from the above files
Compare build.incl.xml to the previous version and make appropriate changes, and
	check that the thirdparty *.jar versions agree with the those in thirdparty\

If necessary, update the Apache-ant folder to a new version
Check _ant_setup.bat corresponds to any new JDK and ANT versions


===================================================================
Patching org.scilab.modules.gui.jar
===================================================================

Delete D:\Work\Scilab\Rebuilds\modules\gui\*

Extract from the source-code Zip file
	build.xml and src\ from modules\gui\
and put in D:\Work\Scilab\Rebuilds\modules\gui\


Modify the source files
	ScilabEventListener.java		"fixes mouse & keyboard interactions"
in D:\Work\Scilab\Rebuilds\modules\gui\src\java\org\scilab\modules\gui\events\
	SwingScilabWaitBar.java			"change progress bar to modal"
in D:\Work\Scilab\Rebuilds\src\modules\gui\src\java\org\scilab\modules\gui\bridge\waitbar\


In D:\Work\Scilab\Rebuilds\modules\gui\
run the CMD commands
	..\..\_ant_setup.bat
	ant

The compiled java file
	org.scilab.modules.gui.jar
is located in
	D:\Work\Scilab\Rebuilds\modules\gui\jar


NB: Do "rmdir/s build" and "rmdir/s jar" on subsequent compilations of ant



===================================================================
Patching org.scilab.modules.scinotes.jar
===================================================================

Save D:\Work\Scilab\Rebuilds\modules\scinotes\jar\org.scilab.modules.scinotes.jar

Delete D:\Work\Scilab\Rebuilds\modules\scinotes\*

Extract from the source-code Zip file
	build.xml and src\ from \modules\scinotes\
and put in D:\Work\Scilab\Rebuilds\modules\scinotes\


Save a copy, then moodify the source file:
	SearchWordInFilesAction.java		"change search combo box from FIFO to LIFO"
in D:\Work\Scilab\Rebuilds\modules\scinotes\src\java\org\scilab\modules\scinotes\actions\


In D:\Work\Scilab\Rebuilds\modules\scinotes\
run the CMD commands
	..\..\_ant_setup.bat
	ant

The compiled java file
	org.scilab.modules.scinotes.jar
is located in
	D:\Work\Scilab\Rebuilds\modules\scinotes\jar


NB: Do "rmdir/s build" and "rmdir/s jar" on subsequent compilations of ant


===================================================================
Clean up
===================================================================

Copy the originals of the patched *.java and *.jar files to
	D:\Work\Scilab\Rebuilds\_Changes\*.*.*\orig\

Copy the patched *.java and new *.jar files to 
	D:\Work\Scilab\Rebuilds\_Changes\*.*.*\


Copy
	_ant_setup.bat
	__ReadMe__.txt
	build.qa.xml

to
	D:\Work\Scilab\Rebuilds\_Changes\*.*.*\

In
	D:\Work\Scilab\Rebuilds\
delete
	*.properties
	build.qa.incl.xml
	modules\
	thirdparty\
