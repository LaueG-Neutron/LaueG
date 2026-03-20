My patches to SCILAB:


gui.jar:events\ScilabEventListener.java		fixed mouse & keyboard interactions
gui.jar:editor\Editor.java			removed right-click "editor" menu on figures

renderer.jar:JoGLView\interaction\DragZoomRotateInteraction.java
						removed rotate axes on right-drag


===============================================================

Download and unzip SCILAB source files

Move modules\ for the zip file to D:\Work\Scilab\Rebuilds\src\

Move from the zip file:
	build.incl.xml
	build.qa.incl.xml
	scilab.properties.in
	scilab-lib.properties.in
	scilab-lib-doc.properties.in
to D:\Work\Scilab\Rebuilds

Remove the .in suffix from the above files

Edit build.incl.xml by changing the following lines:
        <pathelement location="${flexdock.jar}"/>
        <pathelement location="${gluegen2.jar}"/>
        <pathelement location="${jrosetta-API.jar}"/>
        <pathelement location="${jrosetta-engine.jar}"/>
        <pathelement location="${jogl2.jar}"/>
        <pathelement location="${jhall.jar}"/>
        <pathelement location="${scirenderer.jar}"/>
        <pathelement location="${jlatexmath.jar}"/>
                    <attribute name="Specification-Version" value="${SCIVERSION}"
to:
        <pathelement location="${thirdparty.dir}/flexdock-1.2.4.jar"/>
        <pathelement location="${thirdparty.dir}/gluegen2-rt.jar"/>
        <pathelement location="${thirdparty.dir}/jrosetta-API.jar"/>
        <pathelement location="${thirdparty.dir}/jrosetta-engine.jar"/>
        <pathelement location="${thirdparty.dir}/jogl2.jar"/>
        <pathelement location="${thirdparty.dir}/jhall.jar"/>
        <pathelement location="${modules.dir}/scirenderer/jar/scirenderer.jar"/>
        <pathelement location="${thirdparty.dir}/jlatexmath-1.0.3.jar"/>
                    <attribute name="Specification-Version" value="scilab-5.5.2 for LaueG"/>
or the equivalent file names for new file versions.


Copy C:\Program Files\scilab-5.5.2\modules to D:\Work\Scilab\Rebuilds\modules


Download and unzip SCILAB prerequirements files

Move java\ and thirdparty\ folders to D:\Work\Scilab\Rebuilds\

Copy the following to ant_setup.bat:
	set ANT_HOME=D:\Work\Scilab\Rebuilds\java\ant
	set JAVA_HOME=D:\Work\Scilab\Rebuilds\java\jdk


===================================================================
  To build "renderer.jar" 

Delete D:\Work\Scilab\Rebuilds\modules\renderer\*

Copy
	D:\Work\Scilab\Rebuilds\src\modules\renderer\build.xml
	D:\Work\Scilab\Rebuilds\src\modules\renderer\src
to
	D:\Work\Scilab\Rebuilds\modules\renderer\

*** Delete any build\ or jar\ directories from previous builds ***

Change *.java files as required

From a CMD:
	cd D:\Work\Scilab\Rebuilds\modules\renderer\
	..\..\ant_setup.bat
	ant

Which will create org.scilab.modules.renderer.jar in modules\renderer\jar\


===================================================================
  To build "gui.jar" 

Same as for renderer.jar, just change the file/folder names


===================================================================

The build.incl.xml & build.qa.incl.xml files are for Rebuilds\
The build.xml is for Rebuilds\modules\scirenderer\


