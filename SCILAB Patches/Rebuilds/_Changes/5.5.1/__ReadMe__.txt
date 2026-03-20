My patches to SCILAB:


scirender.jar:\jogl\texture\JoGLTextureManager.java
						fixed rendering problem with 'marks'

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
                    <attribute name="Specification-Version" value="scilab-5.5.1 for LaueG"/>
or the equivalent file names for new file versions.


Copy C:\Program Files\scilab-5.5.1\modules to D:\Work\Scilab\Rebuilds\modules


Download and unzip SCILAB prerequirements files

Move java\ and thirdparty\ folders to D:\Work\Scilab\Rebuilds\

Copy the following to ant_setup.bat:
	set ANT_HOME=D:\Work\Scilab\Rebuilds\java\ant
	set JAVA_HOME=D:\Work\Scilab\Rebuilds\java\jdk
	set PATH=%PATH%;%ANT_HOME%\bin
	set CLASSPATH=



==================================================
To build "scirenderer.jar" (special case):
==================================================

Scirenderer is a special case as it has the third-party locations defined in
the local "build.xml" instead of "D:\Work\Scilab\Rebuilds\build.incl.xml".


Delete D:\Work\Scilab\Rebuilds\modules\scirenderer\*

Copy	D:\Work\Scilab\Rebuilds\src\modules\scirenderer\build.xml
	D:\Work\Scilab\Rebuilds\src\modules\scirenderer\src
to
	D:\Work\Scilab\Rebuilds\modules\scirenderer\

*** Delete any build\ or jar\ directories from previous builds ***

Edit build.incl.xml by changing the following lines:
        <pathelement location="${jogl2.jar}"/>
        <pathelement location="${gluegen2.jar}"/>
        <pathelement location="${jlatexmath.jar}"/>
to:
        <pathelement location="../../thirdparty/jogl2.jar"/>
        <pathelement location="../../thirdparty/gluegen2-rt.jar"/>
        <pathelement location="../../thirdparty/jlatexmath.jar"/>

Change *.java files as required

From a CMD:
	cd D:\Work\Scilab\Rebuilds\modules\scirenderer\
	..\..\ant_setup.bat
	ant

Which will create "scirenderer-${version}.jar" in modules\scirenderer\jar\.
Rename this file to "scirenderer.jar".



==================================================
To build "renderer.jar" (and other .jar files):
==================================================

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
