/*
   Patch submitted by Ross Piltz for SCILAB 5.5.1
	Tested on Windows XP, Vista, 7, 8, Linux, and Mac OS X.
	Tested using US, UK, and German keyboards.

/*
Several problems with the original code:
 I assume that eventTranslator.setClickAction() pushes actions onto
 some type of stack. The problem is when you have events from the
 keyboard, mouse-click, and mouse-motion then the actions become
 confused and connecting key-down with a key-up isn't possible.
 It also appears that keyEvent.getKeyChar() & keyEvent.getKeyCode() have
 been used as if they were interchangeable.
*/

/*
 * Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 * Copyright (C) 2008-2008 - INRIA - Bruno JOFRET
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt
 *
 */
package org.scilab.modules.gui.events;

import static org.scilab.modules.graphic_objects.graphicObject.GraphicObjectProperties.__GO_ID__;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import org.scilab.modules.action_binding.InterpreterManagement;
import org.scilab.modules.graphic_objects.graphicController.GraphicController;
import org.scilab.modules.gui.utils.SciTranslator;
/*
 * This class is to manage scilab callback through seteventhandler
 * it means call a dedicated scilab function like this :
 * function my_eventhandler(windowsId, mouse X, mouse Y, mouse Button)
 */
public class ScilabEventListener implements KeyListener, MouseListener, MouseMotionListener {

    private String callback;
    private Integer windowsUID;
    private int mouseX = 0;
    private int mouseY = 0;
    private SciTranslator eventTranslator = new SciTranslator();
    private boolean freedom = true;

/* ROP >>> */
/* Assuming we start with the canvas in focus fixes one of my bugs. */
/* This changes seems too simplistic for me, but it work */
/*    private boolean inCanvas = false; */
	private boolean inCanvas = true;
/* <<< ROP */

    private boolean useHandle = true;

    public ScilabEventListener(String callback, Integer windowsUID) {
        eventTranslator.setClickAction(SciTranslator.UNMANAGED);
        this.callback = callback;
        this.windowsUID	= windowsUID;
    }

    // Remove this constructor
    // once event_handler call are unified using handle 
    public ScilabEventListener(String callback, Integer windowsUID, boolean useHandle) {
        eventTranslator.setClickAction(SciTranslator.UNMANAGED);
        this.callback = callback;
        this.windowsUID = windowsUID;
        this.useHandle = useHandle;
    }

/* ROP >>> */
/* Allows specific events "eventAction" to be reported when useHandle is false. */
/* Otherwise, eventTranslator.getClickAction() is reported as in the original code. */
/*    private void callScilab() { */
	private void callScilab(int eventAction) {
/* <<< ROP */

        // @FIXME : choose to send it to scilab or to display it
        //
        if (useHandle) { 
            InterpreterManagement.requestScilabExec(callback + "(getcallbackobject(" + windowsUID + ")," + mouseX + ',' + mouseY + ',' + eventTranslator.getClickAction() + ')');
        } else {
            int windowsId = (Integer) GraphicController.getController().getProperty(windowsUID, __GO_ID__);

/* ROP >>> */
/* I assume this is where Scilab sends the action to the user-defined event-handler. */
/*            InterpreterManagement.requestScilabExec(callback + '(' + windowsId + ',' + mouseX + ',' + mouseY + ',' + eventTranslator.getClickAction() + ')'); */
		InterpreterManagement.requestScilabExec(callback+'('+windowsId+','+mouseX+','+mouseY+','+eventAction+')');
/* <<< ROP */

        }
        //
        //System.out.println("call " + callback+'('+windowsId+','+mouseX+','+mouseY+','+eventTranslator.getClickAction()+')');
    }

/* ROP >>> */
/* Routine no longer needed */
/*    private void invokeScilab() {
        // @FIXME : choose to send it to scilab or to display it
        //
        if (useHandle) { 
            InterpreterManagement.requestScilabExec(callback + "(getcallbackobject(" + windowsUID + ")," + mouseX + ',' + mouseY + ',' + eventTranslator.javaClick2Scilab() + ')');
        } else {
            int windowsId = (Integer) GraphicController.getController().getProperty(windowsUID, __GO_ID__);
            InterpreterManagement.requestScilabExec(callback + '(' + windowsId + ',' + mouseX + ',' + mouseY + ',' + eventTranslator.javaClick2Scilab() + ')');
        }
        //
        //System.out.println("invoke " + callback+'('+windowsId+','+mouseX+','+mouseY+','+eventTranslator.javaClick2Scilab()+')');
    } */
/* <<< ROP */

    public void keyPressed(KeyEvent keyEvent) {
        if (inCanvas) {

/* ROP >>> */
/* I don't know enough about Java to understand what this code is doing. */
/* I found it unnecessary for my solution, do I deleted it. */
/*            if (Character.isJavaIdentifierStart(keyEvent.getKeyChar())) {
                eventTranslator.setClickAction(SciTranslator.javaKey2Scilab(keyEvent.getKeyChar(), keyEvent.isControlDown()));
                callScilab();
            } else { */
/* <<< ROP */

			int keyChar;

/* ROP >>> */
/* Don't use KeyCode and then change the case if SHIFT is down. */
/* Instead use KeyChar to get the UNICODE character. This also fixes a bug */
/* when CAPS-LOCK is engaged. This bug has not been fixed for xGetMouse(). */
/* Otherwise, "ibut" returned by the event-handler and xGetMouse() are identical. */
/*                if (keyEvent.isShiftDown()) {
                    keyChar = keyEvent.getKeyCode();
                } else {
                    keyChar = Character.toLowerCase(keyEvent.getKeyCode());
                    callScilab(); */
			keyChar = keyEvent.getKeyChar();
/* <<< ROP */

/* ROP >>> */
/* Don't process non-character key presses such as <F4>. These could be
/* obtained using keyEvent.getKeyCode() but Scilab doesn't use them anyway. */ 
			if (keyChar != KeyEvent.CHAR_UNDEFINED) {
/* <<< ROP */

/* ROP >>> */
/* If CTRL pressed, convert any UNICODE {^A to ^Z} to {A - Z} */
			if (keyEvent.isControlDown()) {
/* I assume there is a more elegant way to do this */
				if( (keyChar >= 1) && (keyChar <= 26) ) keyChar=keyChar+64;
/* <<< ROP */

             }

/* ROP >>> */
/* Original registered a key-click action then used callScilab() to service it. */
/* It also used SciTranslator to apparently add 1000 to keyChar when CTRL is pressed. */
/* I have been explicit about adding 1000 to keyChar and have completely avoided
/* SciTranslator. I also prevented TAB, ENTER, etc. presses from being served. */
/*                eventTranslator.setClickAction(SciTranslator.javaKey2Scilab(keyChar,
                                               keyEvent.isControlDown())); */
/* Only process printable characters, not TAB, ENTER, etc. */
/* This must occur after we have converted {^A to ^Z} to {A - Z} */
/* Also not elegant */
			if (keyChar >= 32) {
/* If CTRL is pressed, add 1000 to key character */
				if (keyEvent.isControlDown()) keyChar=keyChar+1000;
/* Explicitly process keyChar using callScilab() */
					callScilab(keyChar);
	            }
            }
        }
/* <<< ROP */

    }

/* ROP >>> */
/* The original code ingores the argument. From my memory it is assumed */
/* the key char can be obtained from the last click-action pushed onto */
/* the stack. */
/*    public void keyReleased(KeyEvent arg0) { */
	public void keyReleased(KeyEvent keyEvent) {
/* <<< ROP */

/* ROP >>> */
/* It appears the original code uses SciTranslator.UNMANAGED as a place-marker */
/* in the actions stack. I'm not sure what the code is actually doing. */
/* You should note that I do not test for inCanvas. This allows a key-up */
/* to be processed correctly if we lose focus. */
/*        if (inCanvas && eventTranslator.getClickAction() != SciTranslator.UNMANAGED) { */
/* <<< ROP */

/* ROP >>> */
/* From memory the original code messes up CTRL characters as it trys to
/* use values of -1032 instead of 968. Not sure, long time ago. */
/*            eventTranslator.setClickAction(-eventTranslator.getClickAction());
            callScilab(); */
/* Calculate the keyChar value using the samme method as in keyPressed(). */
		int keyChar;
		keyChar = keyEvent.getKeyChar();
		int keyChar = keyEvent.getKeyChar();
		if (keyChar == 65535) return;
		if (keyEvent.isControlDown()) {
			if( (keyChar >= 1) && (keyChar <= 26) ) keyChar=keyChar+64;
        }
		if (keyChar < 32) return;
/* Subtract 1000 from keyChar as we use callScilab() with -keyChar */
		if (keyEvent.isControlDown()) keyChar=keyChar-1000;
		callScilab(-keyChar);
/* <<< ROP */

/* ROP >>> */
/* Remove unnecessary place-marker */
/*            eventTranslator.setClickAction(SciTranslator.UNMANAGED); */

    }

    public void keyTyped(KeyEvent arg0) {
        // Do nothing !!!
    }

    public void mouseClicked(MouseEvent arg0) {
        mouseX = arg0.getX();
        mouseY = arg0.getY();
        if (arg0.getClickCount() == 1) {
/*	Set click action to "single click"	*/
            eventTranslator.setClickAction(
                SciTranslator.javaButton2Scilab(
                    arg0.getButton(),
                    SciTranslator.CLICKED,
                    arg0.isControlDown()));
		} else {
/*	Set click action to "double clicked"	*/
            /* Means mouseEvent.getClickCount() >= 2 */
            eventTranslator.setClickAction(
                SciTranslator.javaButton2Scilab(
                    arg0.getButton(),
                    SciTranslator.DCLICKED,
                    arg0.isControlDown()));

/* ROP >>> */
/* I really don't understand Java threading. However, I found by */
/* trial-and-error that this must be removed. */
/*            // To unlock javaClick2Scilab done in launchfilter
            synchronized (eventTranslator) {
                eventTranslator.notify();
            } */
/* <<< ROP*/

        }
    }

    public void mouseEntered(MouseEvent arg0) {
        inCanvas = true;
    }

    public void mouseExited(MouseEvent arg0) {
        inCanvas = false;
    }

    public void mousePressed(MouseEvent arg0) {
        if (this.freedom) {
            this.freedom = false;
            mouseX = arg0.getX();
            mouseY = arg0.getY();
/*	Set click action to "button pressed"	*/
            eventTranslator.setClickAction(
                SciTranslator.javaButton2Scilab(
                    arg0.getButton(),
                    SciTranslator.PRESSED,
                    arg0.isControlDown()));
            Thread launchMe = new Thread() {
                public void run() {

/* ROP >>> */
/* I really don't understand Java threading. However, I found by */
/* trial-and-error that this makes a 300 msec delay which */
/* distinguishes mouse-clicks from mouse-presses. */
/*                     invokeScilab(); */
					eventTranslator.javaClick2Scilab();	/* just use as a 300 msec delay */
					callScilab(eventTranslator.getClickAction());
/* <<< ROP*/

                    freedom = true;

/* ROP >>> */
/* Remove unnecessary place-marker
/*                    eventTranslator.setClickAction(SciTranslator.UNMANAGED); */
/* <<< ROP */

                }
            };
            launchMe.start();
        }
    }

    public void mouseReleased(MouseEvent arg0) {

/* ROP >>> */
/* More Java I am not sure of, but I found a way using trial-and-error */
/*        if (eventTranslator.getClickAction() == SciTranslator.UNMANAGED || 
                eventTranslator.getClickAction() == SciTranslator.SCIMOVED) { */
		if (this.freedom) {
/* <<< ROP */

			mouseX = arg0.getX();
			mouseY = arg0.getY();

/* ROP >>> */
/* Again explicitly send action to callScilab() instead of setClickAction() */
/*            eventTranslator.setClickAction(
                SciTranslator.javaButton2Scilab(arg0.getButton(),
                                                SciTranslator.RELEASED,
                                                arg0.isControlDown())); */
			callScilab(
                SciTranslator.javaButton2Scilab(arg0.getButton(),
                                                SciTranslator.RELEASED,
                                                arg0.isControlDown()));
/* <<< ROP */

        }
    }

    public void mouseDragged(MouseEvent arg0) {

/* ROP >>> */
/* Looking at this code again, it appears I am allowing mouse-dragging using */
/* any mouse button pressed, not just the left-button. This should probably */
/* be fixed. */
/*        if (eventTranslator.getClickAction() == SciTranslator.javaButton2Scilab(MouseEvent.BUTTON1, SciTranslator.PRESSED, false)) {
/* <<< ROP */

            mouseX = arg0.getX();
            mouseY = arg0.getY();

/* ROP >>> */
/* Use callScilab() directly instead of via setClickAction, and remove */
/* the place-marker */
/*            callScilab();
            freedom = true;
            eventTranslator.setClickAction(SciTranslator.SCIMOVED);
        } else {
            eventTranslator.setClickAction(SciTranslator.SCIMOVED);
            mouseX = arg0.getX();
            mouseY = arg0.getY();
            callScilab();
        }
        eventTranslator.setClickAction(SciTranslator.UNMANAGED); */
		callScilab(SciTranslator.SCIMOVED);
/* <<< ROP */

        }

    public void mouseMoved(MouseEvent arg0) {

/* ROP >>> */
/* Use callScilab() directly instead of via setClickAction, and remove */
/* the place-marker */
/*        eventTranslator.setClickAction(SciTranslator.SCIMOVED);
        mouseX = arg0.getX();
        mouseY = arg0.getY();
/*        callScilab();
        eventTranslator.setClickAction(SciTranslator.UNMANAGED); */
		callScilab(SciTranslator.SCIMOVED);
/* <<< ROP */

    }

}
