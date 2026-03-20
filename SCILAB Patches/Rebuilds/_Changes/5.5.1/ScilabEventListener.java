/*
   Patch submitted by Ross Piltz for SCILAB 5.5.1
	Tested on Windows XP, Vista, 7, 8, Linux, and Mac OS X.
	Tested using US, UK, and German keyboards.
*/

/*
 * Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 * Copyright (C) 2008-2008 - INRIA - Bruno JOFRET
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
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
	private boolean inCanvas = true;	/* assume we have focus at start */

    public ScilabEventListener(String callback, Integer windowsUID) {
		eventTranslator.setClickAction(SciTranslator.UNMANAGED);
		this.callback = callback;
		this.windowsUID = windowsUID;
	}

	private void callScilab(int clickAction) {
		int windowsId = (Integer) GraphicController.getController().getProperty(windowsUID, __GO_ID__);
		InterpreterManagement.requestScilabExec(callback+'('+windowsId+','+mouseX+','+mouseY+
									','+clickAction+')');
	}


	public void keyPressed(KeyEvent keyEvent) {
		int keyChar = keyEvent.getKeyChar();
		if (keyChar == 65535) return;		/* ignore special keys */
		if (!inCanvas) return;			/* ignore if lost focus */ 
/* If CTRL pressed, convert keyChar for ^A to ^Z back to A - Z */
		if (keyEvent.isControlDown()) {
			if( (keyChar >= 1) && (keyChar <= 26) ) keyChar=keyChar+64;
		}
		if (keyChar < 32) return;		/* ignore TAB, ENTER, etc. */
		if (keyEvent.isControlDown()) keyChar=keyChar+1000;
		callScilab(keyChar);
	}

	public void keyReleased(KeyEvent keyEvent) {
		int keyChar = keyEvent.getKeyChar();
		if (keyChar == 65535) return;		/* ignore special keys */
/* If CTRL pressed, convert keyChar for ^A to ^Z back to A - Z */
		if (keyEvent.isControlDown()) {
			if( (keyChar >= 1) && (keyChar <= 26) ) keyChar=keyChar+64;
		}
		if (keyChar < 32) return;		/* ignore TAB, ENTER, etc. */
		if (keyEvent.isControlDown()) keyChar=keyChar-1000;
		callScilab(-keyChar);
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
		}
		else {
/*	Set click action to "double clicked"	*/
			/* Means mouseEvent.getClickCount() >= 2 */
			eventTranslator.setClickAction(
					SciTranslator.javaButton2Scilab(
							arg0.getButton(),
							SciTranslator.DCLICKED,
							arg0.isControlDown()));
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
					eventTranslator.javaClick2Scilab();	/* just use as a 300 msec delay */
					callScilab(eventTranslator.getClickAction());
					freedom = true;
				}
			};
			launchMe.start();
		}
	}

	public void mouseReleased(MouseEvent arg0) {
		if (this.freedom) {
			mouseX = arg0.getX();
			mouseY = arg0.getY();
			callScilab(
				SciTranslator.javaButton2Scilab(arg0.getButton(),
								SciTranslator.RELEASED,
								arg0.isControlDown()));
		}
	}

	public void mouseDragged(MouseEvent arg0) {
		mouseX = arg0.getX();
		mouseY = arg0.getY();
		callScilab(SciTranslator.SCIMOVED);
	}

	public void mouseMoved(MouseEvent arg0) {
		mouseX = arg0.getX();
		mouseY = arg0.getY();
		callScilab(SciTranslator.SCIMOVED);
	}



}
