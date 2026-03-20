/*
 * Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 * Copyright (C) 2009-2012 - DIGITEO - Pierre Lando
 * Copyright (C) 2013 - Scilab Enterprises - Calixte DENIZET
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt
 */
package org.scilab.modules.renderer.JoGLView.interaction;

import java.awt.Component;
import java.awt.Cursor;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import org.scilab.modules.commons.OS;
import org.scilab.modules.graphic_objects.axes.Axes;
import org.scilab.modules.graphic_objects.graphicController.GraphicController;
import org.scilab.modules.graphic_objects.graphicObject.GraphicObjectProperties;
import org.scilab.modules.renderer.JoGLView.DrawerVisitor;
import org.scilab.modules.renderer.JoGLView.util.ScaleUtils;

/**
 * This class manage figure interaction.
 *
 * @author Pierre Lando
 */
public class DragZoomRotateInteraction extends FigureInteraction {

    private static final int XY_TRANSLATION_MODIFIER = MouseEvent.BUTTON1_MASK;
    private static final int Z_TRANSLATION_MODIFIER = MouseEvent.BUTTON1_MASK | MouseEvent.ALT_MASK;
    private static final int ROTATION_MODIFIER = MouseEvent.BUTTON3_MASK;
    private static final int MACOSX_ROTATION_MODIFIER = MouseEvent.BUTTON1_MASK | MouseEvent.CTRL_MASK;

    /**
     * The box size is multiply by this value.
     */
    private static final double ZOOM_FACTOR = 1.02;

    private final MouseListener mouseListener;
    private final MouseWheelListener mouseWheelListener;
    private final MouseMotionListener mouseMotionListener;


    private Axes currentAxes;


    /**
     * Default constructor.
     * @param drawerVisitor parent drawer visitor.
     */
    public DragZoomRotateInteraction(DrawerVisitor drawerVisitor) {
        super(drawerVisitor);
        mouseMotionListener = new FigureMouseMotionListener();
        mouseWheelListener = new FigureMouseWheelListener();
        mouseListener = new FigureMouseListener();
    }

    @Override
    protected void changeEnable(boolean isEnable) {
        Component component = getDrawerVisitor().getComponent();
        if (component != null) {
            if (isEnable) {
                component.addMouseListener(mouseListener);
                component.addMouseWheelListener(mouseWheelListener);
            } else {
                component.removeMouseListener(mouseListener);
                component.removeMouseMotionListener(mouseMotionListener);
                component.removeMouseWheelListener(mouseWheelListener);
            }
        }
    }

    public void setTranslationEnable(boolean status) {
    }

    /**
     * This {@see MouseListner} activate the {@see MouseMotionListener} when at least
     * one button is pressed.
     * The event is saved in {@see previousEvent}
     */
    private class FigureMouseListener extends MouseAdapter implements MouseListener {

        private int pressedButtons = 0;

        @Override
        public void mousePressed(MouseEvent e) {
            if (pressedButtons == 0) {
                if (currentAxes == null) {
                    currentAxes = getUnderlyingAxes(e.getPoint());
                    if (currentAxes != null) {
                        getDrawerVisitor().getComponent().addMouseMotionListener(mouseMotionListener);
                    }
                }
            }
            pressedButtons++;
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            if (pressedButtons > 0) {
                pressedButtons--;
            }

            if (pressedButtons == 0) {
                getDrawerVisitor().getComponent().removeMouseMotionListener(mouseMotionListener);
                currentAxes = null;
            }
        }
    }

    /**
     * This {@see MouseWheelListener} manage zoom/un-zoom on the figure.
     */
    private class FigureMouseWheelListener implements MouseWheelListener {

        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
            Axes axes = getUnderlyingAxes(e.getPoint());
            if (axes != null) {
                double scale = Math.pow(ZOOM_FACTOR, e.getUnitsToScroll());
                Double[] bounds = axes.getDisplayedBounds();
                double[][] factors = axes.getScaleTranslateFactors();

                double xDelta = (bounds[1] - bounds[0]) / 2;
                double xMiddle = (bounds[1] + bounds[0]) / 2;
                bounds[0] = xMiddle - xDelta * scale;
                bounds[1] = xMiddle + xDelta * scale;

                double yDelta = (bounds[3] - bounds[2]) / 2;
                double yMiddle = (bounds[3] + bounds[2]) / 2;
                bounds[2] = yMiddle - yDelta * scale;
                bounds[3] = yMiddle + yDelta * scale;

                double zDelta = (bounds[5] - bounds[4]) / 2;
                double zMiddle = (bounds[5] + bounds[4]) / 2;
                bounds[4] = zMiddle - zDelta * scale;
                bounds[5] = zMiddle + zDelta * scale;

                bounds[0] = bounds[0] * factors[0][0] + factors[1][0];
                bounds[1] = bounds[1] * factors[0][0] + factors[1][0];
                bounds[2] = bounds[2] * factors[0][1] + factors[1][1];
                bounds[3] = bounds[3] * factors[0][1] + factors[1][1];
                bounds[4] = bounds[4] * factors[0][2] + factors[1][2];
                bounds[5] = bounds[5] * factors[0][2] + factors[1][2];

                Boolean zoomed = tightZoomBounds(axes, bounds);

                bounds[0] = (bounds[0] - factors[1][0]) / factors[0][0];
                bounds[1] = (bounds[1] - factors[1][0]) / factors[0][0];
                bounds[2] = (bounds[2] - factors[1][1]) / factors[0][1];
                bounds[3] = (bounds[3] - factors[1][1]) / factors[0][1];
                bounds[4] = (bounds[4] - factors[1][2]) / factors[0][2];
                bounds[5] = (bounds[5] - factors[1][2]) / factors[0][2];

                boolean[] logFlags = { axes.getXAxisLogFlag(), axes.getYAxisLogFlag(), axes.getZAxisLogFlag()};
                ScaleUtils.applyInverseLogScaleToBounds(bounds, logFlags);

                GraphicController.getController().setProperty(axes.getIdentifier(), GraphicObjectProperties.__GO_ZOOM_BOX__, bounds);
                GraphicController.getController().setProperty(axes.getIdentifier(), GraphicObjectProperties.__GO_ZOOM_ENABLED__, zoomed);

            }
        }
    }

    private static void applyUnlog(Double[] bounds, Axes axes) {
        if (axes.getXAxisLogFlag()) {
            bounds[0] = Math.pow(10, bounds[0]);
            bounds[1] = Math.pow(10, bounds[1]);
        }

        if (axes.getYAxisLogFlag()) {
            bounds[2] = Math.pow(10, bounds[2]);
            bounds[3] = Math.pow(10, bounds[3]);
        }

        if (axes.getZAxisLogFlag()) {
            bounds[4] = Math.pow(10, bounds[4]);
            bounds[5] = Math.pow(10, bounds[5]);
        }
    }

    /**
     * This {@see MouseMotionListener} manage rotation and translation on the figure.
     */
    private class FigureMouseMotionListener extends MouseMotionAdapter implements MouseMotionListener {

    }
}
