/*
Copyright 2008-2011 Gephi
Authors : Yudi Xue <yudi.xue@usask.ca>, Mathieu Bastian
Website : http://www.gephi.org

This file is part of Gephi.

DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.

Copyright 2011 Gephi Consortium. All rights reserved.

The contents of this file are subject to the terms of either the GNU
General Public License Version 3 only ("GPL") or the Common
Development and Distribution License("CDDL") (collectively, the
"License"). You may not use this file except in compliance with the
License. You can obtain a copy of the License at
http://gephi.org/about/legal/license-notice/
or /cddl-1.0.txt and /gpl-3.0.txt. See the License for the
specific language governing permissions and limitations under the
License.  When distributing the software, include this License Header
Notice in each file and include the License files at
/cddl-1.0.txt and /gpl-3.0.txt. If applicable, add the following below the
License Header, with the fields enclosed by brackets [] replaced by
your own identifying information:
"Portions Copyrighted [year] [name of copyright owner]"

If you wish your version of this file to be governed by only the CDDL
or only the GPL Version 3, indicate your decision by adding
"[Contributor] elects to include this software in this distribution
under the [CDDL or GPL Version 3] license." If you do not indicate a
single choice of license, a recipient has the option to distribute
your version of this file under either the CDDL, the GPL Version 3 or
to extend the choice of license to its licensees as provided above.
However, if you add GPL Version 3 code and therefore, elected the GPL
Version 3 license, then the option applies only if the new code is
made subject to such option by the copyright holder.

Contributor(s):

Portions Copyrighted 2011 Gephi Consortium.
 */
package org.frobic.colorednodes.renderers;

import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfGState;
import java.awt.Color;
import java.util.ArrayList;
import org.frobic.colorednodes.items.NodeItem;
import org.gephi.graph.api.Node;
import org.gephi.preview.api.*;
import org.gephi.preview.spi.Renderer;
import org.gephi.preview.types.DependantColor;
import org.openide.util.NbBundle;
import org.openide.util.lookup.ServiceProvider;
import org.w3c.dom.Element;
import processing.core.PGraphics;

/**
 *
 * @author Yudi Xue, Mathieu Bastian
 */
@ServiceProvider(service = Renderer.class, position = 300)
public class NewNodeRenderer extends org.gephi.preview.plugin.renderers.NodeRenderer {

//Default values
    @Override public void preProcess(PreviewModel previewModel) {
    }

    @Override public void render(Item item, RenderTarget target, PreviewProperties properties) {
        
        if (target instanceof ProcessingTarget) {
            renderProcessing(item, (ProcessingTarget) target, properties);
        } else if (target instanceof SVGTarget) {
            renderSVG(item, (SVGTarget) target, properties);
        } else if (target instanceof PDFTarget) {
            renderPDF(item, (PDFTarget) target, properties);
        }
    }

    @Override public void renderProcessing(Item item, ProcessingTarget target, PreviewProperties properties) {
        //Params
        Float x = item.getData(NodeItem.X);
        Float y = item.getData(NodeItem.Y);
		Float angle = item.getData(NodeItem.ANGLE) ;
        Float size = item.getData(NodeItem.SIZE);
        Color color = item.getData(NodeItem.COLOR);
		Integer nbcolors = item.getData(NodeItem.NBCOLOR);
		Color[] colors = item.getData(NodeItem.COLORS);
        Color borderColor = ((DependantColor) properties.getValue(PreviewProperty.NODE_BORDER_COLOR)).getColor(color);
        float borderSize = properties.getFloatValue(PreviewProperty.NODE_BORDER_WIDTH);
        int alpha = (int) ((properties.getFloatValue(PreviewProperty.NODE_OPACITY) / 100f) * 255f);
        if (alpha > 255) {
            alpha = 255;
        }

        //Graphics
        PGraphics graphics = target.getGraphics();
		
		if (borderSize > 0) {
			graphics.stroke(borderColor.getRed(), borderColor.getGreen(), borderColor.getBlue(), alpha);
			graphics.strokeWeight(borderSize);
		} else {
			graphics.noStroke();
		}

		for (int i = 0 ; i < nbcolors ; i++) {
			graphics.fill(colors[i].getRed(), colors[i].getGreen(), colors[i].getBlue(), alpha);
			if (size >= 0.5) {
				graphics.arc(x,y,size,size,6.2832f*(nbcolors-i-1)/nbcolors-angle,6.2832f*(nbcolors-i)/nbcolors-angle) ;
			}
		}
		
    }

    @Override public void renderSVG(Item item, SVGTarget target, PreviewProperties properties) {
        Node node = (Node) item.getSource();
        //Params
        Float x = item.getData(NodeItem.X);
        Float y = item.getData(NodeItem.Y);
        Float size = item.getData(NodeItem.SIZE);
        size /= 2f;
        Color color = item.getData(NodeItem.COLOR);
		Integer nbcolors = item.getData(NodeItem.NBCOLOR);
		Color[] colors = item.getData(NodeItem.COLORS);
        Color borderColor = ((DependantColor) properties.getValue(PreviewProperty.NODE_BORDER_COLOR)).getColor(color);
        float borderSize = properties.getFloatValue(PreviewProperty.NODE_BORDER_WIDTH);
        float alpha = properties.getIntValue(PreviewProperty.NODE_OPACITY) / 100f;
        if (alpha > 1) {
            alpha = 1;
        }

        Element nodeElem = target.createElement("circle");
        nodeElem.setAttribute("class", node.getNodeData().getId());
        nodeElem.setAttribute("cx", x.toString());
        nodeElem.setAttribute("cy", y.toString());
        nodeElem.setAttribute("r", size.toString());
        nodeElem.setAttribute("fill", target.toHexString(color));
        nodeElem.setAttribute("fill-opacity", "" + alpha);
        if (borderSize > 0) {
            nodeElem.setAttribute("stroke", target.toHexString(borderColor));
            nodeElem.setAttribute("stroke-width", new Float(borderSize * target.getScaleRatio()).toString());
            nodeElem.setAttribute("stroke-opacity", "" + alpha);
        }
        target.getTopElement(SVGTarget.TOP_NODES).appendChild(nodeElem);
    }

    @Override public void renderPDF(Item item, PDFTarget target, PreviewProperties properties) {
        Float x = item.getData(NodeItem.X);
        Float y = item.getData(NodeItem.Y);
        Float size = item.getData(NodeItem.SIZE);
		Float angle = item.getData(NodeItem.ANGLE) ;
        size /= 2f;
        Color color = item.getData(NodeItem.COLOR);
		Integer nbcolors = item.getData(NodeItem.NBCOLOR);
		Color[] colors = item.getData(NodeItem.COLORS);
        Color borderColor = ((DependantColor) properties.getValue(PreviewProperty.NODE_BORDER_COLOR)).getColor(color);
        float borderSize = properties.getFloatValue(PreviewProperty.NODE_BORDER_WIDTH);
        float alpha = properties.getIntValue(PreviewProperty.NODE_OPACITY) / 100f;

        PdfContentByte cb = target.getContentByte();
        cb.setRGBColorStroke(borderColor.getRed(), borderColor.getGreen(), borderColor.getBlue());
        cb.setLineWidth(borderSize);
        cb.setRGBColorFill(color.getRed(), color.getGreen(), color.getBlue());
        if (alpha < 1f) {
            cb.saveState();
            PdfGState gState = new PdfGState();
            gState.setFillOpacity(alpha);
            gState.setStrokeOpacity(alpha);
            cb.setGState(gState);
        }
		for (int i = 0 ; i < nbcolors ; i++) {
			cb.setRGBColorFill(colors[nbcolors-i-1].getRed(), colors[nbcolors-i-1].getGreen(), colors[nbcolors-i-1].getBlue());
			if (size >= 0.5) {
				cb.newPath() ;
				ArrayList ar = cb.bezierArc(x-size,-y+size,x+size,0-size-y,360f*(nbcolors-i-1)/nbcolors+(360f/6.28f)*angle,360f*(1)/nbcolors);
				//cb.arc(x-size,-y+size,x+size,0-size-y,360f*(nbcolors-i-1)/nbcolors+angle,360f*(1)/nbcolors) ;
				cb.moveTo(x,-y) ;
				float pt[] = (float [])ar.get(0);
				cb.moveTo(pt[0], pt[1]);
				for (int k = 0; k < ar.size(); ++k) {
					pt = (float [])ar.get(k);
					cb.curveTo(pt[2], pt[3], pt[4], pt[5], pt[6], pt[7]);
				}
				cb.lineTo(x, -y);
				//strokeAndFill();
				//				cb.ClosePathFillStroke();
			if (borderSize > 0) {
				cb.fill();
			} else {
				cb.fill();
			}
				
				if (borderSize > 0) {
					cb.circle(x, -y, size);
					cb.stroke();
				}
			}
		}
        if (alpha < 1f) {
            cb.restoreState();
        }
    }

    @Override public PreviewProperty[] getProperties() {
        return new PreviewProperty[]{
                    PreviewProperty.createProperty(this, PreviewProperty.NODE_BORDER_WIDTH, Float.class,
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.borderWidth.displayName"),
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.borderWidth.description"),
                    PreviewProperty.CATEGORY_NODES).setValue(defaultBorderWidth),
                    PreviewProperty.createProperty(this, PreviewProperty.NODE_BORDER_COLOR, DependantColor.class,
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.borderColor.displayName"),
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.borderColor.description"),
                    PreviewProperty.CATEGORY_NODES).setValue(defaultBorderColor),
                    PreviewProperty.createProperty(this, PreviewProperty.NODE_OPACITY, Float.class,
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.opacity.displayName"),
                    NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.property.opacity.description"),
                    PreviewProperty.CATEGORY_NODES).setValue(defaultOpacity)};
    }

    @Override public boolean isRendererForitem(Item item, PreviewProperties properties) {
        if (item instanceof NodeItem) {
            return true;
        }
        return false;
    }

    public String getName() {
        return NbBundle.getMessage(NewNodeRenderer.class, "NewNodeRenderer.name");
    }
	
}
