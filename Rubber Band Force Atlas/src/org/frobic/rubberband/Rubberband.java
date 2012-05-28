/*
Copyright 2008-2011 Gephi
Authors : Mathieu Jacomy <mathieu.jacomy@gmail.com>
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
package org.frobic.rubberband;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.gephi.data.attributes.type.TimeInterval;
import org.gephi.data.attributes.type.StringList ;
import org.gephi.dynamic.DynamicUtilities;
import org.gephi.dynamic.api.DynamicController;
import org.gephi.dynamic.api.DynamicModel;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.HierarchicalGraph;
import org.gephi.graph.api.Node;
import org.gephi.graph.api.NodeData;
import org.gephi.graph.api.EdgeData;
import org.frobic.rubberband.ForceFactory.AttractionForce;
import org.frobic.rubberband.ForceFactory.RepulsionForce;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutProperty;
import org.gephi.project.api.Workspace;
import org.openide.util.Exceptions;
import org.openide.util.Lookup;
import org.openide.util.NbBundle;
import java.lang.String;
import java.util.AbstractList;
import java.util.ArrayList;

/**
 * ForceAtlas 2 Layout, manages each step of the computations.
 * @author Mathieu Jacomy
 */
public class Rubberband implements Layout {

    private GraphModel graphModel;
    private HierarchicalGraph graph;
    private final RubberbandBuilder layoutBuilder;
    private DynamicModel dynamicModel;
    private double edgeWeightInfluence;
    private double jitterTolerance;
    private double scalingRatio;
    private double gravity;
    private double factorf;
	private double alpha;

    private double speed;
    private boolean adjustSizes;
	private boolean disableString;
    private boolean barnesHutOptimize;
    private double barnesHutTheta;
    private boolean linLogMode;
    private boolean strongGravityMode;
    private int threadCount;
    private int currentThreadCount;
    private Region rootRegion;
    double outboundAttCompensation = 1;
	ArrayList<ArrayList<Integer>> Comm = new ArrayList<ArrayList<Integer>>();
	boolean Rubberband_First = true;
	int commBegin = -1 ;
    //Dynamic Weight
    private TimeInterval timeInterval;
    private ExecutorService pool;

    public Rubberband(RubberbandBuilder layoutBuilder) {
        this.layoutBuilder = layoutBuilder;
        this.threadCount = Math.min(4, Math.max(1, Runtime.getRuntime().availableProcessors() - 1));
    }

    @Override
    public void initAlgo() {
		
        speed = 1.;

        graph = graphModel.getHierarchicalGraphVisible();
        this.timeInterval = DynamicUtilities.getVisibleInterval(dynamicModel);	
		
        Node[] nodes = graph.getNodes().toArray();
		Edge[] edges = graph.getEdgesAndMetaEdges().toArray();
		
		for (Node n : nodes) {
			if (n.getNodeData().getAttributes().getValue("Comm") != null) {
				StringList temp = (StringList) n.getNodeData().getAttributes().getValue("Comm") ;
				int[] Communautes ;
				Communautes = new int[temp.size()] ;
				for (int i = 0 ; i < temp.size() ; i++) {
					Communautes[i] = Integer.parseInt(temp.getItem(i)) ; 
				}
			}
		}
		int NbCommunautes = 0 ;

        // Initialise layout data
        for (Node n : nodes) {
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof RubberbandLayoutData)) {
                RubberbandLayoutData nLayout = new RubberbandLayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            RubberbandLayoutData nLayout = nData.getLayoutData();
            nLayout.mass = 1 + graph.getDegree(n);
            nLayout.old_dx = 0;
            nLayout.old_dy = 0;
            nLayout.dx = 0;
            nLayout.dy = 0;
			if (n.getNodeData().getAttributes().getValue("Comm") != null) {
				StringList temp = (StringList) n.getNodeData().getAttributes().getValue("Comm") ;
				for (int _i = 0 ; _i < temp.size() ; _i++) {
					int c = Integer.parseInt(temp.getItem(_i)) ;
					if (c > NbCommunautes) {
						NbCommunautes = c ;
					}
				}
			}
        }
		NbCommunautes = NbCommunautes + 1 ;
		Comm.clear() ;
		for (int i = 0 ; i < NbCommunautes ; i++) {
			Comm.add(new ArrayList<Integer>()) ;
		}
		
		int i_n = 0 ;
		
		for (Node n : graph.getNodes()) {
			if (n.getNodeData().getAttributes().getValue("Comm") != null) {
				StringList temp = (StringList) n.getNodeData().getAttributes().getValue("Comm") ;
				for (int _i = 0 ; _i < temp.size() ; _i++) {
					int c = Integer.parseInt(temp.getItem(_i)) ;
					Comm.get(c).add(n.getId());
				}
				i_n++ ;
			}
		}
		
		
		if (Rubberband_First) {
			
			Rubberband_First = false ;
			for (int i = 0 ; i < NbCommunautes ; i++) {
				Node test = graphModel.factory().newNode("Comm"+i) ;
				test.getNodeData().setLabel("Rubberband_Virtual_Node");
				test.getNodeData().setSize(0.1f) ;
				test.getNodeData().setColor(1f,1f,1f) ;
				RubberbandLayoutData nLayout = new RubberbandLayoutData();
				nLayout.mass = 0 ;
				nLayout.old_dx = 0;
				nLayout.old_dy = 0;
				nLayout.dx = 0;
				nLayout.dy = 0;
				test.getNodeData().setLayoutData(nLayout) ;
				graph.addNode(test) ;
				if (commBegin == -1) {
					commBegin = test.getId() ;
				}
				//System.err.println(test.getId()) ;
			}
			
			
			
			
			for (Node n : nodes) {
				if (n.getNodeData().getAttributes().getValue("Comm") != null) {
					StringList temp = (StringList) n.getNodeData().getAttributes().getValue("Comm") ;
					for (int _i = 0 ; _i < temp.size() ; _i++) {
						int c = Integer.parseInt(temp.getItem(_i)) ;
						Edge ficelle = graphModel.factory().newEdge(graph.getNode("Comm"+c),n) ;
						ficelle.getEdgeData().setColor(1f,1f,1f) ;
						ficelle.getEdgeData().setLabel("ficelle");
						graph.addEdge(ficelle) ;
					}
				}
			}
			
			
			float[] cohecomm = new float[NbCommunautes] ;
			for (int i = 0 ; i < NbCommunautes ; i++) {
				double cohe = cohesion(i) ;
				if (cohe < 0.0001) {
					cohe = 0.0001 ;
				}
				cohecomm[i] = (float) (Math.sqrt(Comm.get(i).size())*2.*(1.-Math.log(cohe))) ;
			}
			
			for (Edge e : graph.getEdgesAndMetaEdges()) {
				String verif = "ficelle" ;
				EdgeData ed = e.getEdgeData() ;
				String edLab = ed.getLabel() ;
				if (verif.equals(edLab)) {
					String s = e.getSource().getNodeData().getId().substring(4) ;
					int c = Integer.parseInt(s) ;
					
					e.setWeight(cohecomm[c]) ;
					
				}
			}
			
			
			
		}
		

		
		
		
		//double nc = cohesion(0) ;

        pool = Executors.newFixedThreadPool(threadCount);
        currentThreadCount = threadCount;
    }

    @Override
    public void goAlgo() {
        // Initialize graph data
        if (graphModel == null) {
            return;
        }
        graph = graphModel.getHierarchicalGraphVisible();
        this.timeInterval = DynamicUtilities.getVisibleInterval(dynamicModel);
			

        Node[] nodes = graph.getNodes().toArray();
        Edge[] edges = graph.getEdgesAndMetaEdges().toArray();
		
		int NbCommunautes = 0 ;
		
		for (Node n : graph.getNodes()) {
			if (n.getNodeData().getAttributes().getValue("Comm") != null) {
				StringList temp = (StringList) n.getNodeData().getAttributes().getValue("Comm") ;
				for (int _i = 0 ; _i < temp.size() ; _i++) {
					int c = Integer.parseInt(temp.getItem(_i)) ;
					if (c > NbCommunautes) {
						NbCommunautes = c ;
					}
				}
			}
		}

        // Initialise layout data
        for (Node n : nodes) {
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof RubberbandLayoutData)) {
                RubberbandLayoutData nLayout = new RubberbandLayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            RubberbandLayoutData nLayout = nData.getLayoutData();
            nLayout.mass = 1 + graph.getDegree(n);
            nLayout.old_dx = nLayout.dx;
            nLayout.old_dy = nLayout.dy;
            nLayout.dx = 0;
            nLayout.dy = 0;
			if (n.getNodeData().getLabel().startsWith("Rubberband")) {
				nLayout.mass = 0.1 ;
			}
        }

        // If Barnes Hut active, initialize root region
        if (isBarnesHutOptimize()) {
            rootRegion = new Region(nodes);
            rootRegion.buildSubRegions();
        }

        // If outboundAttractionDistribution active, compensate.

        // Repulsion (and gravity)
        // NB: Muti-threaded
        RepulsionForce Repulsion = ForceFactory.builder.buildRepulsion(isAdjustSizes(), getScalingRatio());

        int taskCount = 8 * currentThreadCount;  // The threadPool Executor Service will manage the fetching of tasks and threads.
        // We make more tasks than threads because some tasks may need more time to compute.
        ArrayList<Future> threads = new ArrayList();
        for (int t = taskCount; t > 0; t--) {
            int from = (int) Math.floor(nodes.length * (t - 1) / taskCount);
            int to = (int) Math.floor(nodes.length * t / taskCount);
            Future future = pool.submit(new NodesThread(nodes, from, to, isBarnesHutOptimize(), getBarnesHutTheta(), getGravity(), (isStrongGravityMode()) ? (ForceFactory.builder.getStrongGravity(getScalingRatio())) : (Repulsion), getScalingRatio(), rootRegion, Repulsion));
            threads.add(future);
        }
        for (Future future : threads) {
            try {
                future.get();
            } catch (InterruptedException ex) {
                Exceptions.printStackTrace(ex);
            } catch (ExecutionException ex) {
                Exceptions.printStackTrace(ex);
            }
        }

        // Attraction
        AttractionForce Attraction = ForceFactory.builder.buildAttraction(isLinLogMode(), false, isAdjustSizes(), 1);
        if (getEdgeWeightInfluence() == 0) {
            for (Edge e : edges) {
                Attraction.apply(e.getSource(), e.getTarget(), 1);
            }
        } else {
            for (Edge e : edges) {
				String verif = "ficelle" ;
				EdgeData ed = e.getEdgeData() ;
				String edLab = ed.getLabel() ;
				double adjsize = getfactorf() ;
				if (isAdjustSizes()) {
					if (adjsize < 1.4) {
						adjsize = 1.4*adjsize ;
					}
				}
				if (verif.equals(edLab)) {
					Node n = e.getSource() ;
					Node n2 = e.getTarget() ;
					NodeData nData = n.getNodeData();
					RubberbandLayoutData nLayout = nData.getLayoutData();
					NodeData nData2 = n2.getNodeData();
					RubberbandLayoutData nLayout2 = nData2.getLayoutData();
					if ((nData.x()-nData2.x())*(nData.x()-nData2.x()) + (nData.y()-nData2.y())*(nData.y()-nData2.y()) > getWeight(e)*getWeight(e)*adjsize*adjsize) {
						double coeff = Math.sqrt((nData.x()-nData2.x())*(nData.x()-nData2.x()) + (nData.y()-nData2.y())*(nData.y()-nData2.y())) - getWeight(e)*adjsize ;
						coeff = coeff/getWeight(e) ;
						if (coeff > 75) {
							coeff = 75+Math.sqrt(coeff-75) ;
							if (coeff > 100) {
								coeff = 100 ;
							}
						}
						if (isDisableString()) {
							
						}
						else {
							if (isAdjustSizes()) {
								Attraction.apply(e.getSource(), e.getTarget(), (5.*(1+coeff)));
							}
							else {
								Attraction.apply(e.getSource(), e.getTarget(), (10.*(1+coeff)));
							}
						}
					}
				}
				else {
					double a = getalpha() ;
					if(getEdgeWeightInfluence() == 1) {
						Attraction.apply(e.getSource(), e.getTarget(),1 - a + a * getWeight(e));
					}
					else {
						Attraction.apply(e.getSource(), e.getTarget(), Math.pow(1 - a + a * getWeight(e),getEdgeWeightInfluence()));
					}
				}
            }
        }

        // Auto adjust speed
        double totalSwinging = 0d;  // How much irregular movement
        double totalEffectiveTraction = 0d;  // Hom much useful movement
        for (Node n : nodes) {
            NodeData nData = n.getNodeData();
            RubberbandLayoutData nLayout = nData.getLayoutData();
            if (!nData.isFixed()) {
                double swinging = Math.sqrt(Math.pow(nLayout.old_dx - nLayout.dx, 2) + Math.pow(nLayout.old_dy - nLayout.dy, 2));
                totalSwinging += nLayout.mass * swinging;   // If the node has a burst change of direction, then it's not converging.
                totalEffectiveTraction += nLayout.mass * 0.5 * Math.sqrt(Math.pow(nLayout.old_dx + nLayout.dx, 2) + Math.pow(nLayout.old_dy + nLayout.dy, 2));
            }
        }
        // We want that swingingMovement < tolerance * convergenceMovement
        double targetSpeed = getJitterTolerance() * getJitterTolerance() * totalEffectiveTraction / totalSwinging;

        // But the speed shoudn't rise too much too quickly, since it would make the convergence drop dramatically.
        double maxRise = 0.4;   // Max rise: 50%
        speed = speed + Math.min(targetSpeed - speed, maxRise * speed);

        // Apply forces
        if (isAdjustSizes()) {
            // If nodes overlap prevention is active, it's not possible to trust the swinging mesure.
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                RubberbandLayoutData nLayout = nData.getLayoutData();
                if (!nData.isFixed()) {

                    // Adaptive auto-speed: the speed of each node is lowered
                    // when the node swings.
                    double swinging = Math.sqrt((nLayout.old_dx - nLayout.dx) * (nLayout.old_dx - nLayout.dx) + (nLayout.old_dy - nLayout.dy) * (nLayout.old_dy - nLayout.dy));
                    double factor = 0.05 * speed / (1f + speed * Math.sqrt(swinging));

                    double df = Math.sqrt(Math.pow(nLayout.dx, 2) + Math.pow(nLayout.dy, 2));
                    factor = Math.min(factor * df, 10.) / df;

                    double x = nData.x() + nLayout.dx * factor;
                    double y = nData.y() + nLayout.dy * factor;

                    nData.setX((float) x);
                    nData.setY((float) y);
                }
            }
        } else {
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                RubberbandLayoutData nLayout = nData.getLayoutData();
                if (!nData.isFixed()) {

                    // Adaptive auto-speed: the speed of each node is lowered
                    // when the node swings.
                    double swinging = Math.sqrt((nLayout.old_dx - nLayout.dx) * (nLayout.old_dx - nLayout.dx) + (nLayout.old_dy - nLayout.dy) * (nLayout.old_dy - nLayout.dy));
                    //double factor = speed / (1f + Math.sqrt(speed * swinging));
                    double factor = speed / (1f + speed * Math.sqrt(swinging));

                    double x = nData.x() + nLayout.dx * factor;
                    double y = nData.y() + nLayout.dy * factor;

                    nData.setX((float) x);
                    nData.setY((float) y);
                }
            }
        }
    }

    @Override
    public boolean canAlgo() {
        return graphModel != null;
    }

    @Override
    public void endAlgo() {
		for (Node n : graph.getNodes()) {
			String verif = "Rubberband_Virtual_Node" ;
			NodeData ed = n.getNodeData() ;
			String edLab = ed.getLabel() ;
			if (verif.equals(edLab)) {
				//graph.writeUnlock() ;
				//graph.removeNode(n) ;
			}
        }
		
		for (Node n : graph.getNodes()) {
            n.getNodeData().setLayoutData(null);
		}
		
		pool.shutdown();


    }

    @Override
    public LayoutProperty[] getProperties() {
        List<LayoutProperty> properties = new ArrayList<LayoutProperty>();
        final String FORCEATLAS2_TUNING = NbBundle.getMessage(getClass(), "Rubberband.tuning");
        final String FORCEATLAS2_BEHAVIOR = NbBundle.getMessage(getClass(), "Rubberband.behavior");
        final String FORCEATLAS2_PERFORMANCE = NbBundle.getMessage(getClass(), "Rubberband.performance");
        final String FORCEATLAS2_THREADS = NbBundle.getMessage(getClass(), "Rubberband.threads");

        try {
            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "Rubberband.scalingRatio.name"),
                    FORCEATLAS2_TUNING,
                    "Rubberband.scalingRatio.name",
                    NbBundle.getMessage(getClass(), "Rubberband.scalingRatio.desc"),
                    "getScalingRatio", "setScalingRatio"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "Rubberband.strongGravityMode.name"),
                    FORCEATLAS2_TUNING,
                    "Rubberband.strongGravityMode.name",
                    NbBundle.getMessage(getClass(), "Rubberband.strongGravityMode.desc"),
                    "isStrongGravityMode", "setStrongGravityMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "Rubberband.gravity.name"),
                    FORCEATLAS2_TUNING,
                    "Rubberband.gravity.name",
                    NbBundle.getMessage(getClass(), "Rubberband.gravity.desc"),
                    "getGravity", "setGravity"));
			
			properties.add(LayoutProperty.createProperty(
														 this, Double.class,
														 NbBundle.getMessage(getClass(), "Rubberband.factorf.name"),
														 FORCEATLAS2_TUNING,
														 "Rubberband.factorf.name",
														 NbBundle.getMessage(getClass(), "Rubberband.factorf.desc"),
														 "getfactorf", "setfactorf"));
			
			properties.add(LayoutProperty.createProperty(
														 this, Double.class,
														 NbBundle.getMessage(getClass(), "Rubberband.alpha.name"),
														 FORCEATLAS2_TUNING,
														 "Rubberband.alpha.name",
														 NbBundle.getMessage(getClass(), "Rubberband.alpha.desc"),
														 "getalpha", "setalpha"));


            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "Rubberband.linLogMode.name"),
                    FORCEATLAS2_BEHAVIOR,
                    "Rubberband.linLogMode.name",
                    NbBundle.getMessage(getClass(), "Rubberband.linLogMode.desc"),
                    "isLinLogMode", "setLinLogMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "Rubberband.adjustSizes.name"),
                    FORCEATLAS2_BEHAVIOR,
                    "Rubberband.adjustSizes.name",
                    NbBundle.getMessage(getClass(), "Rubberband.adjustSizes.desc"),
                    "isAdjustSizes", "setAdjustSizes"));

			properties.add(LayoutProperty.createProperty(
														 this, Boolean.class,
														 NbBundle.getMessage(getClass(), "Rubberband.disableString.name"),
														 FORCEATLAS2_BEHAVIOR,
														 "Rubberband.adjustSizes.name",
														 NbBundle.getMessage(getClass(), "Rubberband.disableString.desc"),
														 "isDisableString", "setDisableString"));

			
            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "Rubberband.edgeWeightInfluence.name"),
                    FORCEATLAS2_BEHAVIOR,
                    "Rubberband.edgeWeightInfluence.name",
                    NbBundle.getMessage(getClass(), "Rubberband.edgeWeightInfluence.desc"),
                    "getEdgeWeightInfluence", "setEdgeWeightInfluence"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "Rubberband.jitterTolerance.name"),
                    FORCEATLAS2_PERFORMANCE,
                    "Rubberband.jitterTolerance.name",
                    NbBundle.getMessage(getClass(), "Rubberband.jitterTolerance.desc"),
                    "getJitterTolerance", "setJitterTolerance"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "Rubberband.barnesHutOptimization.name"),
                    FORCEATLAS2_PERFORMANCE,
                    "Rubberband.barnesHutOptimization.name",
                    NbBundle.getMessage(getClass(), "Rubberband.barnesHutOptimization.desc"),
                    "isBarnesHutOptimize", "setBarnesHutOptimize"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "Rubberband.barnesHutTheta.name"),
                    FORCEATLAS2_PERFORMANCE,
                    "Rubberband.barnesHutTheta.name",
                    NbBundle.getMessage(getClass(), "Rubberband.barnesHutTheta.desc"),
                    "getBarnesHutTheta", "setBarnesHutTheta"));

            properties.add(LayoutProperty.createProperty(
                    this, Integer.class,
                    NbBundle.getMessage(getClass(), "Rubberband.threads.name"),
                    FORCEATLAS2_THREADS,
                    "Rubberband.threads.name",
                    NbBundle.getMessage(getClass(), "Rubberband.threads.desc"),
                    "getThreadsCount", "setThreadsCount"));

        } catch (Exception e) {
            e.printStackTrace();
        }

        return properties.toArray(new LayoutProperty[0]);
    }

    @Override
    public void resetPropertiesValues() {
        int nodesCount = 0;

        if (graphModel != null) {
            nodesCount = graphModel.getHierarchicalGraphVisible().getNodeCount();
        }

        // Tuning
        if (nodesCount >= 100) {
            setScalingRatio(2.0);
        } else {
            setScalingRatio(10.0);
        }
        setStrongGravityMode(false);
        setGravity(1.);
		setfactorf(1.) ;

		setalpha(0.75) ;
		
        // Behavior
        setLinLogMode(false);
        setAdjustSizes(false);
		setDisableString(false);
        setEdgeWeightInfluence(1.);

        // Performance
        if (nodesCount >= 50000) {
            setJitterTolerance(10d);
        } else if (nodesCount >= 5000) {
            setJitterTolerance(1d);
        } else {
            setJitterTolerance(0.1d);
        }
        if (nodesCount >= 1000) {
            setBarnesHutOptimize(true);
        } else {
            setBarnesHutOptimize(false);
        }
        setBarnesHutTheta(1.2);
        setThreadsCount(16);
    }

    @Override
    public LayoutBuilder getBuilder() {
        return layoutBuilder;
    }

    @Override
    public void setGraphModel(GraphModel graphModel) {
        this.graphModel = graphModel;
        Workspace workspace = graphModel.getWorkspace();
        DynamicController dynamicController = Lookup.getDefault().lookup(DynamicController.class);
        if (dynamicController != null && workspace != null) {
            dynamicModel = dynamicController.getModel(workspace);
        }
        // Trick: reset here to take the profile of the graph in account for default values
        resetPropertiesValues();
    }

    public Double getBarnesHutTheta() {
        return barnesHutTheta;
    }

    public void setBarnesHutTheta(Double barnesHutTheta) {
        this.barnesHutTheta = barnesHutTheta;
    }

    public Double getEdgeWeightInfluence() {
        return edgeWeightInfluence;
    }

    public void setEdgeWeightInfluence(Double edgeWeightInfluence) {
        this.edgeWeightInfluence = edgeWeightInfluence;
    }

    public Double getJitterTolerance() {
        return jitterTolerance;
    }

    public void setJitterTolerance(Double jitterTolerance) {
        this.jitterTolerance = jitterTolerance;
    }

    public Boolean isLinLogMode() {
        return linLogMode;
    }

    public void setLinLogMode(Boolean linLogMode) {
        this.linLogMode = linLogMode;
    }

    public Double getScalingRatio() {
        return scalingRatio;
    }

    public void setScalingRatio(Double scalingRatio) {
        this.scalingRatio = scalingRatio;
    }

    public Boolean isStrongGravityMode() {
        return strongGravityMode;
    }

    public void setStrongGravityMode(Boolean strongGravityMode) {
        this.strongGravityMode = strongGravityMode;
    }

    public Double getGravity() {
        return gravity;
    }
	
	public Double getfactorf() {
        return factorf;
    }
	
	public Double getalpha() {
        return alpha;
    }

    public void setGravity(Double gravity) {
        this.gravity = gravity;
    }
	
	public void setfactorf(Double factorf) {
        this.factorf = factorf;
    }
	
	public void setalpha(Double alpha) {
        this.alpha = alpha;
    }

    public Integer getThreadsCount() {
        return threadCount;
    }

    public void setThreadsCount(Integer threadCount) {
        if (threadCount < 1) {
            setThreadsCount(1);
        } else {
            this.threadCount = threadCount;
        }
    }

    public Boolean isAdjustSizes() {
        return adjustSizes;
    }
	
	public Boolean isDisableString() {
        return disableString;
    }

    public void setAdjustSizes(Boolean adjustSizes) {
        this.adjustSizes = adjustSizes;
    }
	
	public void setDisableString(Boolean disableString) {
        this.disableString = disableString;
    }

    public Boolean isBarnesHutOptimize() {
        return barnesHutOptimize;
    }

    public void setBarnesHutOptimize(Boolean barnesHutOptimize) {
        this.barnesHutOptimize = barnesHutOptimize;
    }

    private float getWeight(Edge edge) {
        if (timeInterval != null) {
            return edge.getWeight(timeInterval.getLow(), timeInterval.getHigh());
        } else {
            return edge.getWeight();
        }
    }
	
	public boolean isneigh(int n1, int n2) {
		for (Node n : graph.getNeighbors(graph.getNode(n1))) {
			if (n.getId() == n2) {
				return true ;
			}
		}
		return false ;
	}
	
	public boolean incomm(int node, int comm) {
		StringList temp = (StringList) graph.getNode(node).getNodeData().getAttributes().getValue("Comm") ;
		if (temp != null) {
			for (int _i = 0 ; _i < temp.size() ; _i++) {
				int c = Integer.parseInt(temp.getItem(_i)) ;
				if (c == comm) {
					return true ;
				}
			}
		}
		return false ;
	}
	
	public double cohesion(int nc) {
		int Deltai = 0 ;
		int Deltao = 0 ;
		for (int i = 0 ; i < Comm.get(nc).size() ; i++) {
			int node = Comm.get(nc).get(i) ;
			for (int j = i+1 ; j < Comm.get(nc).size() ; j++) {
				int nodeb = Comm.get(nc).get(j) ;
				if (isneigh(node,nodeb)) {
					Node[] cohenode = graph.getNeighbors(graph.getNode(node)).toArray() ; 
					for (int k = 0 ; k < cohenode.length ; k++) {
						int nodec = cohenode[k].getId() ;
						if (isneigh(nodeb,nodec)) {
							if(incomm(nodec,nc)) {
								if (nodec > nodeb) {
									Deltai++ ;
								}
							}
							else {
								Deltao++ ;
							}
						}
					}
				}
			}
		}
		double density = (double)Deltai*6./( (double) (Comm.get(nc).size() * (Comm.get(nc).size() - 1) * (Comm.get(nc).size() - 2)) ) ;
		double isolation ;
		if (Deltai+Deltao > 0) {
			isolation = (double)Deltai / ((double) Deltai + (double) Deltao) ;
		}
		else {
			isolation = 0.0 ;
		}
		return density*isolation ;
		
	}
}
