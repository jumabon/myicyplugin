package plugins.lagache.sodasuite;

import java.util.ArrayList;
import icy.gui.dialog.MessageDialog;
import icy.plugin.abstract_.Plugin;
import icy.roi.ROI;
import icy.roi.ROI3D;
import icy.type.point.Point3D;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.VarDoubleArrayNative;
import plugins.adufour.vars.lang.VarROIArray;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import plugins.kernel.roi.roi3d.ROI3DArea;


	public class LocalizationsToPoint2D extends Plugin implements Block {
		//en entrée, on prend les coordonnés x,y
		VarDoubleArrayNative x1 = new VarDoubleArrayNative("Positions 1 - x", null);
		VarDoubleArrayNative y1 = new VarDoubleArrayNative("Positions 1 - y", null);		
		//en sortie, on crée une list de ROI
		VarROIArray positions_1 = new VarROIArray("Positions 1 (ROIs)",null);
		
	@Override
	public void run() {
	if (x1.getValue().length != y1.getValue().length)
	{MessageDialog.showDialog("x and y arrays have not the same dimension");return;}
	
	ArrayList<ROI2DPoint> liste = new ArrayList<ROI2DPoint>();
	for (int i=0;i<x1.getValue().length;i++)
	{
		ROI2DPoint pt = new ROI2DPoint(x1.getValue()[i], y1.getValue()[i]);
		pt.setT(0);
		liste.add(pt);
	}
	ROI2DPoint[] tab = new ROI2DPoint[liste.size()];
	for (int i=0;i<liste.size();i++)
	{tab[i]=liste.get(i);}
	positions_1.setValue(tab);
	}
		
	
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add("Positions (x)",x1);
		inputMap.add("Positions (y)",y1);		
	}

	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add("Positions (ROIs)",positions_1);
	}

}
