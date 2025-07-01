package plugins.lagache.sodasuite;

import java.util.ArrayList;

import icy.sequence.Sequence;
import icy.type.point.Point3D;


public class Window3D {
	public ArrayList<Point3D> detectionlist;
	public Window3D(ArrayList<Point3D> liste){
		this.detectionlist = liste;}
	
	public void addToWindow(Point3D p){
		this.detectionlist.add(p);}

	
	public static Window3D[][][] window_tab( ArrayList<Point3D > detections, Sequence sequence, double r_max_xy, double r_max_z) {
		Window3D[][][] imageWindows = new Window3D[(int)(sequence.getWidth()/r_max_xy)+1][(int)(sequence.getHeight()/r_max_xy)+1][(int)(sequence.getSizeZ()/r_max_z)+1];
		for (int i=0;i<imageWindows.length;i++)
		{for (int j=0;j<imageWindows[i].length;j++){
			for (int k=0;k<imageWindows[i][j].length;k++){
			ArrayList<Point3D> liste = new ArrayList<Point3D>();
			imageWindows[i][j][k] = new Window3D(liste);}}}
		//ensuite on parcourt la liste des detections
		for (Point3D d:detections){
			//on récupère les indices de la sous-fenêtre contenant la detection
			int i = Math.max(0,(int)(d.getX()/r_max_xy));int j = Math.max(0,(int)(d.getY()/r_max_xy));int k = Math.max(0,(int)(d.getZ()/r_max_z));
			if ((i<(int)(sequence.getWidth()/r_max_xy)+1)&&(j<(int)(sequence.getHeight()/r_max_xy)+1)&&(k<(int)(sequence.getSizeZ()/r_max_z)+1)){
			imageWindows[i][j][k].addToWindow(d);}
		}
		return imageWindows;	
		}
		
}
