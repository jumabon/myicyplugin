package plugins.lagache.sodasuite;

import java.util.ArrayList;

import icy.sequence.Sequence;
import icy.type.point.Point3D;


public class Window2D {
	public ArrayList<Point3D> detectionlist;
	public Window2D(ArrayList<Point3D> liste){
		this.detectionlist = liste;}
	
	public void addToWindow(Point3D p){
		this.detectionlist.add(p);}

	
	public static Window2D[][] window_tab( ArrayList<Point3D > detections, Sequence sequence, double r_max) {
		Window2D[][] imageWindows = new Window2D[(int)(sequence.getWidth()/r_max)+1][(int)(sequence.getHeight()/r_max)+1];
		for (int i=0;i<imageWindows.length;i++)
		{for (int j=0;j<imageWindows[i].length;j++){
			ArrayList<Point3D> liste = new ArrayList<Point3D>();
			imageWindows[i][j] = new Window2D(liste);}}
		//ensuite on parcourt la liste des detections
		for (Point3D d:detections){
			//on récupère les indices de la sous-fenêtre contenant la detection
			int i = (int)(d.getX()/r_max);int j = (int)(d.getY()/r_max);			
			imageWindows[i][j].addToWindow(d);
		}
		return imageWindows;	
		}
		
}
