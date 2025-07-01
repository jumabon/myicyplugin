package plugins.lagache.sodasuite;

import java.awt.geom.Point2D;
import java.util.ArrayList;


import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROI3D;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import icy.type.point.Point5D;
import plugins.kernel.roi.descriptor.measure.ROIMassCenterDescriptorsPlugin;

// Colocalisation with Ripley function K
// Significant 

public class Methods_distance{

	static double[] results = null;
	static int N_h;
	// pour le test stat
	static double mindist;
	static double maxdist;
	static ArrayList<Double> distance_fit = new ArrayList<Double>();
	static ArrayList<ROI> list_roi = new ArrayList<ROI>();
	static double p_value = 0.0;
	static double log_p_value = 0.0;

	public static double[] distance(ROI[] detections1,ROI[]detections2,Sequence sequence,Sequence sequence2, ArrayList<ROI> roi_t,int t,double volume) throws InterruptedException {
		double[] p=new double[3];		
		if (roi_t==null){			
			return p;}
		
		ArrayList<ROI> detection = detectionsInRoi(detections1, roi_t);
		ArrayList<ROI> detection2 = detectionsInRoi(detections2, roi_t);		
		int nbdetb = detection2.size();
		
		ArrayList<ROI> liste1 = new ArrayList<ROI>();
		for (ROI r1:detection){
			if ((r1.getPosition5D().getT()==t)||(r1.getBounds5D().isInfiniteT())){
			liste1.add(r1);}}
		ROI union1=null;
		if (liste1.isEmpty()){}
		else
		{union1 = ROIUtil.getUnion(liste1);}
		
		if (union1==null)
		{return p;} 
		else {
		double vol1 = union1.getNumberOfPoints();
		double p1 = vol1/volume;
		double sigma1 = Math.sqrt(p1*(1-p1)/(double)nbdetb);
		
		int compteur1 = 0;
		for (ROI r2:detection2){				
			if (r2.getDimension()==3)
			{				
				Point3D pos=Ripley3D.getIntensityCenter((ROI3D)(r2),sequence2);
				if (union1.contains(pos.getX(), pos.getY(), pos.getZ(), union1.getPosition5D().getT(), union1.getPosition5D().getC()))
					compteur1++;
			}
			else
			{
			Point3D pos = Ripley2D.getIntensityCenter((ROI2D)(r2),sequence2);		
			if (union1.contains(pos.getX(), pos.getY(), union1.getPosition5D().getZ(), union1.getPosition5D().getT(), union1.getPosition5D().getC()))
				compteur1++;
			}
		}
		p[0] = (double)(compteur1)/(double)nbdetb;
		//calcul de la pvalue
		//calcul du coeff centr� r�duit
		double R_tilde=(p[0]-p1)/sigma1;
		double pvalue = 0.5*(1-ErrorFunction.erf(R_tilde/Math.sqrt(2)));
		if (R_tilde<4.0){			 
			p[1] = Math.log10(pvalue);			
		}
		else
		{double x = R_tilde/Math.sqrt(2);
		p[1]=-Math.pow(x, 2)/Math.log(10)-Math.log10(2*Math.sqrt(Math.PI)*x);}
		p[2]=pvalue;		
		return p;
		}
	}
	
	public static double[] soda(ROI[] detections1,ROI[]detections2,Sequence sequence, Sequence sequence2, double rmax,double step,ArrayList<ROI> roi_t,int t,double volume) throws InterruptedException {
		double[] retour = new double[4];
		// d�finition du ratio Z/X
				double ratio_zx = 1.0;
				if (sequence.getSizeZ() > 1) {
					ratio_zx = sequence.getPixelSizeZ() / sequence.getPixelSizeX();
				}
				// initialisation boucle en temps
				// temps total de la sequence				
				maxdist=rmax;
				betaCorrection(maxdist / 10, 100);
				/////////////////////////////////////////////
				///////////////////////////////////////////////////
				// et on filtre les d�tections appartenant effectivement � la ROI (t)								
				
				ArrayList<ROI> detection = detectionsInRoi(detections1, roi_t);
				ArrayList<ROI> detection2 = detectionsInRoi(detections2, roi_t);
	
				int nbdeta = detection.size();
				int nbdetb = detection2.size();
				// construction du distance_tab pour calculer la p_value de coloc
				double coeff_radius = 1;
				double step_min = step;// 1.0;//maxdist / 10;

				distance_fit.clear();				
				distance_fit.add(0.);
				double temp = 0.;				
					while (temp + step_min <= maxdist) {
						temp +=step_min;
						distance_fit.add(temp);
					}
				int N_fit = distance_fit.size();
				if (N_fit == 1) {
					distance_fit.add(maxdist);
					N_fit = distance_fit.size();
				}
				////////////////////////////////////////////////////////////////////////////////
			
				// stockage de l'histogram des distances de coloc
				double[] proba_dist = new double[N_fit - 1];				
				double window_size = maxdist;
				double max_z = maxdist/ratio_zx;
				
				// on crée des listes de points 2D ou 3D avec les positions (centres
				// d'intensité ou de masse en 3D) des spots			
				ArrayList<Point3D> positions_1_2D = new ArrayList<Point3D>();
				ArrayList<Point3D> positions_2_2D = new ArrayList<Point3D>();
				//on crée aussi une liste potentielle de points 3d
				ArrayList<Point3D> positions_1_3D = new ArrayList<Point3D>();
				ArrayList<Point3D> positions_2_3D = new ArrayList<Point3D>();
				//sous-listes des detections pour accelerer le temps de calcul
				Window2D[][] imageWindows_1_2D = new Window2D[(int) (sequence.getWidth()/ window_size)][(int) (sequence.getHeight() / window_size)];
				Window2D[][] imageWindows_2_2D = new Window2D[(int) (sequence.getWidth()/ window_size)][(int) (sequence.getHeight() / window_size)];
				Window3D[][][] imageWindows_1_3D = new Window3D[(int) (sequence.getWidth() / window_size)][(int) (sequence.getHeight()/ window_size)][(int) (sequence.getSizeZ() / max_z)];
				Window3D[][][] imageWindows_2_3D = new Window3D[(int) (sequence.getWidth() / window_size)][(int) (sequence.getHeight()/ window_size)][(int) (sequence.getSizeZ() / max_z)];
				/////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
				if (sequence.getSizeZ()>1){
					// on boucle sur detection puis sur detection2
					for (ROI r : detection) {
						if (r.getDimension()==3)
							{Point3D pos = Ripley3D.getIntensityCenter((ROI3D)r,sequence);positions_1_3D.add(pos);}
						else
							{Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, sequence);positions_1_3D.add(pos);}							
							}						
					for (ROI r : detection2) {
						if (r.getDimension()==3)
						{Point3D pos = Ripley3D.getIntensityCenter((ROI3D)r,sequence2);positions_2_3D.add(pos);}
					else
						{Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, sequence2);positions_2_3D.add(pos);}							
						}	
					// on crée ensuite toutes les sous-listes de detections
					imageWindows_1_3D = Window3D.window_tab(positions_1_3D, sequence, maxdist,max_z);				
					imageWindows_2_3D = Window3D.window_tab(positions_2_3D, sequence, maxdist,max_z);}				
					else{				
				for (ROI r : detection) {				
						Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, sequence);
						positions_1_2D.add(pos);}
				for (ROI r : detection2) {				
						Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, sequence2);
						positions_2_2D.add(pos);}											
				// on crée ensuite toutes les sous-listes de detections							
				imageWindows_1_2D = Window2D.window_tab(positions_1_2D, sequence, window_size);			
				imageWindows_2_2D = Window2D.window_tab(positions_2_2D, sequence2, window_size);
				}
				// calcul des parametres globaux: aire totale des rois et nb de detections
				// faut-refaire une fonction calcul de volume
				volume = ROIUtil.getUnion(roi_t).getNumberOfPoints() * ratio_zx;
				if (ROIUtil.getUnion(roi_t).getBounds5D().isInfiniteZ()) {
					volume = volume * sequence.getSizeZ() * ratio_zx;
				}				
				nbdeta = detection.size();
				nbdetb = detection2.size();
				String dataSetName = "unknown sequence";				
				//déclaration des tableaux de p-value
				double[] maximum = new double[N_fit-1];
				double[] pvalue = new double[N_fit-1];
				double[] log_pvalue = new double[N_fit-1];
				
				// esperance du nb de coloc par anneau
				double[] col_ev = new double[N_fit - 1];
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					double[] delta_K = new double[N_fit - 1];
					double[][] Ktemp = new double[N_fit - 1][3];
					////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////

					if (sequence.getSizeZ() > 1) {					
						Ktemp = Ripley3D.correlation_new(roi_t, imageWindows_1_3D,imageWindows_2_3D, volume, N_fit, distance_fit, sequence,nbdeta, nbdetb, ratio_zx);} 
					else {
						Ktemp = Ripley2D.correlation_new(roi_t, imageWindows_1_2D, imageWindows_2_2D, volume, N_fit, distance_fit, nbdeta, nbdetb);}

					for (int k = 0; k < N_fit - 1; k++) {
						delta_K[k] = Ktemp[k][0];
					}
					
					if (sequence.getSizeZ()>1)
					{
						non_parametric_object.main3D_corr(N_fit,sequence,roi_t,positions_1_3D ,positions_2_3D, imageWindows_1_3D,imageWindows_2_3D,distance_fit,N_h,ratio_zx,volume,results, delta_K, Ktemp,proba_dist,maximum,pvalue);					
					}
					// calcul de l'esperance du nb de couples par intervalle et du
					// nb total de couples
					else
					{
						non_parametric_object.main2D_corr(N_fit, sequence, roi_t, positions_1_2D, positions_2_2D, imageWindows_1_2D,imageWindows_2_2D, distance_fit, N_h, volume, results, delta_K, Ktemp, proba_dist,maximum,pvalue);
					}
					int total_coloc_temp =0;	
					for (int j = 0; j < N_fit - 1; j++) {
						col_ev[j] = proba_dist[j] * (nbdetb * nbdeta / volume) * delta_K[j];
						total_coloc_temp += col_ev[j];
						if (maximum[j]<4.){
							log_pvalue[j] = Math.log10(pvalue[j]);}
							else {
								double x = maximum[j] / Math.sqrt(2.);
								double y = (N_fit - 1) / (2 * Math.sqrt(Math.PI) * x);
								double exponent = (Math.log(y) - Math.pow(x, 2)) / Math.log(10);
								log_pvalue[j] = exponent;
							}
					}
				
						retour[0]=(double)total_coloc_temp/((double)nbdetb);
						
						retour[1]=log_pvalue[N_fit-2];
						//calcul de la distance moyenne de coloc
						double dist_moy=0;
						ArrayList<apparatedSpots> liste_app_detect = apparatedSpots.appDetectConstruction(detections1, detections2,proba_dist, distance_fit, sequence, sequence2);
						for (int j=0;j<N_fit-1;j++){
							dist_moy+=apparatedSpots.distance_moyenne(liste_app_detect, distance_fit.get(j+1))*col_ev[j];											
						}												
						retour[2] = dist_moy/total_coloc_temp;
						retour[3] = pvalue[N_fit-2];
						return retour;
					}
	

		private static ArrayList<ROI> detectionsInRoi(ROI[] detec, ArrayList<ROI> roi) throws InterruptedException {
				ArrayList<ROI> ROIDetection = new ArrayList<ROI>();
		ArrayList<ROI> TestList = new ArrayList<ROI>();

		ROI union = ROIUtil.getUnion(roi);
		// fill hashMap
		for (ROI r : detec) {
			boolean b = false;
			for (ROI rt : TestList) {
				if (r == rt) {
					b = true;
				}
			}
			if (b == false) // si le spot n'a pas encore �t� pris en compte
			{
				TestList.add(r);
				if (union.contains(r.getPosition5D().getX(), r.getPosition5D().getY(), union.getPosition5D().getZ(), r.getPosition5D().getT(), union.getPosition5D().getC())) {
					ROIDetection.add(r);
				}
			}
		}
		return (ROIDetection);
	}
		private static void betaCorrection(double pas, int nbN) {

			Double valN = 1 / (1 / (double) nbN);
			N_h = valN.intValue() + 1; // / ATTENTION AUX BORNES

			double[] alpha = new double[N_h + 1];
			results = new double[N_h + 1];

			for (int i = 0; i < results.length; i++) {
				results[i] = 0;
				alpha[i] = 0.0d;

			}

			for (int i = 1; i < results.length; i++) {
				alpha[i] = (i / (double) N_h);
			}
			for (int i = 0; i < results.length; i++) {
				double h, j;
				for (h = alpha[i] + (pas), j = 2; h <= 1; j++) {
					results[i] = results[i] + h * pas / (1 - 1 / Math.PI * Math.acos(alpha[i] / h));
					h = alpha[i] + (j * pas);
				}

				results[i] = results[i] * 2 + alpha[i] * alpha[i];
			}

		}
}