package plugins.lagache.sodasuite;

import icy.image.IcyBufferedImage;
import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import icy.type.point.Point5D;
import plugins.adufour.vars.lang.VarSequence;

import java.awt.*;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Ripley2D {
	// Calcul de la fonction de correlation (Ripley)
	public static double[][] correlation(ROI roi, ArrayList<Point2D> spots, ArrayList<Point2D> spots2, double volume,
			int N, ArrayList<Double> distance, VarSequence sequence1, VarSequence sequence2) throws InterruptedException {

		// ici on peut caster en ROI2D
		BooleanMask2D ma = roi.getBooleanMask2D(0,(int) roi.getPosition5D().getT(), (int) roi.getPosition5D().getC(), true);
		Point[] tab_point = ma.getContourPoints();
		ArrayList<Point> polyg = new ArrayList<Point>();
		for (int i = 0; i < tab_point.length; i++)
			polyg.add(tab_point[i]);

		double result[][] = new double[N - 1][3];
		// results[][1]=K et results[][2]=moyenne distances results[][3]=moyenne
		// distances^2
		double delta_K[] = new double[N - 1];
		double distances_moyennes[] = new double[N - 1];
		double distances_2_moyennes[] = new double[N - 1];
		double compteur[] = new double[N - 1];

		// il faut definir le nb local (par ROI) de detections afin de
		// construire les boucles
		int nbdeta = spots.size();
		int nbdetb = spots2.size();
		double x_a, y_a, x_b, y_b;
		for (int p = 0; p < nbdeta; p++) {
			double weight = 1;
			// distance du point a au bord de la ROI (polygon)
			// Point pt1 = new Point(spots.get(p).mass_center.x,
			// spots.get(p).mass_center.y);
			// x_a = spots.get(p).getPosition5D().getX();
			// y_a = spots.get(p).getPosition5D().getY();
			x_a = spots.get(p).getX();
			y_a = spots.get(p).getY();
			double d = distance2Polygon(x_a, y_a, polyg);
			for (int p2 = 0; p2 < nbdetb; p2++) {
				x_b = spots2.get(p2).getX();
				y_b = spots2.get(p2).getY();
				double temp = Math.sqrt(Math.pow(x_a - x_b, 2) + Math.pow(y_a - y_b, 2));
				// Calcul du poids (Ripley)
				weight = 1;
				if (temp > d) {
					weight = 1 - (Math.acos(d / temp)) / Math.PI;
				}

				// calcul direct de la fonction de correlation de paire
				// delta K
				for (int l = 1; l < N; l++) {
					if ((temp < distance.get(l)) & (temp > distance.get(0))) {
						delta_K[l - 1] += (1 / weight) * volume / (nbdeta * nbdetb);
						compteur[l - 1] += 1;
						distances_moyennes[l - 1] += temp;
						distances_2_moyennes[l - 1] += Math.pow(temp, 2);
						break;
					}
				}
			}
		}

		for (int l = 0; l < N - 1; l++) {
			result[l][0] = delta_K[l];
			if (compteur[l] > 0) {
				result[l][1] = distances_moyennes[l] / compteur[l];
				result[l][2] = distances_2_moyennes[l] / compteur[l];
			} else {
				result[l][1] = 0;
				result[l][2] = 0;
			}
		}
		return result;
	}

	public static double[][] correlation_new(ArrayList<ROI> roi, Window2D[][] imageWindows_1, Window2D[][] imageWindows_2,
			double volume, int N, ArrayList<Double> distance, int nbdeta,
			int nbdetb) throws InterruptedException {

		// ici on peut caster en ROI2D
		BooleanMask2D ma = ROIUtil.getUnion(roi).getBooleanMask2D(0,(int)(ROIUtil.getUnion(roi).getPosition5D().getT()) , -1, true);
		Point[] tab_point = ma.getContourPoints();
		ArrayList<Point> polyg = new ArrayList<Point>();
		for (int i = 0; i < tab_point.length; i++)
			polyg.add(tab_point[i]);

		double result[][] = new double[N - 1][3];
		// results[][1]=K et results[][2]=moyenne distances results[][3]=moyenne
		// distances^2
		double delta_K[] = new double[N - 1];
		double distances_moyennes[] = new double[N - 1];
		double distances_2_moyennes[] = new double[N - 1];
		double compteur[] = new double[N - 1];

		// il faut definir le nb local (par ROI) de detections afin de
		// construire les boucles
		// on commence par boucler sur les window des detections 1
		for (int i = 0; i < imageWindows_1.length; i++) {
			for (int j = 0; j < imageWindows_1[i].length; j++) {
				for (Point3D pos_a : imageWindows_1[i][j].detectionlist) {
					// distance du point a au bord de la ROI (polygon)
					double d = distance2Polygon(pos_a.getX(), pos_a.getY(), polyg);
					// on boucle ensuite sur les points 2 autour de la Window
					// i,j
					for (int k = Math.max(i - 1, 0); k <= Math.min(i + 1, imageWindows_1.length - 1); k++) {
						for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1, imageWindows_1[i].length - 1); l++) {

							for (Point3D pos_b : imageWindows_2[k][l].detectionlist) {
								double temp = Math.sqrt(Math.pow(pos_a.getX() - pos_b.getX(), 2)
										+ Math.pow(pos_a.getY() - pos_b.getY(), 2));
								// Calcul du poids (Ripley)
								double weight = 1;
								if (temp > d) {
									weight = 1 - (Math.acos(d / temp)) / Math.PI;
								}
								// calcul direct de la fonction de correlation
								// de paire
								// delta K
								for (int m = 1; m < N; m++) {
									if ((temp < distance.get(m)) & (temp >= distance.get(0))) {
										delta_K[m - 1] += (1 / weight) * volume / (nbdeta * nbdetb);
										compteur[m - 1] += 1;
										distances_moyennes[m - 1] += temp;
										distances_2_moyennes[m - 1] += Math.pow(temp, 2);
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		for (int l = 0; l < N - 1; l++) {
			result[l][0] = delta_K[l];
			if (compteur[l] > 0) {
				result[l][1] = distances_moyennes[l] / compteur[l];
				result[l][2] = distances_2_moyennes[l] / compteur[l];
			} else {
				result[l][1] = 0;
				result[l][2] = 0;
			}
		}
		return result;
	}

	public static double[] variance_theo_delta_new(ArrayList<ROI> roi_t, Window2D[][] imageWindows_1, Window2D[][] imageWindows_2, double area,
			int nb_a, int nb_b, int N, ArrayList<Double> distance, int N_h, double[] results) throws InterruptedException {

		// la premiere colonne de result[N][1] contient la variance de K, la
		// deuxieme colonne contient quand a elle la variance de delta_K
		double result[] = new double[N - 1];
		// on peut caster en ROI2D
		ROI union = ROIUtil.getUnion(roi_t);
		BooleanMask2D ma = union.getBooleanMask2D(0,(int)(union.getPosition5D().getT()), -1, true);
		Point[] liste_p = ma.getContourPoints();
		ArrayList<Point> polyg = new ArrayList<Point>();
		for (int i = 0; i < liste_p.length; i++)
			polyg.add(liste_p[i]);

		// il faut definir le nb local (par ROI) de detections afin de
		// construire les boucles
		
		for (int k = 1; k < N; k++) {
			double distancek_1 = distance.get(k-1);
			double distancek = distance.get(k);
			
			double d2 = Math.pow(distancek_1, 2);
			double d2bis = Math.pow(distancek, 2);
			double e1 = Math.PI * Math.pow(distance.get(k - 1), 2);
			double e2 = Math.PI * Math.pow(distance.get(k), 2);
			

			double temp_A1 = 0;
			double temp_A2 = 0;
			double temp_A3 = 0;

			double sum_h_a = 0;
			double sum_h_a_bis = 0;

			double distance_ij;

			// il faut definir le nb local (par ROI) de detections afin de
			// construire les boucles
			// on commence par boucler sur les window des detections 1
			for (int i = 0; i < imageWindows_1.length; i++) {
				for (int j = 0; j < imageWindows_1[i].length; j++) {
					for (Point3D pos_a : imageWindows_1[i][j].detectionlist) {
						// distance du point a au bord de la ROI (polygon)
						double dist = distance2Polygon(pos_a.getX(), pos_a.getY(), polyg);
						//calcul du terme de correction de bord				
									if ((dist < distance.get(k)) & (dist >= distance.get(0))) {
										sum_h_a += results[(int) Math.ceil(N_h * dist / distancek)];
				} else {
					sum_h_a += 1;
				}

									if (k>1){
									if ((dist < distance.get(k-1)) & (dist > distance.get(0))) {
										sum_h_a_bis += results[(int) Math.ceil(N_h * dist / distancek_1)];
				} else {
					sum_h_a_bis += 1;
				}}

									//calcul du terme d'intéraction
									// on boucle  sur les points 1 autour de la Window 1
									
									for (int m = Math.max(i - 1, 0); m <= Math.min(i + 1, imageWindows_1.length - 1); m++) {
										for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1, imageWindows_1[i].length - 1); l++) {

											for (Point3D pos_b : imageWindows_1[m][l].detectionlist) {
					
					distance_ij = Math.sqrt(Math.pow(pos_a.getX() - pos_b.getX(), 2) + Math.pow(pos_a.getY() - pos_b.getY(), 2));
					if (distance_ij>0){
					if (distance_ij<2*distancek_1)
					{
						temp_A1+=2 * d2 * Math.acos(distance_ij / (2 * distancek_1))
								- 0.5 * distance_ij * Math.sqrt(4 * d2 - Math.pow(distance_ij, 2));
					}
					if (distance_ij<2*distancek)
					{
						temp_A2 += 2 * Math.pow(distancek, 2) * Math.acos(distance_ij / (2 * distancek))
								- 0.5 * distance_ij * Math.sqrt(4 * Math.pow(distancek, 2) - Math.pow(distance_ij, 2));
					}		
										
					
						
						if (distance_ij < distancek_1 + distancek) {
							if (distance_ij + distancek_1 < distancek) {
								temp_A3 += 2*Math.PI * d2;
							} else {
								temp_A3 += 2*(d2
										* Math.acos((Math.pow(distance_ij, 2) + d2 - d2bis)
												/ (2 * distance_ij * distancek_1))
										+ d2bis * Math.acos((Math.pow(distance_ij, 2) + d2bis - d2)
												/ (2 * distance_ij * distancek))								
										- 0.5 * Math.sqrt(((-distance_ij + distancek_1 + distancek)
												* (distance_ij - distancek_1 + distancek)
												* (distance_ij + distancek_1 - distancek)
												* (distance_ij + distancek_1 + distancek))));
							}}}}
											}
											}
									}
					}}
			
					double I2 = (temp_A1 + temp_A2 - temp_A3
							- (Math.pow(e1, 2) / area + Math.pow(e2, 2) / area - 2 * e1 * e2 / area) * (nb_a * (nb_a - 1)))
							* nb_b / area;

					double I1 = (e2 * sum_h_a - e1 * sum_h_a_bis - nb_a * Math.pow(e2 - e1, 2) / area) * nb_b / area;

					result[k - 1] = Math.pow(area / (nb_b * nb_a), 2) * (I1 + I2);
				}
				return result;

	}

	public static double distance2Polygon(double x_point, double y_point, ArrayList<Point> polyg) {
		double dist = Integer.MAX_VALUE;// roiPolyg.getPerimeter();

		int nbc = polyg.size();

		for (int i = 0; i < nbc; i++) {
			Point pt1 = polyg.get(i);
			double disttmp = Math.sqrt(Math.pow(x_point - pt1.x, 2) + Math.pow(y_point - pt1.y, 2));
			dist = Math.min(dist, disttmp);
		}

		return (dist);

	}

	public static Point3D getIntensityCenter(ROI roi, Sequence seq) throws InterruptedException {
		
		IcyBufferedImage image = seq.getImage(0, 0, 0);
		double x = 0, y = 0;

		final BooleanMask2D mask = roi.getBooleanMask2D(0,(int)(roi.getPosition5D().getT()) , -1, true);
		final boolean m[] = mask.mask;
		final int h = mask.bounds.height;
		final int w = mask.bounds.width;

		int off = 0;
		double total_weight = 0;
		final Point5D pos2d = roi.getPosition5D();
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				if (m[off++]) {
					double weight = image.getData((int) (pos2d.getX() + i), (int) (pos2d.getY() + j), 0);
					total_weight += weight;
					x += i * weight;
					y += j * weight;
				}
			}
		}
		if (total_weight==0.0)
		{return new Point3D.Double(pos2d.getX(), pos2d.getY(),0);}
		else{
		return new Point3D.Double(pos2d.getX() + (x / total_weight), pos2d.getY() + (y / total_weight),0);}
		}
	

	/*public static Point3D getMassCenter(ROI roi) {
		
		double x = 0, y = 0;

		final BooleanMask2D mask = roi.getBooleanMask2D(0,(int)(roi.getPosition5D().getT()) , -1, true);		
		final boolean m[] = mask.mask;
		final int h = mask.bounds.height;
		final int w = mask.bounds.width;

		int off = 0;
		double total_weight = 0;
		final Point5D pos2d = roi.getPosition5D();
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				if (m[off++]) {					
					total_weight += 1;
					x += i ;
					y += j ;
				}
			}
		}
		return new Point3D.Double(pos2d.getX() + (x / total_weight), pos2d.getY() + (y / total_weight),0);}	*/

	public static List safe(List other) {
		return other == null ? Collections.EMPTY_LIST : other;
	}
}
