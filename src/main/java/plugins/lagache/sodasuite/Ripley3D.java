package plugins.lagache.sodasuite;

import icy.image.IcyBufferedImage;
import icy.roi.*;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import icy.type.point.Point5D;

import java.awt.*;
import java.util.ArrayList;

public class Ripley3D {
	  public static double[][] correlation(int dimZ,ROI roi,ArrayList<Point3D> spots,ArrayList<Point3D> spots2,double volume, int nb_a,int nb_b,int N,ArrayList<Double> distance,double ratio_zx) throws InterruptedException {
	double[][] result =new double[N-1][3];
	if (roi==null){return result;}
	  //results[][1]=K et results[][2]=moyenne distances results[][3]=moyenne distances^2 
	double delta_K[] = new double[N-1]; double distances_moyennes[]=new double[N-1]; double distances_2_moyennes[]=new double[N-1]; double compteur[]=new double[N-1];
	  
	  //on fait la boucle sur toutes les rois pour incr�menter le calcul de la fonction de Ripley 
	Point5D pt = roi.getPosition5D();
	  ArrayList<Point3D.Integer> polyg = new ArrayList<Point3D.Integer>(); 
	  //ilfaut extruder si n�cessaire la roi 
	  if (roi.getDimension()==2){ 
		  for (int z=0;z<dimZ;z++){ 
			  BooleanMask2D ma=roi.getBooleanMask2D(z,(int) pt.getT(),(int)pt.getC(), true); 
			  Point[] tab_point = ma.getContourPoints(); 
			  for (int i=0;i<tab_point.length;i++){ Point3D.Integer ptInteger = new Point3D.Integer(); ptInteger.setX((double)tab_point[i].x);ptInteger.setY((double)tab_point[i].y); ptInteger.setZ((double)z); polyg.add(ptInteger); } }	  
	  } else { BooleanMask3D ma = ((ROI3D)roi).getBooleanMask(true);
	  Point3D.Integer[] tab_point = ma.getContourPoints(); for (int i=0;i<tab_point.length;i++) polyg.add(tab_point[i]); }
	  
	  //il faut definir le nb local (par ROI) de detections afin de construire les boucles 
	  int nbdeta=spots.size(); int nbdetb=spots2.size(); double x_a,y_a,x_b,y_b,z_a,z_b; for (int p = 0; p < nbdeta; p++) { double weight = 1; 
	  // distance du point a au bord de la ROI (polygon) //	  
	  x_a=spots.get(p).getX();y_a=spots.get(p).getY();z_a=spots.get(p).getZ();	  
	  double d = distance2Polygon(x_a,y_a,z_a, polyg,ratio_zx);
	  
	  for (int p2 = 0; p2 < nbdetb; p2++) {
	  x_b=spots2.get(p2).getX();y_b=spots2.get(p2).getY();z_b=spots2.get(p2).getZ(); 
	  double temp = Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b,2)+Math.pow((z_a-z_b)*ratio_zx, 2)); 
	  // Calcul du poids (Ripley) 
	  weight =1; if( temp>d){ weight = 0.5+(double)(d/(2*temp)); } //calcul direct de
	  //la fonction de correlation de paire delta K 
	  for (int l = 1; l < N; l++) {
	  if ((temp < distance.get(l))&(temp > distance.get(0))) { delta_K[l-1]+=(1
	  / weight) * volume / (nb_a * nb_b); compteur[l-1]+=1;
	  distances_moyennes[l-1]+=temp;
	  distances_2_moyennes[l-1]+=Math.pow(temp,2); break; } } } }
	  for (int l=0; l<N-1; l++){ result[l][0]=delta_K[l]; if (compteur[l]>0){
	  result[l][1]=distances_moyennes[l]/compteur[l];
	  result[l][2]=distances_2_moyennes[l]/compteur[l];} else
	  {result[l][1]=0;result[l][2]=0;} } return result; }
	 
	public static double[][] correlation_new(ArrayList<ROI> roi, Window3D[][][] imageWindows_1,
			Window3D[][][] imageWindows_2, double volume, int N, ArrayList<Double> distance, Sequence sequence,
			int nbdeta, int nbdetb, double ratio_zx) throws InterruptedException {
		
		// ici on peut caster en ROI2D
		ArrayList<Point3D> liste_point=new ArrayList<Point3D>();
		for (int h=0;h<sequence.getSizeZ();h++){			
		BooleanMask2D ma = ROIUtil.getUnion(roi).getBooleanMask2D(h, (int)(ROIUtil.getUnion(roi).getPosition5D().getT()),(int)(ROIUtil.getUnion(roi).getPosition5D().getC()), true);
		Point[] tab_point = ma.getContourPoints();
		for (Point pt:tab_point)
		{Point3D pt3d = new Point3D.Double(pt.x, pt.y, h);liste_point.add(pt3d);}
		}
		Point3D[] tab_points=new Point3D[liste_point.size()];
		for (int i=0;i<liste_point.size();i++)
		{tab_points[i]=liste_point.get(i);}
		

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
				for (int z = 0; z < imageWindows_1[i][j].length; z++) {
					for (Point3D pos_a : imageWindows_1[i][j][z].detectionlist) {
						// distance du point a au bord de la ROI (polygon)
						double d = distance2Polygon3D(pos_a.getX(), pos_a.getY(), pos_a.getZ(), tab_points, ratio_zx);
						// on boucle ensuite sur les points 2 autour de la
						// Window
						// i,j
						for (int k = Math.max(i - 1, 0); k <= Math.min(i + 1, imageWindows_1.length - 1); k++) {
							for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1, imageWindows_1[i].length - 1); l++) {
								for (int m = Math.max(z - 1, 0); m <= Math.min(z + 1,
										imageWindows_1[i][j].length - 1); m++) {
									for (Point3D pos_b : imageWindows_2[k][l][m].detectionlist) {
										double temp = Math.sqrt(Math.pow(pos_a.getX() - pos_b.getX(), 2)
												+ Math.pow(pos_a.getY() - pos_b.getY(), 2)
												+ Math.pow((pos_a.getZ() - pos_b.getZ()) * ratio_zx, 2));
										// Calcul du poids (Ripley)
										double weight = 1;
										if (temp > d) {
											weight = 0.5 + (double) (d / (2 * temp));
										}
										// calcul direct de la fonction de
										// correlation
										// de paire
										// delta K
										for (int n = 1; n < N; n++) {
											if ((temp < distance.get(n)) & (temp > distance.get(0))) {
												delta_K[n - 1] += (1 / weight) * volume / (nbdeta * nbdetb);
												compteur[n - 1] += 1;
												distances_moyennes[n - 1] += temp;
												distances_2_moyennes[n - 1] += Math.pow(temp, 2);
												break;
											}
										}
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


	public static double[] variance_theo_delta_new(ArrayList<ROI> roi, Window3D[][][] imageWindows_1,
			Window3D[][][] imageWindows_2, double volume, int nb_a, int nb_b, int N, ArrayList<Double> distance,
			 Sequence sequence, double ratio_zx) throws InterruptedException {

		// la premiere colonne de result[N][1] contient la variance de K, la
		// deuxieme colonne contient quand a elle la variance de delta_K
		double result[] = new double[N - 1];		
		// ici on peut caster en ROI2D
		ArrayList<Point3D> liste_point=new ArrayList<Point3D>();
		for (int h=0;h<sequence.getSizeZ();h++){			
		BooleanMask2D ma = ROIUtil.getUnion(roi).getBooleanMask2D(h, (int)(ROIUtil.getUnion(roi).getPosition5D().getT()),(int)(ROIUtil.getUnion(roi).getPosition5D().getC()), true);
		Point[] tab_point = ma.getContourPoints();
		for (Point pt:tab_point)
		{Point3D pt3d = new Point3D.Double(pt.x, pt.y, h);liste_point.add(pt3d);}
		}
		Point3D[] tab_points=new Point3D[liste_point.size()];
		for (int i=0;i<liste_point.size();i++)
		{tab_points[i]=liste_point.get(i);}
		
		// il faut definir le nb local (par ROI) de detections afin de
		// construire les boucles

		for (int k = 1; k < N; k++) {
			double distancek_1 = distance.get(k - 1);
			double distancek = distance.get(k);

			
			double d3 = Math.pow(distancek_1, 3);
			double e1 = 4 / 3 * Math.PI * Math.pow(distance.get(k - 1), 3);
			double e2 = 4 / 3 * Math.PI * Math.pow(distance.get(k), 3);

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
					for (int z = 0; z < imageWindows_1[i][j].length; z++) {
						for (Point3D pos_a : imageWindows_1[i][j][z].detectionlist) {
							// distance du point a au bord de la ROI (polygon)
							double dist = distance2Polygon3D(pos_a.getX(), pos_a.getY(), pos_a.getZ(), tab_points,
									ratio_zx);
							// calcul du terme de correction de bord
							if ((dist < distance.get(k)) & (dist > distance.get(0))) {
								double alpha_ = (double) (dist / distance.get(k));
								sum_h_a += (-4) * Math.pow(alpha_, 3) + 6 * Math.pow(alpha_, 2) - 3 * alpha_ + 2 + 6
										* Math.pow(alpha_, 3) * (alpha_ * Math.log(2 * alpha_) - Math.log(1 + alpha_));
							} else {
								sum_h_a += 1.0;
							}
							if (k > 1) {
								if ((dist < distance.get(k - 1)) & (dist > distance.get(0))) {
									double alpha_ = (double) (dist / distance.get(k - 1));
									sum_h_a_bis += (-4) * Math.pow(alpha_, 3) + 6 * Math.pow(alpha_, 2) - 3 * alpha_ + 2
											+ 6 * Math.pow(alpha_, 3)
													* (alpha_ * Math.log(2 * alpha_) - Math.log(1 + alpha_));
								} else {
									sum_h_a_bis += 1.0;
								}
							}

							// calcul du terme d'intéraction
							// on boucle sur les points 1 autour de la Window 1
							for (int m = Math.max(i - 1, 0); m <= Math.min(i + 1, imageWindows_1.length - 1); m++) {
								for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1,
										imageWindows_1[i].length - 1); l++) {
									for (int n = Math.max(z - 1, 0); n <= Math.min(z + 1,
											imageWindows_1[i][j].length - 1); n++) {
										for (Point3D pos_b : imageWindows_1[m][l][n].detectionlist) {
											distance_ij = Math.sqrt(Math.pow(pos_a.getX() - pos_b.getX(), 2)
													+ Math.pow(pos_a.getY() - pos_b.getY(), 2)
													+ Math.pow((pos_a.getZ() - pos_b.getZ()) * ratio_zx, 2));
											if (distance_ij > 0) {
												if (distance_ij < 2 * distancek_1) {
													temp_A1 += (Math.PI / 12) * (4 * distancek_1 + distance_ij)
															* Math.pow(2 * distancek_1 - distance_ij, 2);
												}
												if (distance_ij < 2 * distancek) {
													temp_A2 += (Math.PI / 12) * (4 * distancek + distance_ij)
															* Math.pow(2 * distancek - distance_ij, 2);
												}

												if (distance_ij < distancek_1 + distancek) {
													if (distance_ij + distancek_1 < distancek) {
														temp_A3 += 2 * 4 / 3 * Math.PI * d3;
													} else {
														temp_A3 += (Math.PI/(12*distance_ij))*Math.pow(distancek+distancek_1-distance_ij,2)*(Math.pow(distance_ij, 2)+2*distance_ij*distancek_1-3*Math.pow(distancek_1, 2)+2*distance_ij*distancek+6*distancek*distancek_1-3*Math.pow(distancek, 2));
														//2 * (Math.PI / (12 * distance_ij))* (Math.pow(distance_ij, 2)+ 2 * distance_ij * (distancek_1 + distancek)- 3 * Math.pow(distancek - distancek_1, 2))* (distancek + distancek_1 - Math.pow(distance_ij, 2));
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			double I2 = (temp_A1 + temp_A2 - temp_A3
					- (Math.pow(e1, 2) / volume + Math.pow(e2, 2) / volume - 2 * e1 * e2 / volume) * (nb_a * (nb_a - 1)))
					* nb_b / volume;

			double I1 = (e2 * sum_h_a - e1 * sum_h_a_bis - nb_a * Math.pow(e2 - e1, 2) / volume) * nb_b / volume;

			result[k - 1] = Math.pow(volume / (nb_b * nb_a), 2) * (I1 + I2);
		}
 		return result;

	}
public static Point3D getIntensityCenter(ROI3D roi, Sequence seq) throws InterruptedException {
		
		IcyBufferedImage image = seq.getImage(0, 0, 0);
		double x = 0, y = 0,z=0;

		final BooleanMask3D mask = roi.getBooleanMask3D(0,(int)(roi.getPosition5D().getT()) , -1, true);		
		final int h = mask.bounds.sizeX;
		final int w = mask.bounds.sizeY;
		final int d = mask.bounds.sizeZ;

		
		int off = 0;
		double total_weight = 0;
		final Point5D pos = roi.getPosition5D();
		for (int k = 0; k < d; k++) {
			for (int j = 0; j < h; j++) {		
				for (int i = 0; i < w; i++) {							
				if (mask.contains((int) (pos.getX() + i), (int) (pos.getY() + j), (int) (pos.getZ() + k))) {
					double weight = image.getData((int) (pos.getX() + i), (int) (pos.getY() + j), (int) (pos.getZ() + k));
					total_weight += weight;
					x += i * weight;
					y += j * weight;
					z += k * weight;
				}
			}
		}
		}
		if (total_weight==0.0)
		{return new Point3D.Double(pos.getX(), pos.getY(),pos.getZ());}
		else{
		return new Point3D.Double(pos.getX() + (x / total_weight), pos.getY() + (y / total_weight),pos.getZ() + (z / total_weight));
		}
		}
		
/*public static Point3D getMassCenter(ROI3D roi) {
		
		double x=0,y=0,z=0;
		final BooleanMask3D mask = roi.getBooleanMask3D(-1, 0, 0, true);
		
		final int w = mask.bounds.sizeX;
		final int h = mask.bounds.sizeY;
		final int d = mask.bounds.sizeZ;
		
		double total_weight = 0;
		final Point3D pos3d = roi.getPosition3D();
		for (int k = 0; k < d; k++) {
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				int xi=(int) (pos3d.getX() + i); int yi=(int) (pos3d.getY() + j); int zi=(int) (pos3d.getZ() + k); 
				if (mask.contains(xi, yi, zi)){
					total_weight+=1;
					x += i;
					y += j;
					z += k;
				}
			}
		}
		}
		if (total_weight==0.0)
		{return roi.getPosition3D();}
		else{
		return new Point3D.Double(pos3d.getX() + (x / total_weight), pos3d.getY() + (y / total_weight),pos3d.getZ() + (z / total_weight));}
		}
*/
	public static double distance2Polygon3D(double x_point, double y_point, double z_point, Point3D[] polyg,
			double ratio_zx) {
		double dist = Integer.MAX_VALUE;// roiPolyg.getPerimeter();

		int nbc = polyg.length;

		for (int i = 0; i < nbc; i++) {
			Point3D pt1 = polyg[i];
			double disttmp = Math.sqrt(Math.pow(x_point - pt1.getX(), 2) + Math.pow(y_point - pt1.getY(), 2)
					+ Math.pow((z_point - pt1.getZ()) * ratio_zx, 2));
			dist = Math.min(dist, disttmp);
		}

		return (dist);

	}

	public static double distance2Polygon(double x_point, double y_point, double z_point,
			ArrayList<Point3D.Integer> polyg, double ratio_zx) {
		double dist = Integer.MAX_VALUE;// roiPolyg.getPerimeter();

		int nbc = polyg.size();

		for (int i = 0; i < nbc; i++) {
			Point3D pt1 = polyg.get(i);
			double disttmp = Math.sqrt(Math.pow(x_point - pt1.getX(), 2) + Math.pow(y_point - pt1.getY(), 2)
					+ Math.pow((z_point - pt1.getZ()) * ratio_zx, 2));
			dist = Math.min(dist, disttmp);
		}

		return (dist);

	}

}
