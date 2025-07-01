package plugins.lagache.sodasuite;

import icy.type.point.Point3D;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;



public class apparatedLocalizations {
		
	Point3D p1;
	Point3D p2;
	double distance;
	double proba;

	public apparatedLocalizations( Point3D p1, Point3D p2, double distance,double proba ) {
		
		this.p1 = p1;
		this.p2=p2;
		this.distance=distance;
		this.proba=proba;
	}

	public static ArrayList<Point3D> coloc1 = new ArrayList<Point3D>();
	public static ArrayList<Point3D> coloc2 = new ArrayList<Point3D>();
	
	public static ArrayList<Point3D> single1 = new ArrayList<Point3D>();
	public static ArrayList<Point3D> single2 = new ArrayList<Point3D>();
	
public static ArrayList<apparatedLocalizations> appDetectConstruction(ArrayList<Point3D> points,ArrayList<Point3D> points2,double[]probas,ArrayList<Double>distances,int dim)
{
	ArrayList<apparatedLocalizations> liste_retour = new ArrayList<apparatedLocalizations>();
	int N = distances.size();
	
	if (points.isEmpty()){}else{
	for (Point3D p:points) {		
		for (Point3D p2:points2) {
		double temp=0;	
		if (dim==2){	
		temp=Math.sqrt(Math.pow(p.getX()-p2.getX(), 2)+Math.pow(p.getY()-p2.getY(), 2));}
		else
		{temp=Math.sqrt(Math.pow(p.getX()-p2.getX(), 2)+Math.pow(p.getY()-p2.getY(), 2)+Math.pow(p.getZ()-p2.getZ(), 2));}
		//calcul de la proba associée
		double proba_temp=0.;
		if (temp<distances.get(N-1))
		{
			int j=N-1;
			while ((temp<distances.get(j))&&(j>0)){j=j-1;}
			proba_temp = probas[j];
		}
		if (proba_temp>0){
		apparatedLocalizations aS = new apparatedLocalizations(p, p2, temp,proba_temp);
		liste_retour.add(aS);}
	}				
	}}
	return liste_retour;
}

//cr�e la sous liste des spots "colocalis�s" par rapport au pourcentage calcul� statistiquement	
public static ArrayList<apparatedLocalizations> appDetectSelect(ArrayList<apparatedLocalizations> liste_app_detect,int ind_max,double min_rad)
{ArrayList<apparatedLocalizations> liste_retour = new ArrayList<apparatedLocalizations>();
//d�termination de la distance "max" accept�e
//cr�ation d'une liste de toutes les distance
ArrayList<Double> distances = new ArrayList<Double>();
for (apparatedLocalizations aS:liste_app_detect)
{
	distances.add(aS.distance);
}
Collections.sort(distances);
//détermination de l'indice min
int indice_min=0;
if (distances.isEmpty()||ind_max==0){}else{
while (distances.get(indice_min)<min_rad)
{indice_min++;}}
//d�termination de l'indice max qui donne la distance max
//int ind_max = (int)(percentage*distances.size());
if (distances.isEmpty()||ind_max==0){} else{
double distance_max=distances.get(ind_max-1+indice_min);
//cr�ation de la sous liste des apparayedSpots dont la distance est < distance_max
for (apparatedLocalizations aS:liste_app_detect)
{
	if ((aS.distance<=distance_max)&(aS.distance>=min_rad))
	{
		liste_retour.add(aS);
	}
}	}
return liste_retour;
}

//remplit un tableau de ROIs correspondant aux spots  qui colocalisent � partir de la s�lection des spots colocalis�s
public static void locColoc( int t, ArrayList<apparatedLocalizations> spotsColoc)
{			
	coloc1.clear();coloc2.clear();	
	for (apparatedLocalizations aS:spotsColoc)
	{			
		coloc1.add(aS.p1);
		coloc2.add(aS.p2);
				
	}
								
}
	//remplit un tableau avec les ROIs single
		public static void locSingle(ArrayList<apparatedLocalizations> spotsColoc,ArrayList<Point3D> detection1,ArrayList<Point3D> detection2)
		{			
			single1.clear();single2.clear();
			int ind1=0;
			while (ind1<detection1.size()){			
			boolean isSingle = true;
			Point3D s1 = detection1.get(ind1);
			for (apparatedLocalizations aS: spotsColoc){
				if (s1.equals(aS.p1)) isSingle=false;}
			if (isSingle)
			{single1.add(s1);}//if (proba1.length==detection1.size()){probaS1.add(proba1[ind1]);}}
			ind1++;
			}				
			
			int ind2 = 0;
			while (ind2<detection2.size()){	
				boolean isSingle = true;
				Point3D s2 = detection2.get(ind2);
				for (apparatedLocalizations aS: spotsColoc){
					if (s2.equals(aS.p2))isSingle=false;}
				if (isSingle)
				{single2.add(s2);}//if (proba2.length==detection2.size()){probaS2.add(proba2[ind2]);}}
			ind2++;	
			}
		}
		
		public static int[] nbSingle(ArrayList<apparatedLocalizations> spotsColoc,ArrayList<Point3D> detection1,ArrayList<Point3D> detection2, double dist_max)
		{
			int[] nb_single = new int[2];
			int ind1=0;
			while (ind1<detection1.size()){			
			boolean isSingle = true;
			Point3D s1 = detection1.get(ind1);
			for (apparatedLocalizations aS: spotsColoc){
				if (s1.equals(aS.p1)&&aS.distance<dist_max) isSingle=false;}
			if (isSingle)
			{nb_single[0]+=1;}//if (proba1.length==detection1.size()){probaS1.add(proba1[ind1]);}}
			ind1++;
			}				
			
			int ind2 = 0;
			while (ind2<detection2.size()){	
				boolean isSingle = true;
				Point3D s2 = detection2.get(ind2);
				for (apparatedLocalizations aS: spotsColoc){
					if (s2.equals(aS.p2)&&aS.distance<dist_max)isSingle=false;}
				if (isSingle)
				{nb_single[1]+=1;}
			ind2++;	
			}
			return nb_single;
		}

//calcul de la distance moyenne
		public static double distance_moyenne(ArrayList<apparatedLocalizations> spotsColoc, double dist_max)
		{
			double distance_moy=0;
			int compteur=0;
			for (apparatedLocalizations aS: spotsColoc){
				if (aS.distance<dist_max) 
					{distance_moy+=aS.distance;compteur+=1;}
		}
			return distance_moy/compteur;
		}


}
		