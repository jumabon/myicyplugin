package plugins.lagache.sodasuite;

import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROI3D;
import icy.sequence.Sequence;
import icy.type.point.Point3D;

import java.util.ArrayList;
import java.util.Collections;


public class apparatedSpots {
	ROI s1;
	ROI s2;
	double distance;
	double proba;

	public apparatedSpots( ROI s1, ROI s2, double distance,double proba ) {
		this.s1 = s1;
		this.s2=s2;
		this.distance=distance;
		this.proba=proba;
	}

	public static ArrayList<ROI> coloc1 = new ArrayList<ROI>();
	public static ArrayList<ROI> coloc2 = new ArrayList<ROI>();
	
	public static ArrayList<ROI> single1 = new ArrayList<ROI>();
	public static ArrayList<ROI> single2 = new ArrayList<ROI>();
	public static ArrayList<Double> probaS1 = new ArrayList<Double>();
	public static ArrayList<Double> probaS2 = new ArrayList<Double>();

public static ArrayList<apparatedSpots> appDetectConstruction(ROI[] spots,ROI[] spots2,double[]probas,ArrayList<Double>distances,Sequence sequence1,Sequence sequence2) throws InterruptedException {
	ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();
	int N = distances.size();
	
	if (spots==null){}else{
	int nbdeta=spots.length;
	int nbdetb=spots2.length;
	double x_a,y_a,z_a,x_b,y_b,z_b;
	for (int p = 0; p < nbdeta; p++) {
		if (spots[p].getDimension()==2){		
		Point3D position_a= Ripley2D.getIntensityCenter((ROI2D)spots[p], sequence1);x_a=position_a.getX();y_a=position_a.getY();z_a=0;}
		else
		{Point3D position_a= Ripley3D.getIntensityCenter((ROI3D)spots[p], sequence1);x_a=position_a.getX();y_a=position_a.getY();z_a=position_a.getZ();}
		for (int p2 = 0; p2 < nbdetb; p2++) {
			if (spots2[p2].getDimension()==2){
			Point3D position_b= Ripley2D.getIntensityCenter((ROI2D)spots2[p2],sequence2);x_b=position_b.getX();y_b=position_b.getY();z_b=0;}
			else {Point3D position_b= Ripley3D.getIntensityCenter((ROI3D)spots2[p2],sequence2);x_b=position_b.getX();y_b=position_b.getY();z_b=position_b.getZ();}
		double temp=Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2)+Math.pow(z_a-z_b, 2));;
		//calcul de la proba associée
		double proba_temp=0.;
		if (temp<distances.get(N-1))
		{
			int j=N-1;
			while ((temp<distances.get(j))&&(j>0)){j=j-1;}
			proba_temp = probas[j];
		}
		if (proba_temp>0){
		apparatedSpots aS = new apparatedSpots(spots[p], spots2[p2], temp,proba_temp);
		liste_retour.add(aS);}
	}				
	}}
	return liste_retour;
}

//cr�e la sous liste des spots "colocalis�s" par rapport au pourcentage calcul� statistiquement	
public static ArrayList<apparatedSpots> appDetectSelect(ArrayList<apparatedSpots> liste_app_detect,int ind_max,double min_rad)
{ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();
//d�termination de la distance "max" accept�e
//cr�ation d'une liste de toutes les distance
ArrayList<Double> distances = new ArrayList<Double>();
for (apparatedSpots aS:liste_app_detect)
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
for (apparatedSpots aS:liste_app_detect)
{
	if ((aS.distance<=distance_max)&(aS.distance>=min_rad))
	{
		liste_retour.add(aS);
	}
}	}
return liste_retour;
}

//remplit un tableau de ROIs correspondant aux spots  qui colocalisent � partir de la s�lection des spots colocalis�s
public static void roiColoc( int t, ArrayList<apparatedSpots> spotsColoc)
{			
	coloc1.clear();coloc2.clear();	
	for (apparatedSpots aS:spotsColoc)
	{			
		coloc1.add(aS.s1);
		coloc2.add(aS.s2);
				
	}
								
}
	//remplit un tableau avec les ROIs single
		public static void roiSingle(int t,ArrayList<apparatedSpots> spotsColoc,ROI[] detection1,ROI[] detection2)
		{			
			single1.clear();single2.clear();
			int ind1=0;
			while (ind1<detection1.length){			
			boolean isSingle = true;
			ROI s1 = detection1[ind1];
			for (apparatedSpots aS: spotsColoc){
				if (s1.equals(aS.s1)) isSingle=false;}
			if (isSingle)
			{single1.add(s1);}//if (proba1.length==detection1.size()){probaS1.add(proba1[ind1]);}}
			ind1++;
			}				
			
			int ind2 = 0;
			while (ind2<detection2.length){	
				boolean isSingle = true;
				ROI s2 = detection2[ind2];
				for (apparatedSpots aS: spotsColoc){
					if (s2.equals(aS.s2))isSingle=false;}
				if (isSingle)
				{single2.add(s2);}//if (proba2.length==detection2.size()){probaS2.add(proba2[ind2]);}}
			ind2++;	
			}
		}
		
		public static int[] nbSingle(int t,ArrayList<apparatedSpots> spotsColoc,ROI[] detection1,ROI[] detection2, double dist_max)
		{
			int[] nb_single = new int[2];
			int ind1=0;
			while (ind1<detection1.length){			
			boolean isSingle = true;
			ROI s1 = detection1[ind1];
			for (apparatedSpots aS: spotsColoc){
				if (s1.equals(aS.s1)&&aS.distance<dist_max) isSingle=false;}
			if (isSingle)
			{nb_single[0]+=1;}//if (proba1.length==detection1.size()){probaS1.add(proba1[ind1]);}}
			ind1++;
			}				
			
			int ind2 = 0;
			while (ind2<detection2.length){	
				boolean isSingle = true;
				ROI s2 = detection2[ind2];
				for (apparatedSpots aS: spotsColoc){
					if (s2.equals(aS.s2)&&aS.distance<dist_max)isSingle=false;}
				if (isSingle)
				{nb_single[1]+=1;}
			ind2++;	
			}
			return nb_single;
		}

//calcul de la distance moyenne
		public static double distance_moyenne(ArrayList<apparatedSpots> spotsColoc, double dist_max)
		{
			double distance_moy=0;
			int compteur=0;
			for (apparatedSpots aS: spotsColoc){
				if (aS.distance<dist_max) 
					{distance_moy+=aS.distance;compteur+=1;}
		}
			return distance_moy/compteur;
		}


}
		