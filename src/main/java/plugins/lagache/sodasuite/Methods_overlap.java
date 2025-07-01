package plugins.lagache.sodasuite;

import java.awt.Point;
import java.util.ArrayList;
import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point5D;

public class Methods_overlap {
	public static double[] Overlap(Sequence seq1, ArrayList<ROI> roi1, ArrayList<ROI> roi2,ArrayList<ROI> roi_t,int nb_mc,int t,double T) throws InterruptedException {
		double[] p = new double[8];
		//p[0]&[1] = coeffs, p[2]&p[3] = mean (simulations), p[4]&[5] = pvalue&log pvalue (coeff1>simus),p[6]&[7] = pvalue&log pvalue (coeff1>simus)
		p[0]=0;p[1]=0;p[2]=0;p[3]=0;p[4]=0;p[5]=0;p[6]=0;p[7]=0;
		if (roi_t.size()==0){			
			return p;}
		//computation of the coefficients
		double[] coeffs = OverlapCoeff(seq1, roi1, roi2,roi_t,t,T);
		p[0] = coeffs[0];p[1] = coeffs[1];
		//computation of the p_value		
		double[] pvalue = new double[6]; 
		pvalue=pvalue_computation_overlap(coeffs[0],coeffs[1], seq1, roi1, roi2, roi_t, nb_mc,t,T);
		p[2]=pvalue[4];
		p[3]=pvalue[5];
		p[4]=pvalue[0];
		p[5]=pvalue[2];
		p[6]=pvalue[1];
		p[7]=pvalue[3];
		return p;		
	}
	public static double[] Manders(Sequence seq1, ArrayList<ROI> roi1, ArrayList<ROI> roi2,ArrayList<ROI> roi_t,int nb_mc,int t) throws InterruptedException {
		double[] p = new double[8];
		//p[0]&[1] = coeffs, p[2]&p[3] = mean (simulations), p[4]&[5] = pvalue&log pvalue (coeff1>simus),p[6]&[7] = pvalue&log pvalue (coeff2>simus)
		p[0]=0;p[1]=0;p[2]=0;p[3]=0;p[4]=0;p[5]=0;p[6]=0;p[7]=0;
		if (roi_t.size()==0){			
			return p;}
		//computation of the coefficients
		double[] coeffs = MandersCoeff(roi1,roi2,roi_t,t);
		p[0] = coeffs[0];p[1] = coeffs[1];
		//computation of the p_value		
		double[] pvalue = new double[6]; 
		pvalue=pvalue_computation(coeffs[0],coeffs[1], seq1, roi1, roi2, roi_t, nb_mc,t);
		p[2]=pvalue[4];
		p[3]=pvalue[5];
		p[4]=pvalue[0];
		p[5]=pvalue[2];
		p[6]=pvalue[1];
		p[7]=pvalue[3];
		return p;		
	}
	public static double[] MandersCoeff(ArrayList<ROI> roi1, ArrayList<ROI> roi2,ArrayList<ROI> roi_t,int t) throws InterruptedException {
		
		double[] p=new double[2];
		p[0]=0;p[1]=0;
		if (roi_t.size()==0){			
			return p;}
		ROI union_t = ROIUtil.getUnion(roi_t);
		ArrayList<ROI> liste1 = new ArrayList<ROI>();
		ArrayList<ROI> liste2 = new ArrayList<ROI>();
		for (ROI r1:roi1){
			if (union_t.contains(r1)){
			liste1.add(r1);}}
		for (ROI r2:roi2){
			if (union_t.contains(r2)){
			liste2.add(r2);}}
		
		ROI union1 = ROIUtil.getUnion(liste1);
		ROI union2 = ROIUtil.getUnion(liste2);
						
		if ((union1==null)||(union2==null))
		{return p;}
		else
		{			
			ROI inter = ROIUtil.getIntersection(union1,union2);
		if (union1.getNumberOfPoints()>0)
			p[0] = inter.getNumberOfPoints()/union1.getNumberOfPoints();
		if (union2.getNumberOfPoints()>0)
			p[1] = inter.getNumberOfPoints()/union2.getNumberOfPoints();		
		return p;}		
	}
public static double[] OverlapCoeff(Sequence seq1, ArrayList<ROI> roi1, ArrayList<ROI> roi2,ArrayList<ROI> roi_t,int t,double T) throws InterruptedException {
	double[] p=new double[2];
	p[0]=0;p[1]=0;
	if (roi_t.size()==0){			
		return p;}
	ROI union_t = ROIUtil.getUnion(roi_t);
	ArrayList<ROI> liste1 = new ArrayList<ROI>();
	ArrayList<ROI> liste2 = new ArrayList<ROI>();
	for (ROI r1:roi1){
		if (union_t.contains(r1)){
		liste1.add(r1);}}
	for (ROI r2:roi2){
		if (union_t.contains(r2)){
		liste2.add(r2);}}
	
	ROI union1 = ROIUtil.getUnion(liste1);
	int compteur_roi = 0;
		for (ROI r2:liste2){
				double volume2 = r2.getNumberOfPoints();
				//calcul intersection avec union 1 bis 
				if (volume2>0){
				double alpha = ROIUtil.getIntersection(r2, union1).getNumberOfPoints()/(volume2);
				if (alpha>T)
					compteur_roi+=1;}
			}		
		p[0] = (double)(compteur_roi)/((double)liste1.size());
		p[1] = (double)(compteur_roi)/((double)liste2.size());
		return p;		
	}
	public static double[] pvalue_computation (double m1,double m2,Sequence seq1,ArrayList<ROI> liste1,ArrayList<ROI> liste2, ArrayList<ROI> cell, int nb_mc,int t) throws InterruptedException {double[] pvalue=new double[6];pvalue[0]=0;pvalue[1]=0;pvalue[2]=0;pvalue[3]=0;pvalue[4]=0;pvalue[5]=0;
	for (int i=1;i<nb_mc;i++){
		ArrayList<ROI> liste_new = randomization(liste2,cell,t);
		double[] coeffs = MandersCoeff(liste1, liste_new,cell,t);
		pvalue[4]+=coeffs[0];
		pvalue[5]+=coeffs[1];
		if (coeffs[0]>m1)
			pvalue[0]+=1;
		if (coeffs[1]>m2)
			pvalue[1]+=1;
	}	
		pvalue[0]=pvalue[0]/nb_mc;
		pvalue[1]=pvalue[1]/nb_mc;
		pvalue[4]=pvalue[4]/nb_mc;
		pvalue[5]=pvalue[5]/nb_mc;
		if (pvalue[0]>0){
		pvalue[2]= Math.log10(pvalue[0]);}
		if (pvalue[1]>0){
			pvalue[3]= Math.log10(pvalue[1]);}		
	return pvalue;
	}
	
	public static double[] pvalue_computation_overlap (double m1,double m2,Sequence seq1,ArrayList<ROI> liste1,ArrayList<ROI> liste2, ArrayList<ROI> cell, int nb_mc,int t,double T) throws InterruptedException {double[] pvalue=new double[6];
	pvalue[0]=0;pvalue[1]=0;pvalue[2]=0;pvalue[3]=0;pvalue[4]=0;pvalue[5]=0;
	for (int i=1;i<nb_mc;i++){
		ArrayList<ROI> liste_new = randomization(liste2,cell,t);
		double[] coeffs = OverlapCoeff(seq1, liste1, liste_new,cell,t,T);
		pvalue[4]+=coeffs[0];
		pvalue[5]+=coeffs[1];
		if (coeffs[0]>m1)
			pvalue[0]+=1;
		if (coeffs[1]>m2)
			pvalue[1]+=1;
	}
	if (pvalue[0]>0){
		pvalue[0] = pvalue[0]/nb_mc;
		pvalue[2] = Math.log10(pvalue[0]);
		}
	if (pvalue[1]>0){
		pvalue[1] = pvalue[1]/nb_mc;
		pvalue[3] = Math.log10(pvalue[1]);
		}
	pvalue[4]=pvalue[4]/nb_mc;
	pvalue[5]=pvalue[5]/nb_mc;
	return pvalue;	
	}
			
	public static ArrayList<ROI> randomization (ArrayList<ROI> liste, ArrayList<ROI> cell,int t) throws InterruptedException {
		ArrayList<ROI> liste_retour = new ArrayList<ROI>();		
		BooleanMask2D mask = ROIUtil.getUnion(cell).getBooleanMask2D(0,t,0,true);
		Point[] mask_pts = mask.getPoints();
		int ind=0;
		while (ind<liste.size()){
			ROI r=liste.get(ind).getCopy();
			int indice = (int) (mask_pts.length*Math.random());
			Point5D newPosition = new Point5D.Double(mask_pts[indice].x, mask_pts[indice].y, 0, t, 0);
			r.setPosition5D(newPosition);
			if (ROIUtil.getUnion(cell).contains(r)){
				liste_retour.add(r);ind++;}			
		}			
		return liste_retour;			
	}
		
}
