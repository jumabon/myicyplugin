package plugins.lagache.sodasuite;


import java.util.ArrayList;
import Jama.*;

import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.type.point.Point3D;

public class non_parametric_object {
	
public static void main3D_corr(int N_fit,Sequence sequence, ArrayList<ROI> roi_t,ArrayList<Point3D> positions1 ,ArrayList<Point3D> positions2, Window3D[][][] imageWindows_1,Window3D[][][] imageWindows_2,ArrayList<Double> distance_fit,int N_h,double ratio_zx,double volume,double[] results, double[] delta_K, double[][] Ktemp,double[] proba_dist,double[] maximum,double[] pvalue) throws InterruptedException {
	int nbdeta = positions1.size();int nbdetb = positions2.size();

	double[] var_delta_new=new double[N_fit-1];

var_delta_new = Ripley3D.variance_theo_delta_new(roi_t, imageWindows_1, imageWindows_2, volume, nbdeta, nbdetb, N_fit, distance_fit, sequence,ratio_zx);	
//calcul du nb de coloc corrigé
//calcul de la matrice des intersection et renormalisation du nombre de couple par anneau
double[][] A=intersection3D_new(imageWindows_1,distance_fit,nbdeta);

    
int N=distance_fit.size();
double T_p=Math.sqrt(2);


//définition des matrices
double[][] Num  = new double[N-1][1];
double [][] mu_tab = new double[N-1][1];
double [][] sigma_T = new double[N-1][1];
double max_delta=-100;
for (int p=0;p<N-1;p++)
{
//nb total de positions dans la couronne [p;p+1]	
Num[p][0]=delta_K[p]*(nbdeta*nbdetb/volume);
//espérance sous l'hypothese nulle de distribution aléatoire
mu_tab[p][0] = 4/3*Math.PI*(Math.pow(distance_fit.get(p+1),3)-Math.pow(distance_fit.get(p),3))*(nbdeta*nbdetb/volume);
//seuil adaptatif: s'actualise avec p
if (p>0){T_p = Math.sqrt(2*Math.log(p+1));}
//sigma_T[p][0]=(nbdeta*nbdetb*T/volume)*Math.sqrt(var_delta_new[p]);
sigma_T[p][0]=(nbdeta*nbdetb*T_p/volume)*Math.sqrt(var_delta_new[p]);
}
Matrix Num_mat=new Matrix(Num);
Matrix mu_tab_mat = new Matrix(mu_tab);
Matrix A_mat=new Matrix(A);
Matrix A_matb = A_mat.times((double)1/nbdeta);
Matrix A_matb_inverse = A_matb.inverse();
Matrix C = A_matb_inverse.times(Num_mat.minus(mu_tab_mat));

//calcul itératif des pvalues
for (int p=0;p<N-1;p++)
{
	if (p>0){T_p = Math.sqrt(2*Math.log(p+1));}
	else{T_p=Math.sqrt(2);}
	
	max_delta=Math.max(max_delta,C.get(p, 0)*T_p/sigma_T[p][0]);
	maximum[p]=max_delta;
	double normcdf=0.5*(1+ErrorFunction.erf(max_delta/(Math.sqrt(2))));
	pvalue[p]=(1-Math.pow(normcdf, p+1));	
}
for (int p=0;p<N-1;p++)
{//thresholding pour ne garder que les colocs C significatives
	if (C.get(p, 0)<sigma_T[p][0])
	C.set(p, 0, 0.);}


	
	//double normcdf=0.5*(1+ErrorFunction.erf(max_delta/Math.sqrt(2)));
	//retour[0]=(1-Math.pow(normcdf, N_fit-1));	
	//retour[4]=max_delta;
	//hard thresholding

	
	double[] total_events=new double[N_fit-1];
	
	for (int i=0;i<N_fit-1;i++)
	{
	    total_events[i]=(nbdetb*nbdeta/volume)*delta_K[i];
	    if (total_events[i]==0){proba_dist[i]=0;}else{
	    proba_dist[i]=C.get(i, 0)/total_events[i];}
	    }	
}
public static void main2D_corr(int N_fit,Sequence sequence, ArrayList<ROI> roi_t,ArrayList<Point3D> positions1 ,ArrayList<Point3D> positions2, Window2D[][] imageWindows_1,Window2D[][] imageWindows_2,ArrayList<Double> distance_fit,int N_h,double volume,double[] results, double[] delta_K, double[][] Ktemp,double[] proba_dist,double[] maximum,double[] pvalue) throws InterruptedException {
	int nbdeta = positions1.size();int nbdetb = positions2.size();
	double[] var_delta_new=new double[N_fit-1];
	var_delta_new = Ripley2D.variance_theo_delta_new(roi_t, imageWindows_1, imageWindows_2, volume, nbdeta, nbdetb, N_fit, distance_fit, N_h, results);	
//calcul du nb de coloc corrigé
//calcul de la matrice des intersection et renormalisation du nombre de couple par anneau
double[][] A=intersection2D_new(imageWindows_1,distance_fit,nbdeta);

    
int N=distance_fit.size();
double T_p=Math.sqrt(2);


//définition des matrices
double[][] Num  = new double[N-1][1];
double [][] mu_tab = new double[N-1][1];
double [][] sigma_T = new double[N-1][1];
double max_delta=-100;
for (int p=0;p<N-1;p++)
{
//nb total de positions dans la couronne [p;p+1]	
Num[p][0]=delta_K[p]*(nbdeta*nbdetb/volume);
//espérance sous l'hypothese nulle de distribution aléatoire
mu_tab[p][0] = Math.PI*(Math.pow(distance_fit.get(p+1),2)-Math.pow(distance_fit.get(p),2))*(nbdeta*nbdetb/volume);
//seuil adaptatif: s'actualise avec p
if (p>0){T_p = Math.sqrt(2*Math.log(p+1));}
//sigma_T[p][0]=(nbdeta*nbdetb*T/volume)*Math.sqrt(var_delta_new[p]);
sigma_T[p][0]=(nbdeta*nbdetb*T_p/volume)*Math.sqrt(var_delta_new[p]);
}
Matrix Num_mat=new Matrix(Num);
Matrix mu_tab_mat = new Matrix(mu_tab);
Matrix A_mat=new Matrix(A);
Matrix A_matb = A_mat.times((double)1/nbdeta);
Matrix A_matb_inverse = A_matb.inverse();
Matrix C = A_matb_inverse.times(Num_mat.minus(mu_tab_mat));

//calcul itératif des pvalues
for (int p=0;p<N-1;p++)
{
	if (p>0){T_p = Math.sqrt(2*Math.log(p+1));}
	else{T_p=Math.sqrt(2);}
	
	max_delta=Math.max(max_delta,C.get(p, 0)*T_p/sigma_T[p][0]);
	maximum[p]=max_delta;
	double normcdf=0.5*(1+ErrorFunction.erf(max_delta/(Math.sqrt(2))));
	pvalue[p]=(1-Math.pow(normcdf, p+1));	
}
for (int p=0;p<N-1;p++)
{//thresholding pour ne garder que les colocs C significatives
	if (C.get(p, 0)<sigma_T[p][0])
	C.set(p, 0, 0.);}


	
	//double normcdf=0.5*(1+ErrorFunction.erf(max_delta/Math.sqrt(2)));
	//retour[0]=(1-Math.pow(normcdf, N_fit-1));	
	//retour[4]=max_delta;
	//hard thresholding

	
	double[] total_events=new double[N_fit-1];
	
	for (int i=0;i<N_fit-1;i++)
	{
	    total_events[i]=(nbdetb*nbdeta/volume)*delta_K[i];
	    if (total_events[i]>0){
	    proba_dist[i]=C.get(i, 0)/total_events[i];}
	    else
	    	proba_dist[i]=0.;
	    }	
}
public static double[][] intersection3D(ArrayList<Point3D> positions1,ArrayList<Double> distance_fit){
int N=distance_fit.size();
int n=positions1.size();
double[][] A = new double[N-1][N-1];
for (int i=0;i<N-1;i++){
	A[i][i]=n;
}
double[][] S = new double[N][N];

for (int i=0;i<N;i++){
	for (int j=0;j<N;j++){
		double m=Math.min(distance_fit.get(i),distance_fit.get(j));
		double M=Math.max(distance_fit.get(i),distance_fit.get(j));
        //on fait ensuite la boucle sur les points
        for (Point3D pta:positions1){
            for (Point3D ptb:positions1){            	            
                double d=Math.sqrt(Math.pow(pta.getX()-ptb.getX(),2)+Math.pow(pta.getY()-ptb.getY(),2)+Math.pow(pta.getZ()-ptb.getZ(),2));
                if (d>0){                                    
                    if (m+d<M)
                        S[i][j]=S[i][j]+4/3*Math.PI*Math.pow(m,3);
                    else
                        if (d<(m+M)){
                        	S[i][j]=S[i][j]+(Math.PI/(12*d))*Math.pow(M+m-d,2)*(Math.pow(d,2)+2*d*m-3*Math.pow(m,2)+2*d*M+6*M*m-3*Math.pow(M,2));}
                }}}}}
                        
for (int i=0;i<N-1;i++){
	for (int j=0;j<N-1;j++){
        double vol=4/3*Math.PI*(Math.pow(distance_fit.get(j+1),3)-Math.pow(distance_fit.get(j),3));
        A[i][j]=A[i][j]+(S[i+1][j+1]+S[i][j]-S[i][j+1]-S[i+1][j])/vol;        
	}}
return A;
}

public static double[][] intersection3D_new(Window3D[][][] imageWindows_1,ArrayList<Double> distance_fit,int nba){
int N=distance_fit.size();
double[][] A = new double[N-1][N-1];
for (int i=0;i<N-1;i++){
	A[i][i]=nba;
}
double[][] S = new double[N][N];
for (int k1=0;k1<N;k1++){
	for (int k2=0;k2<N;k2++){
		double m=Math.min(distance_fit.get(k1),distance_fit.get(k2));
		double M=Math.max(distance_fit.get(k1),distance_fit.get(k2));
		
for (int i = 0; i < imageWindows_1.length; i++) {
	for (int j = 0; j < imageWindows_1[i].length; j++) {
		for (int z = 0; z < imageWindows_1[i][j].length; z++) {
			for (Point3D pta : imageWindows_1[i][j][z].detectionlist) {
				for (int o = Math.max(i - 1, 0); o <= Math.min(i + 1, imageWindows_1.length - 1); o++) {
					for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1,imageWindows_1[i].length - 1); l++) {
						for (int n = Math.max(z - 1, 0); n <= Math.min(z + 1,imageWindows_1[i][j].length - 1); n++) {
							for (Point3D ptb : imageWindows_1[o][l][n].detectionlist) {
								double d=Math.sqrt(Math.pow(pta.getX()-ptb.getX(),2)+Math.pow(pta.getY()-ptb.getY(),2)+Math.pow(pta.getZ()-ptb.getZ(),2));							
                if (d>0){                                    
                    if (m+d<M)
                        S[k1][k2]=S[k1][k2]+4/3*Math.PI*Math.pow(m,3);
                    else
                        if (d<(m+M)){
                        	S[k1][k2]=S[k1][k2]+(Math.PI/(12*d))*Math.pow(M+m-d,2)*(Math.pow(d,2)+2*d*m-3*Math.pow(m,2)+2*d*M+6*M*m-3*Math.pow(M,2));}
                }}}}}}}}}}}
                        
for (int i=0;i<N-1;i++){
	for (int j=0;j<N-1;j++){
        double vol=4/3*Math.PI*(Math.pow(distance_fit.get(j+1),3)-Math.pow(distance_fit.get(j),3));
        A[i][j]=A[i][j]+(S[i+1][j+1]+S[i][j]-S[i][j+1]-S[i+1][j])/vol;        
	}}
return A;}								

public static double[][] intersection2D_new(Window2D[][] imageWindows_1,ArrayList<Double> distance_fit,int nba){
int N=distance_fit.size();
double[][] A = new double[N-1][N-1];
for (int i=0;i<N-1;i++){
	A[i][i]=nba;
}
double[][] S = new double[N][N];
for (int k1=0;k1<N;k1++){
	for (int k2=0;k2<N;k2++){
		double m=Math.min(distance_fit.get(k1),distance_fit.get(k2));
		double M=Math.max(distance_fit.get(k1),distance_fit.get(k2));
		
for (int i = 0; i < imageWindows_1.length; i++) {
	for (int j = 0; j < imageWindows_1[i].length; j++) {		
			for (Point3D pta : imageWindows_1[i][j].detectionlist) {
				for (int o = Math.max(i - 1, 0); o <= Math.min(i + 1, imageWindows_1.length - 1); o++) {
					for (int l = Math.max(j - 1, 0); l <= Math.min(j + 1,imageWindows_1[i].length - 1); l++) {						
							for (Point3D ptb : imageWindows_1[o][l].detectionlist) {
								double d=Math.sqrt(Math.pow(pta.getX()-ptb.getX(),2)+Math.pow(pta.getY()-ptb.getY(),2));							
                if (d>0){                                    
                    if (m+d<M)
                        S[k1][k2]=S[k1][k2]+Math.PI*Math.pow(m,2);
                    else
                        if (d<(m+M)){
                        	double temp1 = Math.pow(m, 2)*Math.acos((Math.pow(d, 2)+Math.pow(m, 2)-Math.pow(M, 2))/(2*d*m));
                        	double temp2 = Math.pow(M, 2)*Math.acos((Math.pow(d, 2)+Math.pow(M, 2)-Math.pow(m, 2))/(2*d*M));
                        	double temp3 = 0.5*Math.sqrt((-d+m+M)*(d+m-M)*(d-m+M)*(d+m+M));
                        	S[k1][k2]=S[k1][k2]+temp1+temp2-temp3;
                }}}}}}}}}}
                        
for (int i=0;i<N-1;i++){
	for (int j=0;j<N-1;j++){
        double vol=Math.PI*(Math.pow(distance_fit.get(j+1),2)-Math.pow(distance_fit.get(j),2));
        A[i][j]=A[i][j]+(S[i+1][j+1]+S[i][j]-S[i][j+1]-S[i+1][j])/vol;        
	}}
return A;
}								
}