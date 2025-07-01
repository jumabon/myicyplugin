package plugins.lagache.sodasuite;

import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import icy.type.rectangle.Rectangle3D;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarSequence;

public class Methods_Correlation {
	
	public static double[] pearson_TCL(Sequence seq1,
			Sequence seq2, int t, ArrayList<ROI> roi_t_liste) throws InterruptedException {
		
		if (roi_t_liste.isEmpty()){
			double[] p=new double[3];p[0]=0;p[1]=0;p[2]=0;
			return p;}
		
		int Z = seq1.getSizeZ();
		
		double[] Pearson = new double[3];
		double mean1 = 0;
		double mean2 = 0;
		double var1 = 0;
		double var2 = 0;
		double m1=0;		
		double sample = 0d;
		for (int z=0;z<Z;z++){
		IcyBufferedImage img1 = seq1.getImage(t, z);		
		IcyBufferedImage img2 = seq2.getImage(t, z);
		ROI roi_t = ROIUtil.getUnion(roi_t_liste);
		Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
			double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(0),
					img1.isSignedDataType());
			double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(0),
					img2.isSignedDataType());

			int sizeX = img1.getSizeX();
			
			int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
			int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
			int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
			int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
			

			for (int x = minX; x <= maxX; x++) {
				for (int y = minY; y <= maxY; y++) {
					if (roi_t.contains(x, y,z,t,-1)) {
						int off = (y * sizeX) + x;

						m1 += tab1[off] * tab2[off];						
						
						mean1 += tab1[off];
						mean2 += tab2[off];
						
						var1 += Math.pow(tab1[off], 2);
						var2 += Math.pow(tab2[off], 2);
						
						sample++;
					}
				}
			}}

			if (sample>0){
			mean1 /= sample;
			mean2 /= sample;
			m1/=sample;
			var1 /= sample;
			var2 /= sample;
			
			//calcul des moments du coeff de Pearson sous l'hypothese nulle (mean=0)
			double sigma1=Math.sqrt(var1-Math.pow(mean1,2));double sigma2=Math.sqrt(var2-Math.pow(mean2,2));			
			double m2_R=1/(sample);
			
			//calcul du coefficient de Pearson			
			Pearson[0] = (m1 - mean1*mean2)/(sigma1*sigma2);
			//calcul du coeff de Pearson centr� r�duit et de son skew
			double R_tilde=Pearson[0]/Math.sqrt(m2_R);
			double pvalue=0.5*(1-ErrorFunction.erf(R_tilde/Math.sqrt(2)));
			if (R_tilde<4.0){
				Pearson[1] = Math.log10(pvalue);
			}
			else
			{double x = R_tilde/Math.sqrt(2);
			Pearson[1]=-Math.pow(x, 2)/Math.log(10)-Math.log10(2*Math.sqrt(Math.PI)*x);}
			Pearson[2]=pvalue;
			}
			
		else {Pearson[0]=0;Pearson[1]=0;Pearson[2]=0;}		

		return Pearson;
	}
	public static double[] pearson_TCL_withC(Sequence seq1,
			Sequence seq2, int c1,int c2,int t, ArrayList<ROI> roi_t_liste) throws InterruptedException {
		
		if (roi_t_liste.isEmpty()){
			double[] p=new double[3];p[0]=0;p[1]=0;p[2]=0;
			return p;}
		
		int Z = seq1.getSizeZ();
		
		double[] Pearson = new double[3];
		double mean1 = 0;
		double mean2 = 0;
		double var1 = 0;
		double var2 = 0;
		double m1=0;		
		double sample = 0d;
		for (int z=0;z<Z;z++){
		IcyBufferedImage img1 = seq1.getImage(t,z,c1);		
		IcyBufferedImage img2 = seq2.getImage(t, z,c2);
		ROI roi_t = ROIUtil.getUnion(roi_t_liste);
		Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
			double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(0),
					img1.isSignedDataType());
			double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(0),
					img2.isSignedDataType());

			int sizeX = img1.getSizeX();
			
			int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
			int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
			int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
			int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
			

			for (int x = minX; x <= maxX; x++) {
				for (int y = minY; y <= maxY; y++) {
					if (roi_t.contains(x, y,z,t,-1)) {
						int off = (y * sizeX) + x;

						m1 += tab1[off] * tab2[off];						
						
						mean1 += tab1[off];
						mean2 += tab2[off];
						
						var1 += Math.pow(tab1[off], 2);
						var2 += Math.pow(tab2[off], 2);
						
						sample++;
					}
				}
			}}

			if (sample>0){
			mean1 /= sample;
			mean2 /= sample;
			m1/=sample;
			var1 /= sample;
			var2 /= sample;
			
			//calcul des moments du coeff de Pearson sous l'hypothese nulle (mean=0)
			double sigma1=Math.sqrt(var1-Math.pow(mean1,2));double sigma2=Math.sqrt(var2-Math.pow(mean2,2));			
			double m2_R=1/(sample);
			
			//calcul du coefficient de Pearson			
			Pearson[0] = (m1 - mean1*mean2)/(sigma1*sigma2);
			//calcul du coeff de Pearson centr� r�duit et de son skew
			double R_tilde=Pearson[0]/Math.sqrt(m2_R);
			double pvalue=0.5*(1-ErrorFunction.erf(R_tilde/Math.sqrt(2)));
			if (R_tilde<4.0){
				Pearson[1] = Math.log10(pvalue);
			}
			else
			{double x = R_tilde/Math.sqrt(2);
			Pearson[1]=-Math.pow(x, 2)/Math.log(10)-Math.log10(2*Math.sqrt(Math.PI)*x);}
			Pearson[2]=pvalue;
			}
		else {Pearson[0]=0;Pearson[1]=0;Pearson[2]=0;}		

		return Pearson;
	}

		public static double[] ICCS_compute(Sequence seq1,
				Sequence seq2,int c1,int c2, int t, ArrayList<ROI> roi_t_liste) throws InterruptedException {
			
			if (roi_t_liste.isEmpty()){
				double[] p=new double[2];p[0]=0;p[1]=1;
				return p;}
			
			int Z = seq1.getSizeZ();
			
			double[] ICCS = new double[2];
			double mean1 = 0;
			double mean2 = 0;
			double corr1 = 0,corr2=0,corr_inter=0;
			double sample = 0d;
			for (int z=0;z<Z;z++){
				IcyBufferedImage img1 = seq1.getImage(t, z,c1);		
				IcyBufferedImage img2 = seq2.getImage(t, z,c2);
				ROI roi_t = ROIUtil.getUnion(roi_t_liste);
				Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
					double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(0),
							img1.isSignedDataType());
					double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(0),
							img2.isSignedDataType());

				int sizeX = img1.getSizeX();
				int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
				int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
				int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
				int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
				
				for (int x = minX; x <= maxX; x++) {
					for (int y = minY; y <= maxY; y++) {
						if (roi_t.contains(x, y,z,t,-1)) {
							int off = (y * sizeX) + x;

							
							mean1 += tab1[off];
							mean2 += tab2[off];
							
							
							sample++;
						}
					}
				}}
			
				if (sample>0){
				mean1 /= sample;
				mean2 /= sample;}
				
				
				sample=0;
				//calcul des cross correlations
				for (int z=0;z<Z;z++){
					IcyBufferedImage img1 = seq1.getImage(t, z,c1);		
					IcyBufferedImage img2 = seq2.getImage(t, z,c2);
					ROI roi_t = ROIUtil.getUnion(roi_t_liste);
					Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
						double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(0),
								img1.isSignedDataType());
						double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(0),
								img2.isSignedDataType());

					int sizeX = img1.getSizeX();
					int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
					int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
					int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
					int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
					
					for (int x = minX; x <= maxX; x++) {
						for (int y = minY; y <= maxY; y++) {
							if (roi_t.contains(x, y,z,t,-1)) {
								int off = (y * sizeX) + x;

								
								corr1 += Math.pow(tab1[off]-mean1,2);
								corr2 += Math.pow(tab2[off]-mean2,2);
								corr_inter += (tab1[off]-mean1)*(tab2[off]-mean2);								
								
								sample++;
							}
						}
					}}

							
				if (sample>0){
				corr1/=sample;corr2/=sample;corr_inter/=sample;}			
				
				//calcul des fractions P1 et P2 de molecules qui interagissent
				if (corr1>0)
				ICCS[0]=corr_inter/corr1;
				
				if (corr2>0)
					ICCS[1]=corr_inter/corr2;				
											
			return ICCS;
		}


private double[] ICQ_compute (Sequence seq1, Sequence seq2, int t,int c1,int c2, ROI roi_t) {

double[] ICQ = new double[2];

if (roi_t==null){
	double[] p=new double[2];p[0]=0;p[1]=1;
	return p;}

int Z = seq1.getSizeZ();

double mean1 = 0;
double mean2 = 0;
double sample = 0d;
for (int z=0;z<Z;z++){
IcyBufferedImage img1 = seq1.getImage(t, z);
IcyBufferedImage img2 = seq2.getImage(t, z);
Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
	double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(c1),
			img1.isSignedDataType());
	double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(c2),
			img2.isSignedDataType());

	int sizeX = img1.getSizeX();
	
	int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
	int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
	int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
	int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
	

	for (int x = minX; x <= maxX; x++) {
		for (int y = minY; y <= maxY; y++) {
			if (roi_t.contains(x, y,z,t,-1)) {
				int off = (y * sizeX) + x;//calcul de mean 1 / mean 2
			mean1 += tab1[off];
			mean2 += tab2[off];
			sample++;
		}
	}
}
if (sample>0){
mean1 /= sample;
mean2 /= sample;}

double icq=0;
//calcul de ICQ
for (int x = minX; x <= maxX; x++) {
	for (int y = minY; y <= maxY; y++) {
		if (roi_t.contains(x, y,z,t,-1)) {
			int off = (y * sizeX) + x;//calcul de mean 1 / mean 2
			if ((tab1[off]-mean1)*(tab2[off]-mean2)>0)
			icq += 1;
		}
	}
}

if (sample>0){
icq/=sample;
icq-=0.5;
ICQ[0]=icq;

//test statistique (sign test) variance sample/4
double sig=0.5/Math.sqrt(sample);
double tem = ErrorFunction.erf(icq/(sig*Math.sqrt(2)));
ICQ[1]=0.5*(1+tem);}
else {ICQ[0]=0;ICQ[1]=0;}
}

return ICQ;
}
}


