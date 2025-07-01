package plugins.lagache.sodasuite;

import java.util.ArrayList;
import org.apache.poi.ss.util.WorkbookUtil;

import icy.file.FileUtil;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.ROI;
import icy.sequence.Sequence;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarDoubleArrayNative;
import plugins.adufour.vars.lang.VarROIArray;

// Colocalisation with Ripley function K
// Significant 

public class TRIOblock extends Plugin implements Block, PluginBundled {

	// EzVarSwimmingObject<DetectionResult> detections = new
	// EzVarSwimmingObject<DetectionResult>("Detections");
	VarROIArray detections1Total = new VarROIArray("List of detections 1 (ROIs)");
	VarROIArray detections2Total = new VarROIArray("List of detections 2 (ROIs)");
	VarROIArray detections3Total = new VarROIArray("List of detections 3 (ROIs)");

	VarROIArray detections1C2 = new VarROIArray("List of detections 1 Coupled with 2 (ROIs)");
	VarROIArray detections2C1 = new VarROIArray("List of detections 2 Coupled with 1  (ROIs)");

	// on fournit aussi la liste des probas et distances associés au couplage
	VarDoubleArrayNative probas12 = new VarDoubleArrayNative("Coupling probabilities with 2", null);
	VarDoubleArrayNative dist12 = new VarDoubleArrayNative("Coupling distances with 2", null);

	VarROIArray detections1C3 = new VarROIArray("List of detections 1 Coupled with 3 (ROIs)");
	VarROIArray detections3C1 = new VarROIArray("List of detections 3 Coupled with 1 (ROIs)");
	// on fournit aussi la liste des probas et distances associés au couplage
	VarDoubleArrayNative probas13 = new VarDoubleArrayNative("Coupling probabilities with 3", null);
	VarDoubleArrayNative dist13 = new VarDoubleArrayNative("Coupling distances with 3", null);

	VarROIArray detections2C3 = new VarROIArray("List of detections 2 Coupled with 3 (ROIs)");
	VarROIArray detections3C2 = new VarROIArray("List of detections 3 Coupled with 2 (ROIs)");
	// on fournit aussi la liste des probas et distances associés au couplage
	VarDoubleArrayNative probas23 = new VarDoubleArrayNative("Coupling probabilities with 3", null);
	VarDoubleArrayNative dist23 = new VarDoubleArrayNative("Coupling distances with 3", null);
	VarBoolean strict_trio = new VarBoolean("Keep strict trio only", false);

	// ROIs (spots) qui colocalisent en sortie
	VarROIArray detections1S = new VarROIArray("Single 1");
	VarROIArray detections2S = new VarROIArray("Single 2");
	VarROIArray detections3S = new VarROIArray("Single 3");

	VarROIArray detections1T = new VarROIArray("Trios: 1 with 2-3");
	VarROIArray detections2T = new VarROIArray("Trios: 2 with 1-3");
	VarROIArray detections3T = new VarROIArray("Trios: 3 with 1-2");
	VarDoubleArrayNative probasT12 = new VarDoubleArrayNative("Trios: proba. 1 with 2", null);
	VarDoubleArrayNative probasT13 = new VarDoubleArrayNative("Trios: proba. 1 with 3", null);
	VarDoubleArrayNative probasT23 = new VarDoubleArrayNative("Trios: proba. 2 with 3", null);
	VarDoubleArrayNative distT12 = new VarDoubleArrayNative("Trios: dist. 1 with 2", null);
	VarDoubleArrayNative distT13 = new VarDoubleArrayNative("Trios: dist. 1 with 3", null);
	VarDoubleArrayNative distT23 = new VarDoubleArrayNative("Trios: dist. 2 with 3", null);

	VarROIArray detections1C2_out = new VarROIArray("Couples: 1 with 2 (only)");
	VarROIArray detections2C1_out = new VarROIArray("Couples: 2 with 1 (only)");
	VarDoubleArrayNative probas1C2_out = new VarDoubleArrayNative("Couples: proba. 1 with 2 (only)", null);
	VarDoubleArrayNative dist1C2_out = new VarDoubleArrayNative("Couples: dist. 1 with 2 (only)", null);

	VarROIArray detections1C3_out = new VarROIArray("Couples: 1 with 3 (only)");
	VarROIArray detections3C1_out = new VarROIArray("Couples: 3 with 1 (only)");
	VarDoubleArrayNative probas1C3_out = new VarDoubleArrayNative("Couples: proba. 1 with 3 (only)", null);
	VarDoubleArrayNative dist1C3_out = new VarDoubleArrayNative("Couples: dist. 1 with 3 (only)", null);

	VarROIArray detections2C3_out = new VarROIArray("Couples: 2 with 3 (only)");
	VarROIArray detections3C2_out = new VarROIArray("Couples: 3 with 2 (only)");
	VarDoubleArrayNative probas2C3_out = new VarDoubleArrayNative("Couples: proba. 2 with 3 (only)", null);
	VarDoubleArrayNative dist2C3_out = new VarDoubleArrayNative("Couples: dist. 2 with 3 (only)", null);

	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add("Detections 1 (Total)", detections1Total);
		inputMap.add("Detections 2 (Total)", detections2Total);
		inputMap.add("Detections 3 (Total)", detections3Total);
		inputMap.add("Couples 1 with 2", detections1C2);
		inputMap.add("Couples 2 with 1", detections2C1);
		inputMap.add("Proba. 1 with 2", probas12);
		inputMap.add("Dist. 1 with 2", dist12);
		inputMap.add("Couples 1 with 3", detections1C3);
		inputMap.add("Couples 3 with 1", detections3C1);
		inputMap.add("Proba. 1 with 3", probas13);
		inputMap.add("Dist. 1 with 3", dist13);
		inputMap.add("Couples 2 with 3", detections2C3);
		inputMap.add("Couples 3 with 2", detections3C2);
		inputMap.add("Proba. 2 with 3", probas23);
		inputMap.add("Dist. 2 with 3", dist23);
		inputMap.add("Keep only strict trios",strict_trio);
	}

	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add("Single 1", detections1S);
		outputMap.add("Single 2", detections2S);
		outputMap.add("Single 3", detections3S);

		outputMap.add("Couples: 1 with 2", detections1C2_out);
		outputMap.add("Couples: 2 with 1", detections2C1_out);
		outputMap.add("Couples: Proba. 1 with 2", probas1C2_out);
		outputMap.add("Couples: Dist. 1 with 2", dist1C2_out);

		outputMap.add("Couples: 1 with 3", detections1C3_out);
		outputMap.add("Couples: 3 with 1", detections3C1_out);
		outputMap.add("Couples: Proba. 1 with 3", probas1C3_out);
		outputMap.add("Couples: Dist. 1 with 3", dist1C3_out);

		outputMap.add("Couples: 2 with 3", detections2C3_out);
		outputMap.add("Couples: 3 with 2", detections3C2_out);
		outputMap.add("Couples: Proba. 2 with 3", probas2C3_out);
		outputMap.add("Couples: Dist. 2 with 3", dist2C3_out);

		outputMap.add("Trios: 1 with 2-3", detections1T);
		outputMap.add("Trios: 2 with 1-3", detections2T);
		outputMap.add("Trios: 3 with 1-2", detections3T);
		outputMap.add("Trios: Proba. 1 with 2", probasT12);
		outputMap.add("Trios: Dist. 1 with 2", distT12);
		outputMap.add("Trios: Proba. 1 with 3", probasT13);
		outputMap.add("Trios: Dist. 1 with 3", distT13);
		outputMap.add("Trios: Proba. 2 with 3", probasT23);
		outputMap.add("Trios: Dist. 2 with 3", distT23);

	}

	// N: nb entre dmin et dmax avec pas donne
	@Override
	public void run() {
		// on cree des listes pour gerer facilement l'ajout/remove de ROI
		ArrayList<ROI> T1 = new ArrayList<ROI>();
		ArrayList<ROI> T2 = new ArrayList<ROI>();
		ArrayList<ROI> T3 = new ArrayList<ROI>();
		ArrayList<Double> probaT12 = new ArrayList<Double>();
		ArrayList<Double> disT12 = new ArrayList<Double>();
		ArrayList<Double> probaT13 = new ArrayList<Double>();
		ArrayList<Double> disT13 = new ArrayList<Double>();
		ArrayList<Double> probaT23 = new ArrayList<Double>();
		ArrayList<Double> disT23 = new ArrayList<Double>();

		ArrayList<ROI> C12 = new ArrayList<ROI>();
		ArrayList<ROI> C21 = new ArrayList<ROI>();
		ArrayList<Double> probaC12 = new ArrayList<Double>();
		ArrayList<Double> disC12 = new ArrayList<Double>();

		ArrayList<ROI> C13 = new ArrayList<ROI>();
		ArrayList<ROI> C31 = new ArrayList<ROI>();
		ArrayList<Double> probaC13 = new ArrayList<Double>();
		ArrayList<Double> disC13 = new ArrayList<Double>();

		ArrayList<ROI> C23 = new ArrayList<ROI>();
		ArrayList<ROI> C32 = new ArrayList<ROI>();
		ArrayList<Double> probaC23 = new ArrayList<Double>();
		ArrayList<Double> disC23 = new ArrayList<Double>();

		// détections des spots single
		ROI[] single1 = singlefinder(detections1Total.getValue(), detections1C2.getValue(), detections1C3.getValue());
		detections1S.setValue(single1);
		ROI[] single2 = singlefinder(detections2Total.getValue(), detections2C1.getValue(), detections2C3.getValue());
		detections2S.setValue(single2);
		ROI[] single3 = singlefinder(detections3Total.getValue(), detections3C1.getValue(), detections3C2.getValue());
		detections3S.setValue(single3);

		// détection des couples only
		couplefinder(C12, C21, probaC12, disC12, detections1C2.getValue(), detections1C3.getValue(),
				detections2C1.getValue(), detections2C3.getValue(), probas12.getValue(), dist12.getValue());
		detections1C2_out.setValue(ListToArray(C12));
		detections2C1_out.setValue(ListToArray(C21));
		probas1C2_out.setValue(ListToArrayD(probaC12));
		dist1C2_out.setValue(ListToArrayD(disC12));

		couplefinder(C13, C31, probaC13, disC13, detections1C3.getValue(), detections1C2.getValue(),
				detections3C1.getValue(), detections3C2.getValue(), probas13.getValue(), dist13.getValue());
		detections1C3_out.setValue(ListToArray(C13));
		detections3C1_out.setValue(ListToArray(C31));
		probas1C3_out.setValue(ListToArrayD(probaC13));
		dist1C3_out.setValue(ListToArrayD(disC13));

		couplefinder(C23, C32, probaC23, disC23, detections2C3.getValue(), detections2C1.getValue(),
				detections3C2.getValue(), detections3C1.getValue(), probas23.getValue(), dist23.getValue());
		detections2C3_out.setValue(ListToArray(C23));
		detections3C2_out.setValue(ListToArray(C32));
		probas2C3_out.setValue(ListToArrayD(probaC23));
		dist2C3_out.setValue(ListToArrayD(disC23));

		// détection des trios
		triofinder(T1, T2, T3, probaT12, probaT13, probaT23, disT12, disT13, disT23, detections1C2.getValue(),
				detections2C1.getValue(), detections1C3.getValue(), detections3C1.getValue(), detections2C3.getValue(),
				detections3C2.getValue(), probas12.getValue(), dist12.getValue(), probas13.getValue(),
				dist13.getValue(), probas23.getValue(), dist23.getValue(),strict_trio.getValue());
		detections1T.setValue(ListToArray(T1));
		detections2T.setValue(ListToArray(T2));
		detections3T.setValue(ListToArray(T3));
		probasT12.setValue(ListToArrayD(probaT12));
		probasT13.setValue(ListToArrayD(probaT13));
		probasT23.setValue(ListToArrayD(probaT23));
		distT12.setValue(ListToArrayD(disT12));
		distT13.setValue(ListToArrayD(disT13));
		distT23.setValue(ListToArrayD(disT23));

	}

	public ROI[] singlefinder(ROI[] detections1Tot, ROI[] detections1C2, ROI[] detections1C3) {
		ArrayList<ROI> retour = new ArrayList<ROI>();
		for (ROI d : detections1Tot) {
			// declaration des variables boolennes
			Boolean single = true;
			// on parcourt les couples avec 2
			for (ROI C2 : detections1C2) {
				if (d == C2) {
					single = false;
					break;
				}
			}
			// on parcourt les couples avec 3
			for (ROI C3 : detections1C3) {
				if (d == C3) {
					single = false;
					break;
				}
			}
			if (single == true)
				retour.add(d);
		}
		ROI[] ret = ListToArray(retour);
		return ret;

	}

	public void couplefinder(ArrayList<ROI> C1, ArrayList<ROI> C2, ArrayList<Double> proba, ArrayList<Double> dist,
			ROI[] C12, ROI[] C13, ROI[] C21, ROI[] C23, double[] p12, double[] d12) {
		// on parcourt ensuite les couples 1-2 et 1-3 pour les ranger..
		int ind12 = 0;

		for (ROI C1_2 : C12) {
			boolean couple_only = true;
			for (ROI C1_3 : C13) {
				if (C1_2 == C1_3) // la detection 1 est avec 2 et 3
				{
					couple_only = false;
					break;
				}
			}
			// il faut ensuite tester que le 2eme spot du couple n'est pas lié à
			// trois
			for (ROI C2_3 : C23) {
				if (C21[ind12] == C2_3) // la detection 1 est avec 2 et 3
				{
					couple_only = false;
					break;
				}
			}
			if (couple_only == true) {
				C1.add(C1_2);
				C2.add(C21[ind12]);
				proba.add(p12[ind12]);
				dist.add(d12[ind12]);
			}
			ind12++;
		}
	}

	public void triofinder(ArrayList<ROI> T1, ArrayList<ROI> T2, ArrayList<ROI> T3, ArrayList<Double> probaT12,
			ArrayList<Double> probaT13, ArrayList<Double> probaT23, ArrayList<Double> disT12, ArrayList<Double> disT13,
			ArrayList<Double> disT23, ROI[] detections1C2, ROI[] detections2C1, ROI[] detections1C3,
			ROI[] detections3C1, ROI[] detections2C3, ROI[] detections3C2, double[] probas12, double[] dist12,
			double[] probas13, double[] dist13, double[] probas23, double[] dist23,boolean strict_trio) {
		// pour les trios tester 1c2 et 1c3 si ok, tester le couple 2c3
		// correspondant et remplir la proba/dist, mettre la proba du lien à 0
		// pas de 2c3
		// il faut tester toutes les 3 combinaisons de 2 liens sur trois
		//si option strict_trio, il faut que les 3 probas soient >0
		if (strict_trio){
			int ind12 = 0;
			// 12 et 13
			for (ROI C1_2 : detections1C2) {
				int ind13 = 0;
				for (ROI C1_3 : detections1C3) {
					int ind23=0;										
					for (ROI C2_3 : detections2C3) {					
					//sur un trio strict, on vérifie que, de plus, la detection 2 est avec 3 (3 liens non nuls)
						if ((C1_2 == C1_3)&(detections2C1[ind12]==C2_3)&(detections3C1[ind13]==detections3C2[ind23])){					
						T1.add(C1_2);
						T2.add(detections2C1[ind12]);
						T3.add(detections3C1[ind13]);
						// on y ajoute les probas distances correspondantes
						probaT12.add(probas12[ind12]);
						disT12.add(dist12[ind12]);
						probaT13.add(probas13[ind13]);
						disT13.add(dist13[ind13]);
						probaT23.add(probas23[ind23]);
						disT23.add(dist23[ind23]);												
						}						
					ind23++;
					}
					ind13++;
				}
				ind12++;
				}
			}			
		else{
		int ind12 = 0;
		// 12 et 13
		for (ROI C1_2 : detections1C2) {
			int ind13 = 0;
			for (ROI C1_3 : detections1C3) {
				if (C1_2 == C1_3) // la detection 1 est avec 2 et 3
				{// on ajoute les spots correspondants
					T1.add(C1_2);
					T2.add(detections2C1[ind12]);
					T3.add(detections3C1[ind13]);
					// on y ajoute les probas distances correspondantes
					probaT12.add(probas12[ind12]);
					disT12.add(dist12[ind12]);
					probaT13.add(probas13[ind13]);
					disT13.add(dist13[ind13]);
					// on explore ensuite la liste des colocs 2C3 pour trouver
					// la proba/dist correspondante
					int ind23 = 0;
					boolean w23 = false;
					for (ROI C2_3 : detections2C3) {
						if ((C2_3 == detections2C1[ind12]) && (detections3C2[ind23] == detections3C1[ind13])) {
							w23 = true;
							disT23.add(dist23[ind23]);
							probaT23.add(probas23[ind23]);
							break;
						}
						ind23++;
					}
					// si on a pas trouvé de coloc, on met la proba à 0 pour ce
					// lien
					if (w23 == false) {
						disT23.add(0.);
						probaT23.add(0.);
					}
				}
				ind13++;
			}
			// on incrémente ensuite l'indice 12
			ind12++;
		}
		// on parcourt ensuite les trios avec les liens 2-1 et 2-3 mais sans 1-3
		int ind21 = 0;
		for (ROI C2_1 : detections2C1) {
			boolean trio_wo13 = false;
			int ind23 = 0;
			for (ROI C2_3 : detections2C3) {
				if (C2_1 == C2_3) // la detection 2 est avec 1 et 3
				{
					trio_wo13 = true;// pour l'instant!
					// on teste ensuite que le lien correspondant 1-3 est de
					// proba nulle, ie que le trio n'a pas encore été compté
					// pour cela,
					int ind13 = 0;
					for (ROI C1_3 : detections1C3) {
						if ((C1_3 == detections1C2[ind21]) && (detections3C1[ind13] == detections3C2[ind23])) {
							trio_wo13 = false;
							break;
						}
						ind13++;
					}
					if (trio_wo13 == true) {
						T1.add(detections1C2[ind21]);
						T2.add(C2_1);
						T3.add(detections3C2[ind23]);
						// on y ajoute les probas distances correspondantes
						probaT12.add(probas12[ind21]);
						disT12.add(dist12[ind21]);
						probaT13.add(0.);// pas de lien
						disT13.add(0.);
						probaT23.add(probas23[ind23]);
						disT23.add(dist23[ind23]);
					}
				}
				// on incrémente ensuite l'indice 32
				ind23++;
			}
			ind21++;
		}
		// reste à parcourir tous les trios 3-1 et 3-2 sans le lien1-2
		int ind31 = 0;
		for (ROI C3_1 : detections3C1) {
			boolean trio_wo12 = false;
			int ind32 = 0;
			for (ROI C3_2 : detections3C2) {
				if (C3_1 == C3_2) // la detection 3 est avec 1 et 2
				{
					trio_wo12 = true;
					// on teste ensuite que le lien correspondant 1-2 est de
					// proba nulle, ie que le trio n'a pas encore été compté
					// pour cela,
					ind12 = 0;
					for (ROI C1_2 : detections1C2) {
						if ((C1_2 == detections1C3[ind31]) && (detections2C1[ind12] == detections2C3[ind32])) {
							trio_wo12 = false;
							break;
						}
						ind12++;
					}
					if (trio_wo12 == true) {
						T1.add(detections1C3[ind31]);
						T2.add(detections2C3[ind32]);
						T3.add(C3_1);
						// on y ajoute les probas distances correspondantes
						probaT12.add(0.);// pas de lien
						disT12.add(0.);
						probaT13.add(probas13[ind31]);
						disT13.add(dist13[ind31]);
						probaT23.add(probas23[ind32]);
						disT23.add(dist23[ind32]);
					}
				}
				// on incrémente ensuite l'indice 32
				ind32++;
			}
			ind31++;
		}}
	}

	public double[] ListToArrayD(ArrayList<Double> liste) {
		double[] tab_double = new double[liste.size()];
		int ind = 0;
		for (Double r : liste) {
			tab_double[ind] = r;
			ind++;
		}
		return tab_double;
	}

	public ROI[] ListToArray(ArrayList<ROI> liste) {
		ROI[] tab_roi = new ROI[liste.size()];
		int ind = 0;
		for (ROI r : liste) {
			tab_roi[ind] = r;
			ind++;
		}
		return tab_roi;
	}

	@Override
	public String getMainPluginClassName() {
		// TODO Auto-generated method stub
		return SODAsuite.class.getName();
	}
}