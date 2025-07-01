package plugins.lagache.sodasuite;

import java.util.ArrayList;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.ROI;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarDouble;
import plugins.adufour.vars.lang.VarDoubleArrayNative;
import plugins.adufour.vars.lang.VarROIArray;
import plugins.adufour.vars.lang.VarSequence;
import plugins.adufour.vars.lang.VarWorkbook;
import plugins.kernel.roi.roi2d.ROI2DRectangle;

// Colocalisation with Ripley function K
// Significant

public class SODAStorm2D extends Plugin implements Block, PluginBundled
{

    // EzVarSwimmingObject<DetectionResult> detections = new
    // EzVarSwimmingObject<DetectionResult>("Detections");
    VarDoubleArrayNative detections1_x = new VarDoubleArrayNative("Detections 1 (x)", null);
    VarDoubleArrayNative detections1_y = new VarDoubleArrayNative("Detections 1 (y)", null);

    VarDoubleArrayNative detections2_x = new VarDoubleArrayNative("Detections 2 (x)", null);
    VarDoubleArrayNative detections2_y = new VarDoubleArrayNative("Detections 2 (y)", null);

    VarDoubleArrayNative detections1C_x = new VarDoubleArrayNative("Coupled Detections 1 (x)", null);
    VarDoubleArrayNative detections1C_y = new VarDoubleArrayNative("Coupled Detections 1 (y)", null);
    VarDoubleArrayNative detections2C_x = new VarDoubleArrayNative("Coupled Detections 2 (x)", null);
    VarDoubleArrayNative detections2C_y = new VarDoubleArrayNative("Coupled Detections 2 (y)", null);
    VarDoubleArrayNative proba12 = new VarDoubleArrayNative("Coupling probabilities", null);
    VarDoubleArrayNative dist12 = new VarDoubleArrayNative("Coupling distances", null);
    VarDoubleArrayNative detections1S_x = new VarDoubleArrayNative("Single Detections 1 (x)", null);
    VarDoubleArrayNative detections1S_y = new VarDoubleArrayNative("Single Detections 1 (y)", null);
    VarDoubleArrayNative detections2S_x = new VarDoubleArrayNative("Single Detections 2 (x)", null);
    VarDoubleArrayNative detections2S_y = new VarDoubleArrayNative("Single Detections 2 (y)", null);

    VarSequence input_sequence = new VarSequence("Sequence", null);
    VarROIArray ROIs = new VarROIArray("ROIs of Analysis (Cell's shape)");
    VarBoolean correction = new VarBoolean("Correction", false);

    VarDouble mindistance_in = new VarDouble("Min. Radius (in pixels)", 0);
    VarDouble maxdistance_in = new VarDouble("Max. Radius (in pixels)", 5);
    VarDouble Step = new VarDouble("Step", 1.0);
    VarBoolean manual = new VarBoolean("Fixed search distance", false);

    VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

    double[] results = null;
    int N_h;
    // pour le test stat
    double maxdist;
    ArrayList<Double> distance_fit = new ArrayList<Double>();
    ArrayList<ROI> list_roi = new ArrayList<ROI>();

    @Override
    public void declareInput(VarList inputMap)
    {

        inputMap.add("Detections 1 (x)", detections1_x);
        inputMap.add("Detections 1 (y)", detections1_y);

        inputMap.add("Detections 2 (x)", detections2_x);
        inputMap.add("Detections 2 (y)", detections2_y);

        inputMap.add("Sequence", input_sequence);
        inputMap.add("ROIs of Analysis (Cells' Shapes)", ROIs);
        inputMap.add("Min. Radius (in x pixels)", mindistance_in);
        inputMap.add("Max. Radius (in x pixels)", maxdistance_in);
        inputMap.add("Min. Step", Step);
        inputMap.add("Manual estimation", manual);
        inputMap.add("Correction", correction);
    }

    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add("Workbook", book);
        outputMap.add("Coupled Detections 1 (x)", detections1C_x);
        outputMap.add("Coupled Detections 1 (y)", detections1C_y);
        outputMap.add("Coupled Detections 2 (x)", detections2C_x);
        outputMap.add("Coupled Detections 2 (y)", detections2C_y);
        outputMap.add("Coupling probabilities", proba12);
        outputMap.add("Coupling disatnces", dist12);
        outputMap.add("Single Detections 1 (x)", detections1S_x);
        outputMap.add("Single Detections 1 (y)", detections1S_y);
        outputMap.add("Single Detections 2 (x)", detections2S_x);
        outputMap.add("Single Detections 2 (y)", detections2S_y);
    }

    // N: nb entre dmin et dmax avec pas donne
    @Override
    public void run()
    {
        maxdist = maxdistance_in.getValue();
        try {
            performAnalysis();
        }
        catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    private void performAnalysis() throws InterruptedException {
        // on construit les 2 listes de Points 3D correspondant aux coordonnées en entrée
        ArrayList<Point3D> detections1 = new ArrayList<Point3D>();
        ArrayList<Point3D> detections2 = new ArrayList<Point3D>();
        for (int i = 0; i < detections1_x.getValue().length; i++)
        {
            Point3D pt = new Point3D.Double(detections1_x.getValue()[i], detections1_y.getValue()[i], 0.0);
            detections1.add(pt);
        }
        for (int i = 0; i < detections2_x.getValue().length; i++)
        {
            Point3D pt = new Point3D.Double(detections2_x.getValue()[i], detections2_y.getValue()[i], 0.0);
            detections2.add(pt);
        }

        // initialisation du workbook
        Workbook wb = new HSSFWorkbook();
        wb.setMissingCellPolicy(Row.MissingCellPolicy.CREATE_NULL_AS_BLANK);
        book.setValue(wb);

        // create the sheet
        String sheetName = "Coloc. Object Analysis";
        Sheet sheet = wb.getSheet(sheetName);
        if (sheet == null)
            sheet = wb.createSheet(sheetName);

        // gestion des rois d'analyse en entr�e

        Sequence sequence = input_sequence.getValue();
        list_roi.clear();
        ROI[] roiArray = ROIs.getValue();
        if (roiArray.length == 0)
        {
            ROI roi = new ROI2DRectangle(sequence.getBounds2D());
            list_roi.add(roi);
        }
        else
        {
            for (ROI roi : roiArray)
            {
                list_roi.add(roi);
            }
        }

        // on teste si les dimension en c sonts compatible
        double c1 = list_roi.get(0).getPosition5D().getC();
        for (ROI r : list_roi)
        {
            double c = r.getPosition5D().getC();
            if ((c != c1) && (c != (-1)))
            {
                new AnnounceFrame("ROI channels are incompatibles");
                return;
            }
        }

        // on teste si les dimensions en temps/z sont incompatibles pour une
        // union
        boolean one_z = false;
        boolean all_z = true;

        for (ROI r : list_roi)
        {
            if (r.getBounds5D().isInfiniteZ() == false)
            {
                all_z = false;
            }
            if (r.getBounds5D().isInfiniteZ())
            {
                one_z = true;
            }
        }
        // gestion de l'exception
        if (one_z == true && all_z == false)
        {
            new AnnounceFrame("Incompatibility in Z dimensions between ROIs");
            return;
        }

        Row header = sheet.getRow(0);
        if (header == null)
        {
            header = sheet.createRow(0);
            header.getCell(0).setCellValue("Sequence Name");
            header.getCell(1).setCellValue("Time");
            header.getCell(2).setCellValue("Nb detections 1");
            header.getCell(3).setCellValue("Nb detections 2");
            header.getCell(4).setCellValue("r 1");
            header.getCell(5).setCellValue("r max");
            header.getCell(6).setCellValue("Maximum of the K function");
            header.getCell(7).setCellValue("p value");
            header.getCell(8).setCellValue("log_10(p value)");
            header.getCell(9).setCellValue("nb single 1");
            header.getCell(10).setCellValue("nb single 2");
            header.getCell(11).setCellValue("nb coupled 1");
            header.getCell(12).setCellValue("Mean number of partners 2");
            header.getCell(13).setCellValue("nb coupled 2");
            header.getCell(14).setCellValue("Mean number of partners 1");
            header.getCell(15).setCellValue("Mean coupling probability");
            header.getCell(16).setCellValue("Coupling index 1");
            header.getCell(17).setCellValue("Coupling index 2");
            header.getCell(18).setCellValue("Mean Coupling distance (pixels)"); // médiane
            // header.getCell(11).setCellValue("Std. of coupling distance");
        }
        // initialisation boucle en temps
        // temps total de la sequence
        betaCorrection(maxdist / 10, 100);
        /////////////////////////////////////////////
        ///////////////////////////////////////////////////
        // on construit emsuite la roi au temps t=0

        ArrayList<Point3D> detection = detectionsInRoi(detections1, list_roi);
        ArrayList<Point3D> detection2 = detectionsInRoi(detections2, list_roi);
        // on trie ensuite les loc avec la methode sparse

        double volume = ROIUtil.getUnion(list_roi).getNumberOfPoints();

        int nbdeta = detection.size();
        int nbdetb = detection2.size();
        // construction du distance_tab pour calculer la p_value de coloc
        double coeff_radius = 1;
        double step_min = Step.getValue();// 1.0;//maxdist / 10;
        distance_fit.clear();

        distance_fit.add(mindistance_in.getValue());
        double temp = distance_fit.get(0);

        while (temp + step_min <= maxdist)
        {
            temp = temp + step_min;
            distance_fit.add(temp);
        }
        int N_fit = distance_fit.size();
        if (N_fit == 1)
        {
            distance_fit.add(maxdist);
            N_fit = distance_fit.size();
            // new AnnounceFrame("Number of detections is insufficient for
            // robust analysis");
        }
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////
        // on ajoute ensuite virtuellement n=5 pas de distance à N_fit pour gérer les colocs pouvant provenir des anneaux supérieurs à rmax
        temp = distance_fit.get(N_fit - 1);
        int tampon = 2;
        for (int a = 1; a <= tampon; a++)
        {
            distance_fit.add(temp + a * step_min);
        }
        // actulaisation du nouveau N_fit (N_fit+tampon)
        N_fit = distance_fit.size();

        // stockage de l'histogram des distances de coloc
        double[] proba_dist = new double[N_fit - 1];
        int ind_row = 0;
        ind_row++;
        Row row = sheet.createRow(ind_row);
        // on crée ensuite toutes les sous-listes de detections
        double window_size = maxdist;
        Window2D[][] imageWindows_1 = new Window2D[(int) (sequence.getWidth()
                / window_size)][(int) (sequence.getHeight() / window_size)];
        imageWindows_1 = Window2D.window_tab(detection, sequence, maxdist);
        Window2D[][] imageWindows_2 = new Window2D[(int) (sequence.getWidth()
                / window_size)][(int) (sequence.getHeight() / window_size)];
        imageWindows_2 = Window2D.window_tab(detection2, sequence, maxdist);
        nbdeta = detection.size();
        nbdetb = detection2.size();

        String dataSetName = "unknown sequence";
        // esperance du nb de coloc par anneau
        double[] col_ev = new double[N_fit - 1];
        double total_coloc_temp = 0.;
        double[] maximum = new double[N_fit - 1];
        double[] pvalue = new double[N_fit - 1];
        double[] log_pvalue = new double[N_fit - 1];
        if (manual.getValue())
        {
            // on remplit les probas de coloc à 1 jusqu'à dist max
            for (int a = 0; a < N_fit - 1; a++)
            {
                proba_dist[a] = 1.0;
            }
            SODAblock.manual_computation2D(detection, detections2, distance_fit, col_ev);

        }
        else
        {
            ////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////
            double[] delta_K = new double[N_fit - 1];
            double[][] Ktemp = new double[N_fit - 1][3];
            Ktemp = Ripley2D.correlation_new(list_roi, imageWindows_1, imageWindows_2, volume, N_fit, distance_fit,
                    nbdeta, nbdetb);

            for (int k = 0; k < N_fit - 1; k++)
            {
                delta_K[k] = Ktemp[k][0];
            }

            non_parametric_object.main2D_corr(N_fit, sequence, list_roi, detection, detection2, imageWindows_1,
                    imageWindows_2, distance_fit, N_h, volume, results, delta_K, Ktemp, proba_dist, maximum, pvalue);
            // calcul de l'esperance du nb de couples par intervalle et du
            // nb total de couples

            for (int j = 0; j < N_fit - 1; j++)
            {
                col_ev[j] = proba_dist[j] * (nbdetb * nbdeta / volume) * delta_K[j];
                total_coloc_temp += col_ev[j];
                // calcul du log de la p_value\
                if (maximum[j] < 4.)
                {
                    log_pvalue[j] = Math.log10(pvalue[j]);
                }
                else
                {
                    double x = maximum[j] / Math.sqrt(2.);
                    double y = (N_fit - 1) / (2 * Math.sqrt(Math.PI) * x);
                    double exponent = (Math.log(y) - Math.pow(x, 2)) / Math.log(10);
                    log_pvalue[j] = exponent;
                }
            }
        }

        ArrayList<apparatedLocalizations> liste_app_detect = apparatedLocalizations.appDetectConstruction(detection,
                detection2, proba_dist, distance_fit, 3);
        double[] c_1_x = new double[liste_app_detect.size()];
        double[] c_1_y = new double[liste_app_detect.size()];
        double[] c_2_x = new double[liste_app_detect.size()];
        double[] c_2_y = new double[liste_app_detect.size()];
        double[] dist = new double[liste_app_detect.size()];
        double[] probas = new double[liste_app_detect.size()];
        int ind = 0;
        for (apparatedLocalizations as : liste_app_detect)
        {
            c_1_x[ind] = as.p1.getX();
            c_1_y[ind] = as.p1.getY();
            c_2_x[ind] = as.p2.getX();
            c_2_y[ind] = as.p2.getY();
            dist[ind] = as.distance;
            probas[ind] = as.proba;
            ind++;
        }
        detections1C_x.setValue(c_1_x);
        detections1C_y.setValue(c_1_y);
        detections2C_x.setValue(c_2_x);
        detections2C_y.setValue(c_2_y);
        proba12.setValue(probas);
        dist12.setValue(dist);
        // calcul des ROIs single
        apparatedLocalizations.locSingle(liste_app_detect, detection, detection2);
        double[] s_1_x = new double[apparatedLocalizations.single1.size()];
        double[] s_1_y = new double[apparatedLocalizations.single1.size()];
        double[] s_2_x = new double[apparatedLocalizations.single2.size()];
        double[] s_2_y = new double[apparatedLocalizations.single2.size()];

        for (int i = 0; i < apparatedLocalizations.single1.size(); i++)
        {
            s_1_x[i] = apparatedLocalizations.single1.get(i).getX();
            s_1_y[i] = apparatedLocalizations.single1.get(i).getY();
        }
        for (int j = 0; j < apparatedLocalizations.single2.size(); j++)
        {
            s_2_x[j] = apparatedLocalizations.single2.get(j).getX();
            s_2_y[j] = apparatedLocalizations.single2.get(j).getY();
        }
        detections1S_x.setValue(s_1_x);
        detections1S_y.setValue(s_1_y);
        detections2S_x.setValue(s_2_x);
        detections2S_y.setValue(s_2_y);

        // export excel des paramètres de couplage
        for (int j = 0; j < N_fit - 1 - tampon; j++)
        {

            row.getCell(0).setCellValue(dataSetName);
            row.getCell(1).setCellValue(0);
            row.getCell(2).setCellValue(nbdeta);
            row.getCell(3).setCellValue(nbdetb);
            row.getCell(4).setCellValue(distance_fit.get(0));
            row.getCell(5).setCellValue(distance_fit.get(j + 1));
            // p_value et log p value pour la derniere ligne
            row.getCell(6).setCellValue(maximum[j]);
            row.getCell(7).setCellValue(pvalue[j]);
            row.getCell(8).setCellValue(log_pvalue[j]);
            // esperance et nb total de couples entre 0 et
            // distance_fit.get(j + 1)
            int nb_total_couples_temp = 0;
            int esp_total_couples_temp = 0;
            for (int i = 0; i < j + 1; i++)
            {
                esp_total_couples_temp += col_ev[i];
                if (proba_dist[i] > 0)
                {
                    nb_total_couples_temp += col_ev[i] / proba_dist[i];
                }
            }
            // nb de single 1 et 2 à cette distance
            int[] nb_single = apparatedLocalizations.nbSingle(liste_app_detect, detection, detection2,
                    distance_fit.get(j + 1));
            row.getCell(9).setCellValue(nb_single[0]);
            row.getCell(10).setCellValue(nb_single[1]);
            // parametres de couplage pour molecules 1
            // nb de molecules 1 couplées
            row.getCell(11).setCellValue(nbdeta - nb_single[0]);
            // nb total de partenaires
            if (nbdeta - nb_single[0] > 0)
            {
                row.getCell(12).setCellValue((double) nb_total_couples_temp / (double) (nbdeta - nb_single[0]));
            }
            else
            {
                row.getCell(12).setCellValue(0.0);
            }
            // nb de molecules 2 couplées
            row.getCell(13).setCellValue(nbdetb - nb_single[1]);
            // nb total de partenaires
            if (nbdetb - nb_single[1] > 0)
            {
                row.getCell(14).setCellValue((double) nb_total_couples_temp / (double) (nbdetb - nb_single[1]));
            }
            else
            {
                row.getCell(14).setCellValue(0.0);
            }
            // proba moyenne de couplage
            if (nb_total_couples_temp > 0)
            {
                row.getCell(15).setCellValue(esp_total_couples_temp / (double) nb_total_couples_temp);
            }
            else
            {
                row.getCell(15).setCellValue(0.0);
            }
            // coupling index 1 et 2
            double index1 = esp_total_couples_temp / (double) nbdeta;
            double index2 = esp_total_couples_temp / (double) nbdetb;
            row.getCell(16).setCellValue(index1);
            row.getCell(17).setCellValue(index2);
            // calcul de la distance moyenne
            double dist_moy = apparatedLocalizations.distance_moyenne(liste_app_detect, distance_fit.get(j + 1));
            row.getCell(18).setCellValue(dist_moy);
            // on crée une nouvelle ligne dans la feuille excel
            ind_row++;
            row = sheet.createRow(ind_row);
        }

    }

    private void betaCorrection(double pas, int nbN)
    {

        Double valN = 1 / (1 / (double) nbN);
        N_h = valN.intValue() + 1; // / ATTENTION AUX BORNES

        double[] alpha = new double[N_h + 1];
        results = new double[N_h + 1];

        for (int i = 0; i < results.length; i++)
        {
            results[i] = 0;
            alpha[i] = 0.0d;

        }

        for (int i = 1; i < results.length; i++)
        {
            alpha[i] = (i / (double) N_h);
        }
        for (int i = 0; i < results.length; i++)
        {
            double h, j;
            for (h = alpha[i] + (pas), j = 2; h <= 1; j++)
            {
                results[i] = results[i] + h * pas / (1 - 1 / Math.PI * Math.acos(alpha[i] / h));
                h = alpha[i] + (j * pas);
            }

            results[i] = results[i] * 2 + alpha[i] * alpha[i];
        }

    }

    private ArrayList<Point3D> detectionsInRoi(ArrayList<Point3D> detections, ArrayList<ROI> roi) throws InterruptedException {

        ArrayList<Point3D> ROIDetection = new ArrayList<Point3D>();
        ROI union = ROIUtil.getUnion(roi);
        for (Point3D pt : detections)
        {
            if (union.contains(pt.getX(), pt.getY(), union.getPosition5D().getZ(), union.getPosition5D().getT(),
                    union.getPosition5D().getC()))
            {
                ROIDetection.add(pt);
            }
        }
        return (ROIDetection);
    }

    @Override
    public String getMainPluginClassName()
    {
        return SODAsuite.class.getName();
    }
}