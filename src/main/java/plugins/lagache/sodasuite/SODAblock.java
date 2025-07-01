package plugins.lagache.sodasuite;

import java.awt.geom.Point2D;
import java.util.ArrayList;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.WorkbookUtil;

import icy.file.FileUtil;
import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.roi.ROI3D;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import icy.type.point.Point5D;
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

public class SODAblock extends Plugin implements Block, PluginBundled
{

    VarROIArray detections1 = new VarROIArray("List of detections 1 (ROIs)");
    VarROIArray detections2 = new VarROIArray("List of detections 2 (ROIs)");

    // ROIs (spots) qui colocalisent en sortie

    VarROIArray detections1C = new VarROIArray("Coupled Detections 1");
    VarROIArray detections2C = new VarROIArray("Coupled Detections 2");

    VarDoubleArrayNative probas12 = new VarDoubleArrayNative("Coupling probabilities", null);
    VarDoubleArrayNative dist12 = new VarDoubleArrayNative("Coupling distances", null);

    VarROIArray detections1S = new VarROIArray("Single Detections 1");
    VarROIArray detections2S = new VarROIArray("Single Detections 2");

    VarROIArray ROIs = new VarROIArray("ROIs of Analysis (Cell's shape)");

    VarSequence input_sequence1 = new VarSequence("Input sequence 1", null);
    VarSequence input_sequence2 = new VarSequence("Input sequence 2", null);

    // VarDouble mindistance_in = new VarDouble("Min. Radius (in pixels)", 0);
    VarDouble maxdistance_in = new VarDouble("Max. Radius (in pixels)", 5);
    VarDouble Step = new VarDouble("Step", 1.0);
    VarBoolean manual = new VarBoolean("Fixed search distance", false);

    VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

    double[] results = null;
    int N_h;
    // pour le test stat
    double mindist;
    double maxdist;
    ArrayList<Double> distance_fit = new ArrayList<Double>();
    ArrayList<ROI> list_roi = new ArrayList<ROI>();
    double p_value = 0.0;
    double log_p_value = 0.0;

    @Override
    public void declareInput(VarList inputMap)
    {

        inputMap.add("Detections 1 (ROIs)", detections1);
        inputMap.add("Detections 2 (ROIs)", detections2);

        inputMap.add("Input Sequence 1", input_sequence1);
        inputMap.add("Input Sequence 2", input_sequence2);
        inputMap.add("ROIs of Analysis (Cells' Shapes)", ROIs);
        inputMap.add("Max. Radius (in x pixels)", maxdistance_in);
        inputMap.add("Min. Step", Step);
        inputMap.add("Manual estimation", manual);
    }

    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add("Couples (1)", detections1C);
        outputMap.add("Couples (2)", detections2C);
        outputMap.add("Coupling Probabilities", probas12);
        outputMap.add("Coupling Distances", dist12);
        outputMap.add("Single Detections 1", detections1S);
        outputMap.add("Single Detections 2", detections2S);
        outputMap.add("Workbook", book);
    }

    // N: nb entre dmin et dmax avec pas donne
    @Override
    public void run()
    {
        apparatedSpots.coloc1.clear();
        apparatedSpots.coloc2.clear();

        maxdist = maxdistance_in.getValue();

        if (maxdist <= 0)
        {
            new AnnounceFrame("Max. radius must be > 0");
            return;
        }
        if (detections1 == null || detections2 == null)
        {
            new AnnounceFrame("Please first select a detection set");
            return;
        }
        if (input_sequence1.getValue() == null)
        {
            new AnnounceFrame("There is no open sequence");
            return;
        }

        try {
            performAnalysis();
        }
        catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    private void performAnalysis() throws InterruptedException {
        // initialisation du workbook
        Workbook wb = new HSSFWorkbook();
        wb.setMissingCellPolicy(Row.MissingCellPolicy.CREATE_NULL_AS_BLANK);
        book.setValue(wb);

        // create the sheet
        String sheetName = "Coloc. Object Analysis";
        Sheet sheet = wb.getSheet(sheetName);
        if (sheet == null)
            sheet = wb.createSheet(sheetName);
        Sequence sequence = input_sequence1.getValue();
        int dim = 2;
        if (sequence.getSizeZ() > 1)
            dim = 3;

        // gestion des rois d'analyse en entr�e

        list_roi.clear();
        ROI[] roiArray = ROIs.getValue();// TODO Ou sont les ROIs???
        if (roiArray.length == 0)
        {
            for (int t = 0; t < input_sequence1.getValue().getSizeT(); t++)
            {
                ROI roi = null;

                ROI2DRectangle r = new ROI2DRectangle(sequence.getBounds2D());
                for (int h = 0; h < sequence.getSizeZ(); h++)
                {
                    r.setZ(h);
                    r.setT(t);
                    roi = r.getUnion(roi);
                }
                list_roi.add(roi);
            }
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

        // create the header row

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
        // d�finition du ratio Z/X
        double ratio_zx = 1.0;
        if (sequence.getSizeZ() > 1)
        {
            ratio_zx = sequence.getPixelSizeZ() / sequence.getPixelSizeX();
        }
        // initialisation boucle en temps
        // temps total de la sequence
        int T = sequence.getSizeT();
        betaCorrection(maxdist / 10, 100);
        /////////////////////////////////////////////

        double step_min = Step.getValue();// 1.0;//maxdist / 10;

        distance_fit.clear();
        /*
         * distance_fit.add((double)mindistance_in.getValue()); double
         * temp=mindistance_in.getValue();
         */
        distance_fit.add(0.);
        double temp = 0.;
        while (temp + step_min <= maxdist)
        {
            temp += step_min;
            distance_fit.add(temp);
        }
        int N_fit = distance_fit.size();
        if (N_fit == 1)
        {
            distance_fit.add(maxdist);
            N_fit = distance_fit.size();
        }
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////

        // stockage de l'histogram des distances de coloc
        double[][] proba_dist = new double[T][N_fit - 1];
        int ind_row = 0;
        // listes pour stocker tous les couples et les ROI single au fur et a mesure
        ArrayList<ROI> c_1 = new ArrayList<ROI>();
        ArrayList<ROI> c_2 = new ArrayList<ROI>();
        ArrayList<ROI> s_1 = new ArrayList<ROI>();
        ArrayList<ROI> s_2 = new ArrayList<ROI>();

        ArrayList<Double> dist = new ArrayList<Double>();
        ArrayList<Double> probas = new ArrayList<Double>();
        // d�marrage boucle en temps
        for (int t = 0; t < T; t += 1)
        {
            ind_row++;
            Row row = sheet.createRow(ind_row);
            // initialisation
            ArrayList<ROI> roi_t = ROI_t(dim, list_roi, t, -1);
            if (roi_t.size() == 0)
            {
            }
            else
            {
                double window_size = maxdist;
                double max_z = maxdist / ratio_zx;
                ROI[] detection = detectionsInRoi(detections1.getValue(), roi_t);
                ROI[] detection2 = detectionsInRoi(detections2.getValue(), roi_t);

                // on crée des listes de points 2D ou 3D avec les positions (centres
                // d'intensité ou de masse en 3D) des spots
                ArrayList<Point3D> positions_1_2D = new ArrayList<Point3D>();
                ArrayList<Point3D> positions_2_2D = new ArrayList<Point3D>();
                // on crée aussi une liste potentielle de points 3d
                ArrayList<Point3D> positions_1_3D = new ArrayList<Point3D>();
                ArrayList<Point3D> positions_2_3D = new ArrayList<Point3D>();
                // sous-listes des detections pour accelerer le temps de calcul
                Window2D[][] imageWindows_1_2D = new Window2D[(int) (sequence.getWidth()
                        / window_size)][(int) (sequence.getHeight() / window_size)];
                Window2D[][] imageWindows_2_2D = new Window2D[(int) (sequence.getWidth()
                        / window_size)][(int) (sequence.getHeight() / window_size)];
                Window3D[][][] imageWindows_1_3D = new Window3D[(int) (sequence.getWidth()
                        / window_size)][(int) (sequence.getHeight() / window_size)][(int) (sequence.getSizeZ()
                                / max_z)];
                Window3D[][][] imageWindows_2_3D = new Window3D[(int) (sequence.getWidth()
                        / window_size)][(int) (sequence.getHeight() / window_size)][(int) (sequence.getSizeZ()
                                / max_z)];
                /////////////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////
                if (dim == 2)
                {
                    // on boucle sur detection puis sur detection2
                    for (ROI r : detection)
                    {
                        Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, input_sequence1.getValue());
                        positions_1_2D.add(pos);
                    }
                    for (ROI r : detection2)
                    {
                        Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, input_sequence2.getValue());
                        positions_2_2D.add(pos);
                    }
                    // on crée ensuite toutes les sous-listes de detections

                    imageWindows_1_2D = Window2D.window_tab(positions_1_2D, input_sequence1.getValue(), window_size);
                    imageWindows_2_2D = Window2D.window_tab(positions_2_2D, input_sequence2.getValue(), window_size);
                }
                else
                {
                    // on boucle sur detection puis sur detection2
                    for (ROI r : detection)
                    {
                        if (r.getDimension() == 3)
                        {
                            Point3D pos = Ripley3D.getIntensityCenter((ROI3D) r, input_sequence1.getValue());
                            positions_1_3D.add(pos);
                        }
                        else
                        {
                            Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, input_sequence1.getValue());
                            positions_1_3D.add(pos);
                        }
                    }
                    for (ROI r : detection2)
                    {
                        if (r.getDimension() == 3)
                        {
                            Point3D pos = Ripley3D.getIntensityCenter((ROI3D) r, input_sequence2.getValue());
                            positions_2_3D.add(pos);
                        }
                        else
                        {
                            Point3D pos = Ripley2D.getIntensityCenter((ROI2D) r, input_sequence2.getValue());
                            positions_2_3D.add(pos);
                        }
                    }
                    // on crée ensuite toutes les sous-listes de detections
                    imageWindows_1_3D = Window3D.window_tab(positions_1_3D, sequence, maxdist, max_z);
                    imageWindows_2_3D = Window3D.window_tab(positions_2_3D, sequence, maxdist, max_z);
                }
                // calcul des parametres globaux: aire totale des rois et nb de detections
                // faut-refaire une fonction calcul de volume
                double volume = ROIUtil.getUnion(roi_t).getNumberOfPoints() * ratio_zx;
                if (ROIUtil.getUnion(roi_t).getBounds5D().isInfiniteZ())
                {
                    volume = volume * sequence.getSizeZ() * ratio_zx;
                }
                int nbdeta = detection.length;
                int nbdetb = detection2.length;
                String dataSetName = "unknown sequence";
                if (input_sequence1.getValue() != null)
                {
                    dataSetName = getDataSetName(input_sequence1.getValue());
                }
                // déclaration des tableaux de p-value
                double[] maximum = new double[N_fit - 1];
                double[] pvalue = new double[N_fit - 1];
                double[] log_pvalue = new double[N_fit - 1];

                // esperance du nb de coloc par anneau
                double[] col_ev = new double[N_fit - 1];
                // nb total de coloc
                if (manual.getValue())
                {
                    // on remplit les probas de coloc à 1 jusqu'à dist max
                    for (int a = 0; a < N_fit - 1; a++)
                    {
                        proba_dist[t][a] = 1.0;
                    }
                    if (dim == 2)
                    {
                        manual_computation2D(positions_1_2D, positions_2_2D, distance_fit, col_ev);
                    }
                    else
                    {
                        manual_computation3D(positions_1_3D, positions_2_3D, distance_fit, col_ev, ratio_zx);
                    }
                }
                else
                {
                    ////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////
                    double[] delta_K = new double[N_fit - 1];
                    double[][] Ktemp = new double[N_fit - 1][3];
                    ////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////

                    if (sequence.getSizeZ() > 1)
                    {
                        Ktemp = Ripley3D.correlation_new(roi_t, imageWindows_1_3D, imageWindows_2_3D, volume, N_fit,
                                distance_fit, input_sequence1.getValue(), nbdeta, nbdetb, ratio_zx);
                    }
                    else
                    {
                        Ktemp = Ripley2D.correlation_new(roi_t, imageWindows_1_2D, imageWindows_2_2D, volume, N_fit,
                                distance_fit, nbdeta, nbdetb);
                    }

                    for (int k = 0; k < N_fit - 1; k++)
                    {
                        delta_K[k] = Ktemp[k][0];
                    }

                    if (dim == 2)
                    {
                        non_parametric_object.main2D_corr(N_fit, sequence, roi_t, positions_1_2D, positions_2_2D,
                                imageWindows_1_2D, imageWindows_2_2D, distance_fit, N_h, volume, results, delta_K,
                                Ktemp, proba_dist[t], maximum, pvalue);
                    }
                    // calcul de l'esperance du nb de couples par intervalle et du
                    // nb total de couples
                    else
                    {
                        non_parametric_object.main3D_corr(N_fit, sequence, roi_t, positions_1_3D, positions_2_3D,
                                imageWindows_1_3D, imageWindows_2_3D, distance_fit, N_h, ratio_zx, volume, results,
                                delta_K, Ktemp, proba_dist[t], maximum, pvalue);
                    }
                    int total_coloc_temp = 0;
                    for (int j = 0; j < N_fit - 1; j++)
                    {
                        col_ev[j] = proba_dist[t][j] * (nbdetb * nbdeta / volume) * delta_K[j];
                        total_coloc_temp += col_ev[j];
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
                ArrayList<apparatedSpots> liste_app_detect = apparatedSpots.appDetectConstruction(detection, detection2,
                        proba_dist[t], distance_fit, input_sequence1.getValue(), input_sequence2.getValue());
                for (apparatedSpots as : liste_app_detect)
                {
                    c_1.add(as.s1);
                    c_2.add(as.s2);
                    dist.add(as.distance);
                    probas.add(as.proba);
                }
                apparatedSpots.roiSingle(t, liste_app_detect, detection, detection2);
                for (ROI s : apparatedSpots.single1)
                {
                    s_1.add(s);
                }
                for (ROI s : apparatedSpots.single2)
                {
                    s_2.add(s);
                }
                // export excel des paramètres de couplage
                for (int j = 0; j < N_fit - 1; j++)
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
                        if (proba_dist[t][i] > 0)
                        {
                            nb_total_couples_temp += col_ev[i] / proba_dist[t][i];
                        }
                    }
                    // nb de single 1 et 2 à cette distance
                    int[] nb_single = apparatedSpots.nbSingle(t, liste_app_detect, detection, detection2,
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
                    double dist_moy = apparatedSpots.distance_moyenne(liste_app_detect, distance_fit.get(j + 1));
                    row.getCell(18).setCellValue(dist_moy);
                    // on crée une nouvelle ligne dans la feuille excel
                    ind_row++;
                    row = sheet.createRow(ind_row);
                }
            }
        }
        ROI[] c_1_tab = new ROI[c_1.size()];
        ROI[] c_2_tab = new ROI[c_2.size()];
        double[] dist_tab = new double[dist.size()];
        double[] probas_tab = new double[probas.size()];
        for (int i = 0; i < c_1.size(); i++)
        {
            c_1_tab[i] = c_1.get(i);
            c_2_tab[i] = c_2.get(i);
            dist_tab[i] = dist.get(i);
            probas_tab[i] = probas.get(i);
        }

        detections1C.setValue(c_1_tab);
        detections2C.setValue(c_2_tab);
        probas12.setValue(probas_tab);
        dist12.setValue(dist_tab);
        // calcul des ROIs single
        ROI[] s_1_tab = new ROI[s_1.size()];
        ROI[] s_2_tab = new ROI[s_2.size()];

        for (int i = 0; i < s_1.size(); i++)
        {
            s_1_tab[i] = s_1.get(i);
        }
        for (int i = 0; i < s_2.size(); i++)
        {
            s_2_tab[i] = s_2.get(i);
        }
        detections1S.setValue(s_1_tab);
        detections2S.setValue(s_2_tab);//

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

    public static ROI[] detectionsInRoi(ROI[] detections, ArrayList<ROI> roi) throws InterruptedException {
        // SEQUENCE A remplacer par w,h

        // create a hashMap with the detections binded to ROI

        if (roi.size() == 0)
        {
            return new ROI[0];
        }
        else
        {
            ArrayList<ROI> ROIDetection = new ArrayList<ROI>();
            ArrayList<ROI> TestList = new ArrayList<ROI>();

            ROI union = ROIUtil.getUnion(roi);
            Point5D pos = union.getPosition5D();
            // fill hashMap
            for (ROI r : detections)
            {
                if (r == null)
                {
                }
                else
                {
                    boolean b = false;
                    for (ROI rt : TestList)
                    {
                        if (r == rt)
                        {
                            b = true;
                        }
                    }
                    if (b == false) // si le spot n'a pas encore �t� pris en compte
                    {
                        TestList.add(r);
                        if (union.contains(r.getPosition5D().getX(), r.getPosition5D().getY(), r.getPosition5D().getZ(),
                                r.getPosition5D().getT(), pos.getC()))
                        {
                            ROIDetection.add(r);
                        }
                    }
                }
            }
            ROI[] ROIDetection_ = new ROI[ROIDetection.size()];
            int ind = 0;
            for (ROI r : ROIDetection)
            {
                ROIDetection_[ind] = r;
                ind++;
            }
            return (ROIDetection_);
        }
    }

    public static ArrayList<ROI> ROI_t(int dim, ArrayList<ROI> list_roi, int t, int c)
    {

        ArrayList<ROI> list_roi_copy = new ArrayList<ROI>();
        for (ROI r : list_roi)
        {
            list_roi_copy.add(r.getCopy());
        }
        // check compatibility en z si dim =3
        if (dim == 3)
        {
            boolean one_z = false;
            boolean all_z = true;
            for (ROI r0 : list_roi_copy)
            {
                if (r0.getBounds5D().isInfiniteZ() == false)
                {
                    all_z = false;
                }
                if (r0.getBounds5D().isInfiniteZ())
                {
                    one_z = true;
                }
            }
            // gestion de l'exception
            if (one_z == true && all_z == false)
            {
                return null;
            }
        }

        ArrayList<ROI> roi_t = new ArrayList<ROI>();

        for (ROI r : list_roi_copy)
        {
            Point5D pt0 = r.getPosition5D();
            if (dim == 2)
            {
                pt0.setZ(-1);
            }
            if ((r.getBounds5D().isInfiniteC()) || (c == -1))
            {
                pt0.setC(c);
            }
            if ((r.getBounds5D().isInfiniteT()) || (t == -1))
            {
                pt0.setT(t);
            }
            r.setPosition5D(pt0);
            if ((pt0.getT() == t) && (pt0.getC() == c))
            {
                roi_t.add(r);
            }
        }
        return roi_t;
    }

    public static void manual_computation2D(ArrayList<Point3D> positions1, ArrayList<Point3D> positions2,
            ArrayList<Double> distances, double[] col_ev)
    {
        for (Point3D s1Pos : positions1)
        {
            for (Point3D s2Pos : positions2)
            {
                double distance_12 = s1Pos.distance(s2Pos);
                for (int i = 1; i < distances.size(); i++)
                {
                    if (distance_12 < distances.get(i))
                    {
                        col_ev[i - 1] += 1;
                        break;
                    }
                }
            }
        }
    }

    public static void manual_computation3D(ArrayList<Point3D> positions1, ArrayList<Point3D> positions2,
            ArrayList<Double> distances, double[] col_ev, double ratio_zx)
    {
        for (Point3D s1Pos : positions1)
        {
            for (Point3D s2Pos : positions2)
            {
                double distance_12 = Math
                        .sqrt(Math.pow(s1Pos.getX() - s2Pos.getX(), 2) + Math.pow(s1Pos.getY() - s2Pos.getY(), 2)
                                + Math.pow((s1Pos.getZ() - s2Pos.getZ()) * ratio_zx, 2));
                for (int i = 1; i < distances.size(); i++)
                {
                    if (distance_12 < distances.get(i))
                    {
                        col_ev[i - 1] += 1;
                        break;
                    }
                }
            }
        }
    }

    private String getDataSetName(Sequence sequence)
    {
        String dataSetName = "";
        // replace the sheet name by the file or sequence name
        dataSetName = FileUtil.getFileName(sequence.getFilename());
        if (dataSetName.isEmpty())
            dataSetName = sequence.getName();

        // make the name "safe"
        return WorkbookUtil.createSafeSheetName(dataSetName);
    }

    @Override
    public String getMainPluginClassName()
    {
        return SODAsuite.class.getName();
    }
}