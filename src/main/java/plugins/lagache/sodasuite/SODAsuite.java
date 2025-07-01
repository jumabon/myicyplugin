package plugins.lagache.sodasuite;

import java.util.ArrayList;

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.ROI;
import icy.roi.ROIUtil;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.*;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

public class SODAsuite extends Plugin implements Block
{

    VarSequence sequence1 = new VarSequence("Sequence 1", null);
    VarSequence sequence2 = new VarSequence("Sequence 2", null);
    VarROIArray cell_mask = new VarROIArray("Cell mask");
    VarROIArray spots1 = new VarROIArray("Spots 1 (ROIS)");
    VarROIArray spots2 = new VarROIArray("Spots 2 (ROIS)");
    VarBoolean PCC = new VarBoolean("PCC (Pearson)", true);
    VarBoolean Manders = new VarBoolean("Manders analysis ", true);
    VarInteger MC = new VarInteger("Number of MC simulations", 100);
    VarBoolean Overlap = new VarBoolean("Overlap (# with overlap>T) ", true);
    VarDouble threshold = new VarDouble("T", 0.5);
    VarBoolean distance = new VarBoolean("Distance based", true);
    VarBoolean soda = new VarBoolean("SODA", true);
    VarDouble rmax = new VarDouble("r max (pixels)", 5.0);
    VarDouble step = new VarDouble("step (pixels)", 1.0);

    ArrayList<ROI> list_roi = new ArrayList<ROI>();
    VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

    @Override
    public void declareInput(VarList inputMap)
    {

        inputMap.add("Sequence 1", sequence1);
        inputMap.add("Sequence 2", sequence2);
        inputMap.add("Cell mask (ROIs)", cell_mask);
        inputMap.add("Spots 1 (ROIs)", spots1);
        inputMap.add("Spots 2 (ROIs)", spots2);
        inputMap.add("PCC (Pearson)", PCC);
        inputMap.add("Manders", Manders);
        inputMap.add("# MC simulations", MC);
        inputMap.add("Overlap (# spots with overlap > T)", Overlap);
        inputMap.add("Overlap threshold (T)", threshold);
        inputMap.add("Distance based (centers of spots 2 within masks 1", distance);
        inputMap.add("SODA", soda);
        inputMap.add("Max. Radius (pixels)", rmax);
        inputMap.add("Step (pixels)", step);
    }

    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add("Workbook", book);
    }

    @Override
    public void run()
    {
        // initialisation du workbook
        Workbook wb = new HSSFWorkbook();
        wb.setMissingCellPolicy(Row.MissingCellPolicy.CREATE_NULL_AS_BLANK);
        book.setValue(wb);

        int row = 0;
        // create the sheet
        String sheetName = "Colocalization Results";
        Sheet sheet = wb.getSheet(sheetName);
        if (sheet == null)
            sheet = wb.createSheet(sheetName);

        Row header = sheet.createRow(0);
        header.getCell(0).setCellValue("t");
        header.getCell(1).setCellValue("PCC");
        header.getCell(2).setCellValue("p value (analytical)");
        header.getCell(3).setCellValue("log p value");
        header.getCell(4).setCellValue("M1");
        header.getCell(5).setCellValue("M2");
        header.getCell(6).setCellValue("# MC simulations");
        header.getCell(7).setCellValue("M1 (Mean of simulations)");
        header.getCell(8).setCellValue("M2 (Mean of simulations)");
        header.getCell(9).setCellValue("p value (M1>simulations");
        header.getCell(10).setCellValue("log p value (M1>simulations");
        header.getCell(11).setCellValue("p value (M2>simulations");
        header.getCell(12).setCellValue("log p value (M2>simulations");
        header.getCell(13).setCellValue("T");
        header.getCell(14).setCellValue("Overlap 1 (#>T)");
        header.getCell(15).setCellValue("Overlap 2 (#>T)");
        header.getCell(16).setCellValue("Overlap 1 (Mean of simulations)");
        header.getCell(17).setCellValue("Overlap 2 (Mean of simulations)");
        header.getCell(18).setCellValue("p value (Overlap 1>simulations");
        header.getCell(19).setCellValue("log p value (Overlap 1>simulations");
        header.getCell(20).setCellValue("p value (Overlap 2>simulations");
        header.getCell(21).setCellValue("log p value (Overlap 2>simulations");
        header.getCell(22).setCellValue("Distance-based");
        header.getCell(23).setCellValue("p value (Analytical)");
        header.getCell(24).setCellValue("log p value");
        header.getCell(25).setCellValue("SODA");
        header.getCell(26).setCellValue("p value (analytical)");
        header.getCell(27).setCellValue("log p value");
        header.getCell(28).setCellValue("Mean coupling distance");
        row++;

        Sequence seq = sequence1.getValue();
        int numT = seq.getSizeT();
        int dim = 2;
        if (seq.getSizeZ() > 1)
        {
            dim = 3;
        }
        // gestion des rois d'analyse en entr�e

        list_roi.clear();
        ROI[] roiArray = cell_mask.getValue();// TODO Ou sont les ROIs???
        if (roiArray.length == 0)
        {
            for (int t = 0; t < seq.getSizeT(); t++)
            {
                ROI roi = null;

                ROI2DRectangle r = new ROI2DRectangle(seq.getBounds2D());
                for (int h = 0; h < seq.getSizeZ(); h++)
                {
                    r.setZ(h);
                    r.setT(t);
                    try {
                        roi = r.getUnion(roi);
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
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
        // définition du ration pixel z/pixel x
        double ratio_zx = 1.0;
        if (seq.getSizeZ() > 1)
        {
            ratio_zx = seq.getPixelSizeZ() / seq.getPixelSizeX();
        }
        // init
        // on fait ensuite une boucle en temps
        for (int t = 0; t < numT; t++)
        {
            // computation of the cell mask at time t
            ArrayList<ROI> roi_t = SODAblock.ROI_t(dim, list_roi, t, -1);

            if (roi_t.size() == 0)
            {
            }
            else
            {
                // calcul du volume de la région d'étude
                double volume_t;
                try {
                    volume_t = ROIUtil.getUnion(roi_t).getNumberOfPoints() * ratio_zx;
                    if (ROIUtil.getUnion(roi_t).getBounds5D().isInfiniteZ())
                    {
                        volume_t = volume_t * seq.getSizeZ() * ratio_zx;
                    }
                }
                catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
                // Correlation computation
                double[] pearson = new double[3];
                double[] manders = new double[8];
                double[] overlap = new double[8];
                double[] dist = new double[3];
                double[] sod = new double[4];
                if (PCC.getValue())
                {
                    try {
                        pearson = Methods_Correlation.pearson_TCL(sequence1.getValue(), sequence2.getValue(), t, roi_t);
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
                // overlap
                // on selectionne les detection) ROI au temps t
                ArrayList<ROI> spots1_t = new ArrayList<ROI>();
                ArrayList<ROI> spots2_t = new ArrayList<ROI>();
                for (ROI r1 : spots1.getValue())
                {

                    if (r1.getBounds5D().isInfiniteZ())
                    {
                        if (dim == 2)
                        {
                            Point5D pos = r1.getPosition5D();
                            pos.setZ(0);
                            r1.setPosition5D(pos);
                        }
                    }
                    if (r1.getPosition5D().getT() == t)
                    {
                        spots1_t.add(r1);
                    }
                    if (r1.getBounds5D().isInfiniteT())
                    {
                        ROI r1bis = r1.getCopy();
                        Point5D pos = r1bis.getPosition5D();
                        pos.setT(t);
                        r1bis.setPosition5D(pos);
                        spots1_t.add(r1bis);
                    }
                }
                for (ROI r2 : spots2.getValue())
                {
                    if (r2.getBounds5D().isInfiniteZ())
                    {
                        if (dim == 2)
                        {
                            Point5D pos = r2.getPosition5D();
                            pos.setZ(0);
                            r2.setPosition5D(pos);
                        }
                    }
                    if (r2.getPosition5D().getT() == t)
                    {
                        spots2_t.add(r2);
                    }
                    if (r2.getBounds5D().isInfiniteT())
                    {
                        ROI r2bis = r2.getCopy();
                        Point5D pos = r2bis.getPosition5D();
                        pos.setT(t);
                        r2bis.setPosition5D(pos);
                        spots2_t.add(r2bis);
                    }
                }
                if (Manders.getValue())
                {
                    try {
                        manders = Methods_overlap.Manders(sequence1.getValue(), spots1_t, spots2_t, roi_t, MC.getValue(),
                                t);
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
                if (Overlap.getValue())
                {
                    try {
                        overlap = Methods_overlap.Overlap(sequence1.getValue(), spots1_t, spots2_t, roi_t, MC.getValue(), t,
                                threshold.getValue());
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
                // Object based
                // from array to vector..
                ROI[] spots_1_vec = new ROI[spots1_t.size()];
                for (int i = 0; i < spots1_t.size(); i++)
                {
                    spots_1_vec[i] = spots1_t.get(i);
                }
                ROI[] spots_2_vec = new ROI[spots2_t.size()];

                for (int i = 0; i < spots2_t.size(); i++)
                {
                    spots_2_vec[i] = spots2_t.get(i);
                }

                if (distance.getValue())
                {
                    try {
                        dist = Methods_distance.distance(spots_1_vec, spots_2_vec, sequence1.getValue(),
                                sequence2.getValue(), roi_t, t, volume_t);
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
                if (soda.getValue())
                {
                    try {
                        sod = Methods_distance.soda(spots_1_vec, spots_2_vec, sequence1.getValue(), sequence2.getValue(),
                                rmax.getValue(), step.getValue(), roi_t, t, volume_t);
                    }
                    catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
                
                Row ro = sheet.createRow(row);
                ro.getCell(0).setCellValue(t);
                ro.getCell(1).setCellValue(pearson[0]);
                ro.getCell(2).setCellValue(pearson[2]);
                ro.getCell(2).setCellValue(pearson[1]);
                ro.getCell(4).setCellValue(manders[0]);
                ro.getCell(5).setCellValue(manders[1]);
                ro.getCell(6).setCellValue(MC.getValue());
                ro.getCell(7).setCellValue(manders[2]);
                ro.getCell(8).setCellValue(manders[3]);
                ro.getCell(9).setCellValue(manders[4]);
                ro.getCell(10).setCellValue(manders[5]);
                ro.getCell(11).setCellValue(manders[6]);
                ro.getCell(12).setCellValue(manders[7]);
                ro.getCell(13).setCellValue(threshold.getValue());
                ro.getCell(14).setCellValue(overlap[0]);
                ro.getCell(15).setCellValue(overlap[1]);
                ro.getCell(16).setCellValue(overlap[2]);
                ro.getCell(17).setCellValue(overlap[3]);
                ro.getCell(18).setCellValue(overlap[4]);
                ro.getCell(19).setCellValue(overlap[5]);
                ro.getCell(20).setCellValue(overlap[6]);
                ro.getCell(21).setCellValue(overlap[7]);
                ro.getCell(22).setCellValue(dist[0]);
                ro.getCell(23).setCellValue(dist[2]);
                ro.getCell(24).setCellValue(dist[1]);
                ro.getCell(25).setCellValue(sod[0]);
                ro.getCell(26).setCellValue(sod[3]);
                ro.getCell(27).setCellValue(sod[1]);
                ro.getCell(28).setCellValue(sod[2]);
                row++;
            }
        }

    }

}
