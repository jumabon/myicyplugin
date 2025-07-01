package plugins.lagache.sodasuite;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.VarDoubleArrayNative;
import plugins.adufour.vars.lang.VarWorkbook;

public class MetaToWorkbook extends Plugin implements Block, PluginBundled
{
    VarStringNative positionX = new VarDoubleArrayNative("positionX", null); 
    VarDoubleArrayNative positionX = new VarDoubleArrayNative("positionX", null);
    VarDoubleArrayNative positionY = new VarDoubleArrayNative("positionY", null);
    VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

    @Override
    public void declareInput(VarList inputMap)
    {

        inputMap.add("positionX", positionX);
        inputMap.add("positionY", positionY);
    }

    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add("Workbook", book);
    }

    @Override
    public String getMainPluginClassName()
    {
        return SODAsuite.class.getName();
    }

    @Override
    public void run()
    {
        // initialisation du workbook
        Workbook wb = new HSSFWorkbook();
        wb.setMissingCellPolicy(Row.MissingCellPolicy.CREATE_NULL_AS_BLANK);
        book.setValue(wb);

        // create the sheet
        String sheetName = "Results";
        String sheetName2 = new String("Results 2");
        Sheet sheet = wb.getSheet(sheetName);
        if (sheet == null)
            sheet = wb.createSheet(sheetName);
        Row header = sheet.getRow(0);
        if (header == null)
        {
            header = sheet.createRow(0);
            header.getCell(0).setCellValue("positionX");
            header.getCell(1).setCellValue("positionY");
        }

        int ind_row = 0;
        ind_row++;
        Row row = sheet.createRow(ind_row);
        // on recupere les tableaux associés aux proba et positionY
        double[] proba = positionX.getValue();
        double[] dist = positionY.getValue();
        int size_max = Math.min(proba.length, 65000);
        for (int j = 0; j < size_max; j++)
        {
            row.getCell(0).setCellValue(proba[j]);
            row.getCell(1).setCellValue(dist[j]);
            // on crée une nouvelle ligne dans la feuille excel
            ind_row++;
            row = sheet.createRow(ind_row);
        }
    }
}
