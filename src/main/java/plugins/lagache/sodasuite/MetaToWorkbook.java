package plugins.lagache.sodasuite;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import icy.main.Icy;
import icy.plugin.PluginLauncher;
import icy.plugin.PluginLoader;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzVarText;
import plugins.adufour.vars.lang.VarDoubleArrayNative;
import plugins.adufour.vars.lang.VarInteger;
import plugins.adufour.vars.lang.VarWorkbook;

public class MetaToWorkbook extends Plugin implements Block, PluginBundled
{
    private EzVarText name;
    VarDoubleArrayNative positionX = new VarDoubleArrayNative("positionX", null);
    VarDoubleArrayNative positionY = new VarDoubleArrayNative("positionY", null);
    VarDoubleArrayNative pixelSizeX = new VarDoubleArrayNative("pixelSizeX", null);
    VarDoubleArrayNative pixelSizeY = new VarDoubleArrayNative("pixelSizeY", null);
    VarInteger sizeX = new VarInteger("sizeX", null);
    VarInteger sizeY = new VarInteger("sizeY", null);
    VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

    public static void main(String[] args)
    {
        Icy.main(args);
        PluginLauncher.start(PluginLoader.getPlugin(MetaToWorkbook.class.getName()));
    }

    @Override
    public void declareInput(VarList inputMap)
    {
    	name = new EzVarText("filename", "", 0);
    	inputMap.add("name", name.getVariable());
        inputMap.add("positionX", positionX);
        inputMap.add("positionY", positionY);
        inputMap.add("pixelSizeX", pixelSizeX);
        inputMap.add("pixelSizeY", pixelSizeY);
        inputMap.add("sizeX", sizeX);
        inputMap.add("sizeY", sizeY);
    }

    @Override
    public void declareOutput(VarList outputMap)
    {
        outputMap.add("Workbook", book);
    }

    @Override
    public String getMainPluginClassName()
    {
        return MetaToWorkbook.class.getName();
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
        Sheet sheet = wb.getSheet(sheetName);
        if (sheet == null)
            sheet = wb.createSheet(sheetName);
        Row header = sheet.getRow(0);
        if (header == null)
        {
            header = sheet.createRow(0);
            header.getCell(0).setCellValue("filename");
            header.getCell(1).setCellValue("positionX");
            header.getCell(2).setCellValue("positionY");
            header.getCell(3).setCellValue("pixelSizeX");
            header.getCell(4).setCellValue("pixelSizeY");
            header.getCell(5).setCellValue("sizeX");
            header.getCell(6).setCellValue("sizeY");
        }
        int i = sheet.getPhysicalNumberOfRows();
        Row row = sheet.createRow(i);
        row.getCell(0).setCellValue(name.getValue());
        row.getCell(1).setCellValue(positionX.getValue()[0]);
        row.getCell(2).setCellValue(positionY.getValue()[0]);
        row.getCell(3).setCellValue(pixelSizeX.getValue()[0]);
        row.getCell(4).setCellValue(pixelSizeY.getValue()[0]);
        row.getCell(5).setCellValue(sizeX.getValue());
        row.getCell(6).setCellValue(sizeY.getValue());
        // on crée une nouvelle ligne dans la feuille excel
    }
}
