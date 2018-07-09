/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import java.io.*;
import java.io.BufferedReader;
import java.io.FileReader;
import static ballista.utility.ConstParam.*;
import java.util.Optional;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonBar.ButtonData;
import javafx.scene.control.ButtonType;
import javafx.scene.image.Image;
import javafx.stage.Stage;

/**
 *
 * @author Kae Tsunematsu
 */
class TopoRead {
    private String FileNameDEMAscii;
    private double [][]data;
    private int headerlinenum;
    


    public TopoRead(String FileNameDEMAscii) {
        this.FileNameDEMAscii = FileNameDEMAscii;//"dem/hakone_part_5m_DEM.asc";//TODO: Read file name from parameter file
        headerlinenum = 0;
        
    }
    public void readCenter(String FileNameCenterPosition){
        try
        {
            //Read Geotiff header part
            FileReader fr2 = new FileReader(FileNameCenterPosition);
            BufferedReader br2 = new BufferedReader(fr2);
 
            String line;
           

            while((line = br2.readLine()) != null) {
                /*--- Reading header lines ---*/
                if(line.length() > 2){
                    String left = line.substring(0,12); String right = line.substring(13,line.length());
                    //System.out.println("left="+left+",  right="+right);
                    if(left.contains("CenterX")){
                        CenterX = Double.parseDouble(right);
                    }else if(left.contains("CenterY")){
                        CenterY = Double.parseDouble(right);
                    }else{
                        System.out.println("Catch the end of file");
                        break;
                    }
                }
            }
          
            br2.close();
            fr2.close();
//            
        }catch(Exception e){
            e.printStackTrace();
            System.err.println("Error of reading center position");
            System.exit(0);   
        }
        //this.data= new double[DATAHEIGHT][DATAWIDTH];
        
            //To show erros. KOKO
            if(DATAWIDTH == 0 || DATAHEIGHT == 0){
                ErrorDialogInReading("Data width or height is not correct.");
            }
            if(GridSize == 0){
                ErrorDialogInReading("Grid Size is not correct. \n"
                        + "Do you have the same grid size in x and y direction?");
            }
            if (XllCorner == 0 || YllCorner==0){
                ErrorDialogInReading("Corner coordinate is not correct. \n"
                        + "Check your DEM range.");        
            }         
            
        
    }
    
    public void readheader(){
        try
        {
            //Read Geotiff header part
            FileReader fr1 = new FileReader(FileNameDEMAscii);
            BufferedReader br1 = new BufferedReader(fr1);
 
            String line;
            int linenum = 0;
            

            while((line = br1.readLine()) != null) {
                
                /*--- Reading header lines ---*/
                String left = line.substring(0,12); String right = line.substring(13,line.length());
                //System.out.println("left="+left+",  right="+right);
                if(left.contains("ncols")){
                    DATAWIDTH = Integer.valueOf(right);
                    System.out.println("DATAWIDTH = "+DATAWIDTH);
                    linenum++;
                }else if(left.contains("nrows")){
                    DATAHEIGHT = Integer.valueOf(right);
                    System.out.println("DATAHEIGHT = "+DATAHEIGHT);
                    linenum++;
                }else if(left.contains("xllcorner")){
                    XllCorner = Double.valueOf(right);
                    System.out.println("XllCorner = "+XllCorner);
                    linenum++;                        
                }else if(left.contains("yllcorner")){
                    YllCorner = Double.valueOf(right);
                    System.out.println("YllCorner = "+YllCorner);
                    linenum++;
                }else if(left.contains("cellsize")){
                    GridSize = Double.valueOf(right);
                    System.out.println("GridSize = "+GridSize);
                    linenum++;
                }else if(left.contains("NODATA_value")){
                    NODATA_VALUE = Integer.valueOf("-9999");
                    System.out.println("NODATA_VALUE = "+NODATA_VALUE);
                    linenum++;
                }else{
                    break;
                }
            }
            headerlinenum = linenum;

            br1.close();
            fr1.close();
//            
        }catch(Exception e){
            e.printStackTrace();
            System.err.println("Error of reading header");
            System.exit(0);   
        }
       
        //Error Dialog when the reading of the header results in some problems.
        
        if(DATAWIDTH ==0 || DATAHEIGHT == 0)
            ErrorDialogInReading("Data width or height is not correct. ");
        if(GridSize == 0){
            ErrorDialogInReading("GridSize is not correct. \n"
                    + "Check your DEM whether it has the same grid size in X and Y direction. ");
        }else{
            NZ = (int)Math.round(calcAZ/GridSize);
            System.out.println("-----------Here is the setting of the NX, NY, NZ ------------------");
            System.out.println("(NX, NY, NZ) ="+DATAWIDTH+"\t"+DATAHEIGHT+"\t"+NZ);
        }
        if(XllCorner==0 || YllCorner == 0)
            ErrorDialogInReading("Lower left coordinate is not correct. ");
        
                 
            
        
    }
    
    
    public double[][] altitude()
    {
//        System.out.println("DATAHEIGHT = "+DATAHEIGHT);
        data = new double[DATAHEIGHT][DATAWIDTH];

        String p = "";
        String line;
        int gridlinenum = 0;
        int linenum = 0;
        try{
            
            //Read Geotiff header part
            FileReader fr1 = new FileReader(FileNameDEMAscii);
            BufferedReader br1 = new BufferedReader(fr1);
//            System.out.println("file is here!"+DATAHEIGHT);
                    //System.out.println("DATAWIDTH="+DATAWIDTH+", DATAHEIGHT="+ DATAHEIGHT);
            while((line = br1.readLine()) != null) { 
                //System.out.println("linenum = "+line);
                if(linenum >= headerlinenum){//If it is after the header lines
                    
                    String [] readtext;
                    readtext = line.split(" ");
//                    System.out.println("num of texts= "+readtext.length);
                    for(int i=1;i<readtext.length; i++){ //Be careful! Ascii file includes one space at the head of line
                        //System.out.print("readtext ="+readtext[i]);
                        //data[gridlinenum][i-1] = Double.parseDouble(readtext[i]);
                        data[DATAHEIGHT-gridlinenum-1][i-1] = Double.parseDouble(readtext[i]);
//                        System.out.println("i= "+i+" data ="+data[gridlinenum][i-1]);
                    }

//                    System.out.println("linenum = "+gridlinenum);
                    gridlinenum++;
                    
                }else{
//                    System.out.println("koko");
                }
                    
                linenum++;
            }
         
            br1.close();
            fr1.close();
 
        }catch(Exception e){
            
        }
        
        return data;
    }
    
    public double getcenterXYaltitude(){
        int nX =  (int)(Math.floor((CenterX- XllCorner)/GridSize));
        int nY =  (int)(Math.floor((CenterY- YllCorner)/GridSize));
        if(nX> DATAWIDTH){
            System.err.println("Center X is out of DEM region +");
            String errorMS = "Center X is out of DEM region +";
            ErrorDialogCenterPosition(errorMS);
        }else if(nX <0){
            System.err.println("Center X is out of DEM region -");
            String errorMS = "Center X is out of DEM region -";
            ErrorDialogCenterPosition(errorMS);
        }
        if(nY> DATAHEIGHT){
            System.err.println("Center Y is out of DEM region +");
            String errorMS = "Center Y is out of DEM region +";
            ErrorDialogCenterPosition(errorMS);
        }else if(nY <0){
            System.err.println("Center Y is out of DEM region -");
            String errorMS = "Center Y is out of DEM region -";
            ErrorDialogCenterPosition(errorMS);
        }
        
        
        System.out.println("GridSize ="+GridSize);
        System.out.println("at Center("+CenterX+", "+CenterY+")");
        System.out.println ("nX = "+nX+", nY ="+nY);
        System.out.println ( ":  altitude is "+data[nY][nX]);
        return data[nY][nX];
    }
        

    public int getDATAWIDTH() {
        return DATAWIDTH;
    }

    public int getDATAHEIGHT() {
        return DATAHEIGHT;
    }

    public double getXllCorner() {
        return XllCorner;
    }

    public double getYllCorner() {
        return YllCorner;
    }

    public double getGridSize() {
        return GridSize;
    }

    public int getNODATA_VALUE() {
        return NODATA_VALUE;
    }
    
    public void ErrorDialogInReading(String ErrorMessage){
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("Error Dialog");
        alert.setHeaderText("Error in reading a DEM file");
        alert.setContentText(ErrorMessage);
        ButtonType buttonTypeCancel = new ButtonType("Ok and Exit System", ButtonData.CANCEL_CLOSE);
        
        alert.getButtonTypes().setAll( buttonTypeCancel);
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        Optional<ButtonType> result = alert.showAndWait();
        if (result.get() == buttonTypeCancel){
            System.exit(-1);
        }
    }
    
    public void ErrorDialogCenterPosition(String errorMS){
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("Error Dialog");
        alert.setHeaderText(errorMS);
        alert.setContentText("Please set the vent center in the DEM region.");
        ButtonType buttonTypeCancel = new ButtonType("Ok and Exit System", ButtonData.CANCEL_CLOSE);
        
        alert.getButtonTypes().setAll( buttonTypeCancel);
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        Optional<ButtonType> result = alert.showAndWait();
        if (result.get() == buttonTypeCancel){
            System.exit(-1);
        }
    }
   
    
}
