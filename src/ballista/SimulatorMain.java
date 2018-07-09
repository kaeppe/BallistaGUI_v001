/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import java.io.*;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringWriter;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Properties;
import java.util.Random;
import static java.lang.Math.*;

import static ballista.utility.ConstParam.*;
import static ballista.utility.FileNameStorage.*;
import static ballista.utility.Output.*;

import javafx.scene.control.Alert;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.image.Image;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.stage.Stage;



/**
 *
 * @author volcano
 */
public class SimulatorMain {
     Alert calculatingDialog;
    
    
    public  SimulatorMain(){
        File inputfile = new File(InputFName);
        /*Read Properties from initial condition file*/
        Properties prop = new Properties();
        try{
            prop.load(new FileReader(InputFName));
        }catch(IOException e){
            e.printStackTrace();
        }
        ReadInit ri = new ReadInit(prop);
                
        /*Find DEM filename*/
        dirDEMName = ri.getString("directory.name");
        System.out.println("DEM directory = "+ dirDEMName);
        String currentDir = System.getProperty("user.dir");
        System.out.println("My current directory is = "+currentDir);
        //TODO: get the path from inputfile place
        //Get Path of Topo files
        String dirDEMPath = inputfile.getParent().replace("INIT", "DEM")+"//"+dirDEMName;
        System.out.println("DEM directory path is "+dirDEMPath);
        File demDir = new File(dirDEMPath);
        if(!demDir.isDirectory()){
            Exception ex = new FileNotFoundException("Could not find a DEM directory ");
            showFileException(ex);
            System.exit(0);
        }
        String[] listFiles = demDir.list();
        String esriAsciiGridFileName = "";
//        String centerPositionFileName = "";
        for(int i =0; i<listFiles.length;i++){
//            System.out.println("list: "+listFiles[i]);
            if(listFiles[i].endsWith(".asc")){
                esriAsciiGridFileName = listFiles[i];
            }
        }
        //Ascii file
        System.out.println("Esri Ascii grid:      "+esriAsciiGridFileName);
       
        
        
        /*--Read Topography---*/
        String demFilePath = dirDEMPath+"//"+esriAsciiGridFileName;
        TopoRead tprd = new TopoRead(demFilePath); 
        tprd.readheader();
        /*- Set altitude data to the memory-*/
        Altitude = tprd.altitude();
        
        /*--Read Center Position from a file--*/
        /*-- read the height of Center position and set to the memory --*/
        CenterZ = tprd.getcenterXYaltitude();
        System.out.println("Center coordinates = "+ CenterX+ ","+CenterY+","+CenterZ+")");
               
        /*--Bursting!--*/ 
        double avgFq = 1.0;
        double sdFq  = 0;
        double maxTime = 1.001;
        
        Random rnd = new Random();
        double time = 0;
        //System.out.println("maxTime ="+maxTime);
        double nb = maxTime/avgFq -1 ; //Number of Burst
        List timeary = new ArrayList();

        if(avgFq == maxTime){
            timeary.add(0.0);
        }else{
            while(time < maxTime){
                time = time + abs(rnd.nextGaussian() * sdFq + avgFq);
                if(time<maxTime){
                    timeary.add(time);
                }
            }
        }
        
        /*--Sort array of tasks with time--*/
        Collections.sort(timeary);
        
        /*-- Output directory path is defined here --*/
        String nameoffile = inputfile.getName();
        outFilePath = inputfile.getParent().replace("INIT", "results//")+nameoffile.replace(".txt", "//");//TODO: calibrate once the code is executable only with execute file
        Path path = Paths.get(outFilePath);
        //if directory exists? if not, make directory
        if (!Files.exists(path)) {
            try {
                Files.createDirectories(path);
            } catch (IOException e) {
                //fail to create directory
                showFileException(e);
                System.exit(0);
            }
        }
        
        //Let's start Simulation
        showNowCalculating();
        Simulation sim = new Simulation(timeary, ri);
        String inputfilename = inputfile.getName();
        sim.iterate(inputfilename);
      
            
//        String filefordeposition = outFilePath+"resdepo.txt";        
//        String filefortrajectory = outFilePath+"restraje.txt";
//        writeCoord1D(filefordeposition, sim.depoparticles);
//        writeTrajectory3D(filefortrajectory, sim.trajectories); //Trun on this line and other trajectory gettin line in Simulation.java, if you want trajectory
        calculatingDialog.show();
        showFinishedInfo();
        System.exit(0);
    }
    
    public void showNowCalculating(){
        calculatingDialog = new Alert(Alert.AlertType.NONE);
	calculatingDialog.setTitle("Ballista ");
	calculatingDialog.setHeaderText(" Ballista is calculating now.");
	String s ="Please wait until the simulation finishes.";
        calculatingDialog.setContentText(s);
        Stage stage = (Stage) calculatingDialog.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        calculatingDialog.show();
       
    }
    
    public void showFinishedInfo(){
        Alert alert = new Alert(AlertType.INFORMATION);
	alert.setTitle("Simulation Finished");
	alert.setHeaderText("Simulation is successfully finished");
	String s ="Please find an output file at:\n"+outFilePath;
	alert.setContentText(s);
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        alert.showAndWait()
            .filter(response -> response == ButtonType.OK)
            .ifPresent(response -> finisheSystem());
    }
    
    public void showFileException(Exception ex){
        Alert alert = new Alert(AlertType.ERROR);
        alert.setTitle("File or Directory");
        alert.setHeaderText("Can't find the direcotry...");
        alert.setContentText("Something wrong with the directory structures.\n "
                + "Did you put the result & DEM directory at the same level with INIT directory?");
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        ex.printStackTrace(pw);
        String exceptionText = sw.toString();

        Label label = new Label("The exception stacktrace was:");

        TextArea textArea = new TextArea(exceptionText);
        textArea.setEditable(false);
        textArea.setWrapText(true);

        textArea.setMaxWidth(Double.MAX_VALUE);
        textArea.setMaxHeight(Double.MAX_VALUE);
        GridPane.setVgrow(textArea, Priority.ALWAYS);
        GridPane.setHgrow(textArea, Priority.ALWAYS);

        GridPane expContent = new GridPane();
        expContent.setMaxWidth(Double.MAX_VALUE);
        expContent.add(label, 0, 0);
        expContent.add(textArea, 0, 1);


        alert.getDialogPane().setExpandableContent(expContent);

        alert.showAndWait();
        
    }

    private void finisheSystem() {
       System.exit(0);
    }
}
