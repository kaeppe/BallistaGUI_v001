/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import static ballista.utility.ConstParam.*;
import static ballista.utility.FileNameStorage.InputFName;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Properties;
import javafx.scene.control.Alert;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.image.Image;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.stage.Stage;


/**
 *
 * @author kts42
 */
public class ParamRead {

    public ParamRead() {
        
    }
   
    public void readandSetInput(){
        System.out.println("Here, Reading method is called");
        Properties prop = new Properties();
        try{
            prop.load(new FileReader(InputFName));
        }catch(IOException e){
            e.printStackTrace();
        }
        
        ReadInit params = new ReadInit(prop);
        try{
        dirDEMName = params.getString("directory.name");
        CenterX = params.getDouble("CenterX");
        CenterY = params.getDouble("CenterY");
        numP = params.getDouble("numOfParticles");
        avgDensity = params.getDouble("density.avg");
        sdDensity = params.getDouble("density.sd");
        avgDiam = params.getDouble("diameter.avg");
        sdDiam = params.getDouble("diameter.sd");
        minDiam = params.getDouble("diameter.min");
        maxDiam = params.getDouble("diameter.max");
        sdDisp = params.getDouble("displacement.sd");
        maxDisp = params.getDouble("displacement.max");
        avgV = params.getDouble("velocity.norm.avg");
        sdV = params.getDouble("velocity.norm.sd");
        ejcDeg = params.getDouble("axis.eject.deg");
        direcBear = 90.0 - params.getDouble("direct.bearing.deg");
        windX = params.getDouble("wind.x");
        windY = params.getDouble("wind.y");
        Cd = params.getDouble("dragCoefficient.cnst");

        }catch(Exception ex){
            showFileException(ex);
            System.err.println("Error of reading input file");
            System.exit(0);
        }
        
        System.out.println("dirDEMName = "+dirDEMName);
    }
    
    /**
     * 
     * @param ex 
     */
    public void showFileException(Exception ex){
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("File or Directory");
        alert.setHeaderText("Can't find the direcotry...");
        alert.setContentText("Something wrong with the directory structures.\n "
                + "Did you make the result directory at the same level with INIT directory?");
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

        // Set expandable Exception into the dialog pane.
        alert.getDialogPane().setExpandableContent(expContent);

        alert.showAndWait();
        
    }
    
    
}
