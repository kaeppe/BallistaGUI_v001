/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import java.io.File;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.layout.GridPane;
import javafx.scene.text.Font;
import javafx.scene.text.TextAlignment;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import static ballista.utility.FileNameStorage.*;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;


/**
 *
 * @author yies-133
 */
public class Main extends Application {
    
    @Override
    public void start(final Stage stage) {
        GridPane grid = new GridPane();
        grid.setAlignment(Pos.CENTER);
        grid.setHgap(10);
        grid.setVgap(10);
        grid.setPadding(new Insets(25, 25, 25, 25));
        
        Label labelTitle = new Label("Ballista");
        labelTitle.setFont(new Font("Century", 32));
        labelTitle.setTextAlignment(TextAlignment.JUSTIFY);
        grid.add(labelTitle, 1, 0); 
        
        Label labelSubTitle = new Label("-3D Ballistic Simulator-");
        labelSubTitle.setFont(new Font("Century", 18));
        labelSubTitle.setTextAlignment(TextAlignment.CENTER);
        grid.add(labelSubTitle, 1,1);

                    
        Label label = new Label("Open an input parameter file");
        label.setFont(new Font("Century",14));
        grid.add(label,1,2);
        

        final FileChooser fileChooser = new FileChooser();
        Button openbtn = new Button();
        openbtn.setText("Open File");
        openbtn.setOnAction(
            new EventHandler<ActionEvent>() {
                @Override
                public void handle(final ActionEvent e) {
                    File file = fileChooser.showOpenDialog(stage);
                    String fpath = file.getAbsolutePath();
                    System.out.println("fpath = "+ fpath);
                    if(file.isFile()&& file.getName().endsWith(".txt")){
                        InputFName =file.getAbsolutePath();
                        ParamRead pr = new ParamRead();
                        ParameterTable tableview = new ParameterTable();
                        pr.readandSetInput();
                        tableview.start(stage);
                    }else{
                        showAlert();

                    }
                }
            }
        );
        
        openbtn.setMaxWidth(Double.MAX_VALUE);
        grid.add(openbtn, 2,2);
        
        Image trajeImage = null;
        try{
            trajeImage = new Image("file:images/LogoBallista.png");
            
        }catch(Exception e){
            System.err.println("image open exception..."+ e);
        }

        ImageView imview = new ImageView();
        imview.setImage(trajeImage);
        grid.add(imview, 1, 3, 2, 1);
        
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        
        Scene scene = new Scene(grid, 500, 600);
        
        stage.setTitle("Ballista");
        stage.setScene(scene);
        stage.setOnCloseRequest(e -> Platform.exit());
        stage.show();
    }
    
    /**
     * 
     */
    public void showAlert(){
        Alert alert = new Alert(AlertType.ERROR);
	alert.setTitle("File Selection Error");
	alert.setHeaderText("File Selection is faile");
	String s ="Please find and select the input file with the extention .txt";
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
	alert.setContentText(s);
	alert.show();
    }
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        launch(args);
    }
    

    
}
