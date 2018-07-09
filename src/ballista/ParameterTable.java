package ballista;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author volcano
 */

import java.io.File;
//import javafx.application.Application;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableColumn.CellEditEvent;
import javafx.scene.control.TableView;
//import javafx.scene.control.TextField;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
//import javafx.scene.layout.VBox;
import javafx.scene.text.Font;
import javafx.stage.Stage;
import static ballista.utility.ConstParam.*;
import static ballista.utility.FileNameStorage.InputFName;
import javafx.scene.layout.GridPane;
 
public class ParameterTable  {

    public ParameterTable() {
       ParamRead pr = new ParamRead();
       pr.readandSetInput();
    }
     
    private TableView<Inputs> table = new TableView<Inputs>();
    private final ObservableList<Inputs> data =
            FXCollections.observableArrayList(
            new Inputs("directory.name", "OntakeDEM", "-"),
            new Inputs("CenterX", "000000", "-"),    
            new Inputs("CenterY", "000000", "-"),
            new Inputs("numOfParticles", "100", "-"),
            new Inputs("density.avg", "2300", "kg/m^3"),
            new Inputs("density.sd", "300", "kg/m^3"),
            new Inputs("diameter.avg", "0.2", "m"),
            new Inputs("diameter.sd ", "0.0", "m"),
            new Inputs("diameter.min ", "0.1", "m"),
            new Inputs("diameter.max ", "0.7", "m"),
            new Inputs("displacement.sd", "5.0", "m"),
            new Inputs("displacement.max", "10.0", "m"),
            new Inputs("displacement.min", "10.0", "m"),
            new Inputs("velocity.norm.avg", "100.0", "m/s"),
            new Inputs("velocity.norm.sd", "0.0", "m/s"),
            new Inputs("axis.eject.deg", "40.0", "degree"),
            new Inputs("direct.bearing.deg", "20.0", "degree"),
            new Inputs("wind.velocity.x", "0.0", "m/s"),
            new Inputs("wind.velocity.y", "0.0", "m/s"),
            new Inputs("dragCoefficient.cnst", "0.8", "-")
            );
    final HBox hb = new HBox();
   

    public void start(Stage stage) {
        Scene scene = new Scene(new Group());
        stage.setTitle("Ballista Input Parameters");
        stage.setWidth(650);
        stage.setHeight(600);
        File inputFile = new File(InputFName);
        String labelStr = "INIT\\ "+inputFile.getName();
        final Label label = new Label(labelStr);
        label.setFont(new Font("Arial", 20));
 
        table.setEditable(false);
        table.setMinHeight(420);
 
        TableColumn paramNameCol = new TableColumn("Parameter Name");
        paramNameCol.setMinWidth(250);
        paramNameCol.setCellValueFactory(
            new PropertyValueFactory<Inputs, String>("paramName"));
        paramNameCol.setCellFactory(TextFieldTableCell.forTableColumn());
        paramNameCol.setOnEditCommit(
            new EventHandler<CellEditEvent<Inputs, String>>() {
                @Override
                public void handle(CellEditEvent<Inputs, String> t) {
//                    ((Inputs) t.getTableView().getItems().get(
//                            t.getTablePosition().getRow())
//                            ).setParamName(t.getNewValue());
                }
            }
        );
 
 
        TableColumn valueCol = new TableColumn("Value");
        valueCol.setMinWidth(200);
        valueCol.setCellValueFactory(
            new PropertyValueFactory<Inputs, String>("values"));
        valueCol.setCellFactory(TextFieldTableCell.forTableColumn());
        valueCol.setOnEditCommit(
            new EventHandler<CellEditEvent<Inputs, String>>() {
                @Override
                public void handle(CellEditEvent<Inputs, String> t) {
//                    ((Inputs) t.getTableView().getItems().get(
//                        t.getTablePosition().getRow())
//                        ).setValues(t.getNewValue());
                }
            }
        );
 
        TableColumn unitCol = new TableColumn("Unit");
        unitCol.setMinWidth(150);
        unitCol.setCellValueFactory(
            new PropertyValueFactory<Inputs, String>("unit"));
        unitCol.setCellFactory(TextFieldTableCell.forTableColumn());
        unitCol.setOnEditCommit(
            new EventHandler<CellEditEvent<Inputs, String>>() {
                @Override
                public void handle(CellEditEvent<Inputs, String> t) {
//                    ((Inputs) t.getTableView().getItems().get(
//                        t.getTablePosition().getRow())
//                        ).setUnit(t.getNewValue());
                }
            }
        );
 
        table.setItems(data);
        table.getColumns().addAll(paramNameCol, valueCol, unitCol);
        data.clear();
        data.add(new Inputs("directory.name", dirDEMName, "-"));
        data.add(new Inputs("CenterX", Double.toString(CenterX), "-"));
        data.add(new Inputs("CenterY", Double.toString(CenterY), "-"));
        data.add(new Inputs("numOfParticles", Double.toString(numP), "-"));
        data.add(new Inputs("density.avg", Double.toString(avgDensity), "kg/m^3"));
        data.add(new Inputs("density.sd", Double.toString(sdDensity), "kg/m^3"));
        data.add(new Inputs("diameter.avg", Double.toString(avgDiam), "m"));
        data.add(new Inputs("diameter.sd ", Double.toString(sdDiam), "m"));
        data.add(new Inputs("diameter.min ", Double.toString(minDiam), "m"));
        data.add(new Inputs("diameter.max ", Double.toString(maxDiam), "m"));
        data.add(new Inputs("displacement.sd", Double.toString(sdDisp), "m"));
        data.add(new Inputs("displacement.max", Double.toString(maxDisp), "m"));
        data.add(new Inputs("velocity.norm.avg", Double.toString(avgV), "m/s"));
        data.add(new Inputs("velocity.norm.sd", Double.toString(sdV), "m/s"));
        data.add(new Inputs("axis.eject.deg", Double.toString(ejcDeg), "degree"));
        data.add(new Inputs("direct.bearing.deg", Double.toString(90.0-direcBear), "degree"));
        data.add(new Inputs("wind.velocity.x", Double.toString(windX),"m/s"));
        data.add(new Inputs("wind.velocity.y", Double.toString(windY),"m/s"));
        data.add(new Inputs("dragCoefficient.cnst", Double.toString(Cd), "-"));
        
        

        final Button gobackButton = new Button("Exit and Edit Parameter File");
        gobackButton.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
               stage.close();
               
            }
        
        });
 
        final Button addButton = new Button("Calculate");
        addButton.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent e) {
                SimulatorMain simMain = new SimulatorMain();
                
            }
        });
        hb.setAlignment(Pos.BOTTOM_CENTER);
        hb.getChildren().addAll( gobackButton, addButton);
        hb.setSpacing(100);
 

        
        GridPane gpane = new GridPane();
      
        gpane.setPadding(new Insets(10, 10, 10, 10));
        gpane.setHgap(5);
        gpane.setVgap(5);
        gpane.add(label, 0, 0, 3, 1);
        gpane.add(table, 0, 1, 3, 5);
        gpane.add(hb,    0, 6, 3, 1);
        
        ((Group) scene.getRoot()).getChildren().addAll(gpane);
 
        stage.setScene(scene);
        stage.show();
    }
    public class Inputs {
 
        private final SimpleStringProperty paramName;
        private final SimpleStringProperty values;
        private final SimpleStringProperty unit;
 
        private Inputs(String pName, String lName, String unitStr) {
            this.paramName = new SimpleStringProperty(pName);
            this.values = new SimpleStringProperty(lName);
            this.unit = new SimpleStringProperty(unitStr);//Please don't update in the program

        }
 
        public String getParamName() {
            return paramName.get();
        }
 
        public void setParamName(String pName){
            paramName.set(pName);
        }
 
        public String getValues() {
            return values.get();
        }
 
        public void setValues(String value) {
            values.set(value);
        }
 
        public String getUnit() {
            return unit.get();
        }
 
        public void setUnit(String unitStr) {
            unit.set(unitStr);
        }
    }
    
}