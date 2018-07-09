package ballista;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.util.ArrayList;
import java.util.List;
import static java.lang.Math.*;

import static ballista.utility.ConstParam.*;
import static ballista.utility.FileNameStorage.outFilePath;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonBar;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Label;
import javafx.scene.control.TextArea;
import javafx.scene.image.Image;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.stage.Stage;




/**
 *
 * @author tsunemat
 */
public class Simulation {
    private List<Particle> particles;
    private List<Particle> bparticles;
    private List btimeary;
    private final double maxTime;
    private final ReadInit ri;
    private Velocity velocity;
    public List[][] depositParticles;
    public List<Particle> depoparticles;
    private double ctime;//current time
    public List<Trajectory> trajectories;
    private PrintWriter outTraje;//Trajectory file writer
    private PrintWriter outDepo;//Deposit file writer
    

    public Simulation( List timeary, ReadInit ri) {
        this.btimeary = timeary;
        this.maxTime = (double) timeary.get((timeary.size()-1)) + Timebuff;
        this.ri = ri;
        particles=makeParticleList();
        bparticles = makeParticleList();
        depoparticles = makeParticleList();
        velocity = new Velocity();
        this.depositParticles = makeParticleField2D();
        trajectories = makeTrajectoryList();
    
    }
    
    public void iterate(String inputfilename) {
        //Initialize PrintWrite and StringWriter for outputting Trajectory and Deposit

        try {
            String trajeFPath = outFilePath+inputfilename.replace(".txt", "_traje.txt");
            System.out.println("trajeFPath ="+ trajeFPath);
            outTraje = name2printWriter(trajeFPath);
            outTraje.println("Time\tId\tX\tY\tZ");
        } catch (IOException ex) {
            Logger.getLogger(Simulation.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try {
            String depoFPath = outFilePath+inputfilename.replace(".txt", "_depo.txt");
            System.out.println("depoFPath ="+ depoFPath);
            outDepo = name2printWriter(depoFPath);
            outDepo.println("X\tY\tZ\tDistance\tMass\tVelocity\tDenisty\tDiameter\tId");
        } catch (IOException ex) {
            Logger.getLogger(Simulation.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        
        
        double timeburst = 0.;
        while(ctime <maxTime){//Time loop
    
            /*---transport all existing particles---*/
            transport(ctime, particles);
            
            /*---check is there any bursting or not?---*/
            /*---This loop continues until the btime is larger than ctime---*/
            while((double)(btimeary.size())>0.0){
                if((ctime+Dt) < (double)btimeary.get(0)){
                    break;
                }
                 timeburst = (double)btimeary.get(0);
                 btimeary.remove(0);
                 double dtb = ctime+Dt - timeburst;
                 //TODO: delete unapplied (supposed to...) bursting
                 //System.out.println("current time="+ctime+", timeburst="+timeburst+", dtb="+dtb);
                 burst(timeburst);

                 particles.addAll(bparticles); 
                 //String filename = "../../burst.txt";
                 //writeCoord1D(filename,bparticles);
                 
                 bparticles.clear();
            }
            
            ctime = ctime+Dt;
            
            //If particle has all removed and already bursting has all finished,
            //Simulation ends
            if(particles.size()<=0 && btimeary.size()<=0) 
                break;
            
        }//End of time loop 
        outTraje.close();
        outDepo.close();
        
        
    }
    
    public void burst(double btime){
        /*Particle draw*/
        GaussianBurst gb = new GaussianBurst(btime);
        bparticles = gb.drawPartList();
        
        for(int i=0;i<bparticles.size();i++){
            Particle bp = (Particle) bparticles.get(i);
        }
    }
    
    /**
     * 
     * @param t: Current Time
     * @param parray : Array of Particles
     */
    public void transport(double t, List<Particle> parray){
    
        for(int ip =0; ip < parray.size(); ip++){
 
            Particle currentP = (Particle) parray.get(ip);
            double px = currentP.getPx();
            double py = currentP.getPy();
            double pz = currentP.getPz();
            
            //Initially, only the velocity at t=1/2 (delta t) is calculated.
            if (t == 0){
                double dz = pz - CenterZ;
                double dist = sqrt(px *px + py*py + dz*dz);
                double[] vel = velocity.getVelocityWind(currentP, (0.5*Dt), dist);
                currentP.setVelocity(vel[0], vel[1], vel[2]);
            }
                       
            int ix = currentP.getIndexX();
            int iy = currentP.getIndexY();
            int iz = currentP.getIndexZ();
            
            double vx = currentP.getVelocity()[0];
            double vy = currentP.getVelocity()[1];
            double vz = currentP.getVelocity()[2];
            
            double storageX = currentP.getStorageX();
            double storageY = currentP.getStorageY();
            double storageZ = currentP.getStorageZ();
            
                                    
            double nextStorageX = storageX + vx*Dt;
            double nextStorageY = storageY + vy*Dt;
            double nextStorageZ = storageZ + vz*Dt;
            
            /*If there is no movement*/
            int nextX = ix;
            int nextY = iy;
            int nextZ = iz;
            
            /*Let's see the movement of particles by the amount of storage */
            //X
            if(abs(nextStorageX)>= GridSize && nextStorageX >= 0){// move to + direction
                nextX = ix + 1;
                nextStorageX = nextStorageX - GridSize;
                storageX = 0.0;
            }else if(abs(nextStorageX)> GridSize && nextStorageX < 0){// move to - direction
                nextX = ix - 1;
                nextStorageX = nextStorageX + GridSize;
                storageX = 0.0;
            }
            //Y
            if(abs(nextStorageY)>= GridSize && nextStorageY >= 0){
                nextY = iy + 1;
                nextStorageY = nextStorageY - GridSize;
                storageY = 0.0;
            }else if(abs(nextStorageY)> GridSize && nextStorageY < 0){
                nextY = iy - 1;
                nextStorageY = nextStorageY + GridSize;
                storageY = 0.0;
            }
            //Z
            if(abs(nextStorageZ)>= GridSize && nextStorageZ >= 0){
                nextZ = iz + 1;
                nextStorageZ = nextStorageZ - GridSize;
                storageZ = 0.0;
            }else if(abs(nextStorageZ)> GridSize && nextStorageZ < 0){
                nextZ = iz - 1;
                nextStorageZ = nextStorageZ + GridSize;
                storageZ = 0.0;
            }
            
            currentP.setIndexX(nextX);
            currentP.setIndexY(nextY);
            currentP.setIndexZ(nextZ);
            
            //System.out.println("(nextX , nextY , nextZ) ="+nextX+"\t"+nextY+"\t"+nextZ);
            px = nextX * GridSize + nextStorageX;
            py = nextY * GridSize + nextStorageY;
            pz = nextZ * GridSize + nextStorageZ;
            //System.out.println("(px , py , pz) ="+px+"\t"+py+"\t"+pz);

            currentP.setPx(px);
            currentP.setPy(py);
            currentP.setPz(pz);
                        
            currentP.setStorageX(nextStorageX);
            currentP.setStorageY(nextStorageY); 
            currentP.setStorageZ(nextStorageZ); 
            
            int nX =  (int)(Math.floor((px+CenterX-XllCorner)/GridSize));
            int nY =  (int)(Math.floor((py+CenterY-YllCorner)/GridSize));
            if(nX > DATAWIDTH){
                String errorsTxt = "Particle exceeds the boundary in x direction + ";
                errorDialog(errorsTxt);
            }else if(nX < 0){
                String errorsTxt = "Particle exceeds the boundary in x direction - ";
                errorDialog(errorsTxt);
            }
            
            if(nY > DATAHEIGHT){
                String errorsTxt = "Particle exceeds the boundary in y direction + ";
                System.out.println("nY = "+nY + " NY= "+DATAHEIGHT);
                errorDialog(errorsTxt);
            }else if(nY < 0){
                String errorsTxt = "Particle exceeds the boundary in y direction - ";
                errorDialog(errorsTxt);
            }
            
            double groundZ = Altitude[nY][nX];
            
            /*--If particle's position is lower than ground position, then this particle is deposited --*/            
            if(groundZ > pz){
                currentP.setPz(groundZ);
                deposit( currentP);//A particle is added to the deposit array
                /*--trajectory data is stored into the trajectory array--*/
                //trajectories.add(new Trajectory(ctime, currentP.getId(), px, py,pz,currentP.collisioncounter));
                
                //the current particle is removed from "particles" array 
                particles.remove(currentP);
            
            /*--If the particle is in the air, its velocity is updated --*/
            }else{
                double dz = pz - CenterZ;
                double dist = sqrt(px *px + py*py + dz*dz);
                //Velocity is calculated with wind velocity
                double[] vel = velocity.getVelocityWind(currentP, Dt, dist);
                currentP.setVelocity(vel[0], vel[1], vel[2]);
                /** A particle is in the air **/
                if(abs(ctime-round(ctime/OutDt) * OutDt) < (1e-8)){
//                    trajectories.add(new Trajectory(ctime, currentP.getId(), px, py,pz,currentP.collisioncounter));
                    double outtime = ctime-1.0;///Why?
                    double x = px + CenterX;
                    double y = py + CenterY;
                    
                    String thistraje = outtime+"\t"+currentP.getId()+"\t"+x+"\t"+ y+ "\t"+ pz;
                    outTraje.println(thistraje);
                }
            }
        }
        
    }//Finish transport
    

    public void deposit( Particle pp){
        double dx = pp.getPx();        double dy = pp.getPy();
        double x = dx+CenterX;        double y = dy+CenterY;        double z = pp.getPz();
        double dist = sqrt(dx*dx + dy*dy);
        double vx = pp.getVelocity()[0];
        double vy = pp.getVelocity()[1];
        double vz = pp.getVelocity()[2];
        double vnorm = sqrt(vx*vx + vy*vy + vz*vz);
        double costh  = (vx*vx + vy*vy)/(sqrt(vx*vx+vy*vy+vz*vz)* sqrt(vx*vx+ vy*vy) );
        double deg = acos(costh)/PI*180;
        
        String outdata = (x+"\t"+y+"\t"+z+"\t"+dist+"\t"+pp.getMass()+"\t"
               +vnorm+"\t"+pp.getDensity()+"\t"+pp.getDiameter()+"\t"+pp.getId());
        outDepo.println(outdata);
        
//        int inx =(int)floor(nextX/GridSize)+DATAWIDTH/2;//delete later
//        int iny =(int)floor(nextY/GridSize)+DATAHEIGHT/2;//delete later
//        depositParticles[iny][inx].add(pp);//delete later
//        depoparticles.add(pp);//delete later
        particles.remove(pp);
    }
    
//    public void deposit(int nextX, int nextY, Particle currentP, double groundZ){
//        int inx =(int)floor(nextX/GridSize)+DATAWIDTH/2;
//        int iny =(int)floor(nextY/GridSize)+DATAHEIGHT/2;
//        depositParticles[iny][inx].add(currentP);
//        depoparticles.add(currentP);
//        particles.remove(currentP);
//        x = CenterX+
//        dist = nextX
//        
////        sb.append(x+"\t"+y+"\t"+z+"\t"+dist+"\t");
////        sb.append(pp.getMass()+"\t"+vnorm+"\t"+pp.getDensity()+"\t");
////        sb.append(pp.getDiameter()+"\t"+pp.getId());
//        String thisdeposition = nextX +"\t" + nextY +"\t" +groundZ +"\t"+ 
//        
//        outDepo.println();
//    }
     
    private List makeParticleList(){
        List pList = new ArrayList();
        return pList;
    }
    
    private List[][] makeParticleField2D() {
        List[][] newField2D = new List[DATAHEIGHT][DATAWIDTH];
        for( int y=0; y<DATAHEIGHT; y++ ) {
            for( int x=0; x<DATAWIDTH; x++ ) {
                newField2D[y][x] = new ArrayList();
            }
        }
        return newField2D;
    }
    
    private List makeTrajectoryList(){
        List tList = new ArrayList();
        return tList;
    }

    public void errorDialog(String errorsTxt){
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle("Error Dialog");
        alert.setHeaderText(errorsTxt);
        alert.setContentText("Please set the larger calculation area.");
        ButtonType buttonTypeCancel = new ButtonType("Ok and Exit System", ButtonBar.ButtonData.CANCEL_CLOSE);
        
        alert.getButtonTypes().setAll( buttonTypeCancel);
        Stage stage = (Stage) alert.getDialogPane().getScene().getWindow();
        stage.getIcons().add(new Image("file:images/IconBallista.png"));
        Optional<ButtonType> result = alert.showAndWait();
        if (result.get() == buttonTypeCancel){
            System.exit(-1);
        }
    }
    
    //For outputting the file
     private static PrintWriter name2printWriter(String filename) throws IOException {
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
        return pw;
    }
  
}
