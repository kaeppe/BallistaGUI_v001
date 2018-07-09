/*
 * This code is for Lagrangian transport for ballistics
 */
package ballista;

import java.util.ArrayList;
import java.util.List;
import static java.lang.Math.*;
import static ballista.utility.ConstParam.*;
import static ballista.utility.Coordinate.*;

/**
 *
 * @author kae
 */
public class Particle {
    private int id;
    private double mass;
    private double diameter; //(m)
    private double density; //kg/m^3
    private double currentTime; //(s)
    
    //Position
    private double px;
    private double py;
    private double pz;
    
    //Indices
    private int indexX;
    private int indexY;
    private int indexZ;
    
    //Velocity
    private double[] velocity;
    
    //Storage
    private double storageX;
    private double storageY;
    private double storageZ;
    
    public int collisioncounter;
    
  
    
    /**
     * This constructor works only when the particle is drawn initially.
     * @param diam
     * @param density
     * @param x
     * @param v
     * @param time 
     */
    public Particle(int id, double diam, double density, double[] x, double[] v,  double time) {
        this.id = id;
        this.diameter = diam;
        this.density = density;
        this.mass = PI*pow(diameter, 3.0)/6.0*density;
        this.px = x[X];
        this.py = x[Y];
        this.pz = x[Z]; //on the ground, z = 0.0
        this.velocity = v;
        this.indexX =(int) (px/GridSize);
        this.indexY =(int) (py/GridSize);
        this.indexZ =(int) (pz/GridSize);
        this.calcStorageX();
        this.calcStorageY();
        this.calcStorageZ();
        this.collisioncounter = 0;
        this.currentTime = time;
        //this.pileTraje = new ArrayList();
        
    }

    public double getCurrentTime() {
        return currentTime;
    }

    public void setCurrentTime(double currentTime) {
        this.currentTime = currentTime;
    }

    public int getId() {
        return id;
    }

    public double getDiameter() {
        return diameter;
    }

    public double getMass() {
        return mass;
    }   
    
    public double getDensity() {
        return density;
    }  
    

    public void setPx(double px) {
        this.px = px;
    }
    public double getPx() {
        return px;
    }

    public void setPy(double py) {
        this.py = py;
    }
    
    public double getPy() {
        return py;
    }
    
    public void setPz(double pz) {
        this.pz = pz;
    }
    
    public double getPz() {
        return pz;
    }

      

    public void setIndexX(int indexX) {
        this.indexX = indexX;
    }
    public int getIndexX() {
        return indexX;
    }

    public int getIndexY() {
        return indexY;
    }

    public void setIndexY(int indexY) {
        
        this.indexY = indexY;
    }
    
    public int getIndexZ() {
        return indexZ;
    }

    public void setIndexZ(int indexZ) {
        this.indexZ = indexZ;
    }

    public double[] getVelocity() {
        return velocity;
    }
    

    public void setVelocity(double vx, double vy, double vz) {
        double [] v = new double [D];
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
        this.velocity = v;
    }

    public double getStorageX() {
        return storageX;
    }
    public double getStorageY() {
        return storageY;
    }
    public double getStorageZ() {
        return storageZ;
    }

    /*If the storages are given manually */
    public void setStorageX(double storageX) {
        this.storageX = storageX;
    }

    public void setStorageY(double storageY) {
        this.storageY = storageY;
    }
    
    public void setStorageZ(double storageZ) {
        this.storageZ = storageZ;
    }
    
    /*** Automatically calculated storage***/
    private void calcStorageX(){
        this.storageX = px%GridSize;
        
    }
    private void calcStorageY(){
        this.storageY = py%GridSize;  
    }
    
    private void calcStorageZ(){
        this.storageZ = pz%GridSize;  
    }
    
    public void incrCollisionCounter(){
        collisioncounter ++;
    }
    
//    public void makeTrajectory(double time, double px,double py, double pz ){
//        new Trajectory(time, px, py, pz);
//        
//    }
}
