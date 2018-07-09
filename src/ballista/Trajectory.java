/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import static ballista.utility.ConstParam.*;

/**
 *
 * @author tsunemat
 */
public class Trajectory {
    private double ctime;
    private int id;
    private double px;
    private double py;
    private double pz;
    private int collnum;

    public Trajectory(double ctime, int id, double px, double py, double pz, int collnum) {
        this.ctime = ctime;
        this.id = id;
        this.px = px;
        this.py = py;
        this.pz = pz;
        this.collnum = collnum;
    }

    public double getCtime() {
        return ctime;
    }

    public int getId() {
        return id;
    }

    public double getPx() {
        return px;
    }

    public double getPy() {
        return py;
    }

    public double getPz() {
        return pz;
    }
    
    public int getCollNum(){
        return collnum;
    }
    
    
    
    
    
}
