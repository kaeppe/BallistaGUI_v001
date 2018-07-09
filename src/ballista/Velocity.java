/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista;

import static ballista.utility.ConstParam.*;
import static ballista.utility.Coordinate.*;
import static java.lang.Math.*;
import java.util.Random;

/**
 *
 * @author tsunemat
 */
public class Velocity {
    private final Random RNG = new Random();
    private Particle currentp;
    
    private int numEq = 3;

    public Velocity() {
        

    }
    
    
    
    public double[] getVelocityWind(Particle pp, double dt, double dist){
        this.currentp = pp;

        double ufx = windX;
        double ufy = windY;
        double ufz = 0.0;

        double[] velocity = vRunge(dt, ufx, ufy, ufz);
//        if(ufx >0)
//            System.out.println(ufx+"\t"+ufz+";\t"+velocity[0]+"\t"+velocity[1]+"\t"+velocity[2]);
        return velocity;
    }
    

    
    public double[] vRunge(double dt, double ufx, double ufy, double ufz){
        int numPara = 3; //ufx, ufy, ufz(Flow Velocities)
        
        /**-- x is unknown parameter. Here it is velocity--**/
        double[] v = new double[numEq];
        double[] gradV = new double[numEq];
        double[] uPara = new double[numPara];
        
        double[] k0 = new double[numEq];
        double[] k1 = new double[numEq];
        double[] k2 = new double[numEq];
        double[] k3 = new double[numEq];

        double[] v1 = new double[numEq];
        double[] v2 = new double[numEq];
        double[] v3 = new double[numEq];
        
        v[X] = currentp.getVelocity()[0];
        v[Y] = currentp.getVelocity()[1];
        v[Z]= currentp.getVelocity()[2];
        uPara[X] = ufx;
        uPara[Y] = ufy;
        uPara[Z] = ufz;
        
        /***** Runge Kutta *****/ /*
        k1 = h * defineFunction(x, y);
        k2 = h * defineFunction(x + h/2, y + k1/2);
        k3 = h * defineFunction(x + h/2, y + k2/2);
        k4 = h * defineFunction(x + h, y + k3);
        */
//		0
        gradV = functions( v, uPara);
        for(int i=0; i<numEq; i++){
                k0[i] = dt*gradV[i];
                v1[i]=v[i] + k0[i]/2.0;
        }

//		1
        gradV = functions( v1, uPara);
        for(int i=0; i<numEq; i++){
                k1[i]=dt*gradV[i];
                v2[i]=v[i] + k1[i]/2.0;
        }


//		2
        gradV = functions( v2, uPara);
        for(int i=0; i<numEq; i++){
                k2[i]=dt*gradV[i];
                v3[i]=v[i]+k2[i];
        }	

//		3
        gradV = functions( v3, uPara);		
        for(int i=0; i<numEq;i++){
                k3[i]=dt*gradV[i];
                v[i]=v[i]+(k0[i] + 2.0*k1[i] + 2.0*k2[i] + k3[i])/6.0;// 加重平均
        }
        
        //System.out.println("(Vx ,Vy, Vz)="+v[X]+"\t"+v[Y]+"\t"+v[Z]);
        return v;
    }
    
    
    public double[] functions( double[] v, double[] uParam){
        double[] dv = new double [numEq];
        
        final double A = PI*pow((currentp.getDiameter()/2.0),2);//Particle surface area
        double beta = A*rhoa*Cd/currentp.getMass()/2;
        
        double ufx = uParam[0];
        double ufy = uParam[1];
        double ufz = uParam[2];
        
        double absV = sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);//Particle Velocity
        double absU = sqrt(ufx*ufx + ufy*ufy + ufz*ufz);//Flow Velocity
        
        dv[X] = -( beta *(v[X]-ufx)*abs(absV-absU)  );
        dv[Y] = -( beta *(v[Y]-ufy)*abs(absV-absU) );
        dv[Z] = -( beta *(v[Z]-ufz)*abs(absV-absU) )-G;
        //System.out.println("t ="+t+", DVx="+dv[X]+", DVy="+dv[Y]+", DVz="+dv[Z]);
       
        return dv;
    }
    
    
    
}
