/*
 * This is for the lagrangian transport of ballistics
 * 
 */
package ballista;

import static java.lang.Math.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static ballista.utility.Coordinate.*;
import static ballista.utility.ConstParam.*;
import java.util.Locale;

/**
 *
 * @author tsunemat
 */

class GaussianBurst {
    private final Random RNG = new Random();
    private final double time;
//    private final ReadInit params;

    private int id = 0;
    
    //public GaussianBurst(double time, ReadInit params) {
    public GaussianBurst(double time) {
        this.time = time;
   
    }
    
    private double drawNormal( double average, double sd ) {
        double res = RNG.nextGaussian()*sd + average;
        while( res < 0 || res > average*2 ) {
            res = RNG.nextGaussian()*sd + average;
        }
        return res;  
    }
    
    public List<Particle> drawPartList() {
        
        List<Particle> particles = new ArrayList<Particle>((int) numP);
        System.out.println();
        for (int i = 0; i < numP; i++) {
            particles.add(drawParticle());
            id++;
        }
        //System.out.println("numP in DrawParticleList="+particles.size());
        return particles;
    }
    
    public double time() {
        return time;
    }

    public boolean involves(Particle p) {
        return false;
    }
    
//    private int drawNumberOfParticles() {
//        final int n = (int) Math.ceil(RNG.nextGaussian() * sdNum + numP);
//        if (n < 0) {
//            return 1;
//        } else {
//            return n;
//        }
//    }
    
    private Particle drawParticle() {
        
        double density = -0.0;
        double diam = -0.0;
        double x[] = new double[3];//position
        double v[] = new double[3];//particle velocity
        double v2[] = new double[3];//particle velocity
        double vF[] = new double[3];//flow velocity
        double vF2[] = new double[3];//flow velocity
        double vNorm = -0.0;
        double theta = -0.0;
        double phi = -0.0;
        double gamma = -0.0;
        double gamma2 = -0.0;
        double vx = -0.0;
        double vy = -0.0;
        double vz = -0.0;
   
     
        while(v2[Z] <= 0.0){
            density = drawNormal(avgDensity, sdDensity);
             /*---- Particle diameter ----*/
             /*--- Gaussian ---*/
            diam = drawNormal(avgDiam,sdDiam);
            while(diam < minDiam){// minimum diameter is 25cm 
                diam = drawNormal(avgDiam,sdDiam);
            }

            
            /*---- Ejection Position ----*/
         
            double distance= maxDisp + 1;
            while(distance > maxDisp){
                //distance = maxDisp + 1;
                for (int i = 0; i < 2; i++) {
                    double random = RNG.nextDouble();
                    boolean flag = RNG.nextBoolean();
                    int sign;
                    if (flag ==true)
                        sign = 1;
                    else{
                        sign = -1;
                    }
                    x[i] = sign * random*maxDisp;
                }
                distance = sqrt(x[X]* x[X] + x[Y]+x[Y]);
                //System.out.println("x = "+x[X]+", y="+x[Y]+", distance ="+distance);

            }
            
            
            //If it is only vent center, it will be only one position
            int nX =  (int)(Math.floor((x[X]+CenterX-XllCorner)/GridSize));
            int nY =  (int)(Math.floor((x[Y]+CenterY-YllCorner)/GridSize));
            
//            System.out.println("x ="+x[X]+", y="+x[Y]);
//            System.out.println("nx ="+nX+", ny="+nY);
            x[Z] = Altitude[nY][nX];
            
         
            /*---  Velocity ---*/
            vNorm = drawNormal(avgV,sdV);

            theta = drawTheta();
            phi = RNG.nextDouble() * 2 * PI;
            gamma = ejcDeg/180 * PI;

            /***  Calculate Particle Velocity   ***/
            vx = vNorm * sin(theta) * cos(phi);
            vy = vNorm * sin(theta) * sin(phi);
            vz = vNorm * cos(theta);

            //Rotation around Y axis : inclined to  eastward (particle)
            v[X] =  cos(gamma) * vx + sin(gamma) * vz;
            v[Y] =  vy;
            v[Z] = -sin(gamma) * vx + cos(gamma) * vz;
            
            //Rotation of around Z axis: directed to northward (particle)
            gamma2 = direcBear/180 * PI;
            v2[X] =  cos(gamma2) * v[X] - sin(gamma2) * v[Y];
            v2[Y] =  sin(gamma2) * v[X] + cos(gamma2) * v[Y];
            v2[Z] =  v[Z];
            
        }
        

        return new Particle(id, diam, density, x, v2, time);
        
    }
    


    private double drawTheta() {
        double avgVAngle = 0.0;
        double sdVAngle = 15.0;
        double theta = PI * (RNG.nextGaussian() * sdVAngle + avgVAngle) / 180;
        while( theta < -PI/2 || theta > PI/2 ) {
           theta = PI * (RNG.nextGaussian() * sdVAngle + avgVAngle) / 180;
        }
        return theta;
    }
    /**
     * 
     * @param avg: mean of Gaussian Profile
     * @param std: standard deviation of Gaussian Profile
     * @param min
     * @param max
     * @return 
     */
    private double gaussianGS(double avg, double std, double min, double max){
        double ramdomDiam =0.0;
        while (ramdomDiam <minDiam || ramdomDiam >maxDiam){
            ramdomDiam = drawNormal(avgDiam,sdDiam);
        }
        return ramdomDiam;
    }
    
    
    
}

    
    
    

    
    
    
