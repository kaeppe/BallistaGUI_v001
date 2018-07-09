/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ballista.utility;

/**
 *
 * @author tsunemat
 */
public class ConstParam {
    public static final int D = 3;
    public static final double Timebuff = 600;
    public static final double calcAX = 5000; //m x=[-250, 250]
    public static final double calcAY = 5000; //m y=[-250, 250]
    public static final double calcAZ = 1000; //m
    
    public static final double  Dt = 0.0001;
    public static final double  G  = 9.81; 
    public static final double OutDt = 0.05; //delta t(time) for output trajectory
    
 
    public static int NY;
    public static int NZ;
    

    
    public static final double rhoa = 0.9;

    
    /*---- Parameters for reading DEM ----*/
    public static int DATAWIDTH;
    public static int DATAHEIGHT;
    public static double XllCorner;
    public static double YllCorner;
    public static double GridSize;
    public static int NODATA_VALUE;
    
    // TODO: It is better to set these parameters from the file
    public static double CenterX ;//= 723588.0;//Ontake case
    public static double CenterY ;//= 3974463.0;//Ontake case
    public static double CenterZ;//Read from DEM file
    public static double [][] Altitude;
    
    /*-----Input Parameters ------*/
    public static String dirDEMName;
    public static double numP;

    public static double avgDensity;
    public static double sdDensity;
    public static double avgDiam;
    public static double sdDiam;
    public static double minDiam;
    public static double maxDiam;
    public static double sdDisp;
    public static double maxDisp;
    public static double avgV;
    public static double sdV;
    public static double ejcDeg;
    public static double direcBear;
    public static double vFNorm;
    public static double windX = 0.0;//WindVelocity in X direction
    public static double windY = 0.0;//WindVelocity in Y direction
    public static double Cd ;//Drag Coefficient. Mean of Alatorre-Ibargüengoitia and Delgado-Granados(2006)


}
