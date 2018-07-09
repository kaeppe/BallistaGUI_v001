package ballista.utility;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import ballista.Particle;
import static java.lang.Math.*;
import ballista.Trajectory;
import static ballista.utility.ConstParam.*;

public class Output {


    /*
     * output particle list
     */
    public static void writeCoord1D(String filename, List  particlelist){
        if(filename == "*burst.txt*")
            System.out.println("Printing burst points distribution");
        try {
            PrintWriter pw = name2printWriter(filename);
           
            StringBuilder sb = new StringBuilder();
            sb.append("X\tY\tZ\tDistance\tMass\tVelocity\tDenisty\tDiameter\tId\n");
            for(int i=0;i<particlelist.size();i++){
                Particle pp = (Particle) particlelist.get(i);
                double deltax =pp.getPx();
                double deltay = pp.getPy();
                double x = deltax + CenterX;
                double y = deltay + CenterY;
                double z = pp.getPz();
                double dist = sqrt(deltax*deltax + deltay*deltay);
                double vx = pp.getVelocity()[0];
                double vy = pp.getVelocity()[1];
                double vz = pp.getVelocity()[2];
                double vnorm = sqrt(vx*vx + vy*vy + vz*vz);
                
                double costh  = (vx*vx + vy*vy)/(sqrt(vx*vx+vy*vy+vz*vz)* sqrt(vx*vx+ vy*vy) );
                double deg = acos(costh)/PI*180;
                sb.append(x+"\t"+y+"\t"+z+"\t"+dist+"\t");
                sb.append(pp.getMass()+"\t"+vnorm+"\t"+pp.getDensity()+"\t");
                sb.append(pp.getDiameter()+"\t"+pp.getId());
                sb.append("\n");
                              
            }
            pw.println(sb.toString());
            if(filename == "burst.txt")
                System.out.println("File for burst distribution is written");
            pw.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public static void writeTrajectory3D(String filename, List  trajectories){
        System.out.println("In Output, trajectory size = "+trajectories.size());
        try {
            PrintWriter pw = name2printWriter(filename);
            StringBuilder sb = new StringBuilder();
            sb.append("Time\tId\tX\tY\tZ\n");
            for(int i=0;i<trajectories.size();i++){
                Trajectory tj = (Trajectory) trajectories.get(i);
                sb.append((tj.getCtime()-1.0)+"\t"+tj.getId()+"\t");
                double x = tj.getPx()+CenterX;
                double y = tj.getPy()+CenterY;
//                sb.append(x+"\t"+y+"\t"+tj.getPz()+"\t"+tj.getCollNum()+"\n");
                sb.append(x+"\t"+y+"\t"+tj.getPz()+"\n");

            }
            pw.println(sb.toString());
            pw.close();
            
        }catch(IOException e) {
            e.printStackTrace();
        }
    }

    private static PrintWriter name2printWriter(String filename) throws IOException {
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
        return pw;
    }
}
