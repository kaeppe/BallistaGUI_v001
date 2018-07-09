/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ballista.math;

/**
 *
 * @author falcone
 */
public class Vector {

    public static double norm( final double[] vec ) {
        double sum = 0;
        for( int i=0; i<vec.length; i++ ) {
            sum += vec[i]*vec[i];
        }
        return Math.sqrt(sum);
    }

    public static double[] add( final double[] vec1, final double[] vec2 ) {
        checkSize(vec1, vec2);
        final double[] res = new double[vec1.length];
        for( int i=0; i < res.length; i++ ) {
            res[i] = vec1[i] + vec2[i];
        }
        return res;
    }

    public static double[] subs( final double[] vec1, final double[] vec2 ) {
        checkSize(vec1, vec2);
        final double[] res = new double[vec1.length];
        for( int i=0; i < res.length; i++ ) {
            res[i] = vec1[i] - vec2[i];
        }
        return res;
    }

    public static double dot( final double[] vec1, final double[] vec2 ) {
        checkSize(vec1, vec2);
        double sum = 0;
        for( int i=0; i<vec1.length; i++ ) {
            sum += vec1[i] * vec2[i];
        }
        return sum;
    }

     public static double[] divide( final double[] vec, final double k ) {
        if ( k== 0) {
            throw new IllegalStateException("Division by zero");
        }
        final double[] res = new double[vec.length];
        for( int i=0; i<vec.length; i++ ) {
            res[i] = vec[i] / k;
        }
        return res;
    }

   public static double[] mult( final double[] vec, final double k ) {
        final double[] res = new double[vec.length];
        for( int i=0; i<vec.length; i++ ) {
            res[i] = vec[i] * k;
        }
        return res;
   }

    private static void checkSize(final double[] vec1, final double[] vec2) throws IllegalArgumentException {
        if (vec1.length != vec2.length) {
            throw new IllegalArgumentException("Vectors have different length");
        }
    }




}
