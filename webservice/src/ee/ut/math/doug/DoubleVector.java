package ee.ut.math.doug;
import java.util.Arrays;

/**
 * 
 */

/**
 * @author Christian P�cher
 *
 * Standard mathematical vector of doubles. Also implements equality of
 * vectors in the usual mathemathical meaning.
 */
public class DoubleVector {
	
	private double[] vector;
	
    /**
	 * Set of vector only through Constructor.
	 * 
	 * @param vector
	 */
	public DoubleVector(double[] vector) {
		this.vector = vector;
	}
	
	public double[] getVector() {
		return vector;
	}
	
	public DoubleVector div(double divisor) {
		double[] result = new double[vector.length];
		for (int i=0; i < vector.length; i++) {
			result[i] = vector[i] / divisor;
		}
		return new DoubleVector(result);
	}
	
	/**
	 * Calculates the euclidean norm (2-Norm, ||x||_2)
	 * @return the norm of this DoubleVector
	 */
	public double norm() {
		double norm = 0;
		for (int i=0; i < vector.length; i++) {
			norm += vector[i] * vector[i];
		}
		norm = Math.sqrt(norm);
		return norm;
	}
	
	/**
	 * @see java.utils.Arrays#equals(double[], double[])
	 */
	public boolean equals(Object o) {
		double[] testVec = (double[]) o;
		boolean result = Arrays.equals(vector, testVec);
		return result;
	}
	
	/**
	 * @see java.lang.Double#hashCode()
	 */
    public int hashCode() {
    	long elementBits;
		long bits = 0;
		for (int i = 0; i < vector.length; i++) {
			elementBits = Double.doubleToLongBits(vector[i]);
			bits = bits ^ elementBits ^ (elementBits >>> 32);
		}
		return (int) bits;
    }

    /**
     * Vector multiplication
     * @param v vector to multiply <code>this</code> with.
     * @return result
     */
	public double timesVector(DoubleVector v) {
		if (v.getVector().length != vector.length)
			throw new ArithmeticException("Vector lengths don't match.");
		double retVal = 0;
		for (int i=0; i < vector.length; i++) {
			retVal += v.getVector()[i] * vector[i];
		}
		return retVal;
	}

	/**
	 * Scalar multiplication
	 * @param s scalar to multiply vector with
	 * @return new vector with result
	 */
	public DoubleVector timesScalar(double s) {
		double[] result = new double[vector.length];
		for(int i=0; i < vector.length; i++) {
			result[i] = vector[i] * s;
		}
		return new DoubleVector(result);
	}

	/**
	 * Subtraction
	 * @param v vector to subtract
	 * @return new vector with result
	 */
	public DoubleVector minus(DoubleVector v) {
		if (v.getVector().length != vector.length)
			throw new ArithmeticException("Vector lengths don't match.");
		double[] result = new double[vector.length];
		for(int i=0; i < vector.length; i++) {
			result[i] = vector[i] - v.getVector()[i];
		}
		return new DoubleVector(result);
	}
}
