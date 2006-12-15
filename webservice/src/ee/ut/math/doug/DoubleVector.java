// DOUG - Domain decomposition On Unstructured Grids
// Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
// Department of Mathematics, University of Bath
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// or contact the authors (University of Tartu, Faculty of Computer Science, Chair
// of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
// mailto:info(at)dougdevel.org)

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
