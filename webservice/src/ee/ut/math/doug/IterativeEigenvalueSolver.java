package ee.ut.math.doug;

import java.io.IOException;


/**
 * 
 */

/**
 * @author poecher
 *
 */
public class IterativeEigenvalueSolver {

	/**
	 * Implemented as in http://www.cs.utk.edu/~dongarra/etemplates/node96.html
	 * 
	 * @param initialGuess
	 * @param shift
	 * @return EigenSpace with one Eigenvalue and at least one Eigenvector
	 */
	public EigenSpace inverseIteration(AssembledMatrix a, DoubleVector initialGuess, 
			double shift, double error) throws IOException {
		
		EigenSpace eSpace;
		DoubleVector v,y;
		double sigma;
		
		y = initialGuess;
		a.plusSigmaTimesIdentity(shift * -1.0);
		a.writeToDisk();
		while (true) {
			v = y.div( y.norm() );
			y = solve(a, v);
			sigma = y.timesVector(v);
			if (y.minus(v.timesScalar(sigma)).norm() <= error * Math.abs(sigma))
				break;
		}
		eSpace = new EigenSpace(shift + 1 / sigma);
		eSpace.addEigenVector(y.div(sigma));
		return eSpace;
	}

	/**
	 * Solves the linea equation system Ax=v.
	 * @param A the matrix
	 * @param v the vector
	 * @return the resultvector x
	 */
	private DoubleVector solve(AssembledMatrix a, DoubleVector v) {
		DoubleVector x = null;
		// TODO Code, rely on DougWSClient
		return x;
	}
}
