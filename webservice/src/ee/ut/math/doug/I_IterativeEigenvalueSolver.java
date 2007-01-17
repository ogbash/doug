package ee.ut.math.doug;

import java.io.IOException;

public interface I_IterativeEigenvalueSolver {

	/**
	 * Implemented as in http://www.cs.utk.edu/~dongarra/etemplates/node96.html
	 * 
	 * @param initialGuess
	 * @param shift
	 * @return EigenSpace with one Eigenvalue and at least one Eigenvector
	 */
	public EigenSpace inverseIteration(AssembledMatrix a,
			DoubleVector initialGuess, double shift, double error)
			throws IOException, DougServiceException;

}