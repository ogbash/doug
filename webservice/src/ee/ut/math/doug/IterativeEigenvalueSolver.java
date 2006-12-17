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

import java.io.IOException;

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
