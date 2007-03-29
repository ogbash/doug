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
public class IterativeEigenvalueSolver implements I_IterativeEigenvalueSolver {

	/* (non-Javadoc)
	 * @see ee.ut.math.doug.I_IterativeEigenvalueSolver#inverseIteration(ee.ut.math.doug.AssembledMatrix, ee.ut.math.doug.DoubleVector, double, double)
	 */
	public EigenSpace inverseIteration(AssembledMatrix a, DoubleVector initialGuess, 
			double shift, double error) throws IOException, DougServiceException {
		EigenSpace eSpace;
		DoubleVector v,y;
		double sigma, left, right;
		int iteration = 0;
		
		y = initialGuess;
		a.plusSigmaTimesIdentity(shift * -1.0);
		while (true) {
			System.out.println("Iteration " + iteration++);
			v = y.div( y.norm() );
			y = solve(a, v);
			sigma = y.timesVector(v);
			left = y.minus(v.timesScalar(sigma)).norm();
			right = error * Math.abs(sigma);
			if (left <= right)
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
	private DoubleVector solve(AssembledMatrix a, DoubleVector v) throws IOException, DougServiceException {
		//TODO: Change DougServiceException into domain internally Exception
		DougWSClient ws = new DougWSClient();
		DoubleVector x = ws.runAssebled(a, v);
		return x;
	}
	
	private static void test() {
		AssembledMatrix a = new AssembledMatrix(9, 3, 3);
		a.addElement(1, 1, 2.23);
		a.addElement(1,2, 1.15);
		a.addElement(1,3, 1.77);
		a.addElement(2,1, -1.15);
		a.addElement(2, 2, 9.25);
		a.addElement(2,3, -2.13);
		a.addElement(3,1, 1.77);
		a.addElement(3,2, 2.13);
		a.addElement(3, 3, 1.56);
		
		DoubleVector initialGuess = new DoubleVector(new double[] {1,0,0});
		double shift = 10;
		double error = 0.1;
		
		IterativeEigenvalueSolver solver = new IterativeEigenvalueSolver();
		
		try {
			EigenSpace es = solver.inverseIteration(a, initialGuess, shift, error);
			System.out.print(es);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		IterativeEigenvalueSolver.test();
//		AssembledMatrix a = AssembledMatrix.readFromDisk(args[0]);
//		DoubleVector initialGuess = DoubleVector.readFromDisk(args[1]);
//		double shift = Double.parseDouble(args[2]);
//		double error = Double.parseDouble(args[3]);
//		
//		IterativeEigenvalueSolver solver = new IterativeEigenvalueSolver();
//		EigenSpace es = solver.inverseIteration(a, initialGuess, shift, error);
//		System.out.print(es);
	}
}
