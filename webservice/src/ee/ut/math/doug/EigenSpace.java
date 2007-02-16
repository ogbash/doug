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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * @author Christian P�cher
 *
 * Not necessary complete Eigenspace: Only some or even none Eigenvectors 
 * may be given.
 */
public class EigenSpace {
	
	private double eigenValue;
	private Set eigenVectors = new HashSet();
	
	public void setEigenValue(double eigenValue) {
		this.eigenValue = eigenValue;
	}
	
	public double getEigenValue() {
		return eigenValue;
	}

	public void addEigenVector(DoubleVector vector) {
		eigenVectors.add(vector);
	}
	
	public DoubleVector firstEigenVector() {
		return (DoubleVector) eigenVectors.iterator().next();
	}
	
	public Set getEigenVectors() {
		return eigenVectors;
	}
	
	public EigenSpace(double eigenvalue) {
		setEigenValue(eigenvalue);
	}
	
	/**
	 * For compability with BeanSerializer only. Can be removed as soon as custom 
	 * serializer is written. 
	 * @depricated
	 */
	public EigenSpace() {
	}
	
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append("Eigenvalue: " + eigenValue + "\n");
		Iterator iter = eigenVectors.iterator();
		while (iter.hasNext()) {
			buf.append("Eigenvector:\n");
			buf.append(iter.next().toString());
		}
		return buf.toString();
	}

	/**
	 * For compability with BeanSerializer only. Can be removed as soon as custom 
	 * serializer is written. 
	 * @depricated
	 */
	public void setEigenVectors(Set eigenVectors) {
		this.eigenVectors = eigenVectors;
	}
}
