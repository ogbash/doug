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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

/**
 * 
 */

/**
 * Implementation of a sparse matrix. Relys on primitives for efficiency.
 * 
 * @author Christian Pöcher
 *
 */
/* Considered making immutable for readability, but efficiency was against it. */
public class AssembledMatrix {

	private int[] indi, indj;
	private double[] val;
	private int nnz, ncols, nrows;
	private int filledUpTo;
	
	private final static String FILE="/home/poecher/dump_matrix.file";  //TODO: put in properties
	
	public AssembledMatrix(int nnz, int ncols, int nrows) {
		indi = new int[nnz];
		indj = new int[nnz];
		val  = new double[nnz];
		this.nnz = nnz;
		this.ncols = ncols;
		this.nrows = nrows;
		filledUpTo = -1;
	}
	
	public AssembledMatrix(int nnz) {
		this(nnz, -1, -1);
	}
	
/* 
 * Needed for immutable implementation, which is not done here
 * because of efficiency concerns. Can be uncommented as soon as hashCode() is
 * implemented.
 */ 	 
//	protected Object clone() {
//		try {
//			AssembledMatrix clone = (AssembledMatrix) super.clone();
//			clone.indi = (int[])indi.clone();
//			clone.indj = (int[])indj.clone();
//			clone.val = (double[])val.clone();
//			return clone;
//		} catch (CloneNotSupportedException ex) {
//			return null; //never invoked
//		}
//	}
//	
//	public boolean equals(Object o) {
//	    if (this == o)
//			return true;
//	    if (o == null)
//			return false;
//		if (!(o instanceof AssembledMatrix))
//			return false;
//		
//		AssembledMatrix m = (AssembledMatrix) o;
//		if ( (this.nnz != m.nnz) || (this.ncols != m.ncols) 
//				|| (this.nrows != m.nrows) )
//			return false;
//		for (int k = 0; k < val.length; k++) {
//			if ((this.indi[k] != m.indi[k]) || (this.indj[k] != m.indj[k])
//					|| (this.val[k] != m.val[k]))
//				return false;
//		}
//		/* still here? then the two is equal! */
//		return true;
//	}
	
	/**
	 * Generates human-readable representation of matrix.
	 * 
	 * @return human-readable String of matrix
	 */
	public String print() {
		StringBuffer buf = new StringBuffer();
		buf.append("nnz = " + nnz + ";  ncols = " + ncols + ";  nrows = " + nrows + "\n");
		for (int k=0; k < val.length; k++) {
			buf.append("i: " + indi[k] + "\tj: " + indj[k] + "\tval: "+ val[k] + "\n");
		}
		return buf.toString();
	}
	
	/**
	 * Outputs String representation in standard DOUG format.
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append(ncols + " " + nnz + "\n");
		for (int k = 0; k<nnz; k++) {
			buf.append(indi[k] + " " +  indj[k] + " " + val[k] + "\n");
		}
		return buf.toString();
	}
	
	/**
	 * no checking for space or double entries
	 * 
	 * @param indi
	 * @param indj
	 * @param val
	 */
	public void addElement(int indi, int indj, double val) {
		filledUpTo++;
		this.indi[filledUpTo] = indi;
		this.indj[filledUpTo] = indj;
		this.val[filledUpTo] = val;
	}
	
	public boolean isFull() {
		return (filledUpTo == val.length - 1);
	}
	
	public void resize(int newNnz) {
		int[] indiTemp = new int[newNnz];
		int[] indjTemp = new int[newNnz];
		double[] valTemp = new double[newNnz];
		if (val.length == newNnz) {
			return;
		} else if (val.length < newNnz) {
			/* Enlarge */
			System.arraycopy(indi, 0, indiTemp, 0, indi.length);
			System.arraycopy(indj, 0, indjTemp, 0, indj.length);
			System.arraycopy(val, 0, valTemp, 0, val.length);
		} else if (val.length > newNnz) {
			/* Shrink */
			System.arraycopy(indi, 0, indiTemp, 0, newNnz);
			System.arraycopy(indj, 0, indjTemp, 0, newNnz);
			System.arraycopy(val, 0, valTemp, 0, newNnz);
		}
		indi = indiTemp;
		indj = indjTemp;
		val = valTemp;
	}
	
	/*
	 * beware of index shifts because of different counting ways in java and fortran.
	 */
	/**
	 * Adds sigma times the identity to the matrix.
	 * That is, A = A + sigma * I 
	 * @param sigma
	 */
	public void plusSigmaTimesIdentity(double sigma) {
		int minn = Math.min(ncols, nrows);
		boolean[] elementAdded = new boolean[minn];
		int counter = 0;
		/* search for existing elements on the diagonals and add */
		for (int k=0; k < nnz; k++) {
			if (indi[k] != indj[k]) continue;
			elementAdded[indi[k]-1] = true;
			val[k] += sigma;
		}
		/* count and resize Matrix */
		for (int k=0; k < elementAdded.length; k++) {
			if (elementAdded[k] == false) {
				counter++;
			}
		}
		if (counter == 0) return; //No need to iterate elementAdded[] without need.
		resize(nnz + counter);
		/* add missing */
		for (int k=0; k < elementAdded.length; k++) {
			if (elementAdded[k] == false) {
				addElement(k+1, k+1, sigma);
			}
		}
	}

	/**
	 * Writes Matrix to file. Uses same format as
	 * DOUG uses for input file.
	 * 
	 * @throws IOException
	 */
	public void writeToDisk(String fname) throws IOException {
		FileWriter writer = new FileWriter(fname, false);
		writer.write(this.toString());
		writer.flush();
		writer.close();
	}
	
	/*
	 * probably more efficient with RegExp instead of Tokenizer,
	 * but harder to understand.
	 */
	/**
	 * Reads a AssembledMatrix in sparse format from a Reader. Assumes same format as
	 * DOUG assumes for input file.
	 * 
	 * @param r A reader containing the matrix data. No need to use a BufferedReader,
	 * as it is buffered by this method.
	 * @return The AssembledMatrix object read from Reader.
	 * @throws IOException If a read operation on the Reader fails.
	 * @throws FileFormatException If the data is different, than expected.
	 */
	public static AssembledMatrix readFromReader(Reader r) throws IOException, DougServiceException {
		AssembledMatrix m;
		String s;
		StringTokenizer toki;
		int nnz, ncols, indi, indj;
		double val;
		BufferedReader buf = new BufferedReader(r);
		try {
			/* header line */
			s = buf.readLine();
			s = s.trim();
			toki = new StringTokenizer(s);
			ncols = Integer.parseInt(toki.nextToken());
			nnz = Integer.parseInt(toki.nextToken());
			m = new AssembledMatrix(nnz, ncols, ncols);
			/* body */
			for (int k=0; k<nnz; k++) {
				s = buf.readLine();
				s = s.trim();
				toki = new StringTokenizer(s);
				indi = Integer.parseInt(toki.nextToken());
				indj = Integer.parseInt(toki.nextToken());
				val = Double.parseDouble(toki.nextToken());
				m.addElement(indi, indj, val);
			}
		} catch (NumberFormatException ex) {
			throw new DougServiceException("A number could not be parsed correctly. Wrong format.", ex);
		} catch (NoSuchElementException ex) {
			throw new DougServiceException("Attemted to read more data than availible. Is nnz not correct?", ex);
		}
		r.close();
		return m;
	}
	
	/**
	 * Reads a AssembledMatrix in sparse format from a Reader. Assumes same format as
	 * DOUG assumes for input file.
	 * 
	 * @return The AssembledMatrix object read from file.
	 * @throws IOException If a read operation on the File fails or the file is not readable.
	 * @throws FileFormatException If the data is different, than expected.
	 * @throws DougServiceException 
	 */
	public static AssembledMatrix readFromDisk() throws IOException, DougServiceException {
		File file = new File(FILE);
		FileReader reader = new FileReader(file);
		AssembledMatrix m = readFromReader(reader);
		return m;
	}
	
	public static void main(String[] args) throws Exception {
		AssembledMatrix a = AssembledMatrix.readFromDisk();
		a.plusSigmaTimesIdentity(9.99);
		System.out.print(a.print());
		a.writeToDisk("wrote.txt");
	}
	
}
