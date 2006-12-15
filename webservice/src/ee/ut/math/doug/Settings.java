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

/*
 * Import by implements good? Don't think so.
 */
// TODO: put into properties? 
/**
 * Settings for DOUG services.
 * 
 * @author Christian Pöcher
 *
 */
public interface Settings {
	public static final String CONTROL_FILE = "DOUG.dat";
	public static final String INFO_FILE = "doug_info.dat";
	public static final String FREEDOM_LISTS_FILE = "doug_element.dat";
	public static final String ELEMENT_RHS_FILE = "doug_system.dat";
	public static final String COORDS_FILE = "doug_coord.dat";
	public static final String FREEMAP_FILE = "doug_freemap.dat";
	public static final String FREEDOM_MASK_FILE = "doug_freemask.dat";
	public static final String SOLUTION_FILE = "solution.file";
	public static final String DUMP_MATRIX_FILE = "assembled.txt";
	public static final String ENDPOINT_ADRESS = "http://localhost:8080/axis/services/doug";
}
