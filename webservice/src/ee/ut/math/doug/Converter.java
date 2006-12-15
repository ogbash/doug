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

import java.io.File;

public class Converter {

	/**
	 * Converts elemental into assembled by using webservice.
	 * 
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		DougWSClient client = new DougWSClient();
		File freedom_lists_file = new File("/home/poecher/srcdata/doug_element.dat");
		File coords_file        = new File("/home/poecher/srcdata/doug_coord.dat");
		File elemmat_rhs_file   = new File("/home/poecher/srcdata/doug_system.dat");
		File freemap_file       = new File("/home/poecher/srcdata/doug_freemap.dat");
		File freedom_mask_file  = new File("/home/poecher/srcdata/dummy");
		File info_file          = new File("/home/poecher/srcdata/doug_info.dat");
		AssembledMatrix matrix = client.elementalToAssembled(freedom_lists_file,
				elemmat_rhs_file, coords_file, freemap_file, freedom_mask_file, 
				info_file);
		System.out.println(matrix);
	}

}
