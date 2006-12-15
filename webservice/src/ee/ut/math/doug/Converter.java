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
