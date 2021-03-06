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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import javax.activation.DataHandler;
import javax.activation.FileDataSource;

import org.apache.commons.io.IOUtils;

/**
 * Webservice, that starts DOUG synchronously and gives back the solution.
 * 
 * @author Christian Poecher
 */
/**
 * @author Christian Pöcher
 *
 */
public class DougService implements I_DougService {

    /* Dir and executable name fixed by producing WS for security reasons. */
    /* TODO: make editable properties file */
    private static final String WORKING_DIR = "/home/poecher/working/";
    private static final String NO_OF_PROCESSES = "1"; //TODO: buggy with #nodes > 1
    private static final String DOUG_MAIN_EXECUTABLE = "DOUG_main";
    private static final String DOUG_AGGR_EXECUTABLE = "DOUG_aggr";
    private static final String[] EXE_CONV = {"mpirun", "-np", "1", DOUG_MAIN_EXECUTABLE, "-q"}; // TODO: better run DOUG without mpirun, but cannot find file. why?
    private static final String[] EXE_MAIN = {"mpirun", "-np", NO_OF_PROCESSES, DOUG_MAIN_EXECUTABLE, "-q"};
    private static final String[] EXE_AGGR = {"mpirun", "-np", NO_OF_PROCESSES, DOUG_AGGR_EXECUTABLE, "-q"};
    private static final int CHUNKSIZE = 4096;

    /**
     * Starts DOUG synchronously and returns stdout and stderr. 
     * 
     * @deprecated Does not handle result correctly, if too large, since it relies
     * on the output to stdout. Use runAssebled or write a new method similar
     * to it.
     * @return String array consisting of stdout and stderr.
     */
    public String[] runDoug() {
        
        Process pro;
        StreamGobbler outGobbler, errGobbler;
        InputStreamReader outIsr, errIsr;
        StringBuffer outBuf, errBuf;
        String[] output = new String[2];
        char[] chunk = new char[CHUNKSIZE];
        int charsRead;
        
        try {
            pro = Runtime.getRuntime().exec(EXE_MAIN, null, new File(WORKING_DIR));
            outGobbler = new StreamGobbler( pro.getInputStream() );
            errGobbler = new StreamGobbler( pro.getErrorStream() );
            outIsr = new InputStreamReader( outGobbler );
            errIsr = new InputStreamReader( errGobbler );
            
            pro.waitFor();
            
            /* the StreamGobblers are full now, put them into a String */
            outBuf = new StringBuffer();
            errBuf = new StringBuffer();
            
            while (true) {
                charsRead = outIsr.read(chunk);
                if (charsRead == -1) 
                	break;
                outBuf.append(chunk, 0, charsRead);
            }
            while (true) {
                charsRead = errIsr.read(chunk);
                if (charsRead == -1) 
                	break;
                errBuf.append(chunk, 0, charsRead);
            }
            
            output[0] = outBuf.toString();
            output[1] = errBuf.toString();
            
            /* control characters < 0x20 are not allowed in XML. */
            /* simply trim Strings, which is removing all leading and trailing chars <= 0x20 */
            output[0] = output[0].trim();
            output[1] = output[1].trim();
            
        } catch (IOException ex) {
            ex.printStackTrace();
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
        return output;
    }

    /* (non-Javadoc)
	 * @see ee.ut.math.doug.IDougService#runAssembled(ee.ut.math.doug.AssembledMatrix, ee.ut.math.doug.DoubleVector, javax.activation.DataHandler)
	 */
    public DoubleVector runAssembled(AssembledMatrix matrix, DoubleVector rhs, DataHandler control_file) {
    	/* save files to working dir */
    	try {
	    	matrix.writeToDisk(WORKING_DIR + Settings.ASSEMBLED_MTX_FILE);
	    	rhs.writeToDisk(WORKING_DIR + Settings.RHS_FILE);
	    	writeFile(control_file, WORKING_DIR + Settings.CONTROL_FILE);
    	} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	}
    	/* run */ 
    	runExecutable(DOUG_AGGR_EXECUTABLE , EXE_AGGR);
    	/* return result */
    	DoubleVector solution = parseSolution();
    	return solution;
    }
    
    /* (non-Javadoc)
	 * @see ee.ut.math.doug.IDougService#runAssembled(javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler)
	 */
	public DataHandler runAssembled(DataHandler matrix, DataHandler rhs, DataHandler control_file) {
		/* save files to working dir */
		try {
			writeFile(matrix, WORKING_DIR + Settings.ASSEMBLED_MTX_FILE);
			writeFile(rhs, WORKING_DIR + Settings.RHS_FILE);
			writeFile(control_file, WORKING_DIR + Settings.CONTROL_FILE);
		} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
		}
    	/* run */ 
    	runExecutable(DOUG_AGGR_EXECUTABLE , EXE_AGGR);
    	/* return result */
    	File solutionFile = new File(WORKING_DIR + Settings.SOLUTION_FILE);
    	DataHandler solution = new DataHandler(new FileDataSource(solutionFile));
    	return solution;
	}
    
    /* (non-Javadoc)
	 * @see ee.ut.math.doug.IDougService#elementalToAssembled(javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler)
	 */
    public AssembledMatrix elementalToAssembled(DataHandler freedom_lists_file,
    		DataHandler elemmat_rhs_file, DataHandler coords_file,
    		DataHandler freemap_file, DataHandler freedom_mask_file,
    		DataHandler info_file, DataHandler control_file) {
    	
    	/* write files */ 
    	writeElementalData(freedom_lists_file, elemmat_rhs_file, coords_file,
    			freemap_file, freedom_mask_file, info_file, control_file);
    	/* run */ 
    	runExecutable(DOUG_MAIN_EXECUTABLE , EXE_CONV);
		/* return result */
    	AssembledMatrix matrix = null;
    	try {
    		matrix = AssembledMatrix.readFromDisk(WORKING_DIR + Settings.DUMP_MATRIX_FILE);
    	} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	} 
    	return matrix;
    }
    
    /* (non-Javadoc)
	 * @see ee.ut.math.doug.IDougService#runElemental(javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler, javax.activation.DataHandler)
	 */
    public DoubleVector runElemental(DataHandler freedom_lists_file,
    		DataHandler elemmat_rhs_file, DataHandler coords_file,
    		DataHandler freemap_file, DataHandler freedom_mask_file,
    		DataHandler info_file, DataHandler control_file) {
    	
    	/* write files */ 
    	writeElementalData(freedom_lists_file, elemmat_rhs_file, coords_file, 
    			freemap_file, freedom_mask_file, info_file, control_file);
    	/* run */ 
    	runExecutable(DOUG_MAIN_EXECUTABLE , EXE_MAIN);
    	/* return result */
    	DoubleVector solution = parseSolution();
    	return solution;
    }

	/**
	 * Checks, if <code>executable</code> is present, then executes 
	 * <code>cmdLine</code>. Blocks until process is finished.
	 * 
	 * @param executable is checked, if it can be read
	 * @param cmdLine is executed
	 */
	private void runExecutable(String executable, String[] cmdLine) {
		/* test if executable is present */
    	File exe = new File(WORKING_DIR + executable);
    	if (exe.canRead() == false) {
    		//TODO: Fault
    		System.out.println("Executable not readable!");
    	}
    	/* run */
    	try {
	    	Process pro = Runtime.getRuntime().exec(cmdLine, null, new File(WORKING_DIR));
			pro.waitFor();
    	} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	} catch (InterruptedException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	}
	}
    
    /**
     * Writes a file in a DataHandler onto the filesystem.
     * @param dh DataHandler containing the file to be written
     * @param filen Filename of the resulting file
     * @throws IOException If either the read from the DataHandler or the write
     * 		to the filesystem fails.
     */
    private void writeFile(DataHandler dh, String filen) throws IOException  {
		InputStream in;
		FileOutputStream out;
		in = dh.getDataSource().getInputStream();
		out = new FileOutputStream(filen);
		IOUtils.copy(in, out);
		out.flush();
		out.close();
	}


	/**
	 * save files to working dir
	 * 
	 * @param freedom_lists_file
	 * @param elemmat_rhs_file
	 * @param coords_file
	 * @param freemap_file
	 * @param freedom_mask_file
	 * @param info_file
	 * @param control_file
	 */
	private void writeElementalData(DataHandler freedom_lists_file, 
			DataHandler elemmat_rhs_file, DataHandler coords_file, 
			DataHandler freemap_file, DataHandler freedom_mask_file, 
			DataHandler info_file, DataHandler control_file) {
    	try  {
    		writeFile(freedom_lists_file, WORKING_DIR + Settings.FREEDOM_LISTS_FILE);
    		writeFile(elemmat_rhs_file, WORKING_DIR + Settings.ELEMENT_RHS_FILE);
    		writeFile(coords_file, WORKING_DIR + Settings.COORDS_FILE);
    		writeFile(freemap_file, WORKING_DIR + Settings.FREEMAP_FILE);
    		writeFile(freedom_mask_file, WORKING_DIR + Settings.FREEDOM_MASK_FILE);
    		writeFile(info_file, WORKING_DIR + Settings.INFO_FILE);
    		writeFile(control_file, WORKING_DIR + Settings.CONTROL_FILE);
    	} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	}
	}
	

	/**
	 * @param solution
	 * @return
	 */
	private DoubleVector parseSolution() {
		DoubleVector solution = null;
		try {
			solution = DoubleVector.readFromDisk(WORKING_DIR + Settings.SOLUTION_FILE);
		} catch (IOException e) {
    		//TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
		}
		return solution;
	}
	
    /*
     * for testing only
     */
    /*
    public static void main(String[] args) {
        String[] result;
        DougService s = new DougService();
        result = s.runDoug();
        System.out.println(result[0]);
        System.out.println(result[1]);
    } 
    */
}
