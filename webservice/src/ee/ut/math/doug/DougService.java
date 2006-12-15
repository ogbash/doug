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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import javax.activation.DataHandler;
import javax.activation.FileDataSource;

import org.apache.commons.io.IOUtils;

/**
 * Webservice, that starts DOUG synchronously and gives back stdout and stderr.
 * 
 * @author Christian Poecher
 */
public class DougService {

    /* Dir and executable name fixed by producing WS for security reasons. */
    /* TODO: make editable properties file */
    private static final String WORKING_DIR = "/home/poecher/working/"; //WebDAV
    private static final String NO_OF_PROCESSES = "4";
    private static final String DOUG_MAIN_EXECUTABLE = "DOUG_main";
    private static final String[] EXE_CONV = {"mpirun", "-np", "1", DOUG_MAIN_EXECUTABLE, "-q"}; // TODO: better run DOUG without mpirun, but cannot find file. why?
    private static final String[] EXE_MAIN = {"mpirun", "-np", NO_OF_PROCESSES, "DOUG_main", "-q"};
    private static final String[] EXE_AGGR = {"mpirun", "-np", NO_OF_PROCESSES, "DOUG_aggr", "-q"};
    private static final int CHUNKSIZE = 4096;

    /**
     * Starts DOUG synchronously and returns stdout and stderr. 
     * 
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
    
    public double[] runDoug1() {
    		
    		Process pro;
    		double[] solution = null;
    		
    		try {
    			pro = Runtime.getRuntime().exec(EXE_MAIN, null, new File(WORKING_DIR));
    			pro.waitFor();
    			solution = parseSolutionFile();
    			
    		} catch (IOException ex) {
             ex.printStackTrace();
	    } catch (InterruptedException ex) {
			ex.printStackTrace();
	    }
    		return solution;
    }
    
    /* 
     * retruns DataHandler and not AssembledMatrix, because of less 
     * interwoven code. making sense of the file is up to the client.
     * service is only providing the data.
     */
    public DataHandler elementalToAssembled(DataHandler freedom_lists_file,
    		DataHandler elemmat_rhs_file, DataHandler coords_file,
    		DataHandler freemap_file, DataHandler freedom_mask_file,
    		DataHandler info_file, DataHandler control_file) {
    	
    	/* save files to working dir */
    	try  {
    		writeFile(freedom_lists_file, WORKING_DIR + Settings.FREEDOM_LISTS_FILE);
    		writeFile(elemmat_rhs_file, WORKING_DIR + Settings.ELEMENT_RHS_FILE);
    		writeFile(coords_file, WORKING_DIR + Settings.COORDS_FILE);
    		writeFile(freemap_file, WORKING_DIR + Settings.FREEMAP_FILE);
    		writeFile(freedom_mask_file, WORKING_DIR + Settings.FREEDOM_MASK_FILE);
    		writeFile(info_file, WORKING_DIR + Settings.INFO_FILE);
    		writeFile(control_file, WORKING_DIR + Settings.CONTROL_FILE);
    	} catch (IOException e) {
//    		TODO: Fault
    		System.out.println(e.getMessage());
    		e.printStackTrace();
    	}
    	/* test if executable is present */
    	File exe = new File(WORKING_DIR + DOUG_MAIN_EXECUTABLE);
    	if (exe.canRead() == false) {
    		//TODO: Fault
    		System.out.println("Executable not readable!");
    	}
    	/* run */
    	try {
	    	Process pro = Runtime.getRuntime().exec(EXE_CONV, null, new File(WORKING_DIR));
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
		/* return result */
    	FileDataSource fds = new FileDataSource(WORKING_DIR + Settings.DUMP_MATRIX_FILE);
    	DataHandler result = new DataHandler(fds);
    	return result;
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

	private double[] parseSolutionFile() throws FileNotFoundException, IOException {
		
		File file = new File(WORKING_DIR + Settings.SOLUTION_FILE);
		FileReader reader = new FileReader(file);
		BufferedReader buf = new BufferedReader(reader);
		String s = buf.readLine();
		s = s.trim();
		int numOfElem = Integer.parseInt(s);
		double[] solution = new double[numOfElem];
		for (int i=0; i<numOfElem; i++) {
			s = buf.readLine();
			s = s.trim();
			solution[i] = Double.parseDouble(s);
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
