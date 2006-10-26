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
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Webservice, that starts DOUG synchronously and gives back stdout and stderr.
 * 
 * @author Christian Poecher
 */
public class DougService {

    /* Dir and executable name fixed by producing WS for security reasons. */
    /* TODO: make editable properties file */
    private static final String WORKING_DIR = "/home/poecher/apache-tomcat-5.5.17/webapps/ROOT"; //WebDAV
    private static final String NO_OF_PROCESSES = "4";
    private static final String[] EXECUTABLE = {"mpirun", "-np", NO_OF_PROCESSES, "DOUG_main"};
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
            pro = Runtime.getRuntime().exec(EXECUTABLE, null, new File(WORKING_DIR));
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
