/**
 * 
 */
package ee.ut.math.doug;

import java.io.File;

/**
 * @author Christian Pöcher
 *
 */
public class AssembledXDRTester {

	/**
	 * @param args matrix rhs solution
	 */
	public static void main(String[] args) throws Exception {
		DougWSClient client = new DougWSClient();
		File matrix   = new File(args[0]);
		File rhs      = new File(args[1]);
		File solution = new File(args[2]);
		client.runAssebled(matrix, rhs, solution);
	}

}
