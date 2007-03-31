/**
 * 
 */
package ee.ut.math.doug;

/**
 * @author Christian Pöcher
 *
 */
public class AssembledTextTester {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		DougWSClient client = new DougWSClient();
		AssembledMatrix m = AssembledMatrix.readFromDisk(args[0]);
		DoubleVector b = DoubleVector.readFromDisk(args[1]);
		DoubleVector x = client.runAssebled(m, b);
		x.writeToDisk(args[2]);
	}

}
