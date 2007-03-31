/**
 * 
 */
package ee.ut.math.doug;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 * @author Christian Pöcher
 *
 */
public class RandomVectorGenerator {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		int n = Integer.parseInt(args[0]);
		StringBuffer buf = new StringBuffer();
		buf.append(n + "\n");
		for (int i=0; i<n; i++) {
			buf.append(Math.random() + "\n");
		}
		File f = new File("/home/poecher/testing_ass/rhs" + n + ".txt");
		FileWriter fw = new FileWriter(f); 
		BufferedWriter bw = new BufferedWriter(fw);
		bw.write(buf.toString());
		bw.close();
	}

}
