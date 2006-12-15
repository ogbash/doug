package ee.ut.math.doug;
import java.util.HashSet;
import java.util.Set;


/**
 * 
 */

/**
 * @author Christian P�cher
 *
 * Not necessary complete Eigenspace: Only some or even none Eigenvectors 
 * may be given.
 */
public class EigenSpace {
	
	private double eigenValue;
	private Set eigenVectors = new HashSet();
	
	public void setEigenValue(double eigenValue) {
		this.eigenValue = eigenValue;
	}
	
	public double getEigenValue() {
		return eigenValue;
	}

	public void addEigenVector(DoubleVector vector) {
		eigenVectors.add(vector);
	}
	
	public DoubleVector firstEigenVector() {
		return (DoubleVector) eigenVectors.iterator().next();
	}
	
	public Set getEigenVectors() {
		return eigenVectors;
	}
	
	public EigenSpace(double eigenvalue) {
		setEigenValue(eigenvalue);
	}
}
