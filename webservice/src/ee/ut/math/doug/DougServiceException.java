/**
 * 
 */
package ee.ut.math.doug;

/**
 * This, or a subclass of this class is thrown, if something 
 * within the application fails. 
 * 
 * @author Christian Pöcher
 */
public class DougServiceException extends Exception {

	/**
	 * 
	 */
	public DougServiceException() {
		super();
	}

	/**
	 * @param message
	 */
	public DougServiceException(String message) {
		super(message);
	}

	/**
	 * @param cause
	 */
	public DougServiceException(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message
	 * @param cause
	 */
	public DougServiceException(String message, Throwable cause) {
		super(message, cause);
	}

}
