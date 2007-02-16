/**
 * 
 */
package ee.ut.math.doug;

import java.io.IOException;

import javax.xml.namespace.QName;
import javax.xml.rpc.ParameterMode;
import javax.xml.rpc.ServiceException;

import org.apache.axis.AxisFault;
import org.apache.axis.client.Call;
import org.apache.axis.client.Service;
import org.apache.axis.encoding.ser.BeanDeserializerFactory;
import org.apache.axis.encoding.ser.BeanSerializerFactory;

/**
 * @author Christian P�cher
 * 
 */
public class IterativeEigenvalueClient {
	public EigenSpace runSolver(AssembledMatrix matrix,
			DoubleVector initialGuess, double shift, double error)
			throws IterativeEVServiceException, IOException {
		/* prepare call */
		Service service = new Service();
		Call call;
		try {
			call = (Call) service.createCall();
		} catch (ServiceException e) {
			throw new IterativeEVServiceException(e.getMessage(), e);
		}
		// TODO: Hardwired for now, change, use Option
		call.setTargetEndpointAddress(Settings.ENDPOINT_ADRESS_EVS);
		call.setOperationName(new QName(Settings.NAMESPACE_ID,
				"inverseIteration")); // This is the target service's method
										// to invoke.

		// Add (de-)serializer for AssebledMatrix.
		QName qnameAssembledMatrix = new QName(Settings.NAMESPACE_ID,
				"AssembledMatrix");
		call.registerTypeMapping(AssembledMatrix.class, qnameAssembledMatrix,
				new BeanSerializerFactory(AssembledMatrix.class,
						qnameAssembledMatrix), new BeanDeserializerFactory(
						AssembledMatrix.class, qnameAssembledMatrix));
		// Add (de-)serializer for DoubleVector.
		QName qnameDoubleVector = new QName(Settings.NAMESPACE_ID,
				"DoubleVector");
		call
				.registerTypeMapping(DoubleVector.class, qnameDoubleVector,
						new BeanSerializerFactory(DoubleVector.class,
								qnameDoubleVector),
						new BeanDeserializerFactory(DoubleVector.class,
								qnameDoubleVector));
		// Add (de-)serializer for EigenSpace.
		QName qnameEigenSpace = new QName(Settings.NAMESPACE_ID, "EigenSpace");
		call.registerTypeMapping(EigenSpace.class, qnameEigenSpace,
				new BeanSerializerFactory(EigenSpace.class, qnameEigenSpace),
				new BeanDeserializerFactory(EigenSpace.class, qnameEigenSpace));

		call.addParameter("a", qnameAssembledMatrix, ParameterMode.IN);
		call.addParameter("initialGuess", qnameDoubleVector, ParameterMode.IN);
		call.addParameter("shift", org.apache.axis.Constants.SOAP_DOUBLE,
				ParameterMode.IN);
		call.addParameter("error", org.apache.axis.Constants.SOAP_DOUBLE,
				ParameterMode.IN);
		call.setReturnType(qnameEigenSpace);

		/* invoke call */
		Object ret = call.invoke(new Object[] { matrix, initialGuess,
				new Double(shift), new Double(error) });

		/* sanity check return data */
		if (null == ret) {
			System.out.println("Received null ");
			throw new AxisFault("", "Received null", null, null);
		}
		if (ret instanceof String) {
			System.out.println("Received problem response from server: " + ret);
			throw new AxisFault("", (String) ret, null, null);
		}
		if (!(ret instanceof EigenSpace)) {
			// The wrong type of object that what was expected.
			System.out.println("Received problem response from server:"
					+ ret.getClass().getName());
			throw new AxisFault("", "Received problem response from server:"
					+ ret.getClass().getName(), null, null);
		}

		/* process return data */
		EigenSpace es = (EigenSpace) ret;
		return es;
	}
}
