import java.net.MalformedURLException;
import java.rmi.RemoteException;

import javax.xml.rpc.ServiceException;

import org.apache.axis.client.Call;
import org.apache.axis.client.Service;
import org.apache.axis.utils.Options;

/**
 * Client for the DOUG webservice
 * 
 * @author Christian Poecher
 */
public class DougWSClient {

    /**
     * Calls Webservice "doug" on specified host and port and puts DOUG's stdout
     * and stderr to the client machine's stdout and stderr.
     * 
     * @param args
     *            optional parameters to specify service URI are: -h <host> -p
     *            <port>
     */
    public static void main(String[] args) {
        try {
            Options options = new Options(args);
            String endpoint = "http://" + options.getHost() + ":"
                    + options.getPort() + "/axis/services/doug";

            System.out.println("Calling Service on endpoint " + endpoint + ".");
            Service service = new Service();
            Call call = (Call) service.createCall();

            call.setTargetEndpointAddress(new java.net.URL(endpoint));
            call.setOperationName("runDoug");
            call.setReturnClass(String[].class);

            /* Call without parameters */
            String[] output = (String[]) call.invoke(new Object[] {});

            System.out.println(output[0]); // stdout -> stdout
            System.err.println(output[1]); // stderr -> stderr
        } catch (MalformedURLException e) {
            System.err.println("Either host or port parameter are malformed.");
            e.printStackTrace();
        } catch (ServiceException e) {
            System.err.println("Error while preparing service object.");
            e.printStackTrace();
        } catch (RemoteException e) {
            System.err
                    .println("Error in server or while parsing result occurred.");
            e.printStackTrace();
        }
    }
}
