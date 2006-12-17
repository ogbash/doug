package ee.ut.math.doug;

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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.rmi.RemoteException;
import java.util.LinkedList;
import java.util.List;

import javax.activation.DataHandler;
import javax.activation.FileDataSource;
import javax.xml.namespace.QName;
import javax.xml.rpc.ParameterMode;
import javax.xml.rpc.ServiceException;

import org.apache.axis.AxisFault;
import org.apache.axis.attachments.PlainTextDataSource;
import org.apache.axis.client.Call;
import org.apache.axis.client.Service;
import org.apache.axis.encoding.ser.JAFDataHandlerDeserializerFactory;
import org.apache.axis.encoding.ser.JAFDataHandlerSerializerFactory;
import org.apache.axis.utils.Options;

/**
 * Client for the DOUG webservice
 * 
 * @author Christian Poecher
 */
public class DougWSClient {

	private String makeControlFileForSolving() {
		StringBuffer buf = new StringBuffer();
		buf.append("solver 2\n");
		buf.append("solve_maxiters 3000\n");
		buf.append("method 1\n");
		buf.append("levels  2\n");
		buf.append("overlap 2\n");
		buf.append("smoothers 2\n");
		buf.append("input_type 2\n");
		buf.append("symmstruct T\n");
		buf.append("symmnumeric T\n");
		buf.append("radius1 2\n");
		buf.append("strong1 0.67e0\n");
		buf.append("minasize1 4\n");
		buf.append("radius2 5\n");
		buf.append("strong2 0.67e0\n");
		buf.append("minasize2 3\n");
		buf.append("matrix_type 1\n");
		buf.append("number_of_blocks 1\n");
		buf.append("initial_guess 2\n");
		buf.append("solve_tolerance 1.0e-6\n");
		buf.append("solution_format 2\n");
		buf.append("solution_file " + Settings.SOLUTION_FILE +"\n");
		buf.append("debug 0\n");
		buf.append("verbose 1\n");
		buf.append("plotting 3\n");
		buf.append("assembled_mtx_file " + Settings.ASSEMBLED_MTX_FILE + "\n");
		return buf.toString();
	}

	/*
	 * Suboptimal design choice: If new control words are added, both DOUG and
	 * client code has to be changed. Still reasonalble, because DOUG will
	 * become a library at some point and control words should then be passed as
	 * arguments to a function call.
	 */
	/**
	 * Makes a String that can be used as control file for converting elemental
	 * input into an AssembledMatrix with DOUG.
	 */
	private String makeControlFileForConverting(boolean freedomListFile,
			boolean elementRHSFile, boolean coordsFile, boolean freemapFile,
			boolean freedomMaskFile) {
		StringBuffer buf = new StringBuffer();
		buf.append("solver 2\n");
		buf.append("method 1\n");
		buf.append("input_type 1\n");
		buf.append("matrix_type 1\n");
		buf.append("info_file " + Settings.INFO_FILE + "\n");
		if (freedomListFile)
			buf.append("freedom_lists_file " + Settings.FREEDOM_LISTS_FILE + "\n");
		if (elementRHSFile)
			buf.append("elemmat_rhs_file " + Settings.ELEMENT_RHS_FILE + "\n");
		if (coordsFile)
			buf.append("coords_file " + Settings.COORDS_FILE + "\n");
		if (freemapFile)
			buf.append("freemap_file " + Settings.FREEMAP_FILE + "\n");
		if (freedomMaskFile)
			buf.append("freedom_mask_file " + Settings.FREEDOM_MASK_FILE + "\n");
		buf.append("number_of_blocks 1\n");
		buf.append("initial_guess 2\n");
		buf.append("solve_tolerance 1.0e-12\n");
		buf.append("debug 0\n");
		buf.append("verbose 10\n");
		buf.append("plotting 0\n");
		buf.append("dump_matrix_only true\n");
		buf.append("dump_matrix_file " + Settings.DUMP_MATRIX_FILE + "\n");
		return buf.toString();
	}

	/**
	 * Runs the assembled version of DOUG to solve the linear equation a x = b.
	 * 
	 * @param a Matrix
	 * @param b right hand side vector
	 * @return x = a^-1 b
	 * @throws IOException
	 * @throws DougServiceException
	 */
	public DoubleVector runAssebled(AssembledMatrix a, DoubleVector b)
		throws IOException, DougServiceException  {
		DataHandler mtxHandler = 
			new DataHandler(new PlainTextDataSource(null, a.toString()));
		DataHandler rhsHandler =
			new DataHandler(new PlainTextDataSource(null, b.toString()));
		
		
		return null;
	}
	
	/**
	 * Converts elemental input into an AssembledMatrix by calling DOUG.
	 * Parameters are local filename. They will be replaced on the server side
	 * by default names.
	 * 
	 * @param freedom_lists_file
	 *            local filename or null, if not needed
	 * @param elemmat_rhs_file
	 *            local filename or null, if not needed
	 * @param coords_file
	 *            local filename or null, if not needed
	 * @param freemap_file
	 *            local filename or null, if not needed
	 * @param freedom_mask_file
	 *            local filename or null, if not needed
	 * @return corresponding assembled matrix object
	 */
	public AssembledMatrix elementalToAssembled(File freedom_lists_file,
			File elemmat_rhs_file, File coords_file, File freemap_file,
			File freedom_mask_file, File info_file) throws IOException,
			 DougServiceException {
		
		/* prepare attachments */
		boolean flf = (freedom_lists_file != null);
		boolean erf = (elemmat_rhs_file != null);
		boolean cf = (coords_file != null);
		boolean ff = (freemap_file != null);
		boolean fmf = (freedom_mask_file != null);

		List files = new LinkedList();  //eats up possible nulls
		files.add(freedom_lists_file);
		files.add(elemmat_rhs_file);
		files.add(coords_file);
		files.add(freemap_file);
		files.add(freedom_mask_file);
		files.add(info_file);
		String control = makeControlFileForConverting(flf, erf, cf, ff, fmf);

		DataHandler[] attachments = new DataHandler[files.size() + 1];
		for (int i = 0; i < files.size(); i++) {
			attachments[i] = new DataHandler(new FileDataSource((File) files.get(i)));
		}
		// XXX keep control file in memory instead of file. absolutly not need
		// to write it.
		attachments[files.size()] = 
			new DataHandler(new PlainTextDataSource(null ,control));
//		attachments[files.size()] = new DataHandler(new FileDataSource(
//				Settings.CONTROL_FILE));

		/* prepare call */
		Service service = new Service();
		Call call;
		try {
			call = (Call) service.createCall();
		} catch (ServiceException e) {
			throw new DougServiceException(e.getMessage(), e);
		}
		// TODO: Hardwired for now, change, use Option
		call.setTargetEndpointAddress(Settings.ENDPOINT_ADRESS); 
		call.setOperationName(new QName("urn:DougService",
				"elementalToAssembled")); // This is the target services
											// method to invoke.
		QName qnameAttachment = new QName("urn:DougService", "DataHandler");
		call.registerTypeMapping(
				attachments[0].getClass(), // Add serializer for attachment.
				qnameAttachment, JAFDataHandlerSerializerFactory.class,
				JAFDataHandlerDeserializerFactory.class);
		call.addParameter("freedom_lists_file", qnameAttachment,
				ParameterMode.IN);
		call.addParameter("elemmat_rhs_file", qnameAttachment,
				ParameterMode.IN);
		call.addParameter("coords_file", qnameAttachment, ParameterMode.IN);
		call.addParameter("freemap_file", qnameAttachment, ParameterMode.IN);
		call.addParameter("freedom_mask_file", qnameAttachment,
				ParameterMode.IN);
		call.addParameter("info_file", qnameAttachment, ParameterMode.IN);
		call.addParameter("control_file", qnameAttachment, ParameterMode.IN);
		call.setReturnType(qnameAttachment);
		// call.setProperty(Call.ATTACHMENT_ENCAPSULATION_FORMAT,
		// 		   Call.ATTACHMENT_ENCAPSULATION_FORMAT_DIME);

		/* invoke call */
		Object ret = call.invoke(attachments);
		
		/* sanity check return data */
		if (null == ret) {
			System.out.println("Received null ");
			throw new AxisFault("", "Received null", null, null);
		}
		if (ret instanceof String) {
			System.out.println("Received problem response from server: " + ret);
			throw new AxisFault("", (String) ret, null, null);
		}
		if (!(ret instanceof DataHandler)) {
			// The wrong type of object that what was expected.
			System.out.println("Received problem response from server:"
					+ ret.getClass().getName());
			throw new AxisFault("", "Received problem response from server:"
					+ ret.getClass().getName(), null, null);
		}

		/* process return data */
		DataHandler rdh = (DataHandler) ret;
		InputStream is = rdh.getDataSource().getInputStream();
		InputStreamReader isr = new InputStreamReader(is);
		AssembledMatrix m = AssembledMatrix.readFromReader(isr);
		return m;
	}

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
			System.err.println("Error in server or while parsing result occurred.");
			e.printStackTrace();
		}
	}
}
