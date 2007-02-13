package ee.ut.math.doug;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class DougWsGui extends JPanel implements ActionListener {

	/**
	 * Validation of user entered data failed.
	 */
	class ValidationException extends Exception {

		public ValidationException() {
			super();
		}

		public ValidationException(String arg0, Throwable arg1) {
			super(arg0, arg1);
		}

		public ValidationException(String arg0) {
			super(arg0);
		}

		public ValidationException(Throwable arg0) {
			super(arg0);
		}
	}
	
	/**
	 *  A file chooser has been canceled.
	 */
	class CanceledException extends Exception {

		public CanceledException() {
			super();
		}

		public CanceledException(String arg0, Throwable arg1) {
			super(arg0, arg1);
		}

		public CanceledException(String arg0) {
			super(arg0);
		}

		public CanceledException(Throwable arg0) {
			super(arg0);
		}
	}
	
	private JTextField txtShift, txtError, txtMatrix, txtGuess;
	
	private DougWsGui() {
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		JLabel lblHeader = new JLabel("Please select data to calulate an Eigenvalue and -vector:");
        add(lblHeader);
        JLabel lblMatrix = new JLabel("Matrix in Assembled Form");
        add(lblMatrix);
        JPanel pMatrix = new JPanel();
        add(pMatrix);
        JLabel lblGuess = new JLabel("Initial Guess of Eigenvector");
        add(lblGuess);
        JPanel pGuess = new JPanel();
        add(pGuess);
        JLabel lblShift = new JLabel("Shift");
        add(lblShift);
        txtShift = new JTextField("1.0");
        add(txtShift);
        JLabel lblError = new JLabel("Error Boundary");
        add(lblError);
        txtError = new JTextField("0.000001");
        add(txtError);
        JButton btnOK = new JButton("OK");
        btnOK.setActionCommand("ok");
        btnOK.addActionListener(this);
        add(btnOK);

        pMatrix.setLayout(new FlowLayout());
        txtMatrix = new JTextField("~/assembled.txt");
        pMatrix.add(txtMatrix);
        JButton btnMatrix = new JButton("Choose file");
        btnMatrix.setActionCommand("matrix");
        btnMatrix.addActionListener(this);
        pMatrix.add(btnMatrix);
        
        pGuess.setLayout(new FlowLayout());
        txtGuess = new JTextField("~/initialGuess.txt");
        pGuess.add(txtGuess);
        JButton btnGuess = new JButton("Choose file");
        btnGuess.setActionCommand("guess");
        btnGuess.addActionListener(this);
        pGuess.add(btnGuess);
	}
	
    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() {
        //Make sure we have nice window decorations.
        JFrame.setDefaultLookAndFeelDecorated(true);

        //Create and set up the window.
        JFrame frame = new JFrame("Calculate Eigenvalue and -vector");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        //Create and set up the content pane.
        JComponent newContentPane = new DougWsGui();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);
        
        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
    public void actionPerformed(ActionEvent e) {
    	if ("matrix".equals(e.getActionCommand())) {
    		try {
    			txtMatrix.setText(openFileChooser());
    		} catch (CanceledException ex) {
    			//noop
    		}
    	} else if ("guess".equals(e.getActionCommand())) {
    		try {
    			txtGuess.setText(openFileChooser());
    		} catch (CanceledException ex) {
    			//noop
    		}
    	} else if ("ok".equals(e.getActionCommand())) {
    		try {
    			validateData();
    			execute();
    		} catch (ValidationException ex) {
    			printErrorMsg(ex.getMessage());
    		}
    	} else {
    		System.err.println("invalid ActionCommand recieved");
    	}
    }
    
    /**
     * Check if data entered makes sense.
     * 
     * @throws ValidationException if and only if data is not valid
     */
    private void validateData() throws ValidationException {
		File matrix = new File(txtMatrix.getText());
		if ( ! matrix.canRead() )
			throw new ValidationException("Cannot read matrix file. Is path correct?");
		File guess = new File(txtGuess.getText());
		if ( ! guess.canRead() )
			throw new ValidationException("Cannot read initial guess file. Is path correct?");
		try {
			Double.parseDouble(txtShift.getText());
		} catch (NumberFormatException e) {
			throw new ValidationException("Cannot parse real value expected for shift.");
		}
		try {
			Double.parseDouble(txtError.getText());
		} catch (NumberFormatException e) {
			throw new ValidationException("Cannot parse real value expected for error boundary.");
		}
	}

	/**
     * Starts the calculation.
     */
    private void execute() {
		
	}

    /**
     * Prints an error message.
     * @param message
     */
	private void printErrorMsg(String message) {
    	JOptionPane.showMessageDialog(this,
    		    message, "Error", JOptionPane.ERROR_MESSAGE);
	}

	/**
     * Opens a modal file chooser 
     * @return the chosen file
     */
    private String openFileChooser() throws CanceledException {
        JFileChooser chooser = new JFileChooser();
        int returnVal = chooser.showOpenDialog(this);
        if(returnVal == JFileChooser.APPROVE_OPTION) {
        	try {
        		return chooser.getSelectedFile().getCanonicalPath();
        	} catch (IOException e) { //weired error occurred, should not happen
        		e.printStackTrace();
        		System.exit(-1);
        		return null; //never reached, just that the compiler is happy
        	}
        } else {
        	throw new CanceledException();
        }
	}

	public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }

}
