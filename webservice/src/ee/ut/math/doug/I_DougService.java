package ee.ut.math.doug;

import javax.activation.DataHandler;

public interface I_DougService {

	public DoubleVector runAssembled(AssembledMatrix matrix, DoubleVector rhs,
			DataHandler control_file);
	
	public DataHandler runAssembled(DataHandler matrix, DataHandler rhs,
			DataHandler controlfile);

	public DoubleVector runElemental(DataHandler freedom_lists_file,
			DataHandler elemmat_rhs_file, DataHandler coords_file,
			DataHandler freemap_file, DataHandler freedom_mask_file,
			DataHandler info_file, DataHandler control_file);
	
	public AssembledMatrix elementalToAssembled(DataHandler freedom_lists_file,
			DataHandler elemmat_rhs_file, DataHandler coords_file,
			DataHandler freemap_file, DataHandler freedom_mask_file,
			DataHandler info_file, DataHandler control_file);

}