import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class PoissonRegression {
	public RealMatrix xMatrix;
	public RealVector yVector;
	public int localRows;
	public RealVector betas;
	public double step;
	public PoissonRegression(double[][] Inx_mat/*data pt->row*/, double[] y_vec, double precision, double instep) {
		double[][] x_mat = new double[Inx_mat.length][Inx_mat[0].length+1];
		for(int i=0;i<Inx_mat.length;i++){
			for(int j=0;j<Inx_mat[0].length;j++){
				x_mat[i][j] = Inx_mat[i][j];
			}
			x_mat[i][x_mat[0].length-1] = 1;
		}
		localRows = x_mat.length;
		step = instep;
		xMatrix = MatrixUtils.createRealMatrix(x_mat);
		yVector = MatrixUtils.createRealVector(y_vec);
		double[] intial_betas = new double[x_mat[0].length];
		betas = MatrixUtils.createRealVector(intial_betas);
		betas = betas.mapAdd(0.0005);
		
		double x_old = 0;
		double x_new = getValue();
		//System.out.printf("x_new: %.5f, betas: %.5f, %f.5, %.5f, Intercept: %.5f\n",x_new, betas.getEntry(0),betas.getEntry(1),betas.getEntry(2),betas.getEntry(3));
		while(Math.abs(x_old-x_new)>precision){
			//RealVector tmp = x_new.copy();
			x_old = x_new;
			betas = betas.subtract(getValueGradient().mapMultiply(step));
			x_new = getValue();
			//System.out.printf("x_new: %.5f, betas: %.5f, %f.5, %.5f, Intercept: %.5f\n",x_new, betas.getEntry(0),betas.getEntry(1),betas.getEntry(2),betas.getEntry(3));
		}
		//betas = x_new;
	}
	public double getValue() {
		double accum = 0;
		for (int row = 0; row < localRows; row++) {
			RealVector x = xMatrix.getRowVector(row);
			double y = yVector.getEntry(row);
			double mu = betas.dotProduct(x);
			accum += y * mu-Math.exp(mu); // omits "
		}
		return -accum;
	}
	public RealVector getValueGradient() {
		RealVector accum = betas.copy().mapMultiply(0);// negative sign on regularization to subtract
		for (int row = 0; row < localRows; row++) {
			RealVector x = xMatrix.getRowVector(row);
			double y = yVector.getEntry(row);
			accum = accum.add(x.mapMultiply(y - Math.exp(betas.dotProduct(x))));
		}
		return accum.mapMultiply(-1);
	}
	/*public double getValue() {
		double accum = 0;
		for(int row=0;row<localRows;row++){
			RealVector x = xMatrix.getRowVector(row);
			double y = yVector.getEntry(row);
			double mu = betas.dotProduct(x);
			
			accum +=  (-Math.exp(mu) + y * mu); //omits "- Gamma.logGamma(1+y)" as it's constant
		}
		return accum - step / 2 * betas.dotProduct(betas);
	}
	public RealVector getValueGradient() {
		RealVector accum = betas.mapMultiply(-step);//negative sign on regularization to subtract
		for (int row=0;row<localRows;row++){
			RealVector x = xMatrix.getRowVector(row);
			double y = yVector.getEntry(row);
			// if these represent population slices, need to multiply by w
			accum = accum.add(x.mapMultiply((y - Math.exp(betas.dotProduct(x)))));
		}
		return accum;
	}*/
}