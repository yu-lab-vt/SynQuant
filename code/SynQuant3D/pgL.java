import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.util.FastMath;


public class pgL implements MultivariateFunction {
	public double [] yest;
	public double [] sigmaest;
	public double [] ci;
	public double [] di;
	public pgL(double[] xMean, double[] xStd, double [] xCi, double[] xDi){
		yest = xMean;
		sigmaest = xStd;
		ci = xCi;
		di = xDi;
	}
	
	public double value(double[] pCur) {
		final double a = pCur[0];
		final double b = pCur[1];
		int M = 10000;           // integration step
		int N = yest.length;   // number of level sets
		BasicMath mBM = new BasicMath();
		double [] yrg = new double[M+1];
		double [] sigma2reg = new double[M+1];
		for(int i=0;i<M+1;i++){
			yrg[i] = (double)i/M;
			sigma2reg[i] = Math.max(yrg[i]*a+b, 0.000001);
		}
		double res = 0;
		for(int ii=0;ii<N;ii++){
		    double[] k = new double[sigma2reg.length];
		    double[] tc = new double[yrg.length];
		    double[] td = new double[sigma2reg.length];
		    double[] t = new double[sigma2reg.length];
		    double[] f_val = new double[sigma2reg.length];
		    for(int jj = 0;jj<sigma2reg.length;jj++){
		    	k[jj] = 1/(2 * Math.PI * Math.sqrt(ci[ii]*di[ii]) * sigma2reg[jj]);     
		    	tc[jj] = (yest[ii]-yrg[jj])*(yest[ii]-yrg[jj])/ci[ii];
		    	td[jj] = (sigmaest[ii]-Math.sqrt(sigma2reg[jj]))*(sigmaest[ii]-Math.sqrt(sigma2reg[jj]))/di[ii];
		    	t[jj] = -1/(2*sigma2reg[jj])*(tc[jj]+td[jj]);
		    	f_val[jj] = k[jj]*Math.exp(t[jj]);
		    	
		    }
		    res = res + trapz(yrg[0],yrg[yrg.length-1],yrg.length-1,f_val);    
		}
		res = -res;
		//System.out.printf("%.3f\n", res);
		return res;
	}

	private double trapz(double a, double b, int N, double[] f) {
		// TODO Auto-generated method stub
	      double h = (b - a) / N;              // step size
	      double sum = 0.5 * (f[0] + f[N]);    // N=f.length-1
	      for (int i = 1; i < N; i++) {
	         sum = sum + f[i];
	      }

	      return sum * h;
	}
}