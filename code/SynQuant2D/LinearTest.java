import java.util.ArrayList;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.apache.commons.math3.util.FastMath;

public class LinearTest {
	public int height, width;
	public int[][] lbSyn;
	public int SynNum;
	public int max_dis=20;//max distance that synapse to dendrite
	public double[] denSynNum;
	public double[] denSynNumLog;
	public double[][] DenFeatures;
	public double[] betas;
	public double[] pvalues;
	public LinearTest(boolean[][] kSynR1, GrowNeurite den){
		height = kSynR1.length;
		width = kSynR1[0].length;
		
		ImageHandling IH = new ImageHandling();
		lbSyn = IH.bwlabel(kSynR1,8);
		SynNum = IH.NextLabel;
		//ArrayList<Integer[]> syn_pos = new ArrayList<Integer[]>();
		int[][] syn_pos = new int[SynNum][4];
		int[] syn_size = new int[SynNum];
		for(int i=0;i<SynNum;i++){
			syn_pos[i][0] = height;
			syn_pos[i][1] = width;
		}
		for(int y=0;y<height;y++){
			for(int x=0;x<width;x++){
				if(lbSyn[y][x]>0){
					if(syn_pos[lbSyn[y][x]-1][0]>y)
						syn_pos[lbSyn[y][x]-1][0]=y;
					if(syn_pos[lbSyn[y][x]-1][1]>x)
						syn_pos[lbSyn[y][x]-1][1]=x;
					if(syn_pos[lbSyn[y][x]-1][2]<y)
						syn_pos[lbSyn[y][x]-1][2]=y;
					if(syn_pos[lbSyn[y][x]-1][3]<x)
						syn_pos[lbSyn[y][x]-1][3]=x;

					syn_size[lbSyn[y][x]-1]++;
				}
			}
		}
		DenFeatures = new double[den.CurNum][3];//dendrite length, scale, intensity
		denSynNum = new double[den.CurNum]; //synapse number on it
		for(int i=0;i<SynNum;i++){
			//syn_pos[i][0]=(int) Math.round((double)syn_pos[i][0]/syn_size[i]);
			//syn_pos[i][1]=(int) Math.round((double)syn_pos[i][1]/syn_size[i]);
			int y = (int) Math.round((double)(syn_pos[i][0]+syn_pos[i][2])/2);
			int x = (int) Math.round((double)(syn_pos[i][1]+syn_pos[i][3])/2);
			int lu_x = Math.max(x-max_dis, 0);
			int lu_y = Math.max(y-max_dis, 0);
			int rd_x = Math.min(x+max_dis, width-1);
			int rd_y = Math.min(y+max_dis, height-1);
			int MinCurLabel = 0;
			double MinCurDis = 10000;
			for(int iy = lu_y;iy<rd_y+1;iy++){
				for(int ix = lu_x;ix<rd_x+1;ix++){
					if(den.CurIm[iy][ix]>0){
						double tmpDis = Math.sqrt((double)((ix-x+1)*(ix-x+1)+(iy-y+1)*(iy-y+1)));
						if(tmpDis<MinCurDis){
							MinCurLabel = den.CurIm[iy][ix];
							MinCurDis = tmpDis;
						}
					}
				}
			}
			if(MinCurLabel>0)
				denSynNum[MinCurLabel-1]+=1;// always binary
		}
		double[] denL = zscore(den.denLength);
		double[] denS = zscore(den.denSize);
		double[] denI = zscore(den.denIntensity);
		for(int i=0;i<den.CurNum;i++){
			DenFeatures[i][0] = denL[i];
			DenFeatures[i][1] = denS[i];
			DenFeatures[i][2] = denI[i];
			//denSynNumLog[i] = Math.log10(denSynNum[i]);//for poisson regression
		}
		//linear regression test
		/*OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
		regression.newSampleData(denSynNum, DenFeatures);
		betas = regression.estimateRegressionParameters();
		double[] standardErrors = regression.estimateRegressionParametersStandardErrors();
		int  residualdf = regression.estimateResiduals().length - betas.length;
		TDistribution tdistribution = new TDistribution(residualdf);
		pvalues = new double[betas.length];
		  //calculate p-value and create coefficient
		  for (int i = 0; i < betas.length; i++)
		  {
		     double tstat = betas[i] / standardErrors[i];
		     pvalues[i] = tdistribution.cumulativeProbability(-FastMath.abs(tstat)) * 2;
		  }*/
		  
		// poisson regression test
		PoissonRegression PR = new PoissonRegression(DenFeatures,denSynNum, 1e-6,0.0001);
		betas = new double[PR.betas.getDimension()];
		for(int i=0;i<betas.length;i++)
			betas[i] = PR.betas.getEntry(i);
		
		//model.fit(regularization);
		for(int i=0;i<den.CurNum;i++){
			DenFeatures[i][0] = den.denLength[i];
			DenFeatures[i][1] = den.denSize[i];
			DenFeatures[i][2] = den.denIntensity[i];
			//denSynNumLog[i] = Math.log10(denSynNum[i]);//for poisson regression
		}
	}
	double[] zscore(double[] invec){
		double [] outvec = new double[invec.length];
		BasicMath BM = new BasicMath();
		double meanVal = BM.sampleMean(invec);
		double var = BM.sampleVar(invec);
		for(int i=0;i<invec.length;i++)
			outvec[i] = (invec[i]-meanVal)/Math.sqrt(var);

		return outvec;
	}
}
