import java.awt.Color;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import ij.*;
import ij.gui.*;
import ij.plugin.*;
import ij.plugin.frame.RoiManager;
import ij.process.*;

/**
 * This plugin 
 * 1. segments synapse/puncta using Probability principled synapse detection
 * algorithm modified from "PPSD: Probability Principled Synapse Detection"
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @version 1.0
 * @date 2016-05-31
 *
 *
 */

public class ppsd3D {
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp; //Current Image
	protected double[][][] G; //2D matrix saving gray scale Image data
	protected double[][][] Gt; //stabilized image
	protected int ntry=1; //decide max number of neighbors considered in zscore calculation
	protected int zSlice; //image width
	protected int width; //image width
	protected int height; //image height
	protected boolean[][][] kSynR1; //bianry synapse Map
	protected int[][][] SynR1Idx; //synapse Map
	protected int nSyn0;
	protected double fdr;//threshold for FDR control, user input
	protected ArrayList<double[]> zscoreList = new ArrayList<double[]>();
	public ppsd3D(short[][] imArray, int inwidth, int inheight, double vox_x, double infdr,boolean fastFlag,int MinSize,int MaxSize) {
		width = inwidth;
		height = inheight;
		zSlice = imArray.length;
		fdr = infdr;
		G = new double[zSlice][height][width];
		Gt = new double[zSlice][height][width];
		//parameter initialization
		paraQ3D q = new paraQ3D(imArray, vox_x,ntry,height, width);
		
		/**if we are debuging it is better to be fast**/
		boolean isDebug = java.lang.management.ManagementFactory.getRuntimeMXBean(). getInputArguments().toString().contains("-agentlib:jdwp");
	
		// put all zstack in a single frame 
		double[][] tmpG = new double[height][width*zSlice];
		for (int zz = 0; zz<zSlice; zz++) {
			int shiftPixels = zz*width;
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					tmpG[i][j+shiftPixels] = (double) imArray[zz][i * width + j] * 255 / q.maxVal;// the data is saved by line
					G[zz][i][j] = (double) imArray[zz][i * width + j] * 255 / q.maxVal;
				}
			}
		}
		//long startTime1=System.nanoTime();     // Testing running time: start
		// use the tmpG to estimate variance and do stablization
		double [][] tmpGt;
		if (!isDebug) {
			tmpGt = getNoiseVST(tmpG,  q);
		}
		else {
			tmpGt = tmpG;
			q.var = 1;
			q.pgEst = new double[2];
			q.pgEst[0] = 1;
			q.pgEst[1]= 1;
		}
		//long endTime1=System.nanoTime();
		//System.out.println("Running time: "+(endTime1-startTime1)/1e9); 
		//System.out.println("Frame Size: "+tmpGt.length+" * "+tmpGt[0].length); 
		//change the tmpGt to Gt and imArray
		for (int zz = 0; zz<zSlice; zz++) {
			int shiftPixels = zz*width;
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					Gt[zz][i][j] = tmpGt[i][j+shiftPixels];// the data is saved by line, what if we use orginal data?
					imArray[zz][i * width + j] = (short) Gt[zz][i][j];
				}
			}
		}
		//q.var = q.var/validSliceCnt;
		//for (int i=0;i<q.pgEst.length;i++)
		//	q.pgEst[i] = q.pgEst[i]/validSliceCnt;
		//System.out.println("Gt Size: "+Gt.length+" * "+Gt[0].length+" * "+Gt[0][0].length); 
		//System.out.println("imArray Size: "+imArray.length+" * "+imArray[0].length); 
		BasicMath mBM = new BasicMath();
		paraP3D p = new paraP3D(fdr,(int)mBM.matrix2DMin(imArray),(int)mBM.matrix2DMax(imArray),MinSize, MaxSize);
		//ppsd start
		boolean [][][] kMask = new boolean[zSlice][height][width];
		for (int i = 0; i < zSlice; i++)
			for (int j = 0; j < height; j++)
				for (int k = 0; k < width; k++)
					kMask[i][j][k] = true;

		if(!fastFlag){/*General version of PPSD*/
			ppsdreal3D ppsd_main = new ppsdreal3D(imArray, G, Gt, kMask, p, q);
			// Post processing: size and zscore
			kSynR1 = ppsd_main.ppsd_post3D(Gt, p, q);
			//nSyn0 = ppsd_main.nSyn0;
			SynR1Idx = ppsd_main.kMap;
			nSyn0 = ppsd_main.nSyn0;
			zscoreList = ppsd_main.zscoreList;
		}
		else{/*Fast version of PPSD(No zscore updating)*/
			// not used
			/*Fastppsdreal ppsd_main = new Fastppsdreal3D(imArray, width,height, p, q);
			kSynR1 = ppsd_main.SynR1;
			ImageHandling IH = new ImageHandling();
			SynR1Idx = IH.bwlabel(kSynR1, 8);
			nSyn0 = IH.NextLabel;*/
			
			//nSyn0 = ppsd_main.nSyn0;
		}
		//long endTime1=System.nanoTime(); 	  // Testing running time: end
		//System.out.println("Running time: "+(endTime1-startTime1)/1e9); 
	}
	//noise stabilization and variance estimation of Foi's paper. 
	//More details can be found: http://www.cs.tut.fi/~foi/sensornoise.html
	public double[][] getNoiseVST(double[][] Org_Im, paraQ3D qq) {
		double[][] double_G = new double[Org_Im.length][Org_Im[0].length];
		for (int i = 0; i < Org_Im.length; i++)
			for (int j = 0; j < Org_Im[0].length; j++)
				double_G[i][j] = Org_Im[i][j] / 255;
		
		//// Poisson Gaussian model fitting
		int edgeTau = 1;
		int nLevels = 600;
		double delta0 = 1.0 / nLevels;

		/* wavelet and scaling functions */
		double[] lp = { 0.035, 0.085, -0.135, -0.460, 0.807, -0.333 };
		double[] hp = { 0.025, -0.06, -0.095, 0.325, 0.571, 0.235 };
		double phi22 = 0.250141019881;// sum(conv2(lp.',hp).^2)
		// discrete wavelet transform
		double[][] zwdet = dwt2(double_G, lp);
		double[][] zwapp = dwt2(double_G, hp);
		
		/* smooting and edge detection */
		// filtering
		double[][] hw = new double[7][7];
		for (int i = 0; i < 7; i++)
			for (int j = 0; j < 7; j++)
				hw[i][j] = 0.0204; // 7*7 filter average filter

		double[][] hl = { { 0.5, 0, 0.5 }, { 0, -2, 0 }, { 0.5, 0, 0.5 } }; // laplasian
		double[][] hs = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };// sobel
		double[][] zsmo = imfilter(zwapp, hw);

		double[][] zwdet_piabs = new double[zwdet.length][zwdet[0].length];
		for (int i = 0; i < zwdet.length; i++)
			for (int j = 0; j < zwdet[0].length; j++)
				zwdet_piabs[i][j] = Math.sqrt(Math.PI / 2) * Math.abs(zwdet[i][j]);

		double[][] s = imfilter(zwdet_piabs, hw);
		double[][] zwappMedian = medfilter(zwapp, 3);

		double[][] zwappL = imfilter(zwappMedian, hl);
		double[][] zwappS = imfilter(zwappL, hs);
		boolean[][] xsmo = new boolean[zwappS.length][zwappS[0].length];
		for (int i = 0; i < zwappS.length; i++)
			for (int j = 0; j < zwappS[0].length; j++)
				xsmo[i][j] = (Math.abs(zwappL[i][j]) + Math.abs(zwappS[i][j])) < edgeTau * s[i][j];
		
		/* level sets */
		ArrayList<ArrayList<Integer[]>> levelSets = new ArrayList<ArrayList<Integer[]>>();
		double[][] zsmo1 = new double[zsmo.length][zsmo[0].length];
		double zsmo1Max = -10;
		for (int i = 0; i < zsmo.length; i++) {
			for (int j = 0; j < zsmo[0].length; j++) {
				if (!xsmo[i][j]) {
					zsmo1[i][j] = -1;
				} else {
					zsmo1[i][j] = zsmo[i][j];
					if (zsmo1[i][j] > zsmo1Max)
						zsmo1Max = zsmo1[i][j];
				}
			}
		}
		int nPair = 0;
		
		for (int i = 0; i < nLevels; i++) {
			double z0 = i * delta0;
			double z1 = z0 + delta0;
			ArrayList<Integer[]> tmpaylist = new ArrayList<Integer[]>();
			if (z0 > zsmo1Max){
				break;
			}
			
			for (int m = 0; m < zsmo1.length; m++) {
				for (int n = 0; n < zsmo1[0].length; n++) {
					if (zsmo1[m][n] > z0 && zsmo1[m][n] < z1) {
						tmpaylist.add(new Integer[]{m,n} );// y x
					}
				}
			}
			if (tmpaylist.size() > 20){
				nPair = nPair + 1;
			}
			levelSets.add(tmpaylist);
		}
		double[] xMean = new double[nPair];
		double[] xStd = new double[nPair];
		double[] xVar = new double[nPair];
		double[] xKai = new double[nPair];
		double[] xCi = new double[nPair];
		double[] xDi = new double[nPair];
		double[] xNi = new double[nPair];
		int n0 = 0;

		for (int i = 0; i < levelSets.size() ; i++) {
			if (levelSets.get(i).size() > 20) {
				// xMean[n0] = mean(zwapp(pix0));
				double sum_x = 0;
				double[] zwdet0 = new double[levelSets.get(i).size()];
				double sum_zwdet0 = 0;
				for (int j = 0; j < levelSets.get(i).size(); j++) {
					// zwapp[(Integer)
					sum_x = sum_x + zwapp[levelSets.get(i).get(j)[0]][levelSets.get(i).get(j)[1]];

					zwdet0[j] = zwdet[levelSets.get(i).get(j)[0]][levelSets.get(i).get(j)[1]];
					sum_zwdet0 = sum_zwdet0 + zwdet0[j];
				}
				xMean[n0] = sum_x / levelSets.get(i).size();
				double zwdet0Mean = sum_zwdet0 / levelSets.get(i).size();
				int ni = zwdet0.length;
				xNi[n0] = ni;
				double kai = 1.0 - 1.0 / 4 / ni - 7.0 / 32 / (ni * ni);
				xKai[n0] = kai;
				double zwdet0_sum2 = 0;
				for (int j = 0; j < zwdet0.length; j++) {
					zwdet0_sum2 = zwdet0_sum2 + (zwdet0[j] - zwdet0Mean) * (zwdet0[j] - zwdet0Mean);
					// xVar[n0] = zwdet0[j] - zwdet0Mean;
				}
				xVar[n0] = zwdet0_sum2 / (ni - 1);
				xStd[n0] = Math.sqrt(zwdet0_sum2 / (ni - 1)) / kai;
				xCi[n0] = phi22 / ni;
				xDi[n0] = (1 - kai * kai) / (kai * kai);
				n0 = n0 + 1;
			}
		}
		// LS Estimation use Apache
		double[][] A = new double[nPair][2];
		double[][] A_Inv = new double[2][nPair];
		for (int i = 0; i < nPair; i++) {
			A[i][0] = xMean[i];
			A_Inv[0][i] = xMean[i];
			A[i][1] = 1;
			A_Inv[1][i] = 1;
		}
		RealMatrix A_real = MatrixUtils.createRealMatrix(A);
		RealMatrix A_Inv_real = MatrixUtils.createRealMatrix(A_Inv);
		RealMatrix AtA = A_Inv_real.multiply(A_real);
		RealVector b = MatrixUtils.createRealVector(xVar);
		RealVector Atb = A_Inv_real.operate(b);
		DecompositionSolver solver = new LUDecomposition(AtA).getSolver();
		RealVector solution = solver.solve(Atb);
		SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
		final pgL pgL_fun = new pgL(xMean, xStd, xCi, xDi);
		// for debugging
		final PointValuePair pEst = optimizer.optimize(new MaxEval(10000), new ObjectiveFunction(pgL_fun), GoalType.MINIMIZE,
				new InitialGuess(new double[] {solution.getEntry(0), solution.getEntry(1)} ), new NelderMeadSimplex(new double[] { 0.0001, 0.0001 }));
		double[] pEst0 = new double[2];
		pEst0[0] = pEst.getPoint()[0];//0.0080;//
		pEst0[1] = pEst.getPoint()[1];//-1.4e-4;//
		//System.out.println("alpha: "+pEst0[0]+" sigma:"+pEst0[1]);
		double alpha = pEst0[0];
		double sigma2 = Math.max(pEst0[1],2e-4); // this is the var of gausssian noise
		double t0 = 2/alpha*Math.sqrt(alpha*0+3/8*alpha*alpha+sigma2);
		double t1 = 2/alpha*Math.sqrt(alpha*1+3/8*alpha*alpha+sigma2);
		
		double[][] outputGt = new double[double_G.length][double_G[0].length];
		for(int i=0;i<double_G.length;i++){
			for(int j=0;j<double_G[0].length;j++){
				outputGt[i][j] = 2/alpha*Math.sqrt(alpha*double_G[i][j]+3/8*alpha*alpha+sigma2);
				outputGt[i][j] = 255*(outputGt[i][j]-t0)/(t1-t0);
			}
		}
		double var1 = 1/((t1-t0)*(t1-t0));
		System.out.println("alpha: "+pEst0[0]+" sigma:"+pEst0[1]+" Variance:"+var1);
		
		if(qq.pgEst==null) {
			qq.pgEst = new double [pEst0.length];
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = pEst0[i];
		}
		else {
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = qq.pgEst[i]+pEst0[i];
		}
		if (Double.isNaN(qq.var))
			qq.var = var1*255*255;
		else
			qq.var = qq.var + var1*255*255;
		
		return outputGt;
		//qq.Nx = double_G[0].length;
		//qq.Ny = double_G.length;
		//qq.ntry = ntry;
		//return qq;
	}

	// 2d filter on images
	public double[][] imfilter(double[][] origIm, double[][] filter) {
		double[][] outIm = new double[origIm.length][origIm[0].length];
		int tmp_h, tmp_w;
		int filter_h = filter.length;
		// int filter_w = filter[0].length;
		int shift = (filter_h - 1) / 2;
		double med_value;
		for (int i = 0; i < origIm.length; i++) {
			for (int j = 0; j < origIm[0].length; j++) {
				outIm[i][j] = 0;
				for (int m = 0; m < filter_h; m++) {
					for (int n = 0; n < filter_h; n++) {
						tmp_h = i - shift + m;
						tmp_w = j - shift + n;
						if (tmp_h >= 0 & tmp_h < origIm.length & tmp_w >= 0 & tmp_w < origIm[0].length) {
							med_value = filter[m][n] * origIm[tmp_h][tmp_w];
						} else {
							med_value = 0;
						}
						outIm[i][j] = outIm[i][j] + med_value;
					}
				}
			}
		}
		return outIm;
	}
	// media filter with a fixed template size
	public double[][] medfilter(double[][] origIm, int templateSize) {

		double vals[];
		int count;
		double[][] outIm = origIm;

		for (int x = (templateSize - 1) / 2; x < origIm[0].length - (templateSize + 1) / 2; x++) {
			for (int y = (templateSize - 1) / 2; y < origIm.length - (templateSize + 1) / 2; y++) {
				vals = new double[templateSize * templateSize];
				count = 0;
				for (int x1 = x - ((templateSize - 1) / 2); x1 < x + ((templateSize + 1) / 2); x1++) {
					for (int y1 = y - ((templateSize - 1) / 2); y1 < y + ((templateSize + 1) / 2); y1++) {
						vals[count] = origIm[y1][x1];
						count++;
					}
				}
				java.util.Arrays.sort(vals);
				double middle = vals[((templateSize * templateSize + 1) / 2)-1];
				outIm[y][x] = middle;
			}
		}

		return outIm;
	}
	// discrete wavelet transform on 2d matrix
	public double[][] dwt2(double[][] origIm, double[] lp) {

		double[][] lp2d_c = new double[lp.length][1];
		double[][] lp2d_r = new double[1][lp.length];

		for (int i = 0; i < lp.length; i++) {
			lp2d_r[0][i] = lp[i];
			lp2d_c[i][0] = lp[i];
		}
		double[][] cA_1 = new double[origIm.length][origIm[0].length];
		cA_1 = Convolution.convolution2DPadded(origIm, origIm.length, origIm[0].length, lp2d_c, lp2d_c.length,
				lp2d_c[0].length);
		double[][] cA_2 = new double[cA_1.length/2][cA_1[0].length];//
		for (int i = 0; i < cA_2.length; i++)
			for (int j = 0; j < cA_2[0].length; j++)
				cA_2[i][j] = cA_1[2*i+1][j];

		double[][] cA_3 = new double[cA_2.length][cA_2[0].length];
		cA_3 = Convolution.convolution2DPadded(cA_2, cA_2.length, cA_2[0].length, lp2d_r, lp2d_r.length,
				lp2d_r[0].length);
		double[][] cA_4 = new double[cA_3.length][cA_3[0].length / 2];// 
		for (int i = 0; i < cA_4.length; i++)
			for (int j = 0; j < cA_4[0].length; j++)
				cA_4[i][j] = cA_3[i][j*2+1];
		return cA_4;
	}
	/*public void DisplayROI(){
		int[][] roi_pt1 = new int[nSyn0][2]; //larger than enough
		long [] roi_size = new long[nSyn0];
		for(int i = 0; i<height;i++){
			for(int j = 0; j<width;j++){
				int tmp = SynR1Idx[i][j];
				if(tmp!=0){
					roi_pt1[tmp-1][0] = i;
					roi_pt1[tmp-1][1] = j;
					roi_size[tmp-1] = roi_size[tmp-1]+1;
				}
			}
		}
		String imtitle = "";
		ImagePlus newimp = NewImage.createByteImage (imtitle, width, height, 1,
				NewImage.FILL_WHITE);
		ImageProcessor impNP = newimp.getProcessor(); 
		for(int i = 0; i<height;i++){
			for(int j = 0; j<width;j++){
				impNP.putPixel(j,i,SynR1Idx[i][j]);
			}
		}
		ImageStack impstack = newimp.getStack();
		// Generate roimanager
		ByteProcessor ip = (ByteProcessor)impstack.getProcessor(1).convertToByte(true);
		RoiManager manager = new RoiManager();
		int[] roi2fiu = new int[nSyn0];
		int roi_cnt = 0;
		double wandVal = 0.01;
		for(int i=0;i<nSyn0;i++){
			roi2fiu[roi_cnt] = i;
			roi_cnt++;
			Wand w = new Wand(ip);
			w.autoOutline(roi_pt1[i][1],roi_pt1[i][0],wandVal,Wand.EIGHT_CONNECTED); 
			if (w.npoints>0) { // we have an roi from the wand... 
				Roi roi = new PolygonRoi(w.xpoints, w.ypoints, w.npoints, Roi.TRACED_ROI);
				imp.setRoi(roi);
				manager.addRoi(roi);
			}
		}
		manager.runCommand("show all with labels");
		manager.setSize(300, 400);
	}*/
}
