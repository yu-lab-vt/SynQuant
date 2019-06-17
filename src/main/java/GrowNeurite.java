import java.util.ArrayList;

import ij.ImagePlus;
import ij.gui.NewImage;
import ij.process.ImageProcessor;

//parameter for dendrite extraction
//Details can be found in Paper: "Automated analysis of neuronal morphology, synapse number and synaptic recruitment"
public class GrowNeurite {
	protected int height;
    protected int width;
    protected boolean[][] neuriteMask;
    protected int min_length = 5;
    protected double xyRes;//one pixel corresponding to such distance in meters
    protected double filterSize = 5e-7;//in meters, corresponding to about 6 pixels in normal resolution;
    protected double maxAddCost = 0.90;// Maximal cost if we still want to add pixel
    protected double connectCostLambda = 0.7;
    protected double[][] rigidity;
    protected double[][][] dirVect;
    protected double lambdaMax;
    protected boolean[][] den_ext;
    protected int[][] CurSkelIm;
    protected int[][] CurIm;
    protected int CurNum;
    protected double[][] dendrite;
    protected int min_intensity = 20;//too dim
    protected double[] denSize;
    protected double[] denIntensity;
    protected double[] denLength;
    protected int maxDisSkel = 20;
    protected int DenMaxLen = 200;
    protected long totallength=0;
    
	public GrowNeurite(short[] imArray, int inwidth, int inheight,double vox_x,boolean displayROI, String RoiTitle){
		width = inwidth;
		height = inheight;
		dendrite = new double[height][width];
		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++)
				dendrite[i][j] = (double) imArray[i * width + j];// the data is saved by line
		
		height = dendrite.length;
		width = dendrite[0].length;
		xyRes = vox_x;
		neuriteMask = new boolean[height][width];
		boolean[][] den_nosiy = new boolean[height][width];
		steerableFilter(dendrite);//get the response of 2d gaussian filters
		
		boolean [][] somaIdx = new boolean[height][width];
		
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				neuriteMask[i][j] = false;
				den_nosiy[i][j] = false;
				somaIdx[i][j] = true;//all pixels could be a start point of dendrite
			}
		}
		
		for(int y=0;y<height;y++){
			for(int x=0;x<width;x++){
				if(!somaIdx[y][x])
					continue;
				ArrayList<Integer[]> sm_part = new ArrayList<Integer[]>();
				sm_part.add(new Integer[] {y,x});
				neuriteMask[y][x] = true;
				ArrayList<Integer[]> pixelQue = new ArrayList<Integer[]>();
				pixelQue.add(new Integer[] {y,x});
				while(!pixelQue.isEmpty()){
					//test if the neighbors belongs to same dendrite, if so start again from the neighbor
					//if not, start from another soma pixel(here from any other pixel that has not been tested)
					ArrayList<Integer[]> neigh = processNeighbourhood(pixelQue.get(0)[1],pixelQue.get(0)[0]);
					pixelQue.remove(0);
					if(!neigh.isEmpty()){
						sm_part.addAll(neigh);
						pixelQue.addAll(neigh);
						for(int i=0;i<neigh.size();i++){
							neuriteMask[neigh.get(i)[0]][neigh.get(i)[1]] = true;
						}
					}
				}
				if(sm_part.size()==1)
					continue;
				for(int i=0;i<sm_part.size();i++){
					den_nosiy[sm_part.get(i)[0]][sm_part.get(i)[1]] = true;
					somaIdx[sm_part.get(i)[0]][sm_part.get(i)[1]] = false;
				}
			}
		}
		
		ImageHandling IH = new ImageHandling();
		den_nosiy = IH.imopen(den_nosiy, 4);
		den_ext = IH.bwareaopen(den_nosiy, 200, 4);
		boolean[][] skel = IH.Skeleton(den_ext);//extract skeleton of dendrite
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(skel[i][j])
					totallength++;
			}
		}
		
		ArrayList<Integer[]> branchendPts = IH.branchendpts(skel);//get branch points and end points
		
		CurSkelIm = CuveSkelImGen(skel,min_length,branchendPts);//cut dendrite skeleton into pieces
		CurIm = CurImGen(CurSkelIm,den_ext,skel);//based on skeleton pieces, cut dendrite into pieces
		ImagePlus outimp = createImage();
		outimp.show();
		outimp.updateAndDraw();
		displayROI = false;//Usually not show dendrite ROI regions
		if(displayROI){
			boolean[][] BWCurveIm = new boolean[CurIm.length][CurIm[0].length];
			for (int i = 0; i < width; i++)
				for (int j = 0; j < height; j++)
					if (CurIm[j][i] > 0)
						BWCurveIm[j][i] = true;

			int[][] denIdx = IH.bwlabel(BWCurveIm, 8);
			IH.DisplayROI(IH.NextLabel, height, width, denIdx, outimp,RoiTitle);
		}
	}
	private int[][] CurImGen(int[][] skelLabel, boolean[][] den, boolean[][] skel) {
		// TODO Auto-generated method stub
		ArrayList<Integer[]> btpts = new ArrayList<Integer[]>();
		for(int y=0;y<skelLabel.length;y++)
			for(int x=0;x<skelLabel[0].length;x++)
				if(skel[y][x])
					btpts.add(new Integer[] {y,x});
		denSize = new double[CurNum];
		int[][] outCurIm = new int[skel.length][skel[0].length];
		for(int y=0;y<den.length;y++){
			for(int x=0;x<den[0].length;x++){
				if(!den[y][x])
					continue;
				int Label = 0;
				double MinDis = 1e6;
				for(int ii=0;ii<btpts.size();ii++){
					int y_d = btpts.get(ii)[0]-y;
					if(y_d>maxDisSkel || y_d<-maxDisSkel)
						continue;
					int x_d = btpts.get(ii)[1]-x;
					if(x_d>maxDisSkel || x_d<-maxDisSkel)
						continue;
					double tmpDis = Math.sqrt(y_d*y_d+x_d*x_d);
					if(tmpDis<MinDis){
						MinDis = tmpDis;
						Label = skelLabel[btpts.get(ii)[0]][btpts.get(ii)[1]];
					}
				}
				if(Label>0){
					outCurIm[y][x] = Label;
					denSize[Label-1]++;
				}
			}
		}
		return outCurIm;
	}
	private int[][] CuveSkelImGen(boolean[][] skel, int MinLen, ArrayList<Integer[]> bePts) {
		int L = skel.length;
		int W = skel[0].length;
		boolean[][] BW = new boolean[L][W];
		for(int i=0;i<L;i++)
			for(int j=0;j<W;j++)
				BW[i][j] = skel[i][j];
		for(int i=0;i<bePts.size();i++){
			BW[bePts.get(i)[0]][bePts.get(i)[1]] = false;
		}
		ImageHandling IH = new ImageHandling();
		BW = IH.bwareaopen(BW, MinLen-1,8);
		int[][] CurLabel = IH.Regbwlabel(BW,8,DenMaxLen);
		CurNum = IH.NextLabel;
		for(int i=0;i<bePts.size();i++){
			int y = bePts.get(i)[0];
			int x = bePts.get(i)[1];
			CurLabel[y][x] = NeighFind(CurLabel,x,y);
		}
		long[] sumIntensity = new long[CurNum];
		int[] pixSize = new int[CurNum];
		for(int i=0;i<L;i++){
			for(int j=0;j<W;j++){
				if(CurLabel[i][j]>0){
					sumIntensity[CurLabel[i][j]-1]+=dendrite[i][j];
					pixSize[CurLabel[i][j]-1]++;
				}
			}
		}
		int[] labeltransfer = new int[CurNum];
		int[] meanIntensity = new int[CurNum];
		int DelCnt = 0;
		for(int i=0;i<CurNum;i++){
			meanIntensity[i] = (int) Math.round((sumIntensity[i]/pixSize[i]));
			if(meanIntensity[i]<min_intensity){
				DelCnt++;
				for(int y=0;y<L;y++){
					for(int x=0;x<W;x++){
						if(CurLabel[y][x]==(i+1)){
							CurLabel[y][x]=0;
						}
					}
				}
			}
			else
				labeltransfer[i] =  i+1-DelCnt;
		}
		CurNum = CurNum-DelCnt;
	    denIntensity = new double[CurNum];
	    denLength = new double[CurNum];
	    for(int i=0;i<CurNum+DelCnt;i++){
	    	if(labeltransfer[i]>0){
	    		denIntensity[labeltransfer[i]-1] = meanIntensity[i];
		    	denLength[labeltransfer[i]-1]=pixSize[i];
		    }
	    }
		for(int y=0;y<L;y++){
			for(int x=0;x<W;x++){
				if(CurLabel[y][x]>0){
					CurLabel[y][x]=labeltransfer[CurLabel[y][x]-1];
				}
			}
		}
		return CurLabel;
	}
	private int NeighFind(int[][] CurLabel, int x, int y){
		int outN = 0;
		double minDis = 100000;
		for(int ni=-1; ni<=1; ni++) {
            for(int nj=-1; nj<=1; nj++) {
				if (x + ni < 0 || y + nj < 0 || x + ni > width - 1 || y + nj > height - 1) {
					continue;
				} else {
					if(ni==0 && nj==0) continue;
					if(CurLabel[y+nj][x+ni]>0){
						double tmpDis = Math.sqrt((double)nj*nj+ni*ni);
						if(tmpDis<minDis){
							minDis = tmpDis;
							outN = CurLabel[y+nj][x+ni];
						}
							
					}
				}
                }
            }
		return outN;
	}
	private ArrayList<Integer[]> processNeighbourhood(int x, int y) {
		ArrayList<Integer[]> nHoodMask = new ArrayList<Integer[]>();
		for (int ni = -1; ni <= 1; ni++) {
			for (int nj = -1; nj <= 1; nj++) {
				if (x + ni < 0 || y + nj < 0 || x + ni > width - 1 || y + nj > height - 1) {
					continue;
				} else {
					if (ni == 0 && nj == 0)
						continue;
					if(!neuriteMask[y+nj][x+ni] & connectCost(x,y,x+ni,y+nj)< maxAddCost){
						nHoodMask.add(new Integer[] {y+nj, x+ni});
					}
				}
			}
		}
		return nHoodMask;
	}
	private double connectCost(int x, int y, int nx, int ny){
		double ccLambda;
		if(rigidity[ny][nx]>0)
			ccLambda = 1-rigidity[ny][nx]/lambdaMax;
		else
			ccLambda = 1;
		
		
		double norm_ab = Math.sqrt((nx-x)*(nx-x)+(ny-y)*(ny-y));
		double [] dirAB = {(ny-y)/norm_ab,(nx-x)/norm_ab};
		
		double[] tmpDirA = {dirVect[y][x][0],dirVect[y][x][1]};
		double norm_a = Math.sqrt(tmpDirA[0]*tmpDirA[0]+tmpDirA[1]*tmpDirA[1]);
		double[] wA = new double[2];
		wA[0] = tmpDirA[0]/norm_a;
		wA[1] = tmpDirA[1]/norm_a;
		double tmp = Math.abs(dirAB[0]*wA[0]+dirAB[1]*wA[1]);
		double ccV = 0.5*(Math.sqrt(1-tmp)+Math.sqrt(1+tmp));
		double cc = connectCostLambda*ccLambda+(1-connectCostLambda)*ccV;
				
		return cc;
	}
	public void steerableFilter(double [][] dendrite){
		//GaussianDerivatives2D
		double [][]  g0x, g0y, g1x, g1y, g2x, g2y;
		double sigma = filterSize/xyRes;
		double W = 2*Math.sqrt(2)*sigma;
		double[] R = new double[(int)(2*W)+1];
		int vecsize = R.length;
		for(int i=0;i<R.length;i++)
			R[i] = i-W;
		double k = 1/(sigma*sigma);
		g0x = new double[1][vecsize];
		g0y = new double[vecsize][1];
		g1x = new double[1][vecsize];
		g1y = new double[vecsize][1];
		g2x = new double[1][vecsize];
		g2y = new double[vecsize][1];
		for(int i=0;i<vecsize;i++){
			double g0 = Math.exp(-k*R[i]*R[i]);
			g0x[0][i] = g0;
			g0y[i][0] = g0;
			double g1 = -2*k*R[i]*g0;
			g1x[0][i] = g1;
			g1y[i][0] = g1;
			double g2 = 2*k*(2*k*R[i]*R[i]-1)*g0;
			g2x[0][i] = g2;
			g2y[i][0] = g2;
		}
		double N = 1;//4/(Math.pow(sigma, 4));
		double [][] Rx2 = Convolution.convolution2DPadded(dendrite,  height,width, g2x,1, vecsize);
		double [][] Rxx = Convolution.convolution2DPadded(Rx2,Rx2.length, Rx2[0].length, g0y, vecsize, 1);
		
		double [][] Rx1 = Convolution.convolution2DPadded(dendrite, height, width, g1x,1, vecsize);
		double [][] Rxy = Convolution.convolution2DPadded(Rx1,Rx1.length, Rx1[0].length, g1y, vecsize, 1);
		
		double [][] Rx0 = Convolution.convolution2DPadded(dendrite, height, width, g0x,1, vecsize);
		double [][] Ryy = Convolution.convolution2DPadded(Rx0,Rx0.length, Rx0[0].length, g2y, vecsize, 1);

		
		
		for(int i=0;i<Rxx.length;i++){
			for(int j=0;j<Rxx[0].length;j++){
				Rxx[i][j] = -N*Rxx[i][j];
				Rxy[i][j] = -N*Rxy[i][j];
				Ryy[i][j] = -N*Ryy[i][j];
			}
			
		}
		/*ImageHandling IH = new ImageHandling();
		IH.IMwrite(Rxx, "src/res/Rxx.png");
		IH.IMwrite(Rxy, "src/res/Rxy.png");
		IH.IMwrite(Ryy, "src/res/Ryy.png");*/
		
		rigidity  = new double [Rxx.length][Rxx[0].length];
		dirVect  = new double [Rxx.length][Rxx[0].length][2];
		lambdaMax = 0;
		for(int i=0;i<Rxx.length;i++){
			for(int j=0;j<Rxx[0].length;j++){
				double tmp1,tmp2;
				tmp1 = Rxx[i][j]+Ryy[i][j];
				tmp2 = Math.sqrt(tmp1*tmp1+Rxy[i][j]*Rxy[i][j]-Rxx[i][j]*Ryy[i][j]);
				rigidity[i][j] = tmp1+tmp2;
				dirVect[i][j][0] = 1;
				dirVect[i][j][1] = (tmp1-tmp2-Rxx[i][j])/Rxy[i][j];
				double tmp = Math.sqrt(dirVect[i][j][0]*dirVect[i][j][0]+dirVect[i][j][1]*dirVect[i][j][1]);
				dirVect[i][j][0] = dirVect[i][j][0]/tmp;
				dirVect[i][j][1] = dirVect[i][j][1]/tmp;
				if(rigidity[i][j]>lambdaMax)
					lambdaMax = rigidity[i][j];
			}
		}
		
	}
	// Color image by combining detection results(white) with original image(dendrite)
	public ImagePlus createImage(){
		ImagePlus outimp = NewImage.createRGBImage("Dendrite Extraction Results", width, height, 1,NewImage.FILL_BLACK);
		ImageProcessor outIP =  outimp.getProcessor();
		
		for (int i = 0; i < width; i++) {
			for(int j = 0;j<height;j++){
			// put channel values in an integer
			if(CurIm[j][i]>0)
				outIP.set(i, j, (((int)200 & 0xff) << 16)
				          + (((int)dendrite[j][i] & 0xff) << 8)
				          +  ((int)0 & 0xff));
			else
				outIP.set(i, j, (((int)0 & 0xff) << 16) 
				          + (((int)dendrite[j][i] & 0xff) << 8)
				          +  ((int)0 & 0xff));

			}
		}
		//outimp = NewImage.createRGBImage(pixels, width, height, true);
		return outimp;
	}
	
}