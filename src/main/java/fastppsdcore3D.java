import java.util.ArrayList;
import java.util.Collections;


public class fastppsdcore3D{
	public int [][][] kMap;//
	public double [][][] zMap;
	public double thrZ;
	protected int nSyn0;
	protected int zSlice; //image width
	protected int width; //image width
	protected int height; //image height
	//protected ComponentTree3D4Fast[] CT;
	protected ComponentTree3D4Fast CT;
	
    public fastppsdcore3D(short[][] imArray, int inWidth, int inHeight, int inZSlice, paraP3D p, paraQ3D q) {
    	// STEP 1. initialization
    	width = inWidth;
    	height = inHeight;
    	zSlice = inZSlice;
    	thrZ = 1000;
    	//CT = new ComponentTree3D4Fast[3];
		byte[] usedN = new byte[imArray[0].length*imArray.length];
		for (int i=0;i<usedN.length;i++)
			usedN[i] = (byte)0;
		// STEP 2. build the component tree
		CT = new ComponentTree3D4Fast(imArray,usedN, p, q);
//		for (int i=0;i<CT.length;i++) {
//			CT[i] = new ComponentTree3D4Fast(imArray,usedN, p, q);
//			for (int j=0;j<usedN.length;j++) {
//				if (CT[i].outputArray[j] > 0) {
//					usedN[j] = (byte)2;  // this pixel will not be used as neighbors anymore
//				}
//			}
//		}
		cumulateZscore(q);
    }
    public void cumulateZscore(paraQ3D q) {
		// STEP 3. start building z-score map
		zMap = new double[zSlice][height][width];
		int nVoxels=width*height*zSlice; // voxels in a 3D image
		int nPixels = width*height; //pixels in a single slice
		int rmder, z, y, x;
		ImageHandling IM = new ImageHandling();
		// Step 3.1. build synapse labels
		// double[][][] fgzscoreMap = new double[zSlice][height][width];
		for (int j = nVoxels-1; j >= 0; j--)
		{
			if ( CT.outputArray[j] > 0) {
				rmder = j % nPixels;
				z = j / nPixels;
				y=rmder/width;
				x=rmder-y*width;
				zMap[z][y][x] = CT.zscore[j];
			}
		}
		kMap = IM.bwlabel3D(zMap, 26);
		nSyn0 = IM.NextLabel;
		// Step 3.2. cumulate the zscore from previous channels: add z-score of pre-synapse to post-synapse if there is overlap
		if (q.NumChannelProcessed > 0){
			if( q._wayCombinePrePost==1) 
			{
				double[] pre_z = new double[nSyn0];
				for(int i=0; i< pre_z.length; i++)
					pre_z[i] = 0;
				for(int i=0; i<kMap.length; i++) {
					for(int j=0; j<kMap[0].length; j++) {
						for (int k=0; k<kMap[0][0].length; k++) {
							if (kMap[i][j][k]>0 ) {
								if(pre_z[kMap[i][j][k]-1] < q.synZscore[q.curTps][i][j][k]) {
									pre_z[kMap[i][j][k]-1] = q.synZscore[q.curTps][i][j][k];
								}
							}
						}
					}
				}
				for (int j = nVoxels-1; j >= 0; j--)
				{
					if ( CT.outputArray[j] > 0) {
						rmder = j % nPixels;
						z = j / nPixels;
						y=rmder/width;
						x=rmder-y*width;
						zMap[z][y][x] = CT.zscore[j] + pre_z[kMap[z][y][x]-1];
						if (thrZ > zMap[z][y][x])
							thrZ = zMap[z][y][x];
					}
				}
			}else
			{
				if(q.ExtendedDistance>0) {
					ExtendPreChannelResult NewPre=new ExtendPreChannelResult(q);
					q.synZscore=NewPre.ExtendedResult;
				}
				
				double[] pre_z = new double[nSyn0];
				for(int i=0; i< pre_z.length; i++)
					pre_z[i] = 0;
				for(int i=0; i<kMap.length; i++) {
					for(int j=0; j<kMap[0].length; j++) {
						for (int k=0; k<kMap[0][0].length; k++) {
							if (kMap[i][j][k]>0) {
								if(q.synZscore[q.curTps][i][j][k]>0) {
									if(pre_z[kMap[i][j][k]-1] < q.synZscore[q.curTps][i][j][k]) {
										pre_z[kMap[i][j][k]-1] = q.synZscore[q.curTps][i][j][k];
									}
								}else {
									zMap[i][j][k] = 0;
									kMap[i][j][k] = 0;
								}
							}
						}
					}
				}
				for (int j = nVoxels-1; j >= 0; j--)
				{
					if ( CT.outputArray[j] > 0) {
						rmder = j % nPixels;
						z = j / nPixels;
						y=rmder/width;
						x=rmder-y*width;
						if (kMap[z][y][x]>0) {
							if (pre_z[kMap[z][y][x]-1]>0)
								zMap[z][y][x] = CT.zscore[j] + pre_z[kMap[z][y][x]-1];
							else {
								zMap[z][y][x] = 0;
								kMap[z][y][x] = 0;
							}
							if (thrZ > zMap[z][y][x])
								thrZ = zMap[z][y][x];
						}
						else {
							CT.outputArray[j] = 0;
						}
					}
				}
			}
		}else 
		{
			for (int j = nVoxels-1; j >= 0; j--)
			{
				if ( CT.outputArray[j] > 0) {
					if (thrZ > CT.zscore[j])
						thrZ = CT.zscore[j];;
				}
			}

		}
    }
    public void fdrControl(paraP3D p, paraQ3D q) {
		BasicMath mBM = new BasicMath();
		// Step 3.3. fdr control. This happens when we are processing the last channel
		if (q._NumChannel == q.NumChannelProcessed + 1 & p.zscore_thres <=0) {
			// step 3.3.1 get all zscores
			double[] z_vec = new double[nSyn0];
			for(int i=0; i< z_vec.length; i++)
				z_vec[i] = 0;
			
			for(int i=0; i<kMap.length; i++) {
				for(int j=0; j<kMap[0].length; j++) {
					for (int k=0; k<kMap[0][0].length; k++) {
						if (kMap[i][j][k]>0) {
							if(z_vec[kMap[i][j][k]-1] < zMap[i][j][k]) {
								if (z_vec[kMap[i][j][k]-1] > 0)
									System.out.println("we give a synapse two different zscores!");
								z_vec[kMap[i][j][k]-1] = zMap[i][j][k];
							}
						}
					}
				}
			}
			// step 3.3.2 start to do fdr control
			int total_len = nSyn0;
			double[] FDR_Con = new double[total_len];
			for(int i=0;i<total_len;i++)
				FDR_Con[i] = p.fdr*((double)(i+1))/total_len;
			if(p.fdr_den){
				for(int i=0;i<total_len;i++){
					FDR_Con[i] = FDR_Con[i]/Math.log(total_len);
				}
			}
			ArrayList<Double> pvalue = new ArrayList<Double>();
			for(int i=0;i<z_vec.length;i++){
				pvalue.add(mBM.zTop(z_vec[i]));
			}
			Collections.sort(pvalue);// Ascend; add Collections.reverseOrder() for descend
			thrZ = 0;
			for(int i=pvalue.size()-1;i>=0;i--){
				if(pvalue.get(i)<=FDR_Con[i]){
					thrZ = mBM.pToZ(pvalue.get(i));
					break;
				}
				nSyn0--;
			}
			// generate zscore Map
			for(int i=0; i<kMap.length; i++) {
				for(int j=0; j<kMap[0].length; j++) {
					for (int k=0; k<kMap[0][0].length; k++) {
						if (kMap[i][j][k]>0) {
							if (z_vec[kMap[i][j][k]-1] < thrZ) {
								zMap[i][j][k] = 0;
								kMap[i][j][k] = 0;
							}
							else {
								zMap[i][j][k] = z_vec[kMap[i][j][k]-1];
							}
						}
					}
				}
			}
		}		
    }
   
   

}