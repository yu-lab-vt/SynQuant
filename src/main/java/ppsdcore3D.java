import java.util.ArrayList;


public class ppsdcore3D{
	protected int [][][] kMap;//
	protected double [][][] zMap;
	protected double [][][] thrMap;
	protected ArrayList<double[]> zscoreList = new ArrayList<double[]>();
	protected double thrZ;
	protected int nSyn0;
	
    public ppsdcore3D(short[][] imArray, double [][][] G,double [][][] Gt, boolean[][][] kMask, paraP3D p, paraQ3D q) {
    	int zSlice = Gt.length;
    	int height = Gt[0].length;
    	int width = Gt[0][0].length;
    	kMap = new int[zSlice][height][width];
    	zMap = new double[zSlice][height][width];
    	thrMap = new double[zSlice][height][width];
    	//initialization
    	boolean doneAll = false;
    	int loopCnt = 0;
    	nSyn0 = 1;
    	BasicMath BM = new BasicMath();
    	long startTime1=System.nanoTime();
    	scanAllThres3D scanMap = new scanAllThres3D(imArray, G,Gt,kMask,p,q);// thresholding images with multi thresholds(50:10:220)
    	
    	/*ImageHandling IM = new ImageHandling();
    	for(int i=0;i<scanMap.kMapx.length;i++) {
    		System.out.println("ROI Z-socre: "+ scanMap.zMapx[i][7][91][77]);//+" * "+scanMap.zMapx[0][0].length);
    		IM.maskDisplay3D(scanMap.kMapx[i], "Threshold #"+i);
    	}*/
    	//byte[] usedN = new byte[imArray.length];
    	//ComponentTree CT_scanMap = new ComponentTree(imArray,usedN, G[0].length, G.length,p, q);
    	long endTime1=System.nanoTime();
    	System.out.println("scanAll3 First Running time: "+(endTime1-startTime1)/1e9); 
    	ArrayList<Double> Iter_PopUplist = new ArrayList<Double>(); //zscores of detected synapses
    	startTime1=System.nanoTime();
    	double bs = 0,sib=0,suc=0;
    	while(!doneAll){
    		loopCnt = loopCnt + 1;
    		ArrayList<Integer[]> idxUpdt = new ArrayList<Integer[]>();
    		long startTime2=System.nanoTime();
    		BestSyn3D bestSynMap = new BestSyn3D(scanMap, p, Iter_PopUplist);//find best region with highest zscore
    		//for debug
    		//if (bestSynMap.kMap0[4][79][24]) {
    		//	System.out.println("ROI");
    		//	IM.maskDisplay3D(bestSynMap.kMap0, "bestSynMap");
    		//}
    		
    		thrZ = bestSynMap.thrZ;
    		//System.out.println("DO we found synapse: "+BM.matrix3DMax(bestSynMap.kMap0)+"Max zScore: "+thrZ);
    		long endTime2=System.nanoTime();
        	bs = bs+(endTime2-startTime2)/1e9;
    		Iter_PopUplist = new ArrayList<Double>(bestSynMap.PopUplist);
    		// further scan within best region -----
    		if(BM.Allfalse3D(bestSynMap.kMap0)){
    			doneAll = true;
    		}
    		else{//start scan, find if there is a region score higher inside the best region
    			startTime2=System.nanoTime();
    			SynInBest3D SynBBest = new SynInBest3D(Gt,bestSynMap,p,q, false/*containFlag==false*/);
    			endTime2=System.nanoTime();
    			sib = sib+(endTime2-startTime2)/1e9;
        		// for debug
        		//if (bestSynMap.kMap0[4][79][24]) {
        		//	System.out.println("ROI");
        		//	IM.maskDisplay3D(SynBBest.kMap2, "SyninBest");
        		//}
    			boolean firstFlag = true; //only need to add zscore to list once
    			if(SynBBest.foundOne){
    				//insert the found region(marked in SynBBest) into kMap, zMap, thrMap
    				for(int i=0;i<SynBBest.kMap2.length;i++){
    					for(int j=0;j<SynBBest.kMap2[0].length;j++){
    						for(int kk=0;kk<SynBBest.kMap2[0][0].length;kk++){
    							if(SynBBest.kMap2[i][j][kk]){
    								kMap[i][j][kk] = nSyn0;
    								zMap[i][j][kk] = SynBBest.zMap2[i][j][kk];
    								thrMap[i][j][kk] = SynBBest.thrMap2[i][j][kk];
    								if(firstFlag) {
    									firstFlag = false;
    									zscoreList.add(new double[] {nSyn0,SynBBest.zMap2[i][j][kk], SynBBest.thrMap2[i][j][kk]});
    								}
    							}
    						}
    					}
    				}
    			}
    			else{//otherwise, insert the original region(marked in bestSynMap) into kMap, zMap, thrMap
    				for(int i=0;i<bestSynMap.kMap0.length;i++){
    					for(int j=0;j<bestSynMap.kMap0[0].length;j++){
    						for(int kk=0;kk<bestSynMap.kMap0[0][0].length;kk++){
    							if(bestSynMap.kMap0[i][j][kk]){
    								kMap[i][j][kk] = nSyn0;
    								zMap[i][j][kk] = bestSynMap.zMap0[i][j][kk];
    								thrMap[i][j][kk] = bestSynMap.thrMap0[i][j][kk];
    								if(firstFlag) {
    									firstFlag = false;
    									zscoreList.add(new double[] {nSyn0,SynBBest.zMap2[i][j][kk], SynBBest.thrMap2[i][j][kk]});
    								}
    							}
    						}
    					}
    				}
    			}
    			nSyn0++;
    			//update kMask: kMask(kMap>0) = 0;
    			//update idxUpdt: idxUpdt = find(kMap2>0);
    			for(int i=0;i<kMask.length;i++){
    				for(int j=0;j<kMask[0].length;j++){
    					for(int kk=0;kk<kMask[0][0].length;kk++){
    						if(kMap[i][j][kk]>0)
    							kMask[i][j][kk] = false;
    						if(SynBBest.kMap2[i][j][kk])
    							idxUpdt.add(new Integer[] {i,j,kk});
    					}
    				}
    			}
    			startTime2=System.nanoTime();
    			scanMap.scanUpdtCrop(G,Gt,kMask,idxUpdt,p,q);
    			if(scanMap.CT.Zscore_Vec.size()==0) {
    				doneAll = true;
    			}
    				
    			endTime2=System.nanoTime();
    			suc = sib+(endTime2-startTime2)/1e9;
    		}
    	}
    	endTime1=System.nanoTime();
    	System.out.println("While Loop Running time: "+(endTime1-startTime1)/1e9+" Loop:"+loopCnt); 
    	System.out.println("BestSyn,"+bs+" SynInBest,"+sib+" scanUpdtCrop,"+suc); 
    	nSyn0--;
    }
    // post processing with size constrain and pvalue constrain
    public boolean [][][] ppsd_post3D(double [][][] Gt, paraP3D p, paraQ3D q){
    	//System.out.println("Analysis----");
    	ImageHandling IH = new ImageHandling();

    	boolean[][][] kMapi = new boolean[kMap.length][kMap[0].length][kMap[0][0].length];
    	
    	double thrZx = thrZ;// <-2? -2:thrZ; //pvalue threshold should not be too small, otherwise FDR control will be not useful
    	int[] roi_size = new int[nSyn0];
    	for(int i=0;i<kMap.length;i++)
    		for(int j=0;j<kMap[0].length;j++)
    			for(int k=0;k<kMap[0][0].length;k++)
	    			if(kMap[i][j][k]!=0)
	    				roi_size[kMap[i][j][k]-1]++;
    	
    	for(int jj=0;jj<nSyn0;jj++){
    		int idxSyn0 = jj+1;
    		int[][] idx = new int[roi_size[jj]][3];
    		int idx_cnt = 0;
        	for(int i=0;i<kMap.length;i++)
        		for(int j=0;j<kMap[0].length;j++)
        			for(int k=0;k<kMap[0][0].length;k++)
        				if(kMap[i][j][k]==idxSyn0){
        					idx[idx_cnt][0] = i;
        					idx[idx_cnt][1] = j;
        					idx[idx_cnt++][2] = k;
        				}

        	if(idx.length<10*q.Pix_per_synapse){
        		boolean[][][] bmask  = new boolean[q.Nz][q.Ny][q.Nx];
        		boolean[][][] tmask  = new boolean[q.Nz][q.Ny][q.Nx];
				for (int i = 0; i < q.Nz; i++) {
					for (int j = 0; j < q.Ny; j++) {
						for (int k = 0; k < q.Nx; k++) {
						bmask[i][j][k] = true;
						tmask[i][j][k] = false;
						}
					}
				}
				double zScore = IH.scanOneSyn3D(Gt,tmask,bmask,idx,p,q);
        		if(zScore>thrZx){
        			for(int i=0;i<idx.length;i++){
        				kMapi[idx[i][0]][idx[i][1]][idx[i][2]] = true;
        			}
        		}
        	}
    	}
    	
    	return kMapi;
    }
   

}