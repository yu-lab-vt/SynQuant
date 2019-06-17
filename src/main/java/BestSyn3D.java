import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.apache.commons.math3.special.Erf;

public class BestSyn3D {
	protected boolean [][][] kMap0;//
	protected double [][][] zMap0;
	protected int[][][] thrMap0;
	protected double thrZ;
	protected ArrayList<Double> PopUplist;// = new ArrayList<Integer>();
	
	public BestSyn3D(scanAllThres3D scanMap, paraP3D p, ArrayList<Double> In_PopUplist){
		//find region with highest zcore and test if it meet the requirement of fdr control
		extractBestSyn_fdr(scanMap, p, In_PopUplist);
	}
	public BestSyn3D(scanAllThres3D scanMap, paraP3D p){
		//find region with highest zscore
		extractBestSyn(scanMap, p);
	}
	public void extractBestSyn(scanAllThres3D scanMap, paraP3D p) {
		kMap0 = new boolean[scanMap.kMapx[0].length][scanMap.kMapx[0][0].length][scanMap.kMapx[0][0][0].length];// temp use
		zMap0 = new double[scanMap.zMapx[0].length][scanMap.zMapx[0][0].length][scanMap.zMapx[0][0][0].length];// temp use
		/* extractOneMap extract the best map*/
		int nThr = scanMap.kMapx.length;
		double[] maxZ = new double[nThr];
		BasicMath mBM = new BasicMath();
		ImageHandling IH = new ImageHandling();
		int nSynPos = 0;
		boolean [][][] zmask = new boolean[scanMap.zMapx[0].length][scanMap.zMapx[0][0].length][scanMap.zMapx[0][0][0].length];
		double C=-1000; //largest zscore
		int I = -10;//map number of the largest zscore
		for(int ii=0;ii<nThr;ii++){
			for(int i=0;i<kMap0.length;i++){
				for(int j=0;j<kMap0[0].length;j++){
					for(int kk=0;kk<kMap0[0][0].length;kk++){
						zmask[i][j][kk] = scanMap.zMapx[ii][i][j][kk]!=0;
						zMap0[i][j][kk] = scanMap.zMapx[ii][i][j][kk];
					}
				}
			}
			int max_kMap0 = mBM.matrix3DMax(scanMap.kMapx[ii]);
			zmask = IH.bwareaopen3D(zmask, p.min_size);
			ArrayList<Double> zMapMM = new ArrayList<Double>();
			for(int i=0;i<zMap0.length;i++){
				for(int j=0;j<zMap0[0].length;j++){
					for(int kk=0;kk<zMap0[0][0].length;kk++){
						if(!zmask[i][j][kk])
							zMap0[i][j][kk] = 0;
						/*unique zscores, this might wrongly remove ones with identical zscore, though minor influence*/
						if(zMap0[i][j][kk]!=0){
							boolean addflag = true;
							for(int k=zMapMM.size()-1;k>=0;k--){
								if(Math.abs(zMap0[i][j][kk]-zMapMM.get(k))<0.0001){
									addflag = false;
									break;
								}	
							}
							if(addflag)
								zMapMM.add(zMap0[i][j][kk]);
						}
					}
				}
			}
			if(zMapMM.size()==0){
				maxZ[ii] = -1000;
				continue;
			}
			maxZ[ii] = Collections.max(zMapMM);
			nSynPos += max_kMap0;
			if(C<maxZ[ii]){
				C = maxZ[ii];
				I = ii;
			}
		}
		thrZ = -1000;// use fdr control, therefore here thrZ is of no use
		/*PopUp the best ROI*/
		if (!PopUpROI(scanMap, I, C,p)) {
			thrZ = 0;
		}
	}
	public void extractBestSyn_fdr(scanAllThres3D scanMap, paraP3D p, ArrayList<Double> In_PopUplist) {
		/*extract the best map*/
		long startTime1=System.nanoTime();
		BasicMath mBM = new BasicMath();
		scanMap.CT.extractBestSyn(p);
		long endTime1=System.nanoTime();
    	System.out.println("extractBestSyn Running time: "+(endTime1-startTime1)/1e9); 
    	
    	startTime1=System.nanoTime();
		double C = scanMap.CT.C;
		int I = scanMap.CT.I;
		ArrayList<Double> Zscore_Vec = scanMap.CT.Zscore_Vec;
		/*FDR Control*/
		ArrayList<Double> pvalue = new ArrayList<Double>();
		for(int i=0;i<Zscore_Vec.size();i++){
			pvalue.add(mBM.zTop(Zscore_Vec.get(i)));
		}
		int total_len = pvalue.size()+In_PopUplist.size();
		double[] FDR_Con = new double[total_len];
		for(int i=0;i<total_len;i++)
			FDR_Con[i] = p.fdr*((double)(i+1))/total_len;
		if(p.fdr_den){
			for(int i=0;i<total_len;i++){
				FDR_Con[i] = FDR_Con[i]/Math.log(total_len);
			}
		}
		pvalue.addAll(In_PopUplist);
		Collections.sort(pvalue);// Ascend; add Collections.reverseOrder() for descend
		double OK_pvalue = 0;
		for(int i=pvalue.size()-1;i>=0;i--){
			if(pvalue.get(i)<=FDR_Con[i]){
				OK_pvalue = pvalue.get(i);
				break;
			}
		}
		double bestP = mBM.zTop(C);
		if(bestP<=OK_pvalue){
			thrZ = -1000;
			In_PopUplist.add(bestP);
		}
		else{
			
			thrZ = 1000;
		}
		endTime1=System.nanoTime();
    	System.out.println("prePare PopUp Running time: "+(endTime1-startTime1)/1e9); 
		/*PopUp the best ROI and save it in PopUplist*/
    	startTime1=System.nanoTime();
		if (!PopUpROI(scanMap, I, C,p)) {
			System.out.println("We did not found particle.");
			thrZ = mBM.pToZ(Collections.max(In_PopUplist));
		}
		PopUplist = new ArrayList<Double>(In_PopUplist);
		endTime1=System.nanoTime();
    	System.out.println("PopUp Running time: "+(endTime1-startTime1)/1e9); 
	}
	
	//PopUp the best ROI, if there is more, choose the larger one
	public boolean PopUpROI(scanAllThres3D scanMap, int I, double C,paraP3D p){
		int Nthr = (p.thr1-p.thr0)/p.thrg +1;
		int[] tt = new int[Nthr];
		for(int i = 0;i<Nthr;i++)
			tt[i] = p.thr0+i*p.thrg; 
		
		BasicMath mBM = new BasicMath();
		ImageHandling IH = new ImageHandling();
		int zSlice = scanMap.kMapx[0].length;
		int height = scanMap.kMapx[0][0].length;
		int width = scanMap.kMapx[0][0][0].length;
		kMap0 = new boolean[zSlice][height][width];//renew
		zMap0 = new double[zSlice][height][width];//
		thrMap0 = new int[zSlice][height][width];	
		if(I<0){
			return false;
		}
		double[][][] zMap0x = new double[zSlice][height][width];
		boolean foundOne = false;
		for(int i=0;i<zSlice;i++)
			for(int j=0;j<height;j++)
				for(int k=0;k<width;k++)
					zMap0x[i][j][k] = scanMap.zMapx[I][i][j][k];
		if(C>=thrZ){
			boolean[][][] K1 = new boolean[zSlice][height][width];
			for(int i=0;i<zSlice;i++)
				for(int j=0;j<height;j++)
					for(int k=0;k<width;k++)
						K1[i][j][k] = Math.abs(zMap0x[i][j][k]-C)<0.001;
			System.out.println("The particle we are detecting exist or not? "+mBM.matrix3DMax(K1));
			boolean[][][] tmask = IH.bwareaopen3D(K1, p.min_size);
			int[][][] tmask_cc = IH.bwlabel3D(tmask,26);
			int ccN = mBM.matrix3DMax(tmask_cc);
			// if more than one synapse share the same z-score, choose largest one
			if (ccN > 1) {
				foundOne = true;
				int[] roi_size = new int[ccN];
				for (int i = 0; i < tmask_cc.length; i++)
					for (int j = 0; j < tmask_cc[0].length; j++)
						for (int k = 0; k < tmask_cc[0][0].length; k++)
							if (tmask_cc[i][j][k] > 0)
								roi_size[tmask_cc[i][j][k] - 1]++;
				int max_size = mBM.vectorMax(roi_size);
				for (int k = 0; k < ccN; k++) {
					if (roi_size[k] == max_size) {
						for (int i = 0; i < tmask_cc.length; i++) {
							for (int j = 0; j < tmask_cc[0].length; j++)
								for (int zz = 0; zz < tmask_cc[0][0].length; zz++){
									if (tmask_cc[i][j][zz] == k + 1) {
										zMap0[i][j][zz] = zMap0x[i][j][zz];
										kMap0[i][j][zz] = true;
										thrMap0[i][j][zz] = tt[I];
									}
								}
						}
						break;
					}

				}
			}
			if (ccN == 1) {
				foundOne = true;
				for (int i = 0; i < tmask_cc.length; i++) {
					for (int j = 0; j < tmask_cc[0].length; j++) {
						for (int k = 0; k < tmask_cc[0][0].length; k++) {
							if (tmask_cc[i][j][k] != 0) {
								zMap0[i][j][k] = zMap0x[i][j][k];
								kMap0[i][j][k] = true;
								thrMap0[i][j][k] = tt[I];
							}
						}
					}
				}
			}
		}
		return foundOne;
	}

}
