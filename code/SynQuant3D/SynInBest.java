import java.util.ArrayList;
import java.util.Collections;

public class SynInBest {
	protected boolean foundOne;
	protected boolean [][] kMap2;
	protected double [][] zMap2;
	protected int [][] thrMap2;
	protected int [][][] kMapx2;//
	protected double [][][] zMapx2;
	
	//find if there is a better synapse inside the region with highest zscore
	public SynInBest(double[][] Gt, BestSyn bestSynMap, ParaP p, ParaQ q) {
		// TODO Auto-generated constructor stub
		boolean doneEle = false;
		foundOne = false;
		kMap2 = new boolean[Gt.length][Gt[0].length];
		zMap2 = new double[Gt.length][Gt[0].length];
		thrMap2 = new int[Gt.length][Gt[0].length];
		BasicMath mBM = new BasicMath();
		for(int i=0;i<Gt.length;i++){
			for(int j=0;j<Gt[0].length;j++){
				kMap2[i][j] = bestSynMap.kMap0[i][j];
				zMap2[i][j] = bestSynMap.zMap0[i][j];
				thrMap2[i][j] = bestSynMap.thrMap0[i][j];
			}
		}
		
		while(!doneEle){
			ArrayList<Integer[]> idx0 = new ArrayList<Integer[]>();
			int thr0 = 0;
			for(int i=0;i<kMap2.length;i++){
				for(int j=0;j<kMap2[0].length;j++){
					if(kMap2[i][j]){
						idx0.add(new Integer[] {i,j});
						thr0+=thrMap2[i][j];	
					}
				}
			}
			thr0 = thr0/idx0.size();
			if(idx0.size() > (p.max_ratio*q.Pix_per_synapse) && thr0<=p.thr1-p.thrg){
				boolean[][] Kmask2 = thresholding(Gt,thr0,kMap2);
				scanAll3 scanAll3Crop = new scanAll3(Gt, Kmask2, p, q,thr0);
				BestSyn bestSynMap2n = new BestSyn(scanAll3Crop, p);
				int MaxthreMap = mBM.matrix2DMax(bestSynMap2n.thrMap0);
				if(mBM.Allfalse(bestSynMap2n.kMap0) || MaxthreMap<=thr0)
					doneEle = true;
				else{
					//System.out.println("-- hit");
					foundOne = true;
					for(int i=0;i<Gt.length;i++){
						for(int j=0;j<Gt[0].length;j++){
							kMap2[i][j] = bestSynMap2n.kMap0[i][j];
							zMap2[i][j] = bestSynMap2n.zMap0[i][j];
							thrMap2[i][j] = bestSynMap2n.thrMap0[i][j];
						}
					}
				}
			}
			else
				doneEle = true;
		}
	}

	public boolean [][] thresholding(double[][] Gt, int thr, boolean [][] kMask){
		boolean [][] K1 = new boolean [Gt.length][Gt[0].length];
		for(int i=0;i<Gt.length;i++)
			for(int j=0;j<Gt[0].length;j++)
				K1[i][j] = (Gt[i][j]>thr) && kMask[i][j];
		return K1;
	}

}
