import java.util.ArrayList;
import java.util.Collections;

public class SynInBest3D {
	protected boolean foundOne;
	protected boolean [][][] kMap2;
	protected double [][][] zMap2;
	protected int [][][] thrMap2;
	protected int [][][][] kMapx2;//
	protected double [][][][] zMapx2;

	//find if there is a better synapse inside the region with highest zscore
	public SynInBest3D(double[][][] Gt, BestSyn3D bestSynMap, paraP3D p, paraQ3D q, boolean containFlag) {
		// TODO Auto-generated constructor stub
		boolean doneEle = false;
		foundOne = false;
		kMap2 = new boolean[Gt.length][Gt[0].length][Gt[0][0].length];
		zMap2 = new double[Gt.length][Gt[0].length][Gt[0][0].length];
		thrMap2 = new int[Gt.length][Gt[0].length][Gt[0][0].length];
		BasicMath mBM = new BasicMath();
		for(int i=0;i<Gt.length;i++){
			for(int j=0;j<Gt[0].length;j++){
				for(int k=0;k<Gt[0][0].length;k++){
					kMap2[i][j][k] = bestSynMap.kMap0[i][j][k];
					zMap2[i][j][k] = bestSynMap.zMap0[i][j][k];
					thrMap2[i][j][k] = bestSynMap.thrMap0[i][j][k];
				}
			}
		}
		/* here for particle we believe it is not needed to detect particles from a larger particle
		 * SO in this 3D version, we add one more parameter "containFlag" for this function
		 * if containFlag==false, we do not search for better particle inside the old one
		 */
		if(containFlag) {
			while(!doneEle){
				ArrayList<Integer[]> idx0 = new ArrayList<Integer[]>();
				int thr0 = 0;
				for(int i=0;i<kMap2.length;i++){
					for(int j=0;j<kMap2[0].length;j++){
						for(int k=0;k<kMap2[0][0].length;k++){
							if(kMap2[i][j][k]){
								idx0.add(new Integer[] {i,j,k});
								thr0+=thrMap2[i][j][k];	
							}
						}
					}
				}
				thr0 = thr0/idx0.size();
				if(idx0.size() > (p.max_ratio*q.Pix_per_synapse) && thr0<=p.thr1-p.thrg){
					boolean[][][] Kmask2 = thresholding3D(Gt,thr0,kMap2);
					if(mBM.matrix3DMax(Kmask2)>=1) {
						//System.out.println("Got a Particle? "+mBM.matrix3DMax(Kmask2)+"Orginal Particle?"+mBM.matrix3DMax(kMap2));
						scanAllThres3D scanAll3Crop = new scanAllThres3D(Gt, Kmask2, p, q,thr0);
						for(int i=0;i<scanAll3Crop.kMapx.length;i++){
							System.out.println("ROI Z-socre: "+ scanAll3Crop.zMapx[i][4][79][24]);//+" * "+scanMap.zMapx[0][0].length);
						}
						BestSyn3D bestSynMap2n = new BestSyn3D(scanAll3Crop, p);
						int MaxthreMap = mBM.matrix3DMax(bestSynMap2n.thrMap0);
						if(mBM.Allfalse3D(bestSynMap2n.kMap0) || MaxthreMap<=thr0)
							doneEle = true;
						else{
							//System.out.println("-- hit");
							foundOne = true;
							for(int i=0;i<Gt.length;i++){
								for(int j=0;j<Gt[0].length;j++){
									for(int kk=0;kk<Gt[0][0].length;kk++){
										kMap2[i][j][kk] = bestSynMap2n.kMap0[i][j][kk];
										zMap2[i][j][kk] = bestSynMap2n.zMap0[i][j][kk];
										thrMap2[i][j][kk] = bestSynMap2n.thrMap0[i][j][kk];
									}
								}
							}
						}
					}
					else {
						doneEle = true;
					}
				}
				else {
					doneEle = true;
				}
			}
		}
	}

	public boolean [][][] thresholding3D(double[][][] Gt, int thr, boolean [][][] kMask){
		boolean [][][] K1 = new boolean [Gt.length][Gt[0].length][Gt[0][0].length];
		for(int i=0;i<Gt.length;i++)
			for(int j=0;j<Gt[0].length;j++)
				for(int k=0;k<Gt[0][0].length;k++)
					K1[i][j][k] = (Gt[i][j][k]>thr) && kMask[i][j][k];
		return K1;
	}

}
