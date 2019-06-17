import java.util.ArrayList;
import java.util.Collections;

public class scanAllThres3D{
	protected int [][][][] kMapx;//regions numbering for each threshold
	protected double [][][][] zMapx;//zscore maps for each threshold
	protected int Nthr; //number of thresholds used
	protected int[] thrrg; // thresholds
	protected ComponentTree3D CT;

	//SCANALL scan all thresholds and build the zscore map of each threshold
	//Parameters are initialization for the first iteration
	//kMask is binary synapse map with background 1 and detected synaspe 0
	public scanAllThres3D(short[][] imArray, double[][][] G, double[][][] Gt, boolean[][][] kMask, paraP3D p, paraQ3D q) {
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		thrrg = new int[Nthr];
		for(int i = 0;i<Nthr;i++)
			thrrg[i] = p.thr0+i*p.thrg;

		byte[] usedN = new byte[imArray[0].length*imArray.length];
		// build the component tree
		CT = new ComponentTree3D(imArray,usedN, p, q);
		//component tree to zscore maps
		CT.cpt2map(p);
		kMapx = CT.kMapx;
		zMapx = CT.zMapx;
		BasicMath BM = new BasicMath();
		//for (int i=0;i<CT.kMapx.length;i++)
		//	System.out.println("Component Tree got #Particles: "+BM.matrix3DMax(CT.kMapx[i])+", Max zScore: "+BM.matrix3DMax(CT.zMapx[i]));

	}
	/**************Something is wrong*******************/
	//SCANALL scan all thresholds and find the best map in a cropped region
	public scanAllThres3D/*scanAll3Crop*/(double[][][] Gt, boolean[][][] Kmask,paraP3D p,paraQ3D q, int thr0) {
		//find the best synapse region
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		thrrg = new int[Nthr];
		for(int i = 0;i<Nthr;i++)
			thrrg[i] = p.thr0+i*p.thrg; 
		ImageHandling IH = new ImageHandling();

		int Nz = Gt.length;
		int Ny = Gt[0].length;
		int Nx = Gt[0][0].length;		
		
		kMapx = new int[Nthr][Nz][Ny][Nx];
		zMapx = new double[Nthr][Nz][Ny][Nx];
		int ofst = 5;
		for(int ii=0;ii<Nthr;ii++){
			int thr = thrrg[ii];
			int [][][] Kx = kMapx[ii];
			double [][][] Zx = zMapx[ii];
			if(!(thr>thr0)){
				continue;
			}
			boolean [][][] K = thresholding3D(Gt, thr,Kmask);
			//crop
			ArrayList<Integer> r = new ArrayList<Integer>(); // rows
			ArrayList<Integer> c = new ArrayList<Integer>(); // columns
			ArrayList<Integer> zd = new ArrayList<Integer>(); // z-direction
			for(int i=0;i<K.length;i++){
				for(int j=0;j<K[0].length;j++){
					for(int kk=0;kk<K[0][0].length;kk++){
						if(K[i][j][kk]){
							zd.add(i);
							r.add(j);
							c.add(kk);
						}
					}
				}
			}

			if(r.size()==0 || c.size()==0 || zd.size()==0){
				continue;
			}

			//Collections.sort(c);
			int rgzd_min = Math.max(0, Collections.min(zd)-ofst);
			int rgzd_max = Math.min(K.length-1, Collections.max(zd)+ofst);
			int rgr_min = Math.max(0, Collections.min(r)-ofst);
			int rgr_max = Math.min(K[0].length-1, Collections.max(r)+ofst);
			int rgc_min = Math.max(0, Collections.min(c)-ofst);
			int rgc_max = Math.min(K[0][0].length-1, Collections.max(c)+ofst);
			
			double[][][] Gt1 = new double[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][][] Kmask1 = new boolean[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][][] K1 = new boolean[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			double[][][] Zx1 = new double[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			int[][][] Kx1 = new int[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];

			for(int i=0;i<Gt1.length;i++){
				for(int j=0;j<Gt1[0].length;j++){
					for(int kk=0;kk<Gt1[0][0].length;kk++){
						Gt1[i][j][kk] = Gt[rgzd_min+i][rgr_min+j][rgc_min+kk];
						Kmask1[i][j][kk] = Kmask[rgzd_min+i][rgr_min+j][rgc_min+kk];
						K1[i][j][kk] = K[rgzd_min+i][rgr_min+j][rgc_min+kk];
						Zx1[i][j][kk] = Zx[rgzd_min+i][rgr_min+j][rgc_min+kk];
						Kx1[i][j][kk] = Kx[rgzd_min+i][rgr_min+j][rgc_min+kk];
					}
				}
			}

			boolean [][][] K1b = IH.imMorpho3D(K1, Kmask1, p.min_size);
			int [][][] K1b_cc = IH.bwlabel3D(K1b,26);
			//IH.maskDisplay3D(K1b_cc, "SynIn");
			int ccN1 = IH.NextLabel;
			if(ccN1==0)
			{
				continue;
			}
			int [] reg_size = new int[ccN1];
			double [] z_score = new double[ccN1];

			for(int i=0;i<K1b_cc.length;i++){
				for(int j=0;j<K1b_cc[0].length;j++){
					for(int k=0;k<K1b_cc[0][0].length;k++){
						if(K1b_cc[i][j][k]!=0){
							reg_size[K1b_cc[i][j][k]-1]++;
						}
					}
				}
			}
			for(int i=0;i<ccN1;i++){
				int[][] smaskIdx = new int[reg_size[i]][3];
				int smaskIdx_cnt = 0;
				for(int z=0;z<K1b_cc.length;z++) {
					for(int y=0;y<K1b_cc[0].length;y++){
						for(int x=0;x<K1b_cc[0][0].length;x++){
							if(K1b_cc[z][y][x]==i+1){
								smaskIdx[smaskIdx_cnt][0] = z;
								smaskIdx[smaskIdx_cnt][1] = y;
								smaskIdx[smaskIdx_cnt++][2] = x;
							}
						}
					}
				}

				if(reg_size[i]>p.max_size || reg_size[i]<p.min_size){
					z_score[i] = 0;
					//continue;
				}
				else{
					z_score[i] = IH.scanOneSyn3D(Gt1,K1b,Kmask1,smaskIdx,p,q);
				}
				for(int j=0;j<smaskIdx.length;j++){
					kMapx[ii][smaskIdx[j][0]+rgzd_min][smaskIdx[j][1]+rgr_min][smaskIdx[j][2]+rgc_min] = j+1;
				}
				if(z_score[i]>0)
					for(int j=0;j<smaskIdx.length;j++)
						zMapx[ii][smaskIdx[j][0]+rgzd_min][smaskIdx[j][1]+rgr_min][smaskIdx[j][2]+rgc_min] = z_score[i];
			}
		}
	}
	/*scanUpdtCrop: 3D version
	 *Once we detect a synapse, the other regions that overlapped with it needs to be updated. 
	 *For these regions, there are two conditions. One is this region is wholly included in the
	 *detected synapse. The other is the detected synapse is wholly included in this region. 
	 *We have Nthr zscore maps(zMapx) and index maps(kMapx). We scan all them to find such kind of regions.
	 *
	 *G is the orginal image
	 *Gt is the image after variance stablization
	 *bmask is the background mask
	 *idxUpdt saves the pixel positions of a detected synapse
	 *
	 */
	public void scanUpdtCrop(double[][][] G, double[][][] Gt, boolean[][][] bmask, ArrayList<Integer[]> idxUpdt, paraP3D p,
			paraQ3D q) {

		//q.ntry = 1 here, but what is the right range of neighboring size?
		ImageHandling IH = new ImageHandling();
		BasicMath mBM = new BasicMath();
		//long startTime2=System.nanoTime();
		//double bs=0;
		for(int ii=0;ii<Nthr;ii++){
			/**Find the region need to be updated**/
			int[][][] Kx = kMapx[ii];
			double [][][] Zx = zMapx[ii];
			int[] u0 = new int [idxUpdt.size()];
			ArrayList<Double> zu0 = new ArrayList<Double>();

			for(int i=0;i<idxUpdt.size();i++){
				u0[i] = Kx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]][idxUpdt.get(i)[2]];
				boolean zu0f = true;
				for (int j=0;j<zu0.size();j++){
					if(Math.abs(zu0.get(j)-Zx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]][idxUpdt.get(i)[2]])<0.0001){
						zu0f = false;
						break;
					}
				}
				if(zu0f)
					zu0.add(Zx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]][idxUpdt.get(i)[2]]);
			}
			int lu0 = mBM.unique(u0).length;
			if(lu0==1 && mBM.vectorMax(u0)==0){//no region,just background
				continue;
			}
			for(int i=0;i<zu0.size();i++){
				for(int j=0;j<CT.Zscore_Vec.size();j++){
					if(CT.levels.get(j)==ii & Math.abs(CT.Zscore_Vec.get(j)-zu0.get(i))<0.0001){
						CT.Zscore_Vec.remove(j);
						CT.levels.remove(j);
						break;
					}
				}
			}

			if(lu0>1){//condtion 1: the region is wholly included in the detected region
				for(int i=0;i<idxUpdt.size();i++){
					Kx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]][idxUpdt.get(i)[2]]=0;
					Zx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]][idxUpdt.get(i)[2]]=0;
					//if()
				}

				continue;
			}

			//condtion 2: the region contains detected region
			/*First thing: cut the region out*/
			int idx0 = (int) Math.round((double)mBM.vectorSum(u0)/u0.length);//there should be only one value in u0
			boolean[][][] maskUpdt = new boolean[Kx.length][Kx[0].length][Kx[0][0].length];
			for(int i=0;i<Kx.length;i++){
				for(int j=0;j<Kx[0].length;j++){
					for(int kk=0;kk<Kx[0][0].length;kk++){
						if(Kx[i][j][kk]==idx0){
							maskUpdt[i][j][kk] = true;
							Kx[i][j][kk] = 0;//
							Zx[i][j][kk] = 0;//
						}
						else
							maskUpdt[i][j][kk] = false;
					}
				}
			}
			if(true) { 
				// for those regions contain this highest z-score region, they are no longer convex
				// for such kind of region, we can safely remove them?
				continue;
			}
			int nSyn0 = mBM.matrix3DMax(Kx);
			int thr = thrrg[ii];
			boolean [][][] K = thresholding3D(Gt, thr,maskUpdt);
			//crop
			ArrayList<Integer> zd = new ArrayList<Integer>();
			ArrayList<Integer> r = new ArrayList<Integer>();
			ArrayList<Integer> c = new ArrayList<Integer>();
			for(int i=0;i<K.length;i++){
				for(int j=0;j<K[0].length;j++){
					for(int kk=0;kk<K[0][0].length;kk++){
						if(K[i][j][kk]){
							zd.add(i);
							r.add(j);
							c.add(kk);
						}
					}
				}
			}
			if(r.size()==0 || c.size()==0 || zd.size()==0){
				continue;
			}

			int ofst = 5;// cut the region with z+10 * width+10 * height+10
			int rgzd_min = Math.max(0, Collections.min(zd)-ofst);
			int rgzd_max = Math.min(K.length-1, Collections.max(zd)+ofst);
			int rgr_min = Math.max(0, Collections.min(r)-ofst);
			int rgr_max = Math.min(K[0].length-1, Collections.max(r)+ofst);
			int rgc_min = Math.max(0, Collections.min(c)-ofst);
			int rgc_max = Math.min(K[0][0].length-1, Collections.max(c)+ofst);
			
			double[][][] Gt1 = new double[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			double[][][] G1 = new double[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][][] bmask1 = new boolean[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][][] K1 = new boolean[rgzd_max-rgzd_min+1][rgr_max-rgr_min+1][rgc_max-rgc_min+1];

			for(int i=0;i<Gt1.length;i++){
				for(int j=0;j<Gt1[0].length;j++){
					for(int kk=0;kk<Gt1[0][0].length;kk++){
					G1[i][j][kk] = G[rgzd_min+i][rgr_min+j][rgc_min+kk];
					Gt1[i][j][kk] = Gt[rgzd_min+i][rgr_min+j][rgc_min+kk];
					bmask1[i][j][kk] = bmask[rgzd_min+i][rgr_min+j][rgc_min+kk];
					K1[i][j][kk] = K[rgzd_min+i][rgr_min+j][rgc_min+kk];
				}
				}
			}
			//The region is saved in Gt1, a croped rectangle from Gt
			boolean [][][] tmask1 = IH.imMorpho3D(K1, bmask1, p.min_size);
			int [][][] tmask1_cc = IH.bwlabel3D(tmask1,26);//saves the indexes of regions for updating
			int ccN1 = IH.NextLabel;// the number of regions need to be updated. !ccN1 is not needed to be 1
			if(ccN1==0)//after morphorlogy, no region left
			{
				continue;
			}
			double [] reg_mean = new double[ccN1];
			int [] reg_size = new int[ccN1];//save the region size
			double [] z_score = new double[ccN1];
			int [][] reg_cor = new int[ccN1][6]; //coordinates of regions. [zmin ymin xmin zmax ymax xmax]
			for(int i=0;i<ccN1;i++){
				reg_cor[i][0] = tmask1_cc.length;
				reg_cor[i][1] = tmask1_cc[0].length;
				reg_cor[i][2] = tmask1_cc[0][0].length;
			}

			for(int i=0;i<tmask1_cc.length;i++){
				for(int j=0;j<tmask1_cc[0].length;j++){
					for(int kk=0;kk<tmask1_cc[0][0].length;kk++){
						if(tmask1_cc[i][j][kk]!=0){
							reg_mean[tmask1_cc[i][j][kk]-1] += G1[i][j][kk];
							reg_size[tmask1_cc[i][j][kk]-1]++;
							if(i<reg_cor[tmask1_cc[i][j][kk]-1][0])
								reg_cor[tmask1_cc[i][j][kk]-1][0]=i;
							if(j<reg_cor[tmask1_cc[i][j][kk]-1][1])
								reg_cor[tmask1_cc[i][j][kk]-1][1]=j;
							if(kk<reg_cor[tmask1_cc[i][j][kk]-1][2])
								reg_cor[tmask1_cc[i][j][kk]-1][2]=kk;

							if(i>reg_cor[tmask1_cc[i][j][kk]-1][3])
								reg_cor[tmask1_cc[i][j][kk]-1][3]=i;
							if(j>reg_cor[tmask1_cc[i][j][kk]-1][4])
								reg_cor[tmask1_cc[i][j][kk]-1][4]=j;
							if(kk>reg_cor[tmask1_cc[i][j][kk]-1][5])
								reg_cor[tmask1_cc[i][j][kk]-1][5]=kk;
						}
					}
				}
			}
			/*second: Scan the regions and update their index and zscores in zMapx and kMapx*/
			for(int i=0;i<ccN1;i++){
				int[][] smaskIdx = new int[reg_size[i]][3];
				int smaskIdx_cnt = 0;
				for(int z=0;z<tmask1_cc.length;z++){
					for(int y=0;y<tmask1_cc[0].length;y++){
						for(int x=0;x<tmask1_cc[0][0].length;x++){
							if(tmask1_cc[z][y][x]==i+1){
								smaskIdx[smaskIdx_cnt][0] = z;
								smaskIdx[smaskIdx_cnt][1] = y;
								smaskIdx[smaskIdx_cnt++][2] = x;
							}
						}
					}
				}
				//test prior knowledge(constrains)
				if(reg_size[i]>p.max_size || reg_size[i]<p.min_size){
					z_score[i] = 0;
				}
				else if(reg_mean[i]/reg_size[i]<p.minIntensity){
					z_score[i] = 0;
				}
				else{
					int LH = reg_cor[i][4]-reg_cor[i][1]+1;
					int LW = reg_cor[i][5]-reg_cor[i][2]+1;
					double ratio = LH>LW? LH/LW: LW/LH;
					if(ratio>p.maxWHratio || reg_size[i]/(double)(LH*LW)<p.minfill){
						z_score[i] = 0;
					}
					else{//if all constrains satisfied, calculate z-score
						//long startTime=System.nanoTime();
						z_score[i] = IH.scanOneSyn3D(Gt1,tmask1,bmask1,smaskIdx,p,q);
						//long endTime=System.nanoTime();
						//bs = bs+(endTime-startTime)/1e9;
					}
				}
				//update the kMapx[ii] and zMapx[ii] with new zscore.
				for(int j=0;j<smaskIdx.length;j++){
					kMapx[ii][smaskIdx[j][0]+rgzd_min][smaskIdx[j][1]+rgr_min][smaskIdx[j][2]+rgc_min] = i+1+nSyn0;
				}
				if(z_score[i]!=0){
					for(int j=0;j<smaskIdx.length;j++)
						zMapx[ii][smaskIdx[j][0]+rgzd_min][smaskIdx[j][1]+rgr_min][smaskIdx[j][2]+rgc_min] = z_score[i];
					CT.Zscore_Vec.add(z_score[i]);
					CT.levels.add(ii);
				}
			}
		}
		//long endTime2=System.nanoTime();
		//System.out.println("ScanSyn timeï¼š"+bs);
		//System.out.println("Total timeï¼š"+(endTime2-startTime2)/1e9);

	}
	//threshold the image with a threshold and mask
	public boolean [][][] thresholding3D(double[][][] Gt, int thr, boolean [][][] kMask){
		boolean [][][] K1 = new boolean [Gt.length][Gt[0].length][Gt[0][0].length];
		for(int i=0;i<Gt.length;i++)
			for(int j=0;j<Gt[0].length;j++)
				for(int k=0;k<Gt[0][0].length;k++)
					K1[i][j][k] = (Gt[i][j][k]>thr) && kMask[i][j][k];
		return K1;
	}

}