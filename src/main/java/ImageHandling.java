import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.plugin.frame.RoiManager;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

//Congchao's functions for 2D image handling
public class ImageHandling {
	protected int NextLabel=0;
	protected int saveImCnt = 0;
	//Image morphology, open operation, fill holes 
	public boolean [][] imMorpho(boolean[][] K1,boolean[][] kMask, int min_size){
		boolean [][] K1a = bwareaopen(K1, min_size,8);
		//IH.IMwrite(K1a, "src/res/scanAll3_a_"+ii+".png");
		boolean [][] K11a = new boolean [K1a.length][K1a[0].length];
		for(int i=1;i<K1a.length-1;i++){
			for(int j=1;j<K1a[0].length-1;j++){
				if(!K1a[i][j]){
					int myInt = (K1a[i+1][j]) ? 1 : 0;
					myInt += (K1a[i-1][j]) ? 1 : 0;
					myInt += (K1a[i][j+1]) ? 1 : 0;
					myInt += (K1a[i][j-1]) ? 1 : 0;
					if(myInt>1)
						K11a[i][j]=true;
					else
						K11a[i][j] = false;
				}
				else
					K11a[i][j] = true;
			}
		}
		for(int i=0;i<K1a.length;i+=K1a.length-1){
			for(int j=0;j<K1a[0].length;j++){
				K11a[i][j] = K1a[i][j];
			}
		}
		for(int i=0;i<K1a.length;i++){
			for(int j=0;j<K1a[0].length;j+=K1a[0].length-1){
				K11a[i][j] = K1a[i][j];
			}
		}
		//IH.IMwrite(K11a, "src/res/scanAll3_a2_"+ii+".png");
		boolean [][] K1b = imfill(K11a);
		//IH.IMwrite(K1b, "src/res/scanAll3_b_"+ii+".png");
		for(int i=0;i<K1b.length;i++)
			for(int j=0;j<K1b[0].length;j++)
				K1b[i][j] = K1b[i][j] && kMask[i][j];
		//int [][] K1b_cc = bwlabel(K1b);
		//IMwrite(K1b,"src/res/Thre_"+saveImCnt+++".png");
		return K1b;
	}
	public boolean [][][] imMorpho3D(boolean[][][] K1,boolean[][][] kMask, int min_size){
		boolean [][][] K1a = K1;
		if (kMask != null)
			K1a = bwareaopen3D(K1, min_size);
		//IH.IMwrite(K1a, "src/res/scanAll3_a_"+ii+".png");
		boolean [][][] K11a = new boolean [K1a.length][K1a[0].length][K1a[0][0].length];
		for(int i=1;i<K1a.length-1;i++){
			for(int j=1;j<K1a[0].length-1;j++){
				for(int kk=1;kk<K1a[0][0].length-1;kk++){
					if(!K1a[i][j][kk]){
						int myInt = (K1a[i+1][j][kk]) ? 1 : 0;
						myInt += (K1a[i-1][j][kk]) ? 1 : 0;
						myInt += (K1a[i][j+1][kk]) ? 1 : 0;
						myInt += (K1a[i][j-1][kk]) ? 1 : 0;
						myInt += (K1a[i][j][kk-1]) ? 1 : 0;
						myInt += (K1a[i][j][kk+1]) ? 1 : 0;

						if(myInt>1)// should we use >=1?
							K11a[i][j][kk]=true;
						else
							K11a[i][j][kk] = false;
					}
					else
						K11a[i][j][kk] = true;
				}
			}
		}
		//for boundaries
		for(int i=0;i<K1a.length;i+=K1a.length-1){
			for(int j=0;j<K1a[0].length;j++){
				for(int kk=0;kk<K1a[0][0].length;kk++){
					K11a[i][j][kk] = K1a[i][j][kk];
				}
			}
		}
		for(int i=0;i<K1a.length;i++){
			for(int j=0;j<K1a[0].length;j+=K1a[0].length-1){
				for(int kk=0;kk<K1a[0][0].length;kk++){
					K11a[i][j][kk] = K1a[i][j][kk];
				}
			}
		}
		for(int i=0;i<K1a.length;i++){
			for(int j=0;j<K1a[0].length;j++){
				for(int kk=0;kk<K1a[0][0].length;kk+=K1a[0][0].length-1){
					K11a[i][j][kk] = K1a[i][j][kk];
				}
			}
		}
		//IH.IMwrite(K11a, "src/res/scanAll3_a2_"+ii+".png");
		boolean [][][] K1b = imfill3D(K11a);
		//IH.IMwrite(K1b, "src/res/scanAll3_b_"+ii+".png");
		//System.out.println("size of Matrix in imMorph3D: "+K1b.length+" "+K1b[0].length+" "+K1b[0][0].length); 
		//System.out.println("size of Matrix in imMorph3D: "+kMask.length+" "+kMask[0].length+" "+kMask[0][0].length); 
		if (kMask == null) {
			for(int i=0;i<K1b.length;i++)
				for(int j=0;j<K1b[0].length;j++)
					for(int kk=0;kk<K1b[0][0].length;kk++)
						K1b[i][j][kk] = K1b[i][j][kk];
		}else {
			for(int i=0;i<K1b.length;i++)
				for(int j=0;j<K1b[0].length;j++)
					for(int kk=0;kk<K1b[0][0].length;kk++)
						K1b[i][j][kk] = K1b[i][j][kk] && kMask[i][j][kk];
		}
		//int [][] K1b_cc = bwlabel(K1b);
		//IMwrite(K1b,"src/res/Thre_"+saveImCnt+++".png");
		return K1b;
	}
	//bwareaopen: the same as the function in matlab
	boolean [][] bwareaopen(boolean[][] BWIM, int min_size,int CONNECT){
		int [][] bwl = bwlabel(BWIM, CONNECT);
		//IH.IMwrite(K1a, "src/res/scanAll3_a_"+ii+".png");
		boolean [][] bw_out = new boolean [BWIM.length][BWIM[0].length];
		for(int i=0; i<bwl.length; i++)
			for(int j=0; j<bwl[0].length; j++)
				bw_out[i][j]=BWIM[i][j];
		//ArrayList<Integer> reg_size = new ArrayList<Integer>();
		//long [] reg_size = new long[1000000];
		int labelnum = 0;//need to check whether region number > 65535
		for(int i=0; i<bwl.length; i++) {
			for(int j=0; j<bwl[0].length; j++) {
				if(bwl[i][j]>labelnum)
					labelnum = bwl[i][j];
			}
		}
		long [] reg_size = new long[labelnum];
		for(int i=0; i<bwl.length; i++) {
			for(int j=0; j<bwl[0].length; j++) {
				if(bwl[i][j]>0)
					reg_size[bwl[i][j]-1]++;
			}
		}
		for(int k = 1;k<labelnum+1;k++){
			if(reg_size[k-1]<min_size){
				for(int i=0; i<bwl.length; i++) {
					for(int j=0; j<bwl[0].length; j++) {
						if(bwl[i][j]==k)
							bw_out[i][j]=false;
					}
				}	
			}
		}
		//bwl = bwlabel(BWIM);
		return bw_out;
	}
	//bwareaopen3D: the same as the function bwareaopen in matlab but for 3D
	boolean [][][] bwareaopen3D(boolean[][][] BWIM, int min_size){
		int [][][] bwl = bwlabel3D(BWIM, 26);
		//IH.IMwrite(K1a, "src/res/scanAll3_a_"+ii+".png");
		boolean [][][] bw_out = new boolean [BWIM.length][BWIM[0].length][BWIM[0][0].length];
		for(int i=0; i<bwl.length; i++)
			for(int j=0; j<bwl[0].length; j++)
				for(int k=0; k<bwl[0][0].length; k++)
					bw_out[i][j][k]=BWIM[i][j][k];
		//ArrayList<Integer> reg_size = new ArrayList<Integer>();
		//long [] reg_size = new long[1000000];
		int labelnum = 0;//need to check whether region number > 65535
		for(int i=0; i<bwl.length; i++) {
			for(int j=0; j<bwl[0].length; j++) 
				for(int k=0; k<bwl[0][0].length; k++){
					if(bwl[i][j][k]>labelnum)
						labelnum = bwl[i][j][k];
				}
		}
		long [] reg_size = new long[labelnum];
		for(int i=0; i<bwl.length; i++) {
			for(int j=0; j<bwl[0].length; j++) 
				for(int k=0; k<bwl[0][0].length; k++){
					if(bwl[i][j][k]>0)
						reg_size[bwl[i][j][k]-1]++;
				}
		}
		for(int k = 1;k<labelnum+1;k++){
			if(reg_size[k-1]<min_size){
				for(int i=0; i<bwl.length; i++) {
					for(int j=0; j<bwl[0].length; j++) {
						for(int zz = 0; zz<bwl[0][0].length; zz++){
							if(bwl[i][j][zz]==k)
								bw_out[i][j][zz]=false;
						}
					}
				}

			}
		}
		//bwl = bwlabel(BWIM);
		return bw_out;
	}
	//the same as function in imfill
	boolean [][] imfill(boolean[][] BWIM_in){
		//1-BWIM_in
		boolean[][] BWIM = new boolean [BWIM_in.length][BWIM_in[0].length];
		for(int i=0;i<BWIM.length;i++){
			for(int j=0;j<BWIM[0].length;j++){
				BWIM[i][j] = !BWIM_in[i][j];
			}
		}

		int [][] bwl = bwlabel(BWIM,8);
		boolean [][] bw_out = new boolean [BWIM.length][BWIM[0].length];
		ArrayList<Integer> free_reg = new ArrayList<Integer>();
		for(int i=0;i<BWIM.length;i+=BWIM.length-1){
			for(int j=0;j<BWIM[0].length;j++){
				bw_out[i][j] = BWIM_in[i][j];
				if(bwl[i][j]==0)
					continue;
				boolean saveflag = true;
				for(int k=0;k<free_reg.size();k++){
					if(bwl[i][j]==free_reg.get(k)){
						saveflag = false;
						break;
					}	
				}
				if(saveflag)
					free_reg.add(bwl[i][j]);
			}
		}
		for(int i=0;i<BWIM.length;i++){
			for(int j=0;j<BWIM[0].length;j+=BWIM[0].length-1){
				bw_out[i][j] = BWIM_in[i][j];
				if(bwl[i][j]==0)
					continue;
				boolean saveflag = true;
				for(int k=0;k<free_reg.size();k++){
					if(bwl[i][j]==free_reg.get(k)){
						saveflag = false;
						break;
					}	
				}
				if(saveflag)
					free_reg.add(bwl[i][j]);
			}
		}
		for(int i=1;i<BWIM.length-1;i++){
			for(int j=1;j<BWIM[0].length-1;j++){

				if(bwl[i][j]==0){
					bw_out[i][j] = true;
					continue;
				}
				boolean removeflag = true;
				for(int k=0;k<free_reg.size();k++){
					if(bwl[i][j]==free_reg.get(k)){
						removeflag = false; //the background is on edge
						break;
					}	
				}
				if(removeflag) //the background is NOT on edge
					bw_out[i][j] = true;
				else
					bw_out[i][j] = false;

			}
		}
		return bw_out;
	}
	//3D function of imfill
	boolean [][][] imfill3D(boolean[][][] BWIM_in){
		//1-BWIM_in
		boolean[][][] BWIM = new boolean [BWIM_in.length][BWIM_in[0].length][BWIM_in[0][0].length];
		for(int i=0;i<BWIM.length;i++){
			for(int j=0;j<BWIM[0].length;j++){
				for(int kk=0;kk<BWIM[0][0].length;kk++){
					BWIM[i][j][kk] = !BWIM_in[i][j][kk];
				}
			}
		}
		// the ROIs in (1-BWIM_in) is the holes we need to fill
		int [][][] bwl = bwlabel3D(BWIM,26);
		boolean [][][] bw_out = new boolean [BWIM.length][BWIM[0].length][BWIM_in[0][0].length];
		ArrayList<Integer> free_reg = new ArrayList<Integer>();
		// check if the roi is on edge
		// y-axis
		for(int i=0;i<BWIM.length;i+=BWIM.length-1){
			for(int j=0;j<BWIM[0].length;j++){
				for(int kk=0;kk<BWIM[0][0].length;kk++){
					bw_out[i][j][kk] = BWIM_in[i][j][kk];
					if(bwl[i][j][kk]==0)
						continue;
					boolean saveflag = true;
					for(int k=0;k<free_reg.size();k++){
						if(bwl[i][j][kk]==free_reg.get(k)){
							saveflag = false;
							break;
						}
					}
					if(saveflag)
						free_reg.add(bwl[i][j][kk]);
				}
			}
		}
		// x-axis
		for(int i=0;i<BWIM.length;i++){
			for(int j=0;j<BWIM[0].length;j+=BWIM[0].length-1){
				for(int kk=0;kk<BWIM[0][0].length;kk++){
					bw_out[i][j][kk] = BWIM_in[i][j][kk];
					if(bwl[i][j][kk]==0)
						continue;
					boolean saveflag = true;
					for(int k=0;k<free_reg.size();k++){
						if(bwl[i][j][kk]==free_reg.get(k)){
							saveflag = false;
							break;
						}
					}
					if(saveflag)
						free_reg.add(bwl[i][j][kk]);
				}
			}
		}
		// z-axis
		for(int i=0;i<BWIM.length;i++){
			for(int j=0;j<BWIM[0].length;j++){
				for(int kk=0;kk<BWIM[0][0].length;kk+=BWIM[0][0].length-1){
					bw_out[i][j][kk] = BWIM_in[i][j][kk];
					if(bwl[i][j][kk]==0)
						continue;
					boolean saveflag = true;
					for(int k=0;k<free_reg.size();k++){
						if(bwl[i][j][kk]==free_reg.get(k)){
							saveflag = false;
							break;
						}
					}
					if(saveflag)
						free_reg.add(bwl[i][j][kk]);
				}
			}
		}
		for(int i=1;i<BWIM.length-1;i++){
			for(int j=1;j<BWIM[0].length-1;j++){
				for(int kk=0;kk<BWIM[0][0].length;kk++){
					if(bwl[i][j][kk]==0){
						bw_out[i][j][kk] = true;
						continue;
					}
					boolean removeflag = true;
					for(int k=0;k<free_reg.size();k++){
						if(bwl[i][j][kk]==free_reg.get(k)){
							removeflag = false; //the background is on edge
							break;
						}
					}
					if(removeflag) //the background is NOT on edge
						bw_out[i][j][kk] = true;
					else
						bw_out[i][j][kk] = false;
				}
			}
		}
		return bw_out;
	}
	//the same as bwlabel in matlab
	int [][] bwlabel(boolean[][] BWIM, int CONNECT){
		int[][] labels = new int[BWIM.length][BWIM[0].length];

		NextLabel = 1;
		//First pass
		if(CONNECT==8){
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]==true & labels[i][j]==0) {
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							//System.out.println(""+loop_cnt);
							//if(loop_cnt==443)
							//	loop_cnt = 443;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							for(int ni=-1; ni<=1; ni++) {
								for(int nj=-1; nj<=1; nj++) {
									if(cur_i+ni<0 || cur_j+nj<0 || cur_i+ni>labels.length-1 || cur_j+nj>labels[0].length-1) {
										continue;
									}
									else {
										if(ni == 0 && nj == 0) continue;
										if(labels[cur_i+ni][cur_j+nj] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+ni][cur_j+nj]){
											roi_pos.add(new Integer[] {cur_i+ni,cur_j+nj});
											labels[cur_i+ni][cur_j+nj] = NextLabel;
										}
									}
								}
							}
							//deeplabel(BWIM, labels, i, j ,NextLabel);
						}
						NextLabel++;
					}
				}
			}
		}
		else{
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]==true & labels[i][j]==0) {
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							//System.out.println(""+loop_cnt);
							//if(loop_cnt==443)
							//	loop_cnt = 443;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							if(cur_i+1<BWIM.length){
								if(labels[cur_i+1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+1][cur_j]){
									roi_pos.add(new Integer[] {cur_i+1,cur_j});
									labels[cur_i+1][cur_j] = NextLabel;
								}
							}
							if(cur_j+1<BWIM[0].length){
								if(labels[cur_i][cur_j+1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j+1]){
									roi_pos.add(new Integer[] {cur_i,cur_j+1});
									labels[cur_i][cur_j+1] = NextLabel;
								}
							}
							if(cur_i-1>0){
								if(labels[cur_i-1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i-1][cur_j]){
									roi_pos.add(new Integer[] {cur_i-1,cur_j});
									labels[cur_i-1][cur_j] = NextLabel;
								}
							}
							if(cur_j-1>0){
								if(labels[cur_i][cur_j-1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j-1]){
									roi_pos.add(new Integer[] {cur_i,cur_j-1});
									labels[cur_i][cur_j-1] = NextLabel;
								}
							}
							//deeplabel(BWIM, labels, i, j ,NextLabel);
						}
						NextLabel++;
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}
	//bwlabel3D for 3D boolean data
	int [][][] bwlabel3D(boolean[][][] BWIM, int CONNECT){
		int[][][] labels = new int[BWIM.length][BWIM[0].length][BWIM[0][0].length];
		// we only consider 26 or 10 neighbors in this version, 10 means 8-xy neighbors and 2 z-neighbors
		NextLabel = 1;
		if (CONNECT == 26) {
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					for (int k=0; k<BWIM[0][0].length; k++) {
						if(BWIM[i][j][k]==true & labels[i][j][k]==0) {
							//Labels of neighbors
							ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
							roi_pos.add(new Integer[] {i,j,k});
							labels[i][j][k]=NextLabel;
							int loop_cnt = 0;
							while(true){
								if(loop_cnt>=roi_pos.size())
									break;
								int cur_i = roi_pos.get(loop_cnt)[0];
								int cur_j = roi_pos.get(loop_cnt)[1];
								int cur_k = roi_pos.get(loop_cnt++)[2];

								int i0 = Math.max(0, cur_i-1);
								int j0 = Math.max(0, cur_j-1);
								int k0 = Math.max(0, cur_k-1);
								int i2 = Math.min(BWIM.length-1, cur_i+1);
								int j2 = Math.min(BWIM[0].length-1, cur_j+1);
								int k2 = Math.min(BWIM[0][0].length-1, cur_k+1);
								for(int ni=i0; ni<=i2; ni++) {
									for(int nj=j0; nj<=j2; nj++) {
										for(int nk=k0;nk<=k2;nk++) {
											if(cur_i==ni & cur_j==nj & cur_k==nk) {
												continue;
											}
											else {
												if(labels[ni][nj][nk] == 0 & BWIM[cur_i][cur_j][cur_k] == BWIM[ni][nj][nk]){
													roi_pos.add(new Integer[] {ni,nj,nk});
													labels[ni][nj][nk] = NextLabel;
												}
											}
										}
									}
								}
							}
							NextLabel++;
						}
					}
				}
			}
		}else {
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					for (int k=0; k<BWIM[0][0].length; k++) {
						if(BWIM[i][j][k]==true & labels[i][j][k]==0) {
							//Labels of neighbors
							ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
							roi_pos.add(new Integer[] {i,j,k});
							labels[i][j][k]=NextLabel;
							int loop_cnt = 0;
							while(true){
								if(loop_cnt>=roi_pos.size())
									break;
								int cur_i = roi_pos.get(loop_cnt)[0];
								int cur_j = roi_pos.get(loop_cnt)[1];
								int cur_k = roi_pos.get(loop_cnt++)[2];

								int i0 = Math.max(0, cur_i-1);
								int j0 = Math.max(0, cur_j-1);
								int k0 = Math.max(0, cur_k-1);
								int i2 = Math.min(BWIM.length-1, cur_i+1);
								int j2 = Math.min(BWIM[0].length-1, cur_j+1);
								int k2 = Math.min(BWIM[0][0].length-1, cur_k+1);
								for(int ni=i0; ni<=i2; ni++) {
									if (ni == cur_i) {
										for(int nj=j0; nj<=j2; nj++) {
											for(int nk=k0;nk<=k2;nk++) {
												if(cur_j==nj & cur_k==nk) {
													continue;
												}
												else {
													if(labels[ni][nj][nk] == 0 & BWIM[cur_i][cur_j][cur_k] == BWIM[ni][nj][nk]){
														roi_pos.add(new Integer[] {ni,nj,nk});
														labels[ni][nj][nk] = NextLabel;
													}
												}
											}
										}
									}
									else {
										if(labels[ni][cur_j][cur_k] == 0 & BWIM[cur_i][cur_j][cur_k] == BWIM[ni][cur_j][cur_k]){
											roi_pos.add(new Integer[] {ni,cur_j,cur_k});
											labels[ni][cur_j][cur_k] = NextLabel;
										}
					
									}
								}
							}
							NextLabel++;
						}
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}

	//bwlabel for double data
	int [][] bwlabel(double[][] BWIM, int CONNECT, ArrayList<Double> zscores){
		int[][] labels = new int[BWIM.length][BWIM[0].length];

		NextLabel = 1;

		//First pass
		if(CONNECT==8){
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]!=0 & labels[i][j]==0) {
						zscores.add(BWIM[i][j]);
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							for(int ni=-1; ni<=1; ni++) {
								for(int nj=-1; nj<=1; nj++) {
									if(cur_i+ni<0 || cur_j+nj<0 || cur_i+ni>labels.length-1 || cur_j+nj>labels[0].length-1) {
										continue;
									}
									else {
										if(ni == 0 && nj == 0) continue;
										if(labels[cur_i+ni][cur_j+nj] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+ni][cur_j+nj]){
											roi_pos.add(new Integer[] {cur_i+ni,cur_j+nj});
											labels[cur_i+ni][cur_j+nj] = NextLabel;
										}
									}
								}
							}
						}
						NextLabel++;
					}
				}
			}
		}
		else{
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]!=0 & labels[i][j]==0) {
						zscores.add(BWIM[i][j]);
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							//System.out.println(""+loop_cnt);
							//if(loop_cnt==443)
							//	loop_cnt = 443;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							if(cur_i+1<BWIM.length){
								if(labels[cur_i+1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+1][cur_j]){
									roi_pos.add(new Integer[] {cur_i+1,cur_j});
									labels[cur_i+1][cur_j] = NextLabel;
								}
							}
							if(cur_j+1<BWIM[0].length){
								if(labels[cur_i][cur_j+1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j+1]){
									roi_pos.add(new Integer[] {cur_i,cur_j+1});
									labels[cur_i][cur_j+1] = NextLabel;
								}
							}
							if(cur_i-1>0){
								if(labels[cur_i-1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i-1][cur_j]){
									roi_pos.add(new Integer[] {cur_i-1,cur_j});
									labels[cur_i-1][cur_j] = NextLabel;
								}
							}
							if(cur_j-1>0){
								if(labels[cur_i][cur_j-1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j-1]){
									roi_pos.add(new Integer[] {cur_i,cur_j-1});
									labels[cur_i][cur_j-1] = NextLabel;
								}
							}
							//deeplabel(BWIM, labels, i, j ,NextLabel);
						}
						NextLabel++;
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}
	//bwlabel3D for 3D double data with z-scores
	int [][][] bwlabel3D(double[][][] BWIM, int CONNECT, ArrayList<Double> zscores){
		int[][][] labels = new int[BWIM.length][BWIM[0].length][BWIM[0][0].length];
		// we only consider 26 neighbors in this version, so CONNET = 26
		NextLabel = 1;
		for(int i=0; i<BWIM.length; i++) {
			for(int j=0; j<BWIM[0].length; j++) {
				for (int k=0; k<BWIM[0][0].length; k++) {
					if(BWIM[i][j][k]!=0 & labels[i][j][k]==0) {
						zscores.add(BWIM[i][j][k]);
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j,k});
						labels[i][j][k]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt)[1];
							int cur_k = roi_pos.get(loop_cnt++)[2];

							int i0 = Math.max(0, cur_i-1);
							int j0 = Math.max(0, cur_j-1);
							int k0 = Math.max(0, cur_k-1);
							int i2 = Math.min(BWIM.length-1, cur_i+1);
							int j2 = Math.min(BWIM[0].length-1, cur_j+1);
							int k2 = Math.min(BWIM[0][0].length-1, cur_k+1);
							for(int ni=i0; ni<=i2; ni++) {
								for(int nj=j0; nj<=j2; nj++) {
									for(int nk=k0;nk<=k2;nk++) {
										if(cur_i==ni & cur_j==nj & cur_k==nk) {
											continue;
										}
										else {
											if(labels[ni][nj][nk] == 0 & BWIM[cur_i][cur_j][cur_k] == BWIM[ni][nj][nk]){
												roi_pos.add(new Integer[] {ni,nj,nk});
												labels[ni][nj][nk] = NextLabel;
											}
										}
									}
								}
							}
						}
						NextLabel++;
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}
	//bwlabel3D for 3D double data with z-scores
	int [][][] bwlabel3D(double[][][] BWIM, int CONNECT){
		int[][][] labels = new int[BWIM.length][BWIM[0].length][BWIM[0][0].length];
		// we only consider 26 neighbors in this version, so CONNET = 26
		NextLabel = 1;
		for(int i=0; i<BWIM.length; i++) {
			for(int j=0; j<BWIM[0].length; j++) {
				for (int k=0; k<BWIM[0][0].length; k++) {
					if(BWIM[i][j][k]!=0 & labels[i][j][k]==0) {
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j,k});
						labels[i][j][k]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(loop_cnt>=roi_pos.size())
								break;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt)[1];
							int cur_k = roi_pos.get(loop_cnt++)[2];

							int i0 = Math.max(0, cur_i-1);
							int j0 = Math.max(0, cur_j-1);
							int k0 = Math.max(0, cur_k-1);
							int i2 = Math.min(BWIM.length-1, cur_i+1);
							int j2 = Math.min(BWIM[0].length-1, cur_j+1);
							int k2 = Math.min(BWIM[0][0].length-1, cur_k+1);
							for(int ni=i0; ni<=i2; ni++) {
								for(int nj=j0; nj<=j2; nj++) {
									for(int nk=k0;nk<=k2;nk++) {
										if(cur_i==ni & cur_j==nj & cur_k==nk) {
											continue;
										}
										else {
											if(labels[ni][nj][nk] == 0 & BWIM[cur_i][cur_j][cur_k] == BWIM[ni][nj][nk]){
												roi_pos.add(new Integer[] {ni,nj,nk});
												labels[ni][nj][nk] = NextLabel;
											}
										}
									}
								}
							}
						}
						NextLabel++;
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}
	void deeplabel(boolean[][] matrix, int[][] labels, int i, int j ,int NextLabel){
		for(int ni=-1; ni<=1; ni++) {
			for(int nj=-1; nj<=1; nj++) {
				if(i+ni<0 || j+nj<0 || i+ni>labels.length-1 || j+nj>labels[0].length-1) {
					continue;
				}
				else {
					if(ni == 0 && nj == 0) continue;
					if(labels[i+ni][j+nj] == 0 & matrix[i][j] == matrix[i+ni][j+nj]){
						labels[i+ni][j+nj] = NextLabel;
						deeplabel(matrix, labels, i+ni, j+nj ,NextLabel);
					}
				}
			}
		}
	}
	//score one region with several constrains in tmask, bmask
	double scanOneSyn(double[][] d/*Gt*/, boolean[][] tmask/*binary map>thr*/, boolean[][] bmask/*detected synapses=0*/, int[][] reg0, ParaP p, ParaQ q) {
		// TODO Auto-generated method stub
		return getOrderStatUnBalanced( d, reg0, bmask, tmask, q.var, q.ntry, p.mu, p.sigma, p.test );
	}
	double getOrderStatUnBalanced(double[][] d, int[][] reg0,boolean[][] bmask, boolean[][] tmask, double qVar,int qNtry,double[][] pMu,double[][] pSigma,String pTest ){

		double [] g1 = new double[reg0.length];
		int nPix = reg0.length;
		boolean[][] tmask1 = new boolean [tmask.length][tmask[0].length];
		for(int i=0;i<tmask.length;i++){
			for(int j=0;j<tmask[0].length;j++){
				tmask1[i][j] = tmask[i][j];
			}
		}
		int[][] smask = new int [tmask.length][tmask[0].length];//synapse candidate mask
		ArrayList<Integer[]>RegCoor = new ArrayList<Integer[]>();
		for(int i=0;i<reg0.length;i++){
			tmask1[reg0[i][0]][reg0[i][1]] = false; //tmask1 is the foreground mask(now remove the candidate first), got by thresholding thr0:thr1
			smask[reg0[i][0]][reg0[i][1]] = 1;
			RegCoor.add(new Integer[] {reg0[i][0], reg0[i][1]});
			g1[i] = d[reg0[i][0]][reg0[i][1]];
		}

		boolean done = false;
		int ntry = 1;
		int nNeibPrev = 0;
		ArrayList<Integer[]>NeibCoor = new ArrayList<Integer[]>();
		double [] g2 = new double[1];
		while(!done){
			for(int i=0;i<RegCoor.size();i++){
				for(int hs = -1;hs<2;hs++){
					for(int vs = -1;vs<2;vs++){
						/*if(hs==0 && vs==0) //the point itself needs to be tested too
							continue;*/ 
						int tmp_y = vs+RegCoor.get(i)[0];
						int tmp_x = hs+RegCoor.get(i)[1];
						if(tmp_y>=tmask.length || tmp_x>=tmask[0].length || tmp_y<0 || tmp_x<0)
							continue;
						if(smask[tmp_y][tmp_x]==0 && bmask[tmp_y][tmp_x] && !tmask1[tmp_y][tmp_x]){
							NeibCoor.add(new Integer[] {tmp_y, tmp_x});
						}
					}
				}
			}
			ArrayList<Integer[]> uniqueNeiCoor = new ArrayList<Integer[]>();
			if(NeibCoor.size()>0)
				uniqueNeiCoor.add(NeibCoor.get(0));
			for(int i=1;i<NeibCoor.size();i++){
				boolean add_flag = true;
				for(int j=0;j<uniqueNeiCoor.size();j++){
					if(uniqueNeiCoor.get(j)[0]==NeibCoor.get(i)[0] && uniqueNeiCoor.get(j)[1]==NeibCoor.get(i)[1]){
						add_flag = false;
						break;
					}
				}
				if(add_flag){
					uniqueNeiCoor.add(NeibCoor.get(i));
				}
			}
			int nNeib = uniqueNeiCoor.size();
			if(nNeib>=nPix || nNeibPrev==nNeib || ntry==qNtry){
				g2 = new double [nNeib];
				java.util.Iterator<Integer[]> it = uniqueNeiCoor.iterator();
				int g2_cnt = 0;
				while(it.hasNext()){
					//System.out.println(""+d[it.next()[0]][it.next()[1]]);
					Integer[] tmp_it = it.next();
					g2[g2_cnt++]=d[tmp_it[0]][tmp_it[1]];
					//
				}
				done = true;
			}
			else {
				java.util.Iterator<Integer[]> it = NeibCoor.iterator();
				while(it.hasNext()){
					RegCoor.add(it.next());
				}
			}
			ntry++;
			nNeibPrev = nNeib;
		}
		//test statistics: g2 save the neighboring pixel coordinates, g1 is the pixels inside ROI
		BasicMath mBM = new BasicMath();
		double t0 = mBM.sampleMean(g1)-mBM.sampleMean(g2);
		int M = g1.length;
		int N = g2.length;
		if(M < 4)
			M = 4;

		if(N<=3 || N < (M/10))
			return 0;
		if(M>100)
			M = 100;
		if(N>100)
			N = 99;
		double mu = pMu[M-1][N-1];
		double sigma = pSigma[M-1][N-1];
		mu = mu*Math.sqrt(qVar);
		sigma = sigma*Math.sqrt(qVar);

		double zScore = (t0-mu)/sigma;
		//t0 = t0/Math.sqrt(qVar);

		return zScore;
	}
	// zscore of a single 3D particle
	double scanOneSyn3D(double[][][] d/*Gt*/, boolean[][][] tmask/*binary map>thr*/, boolean[][][] bmask/*detected synapses=0*/, int[][] reg0, paraP3D p, paraQ3D q) {
		// TODO Auto-generated method stub
		return getOrderStatUnBalanced3D( d, reg0, bmask, tmask, q.var, q.ntry, p.mu, p.sigma, p.test , p.max_size, p.min_size);
	}
	// zscore of a single 3D particle
	private double getOrderStatUnBalanced3D(double[][][] d, int[][] reg0,boolean[][][] bmask, boolean[][][] tmask, double qVar,int qNtry,double[][] pMu,double[][] pSigma,String pTest, int pMax, int pMin) {
		double [] g1 = new double[reg0.length];
		int nPix = reg0.length;
		boolean[][][] tmask1 = new boolean [tmask.length][tmask[0].length][tmask[0][0].length];
		for(int i=0;i<tmask.length;i++){
			for(int j=0;j<tmask[0].length;j++){
				for(int k=0;k<tmask[0][0].length;k++){
					tmask1[i][j][k] = tmask[i][j][k];
				}
			}
		}
		int[][][] smask = new int [tmask.length][tmask[0].length][tmask[0][0].length];//3D particle candidate mask
		ArrayList<Integer[]>RegCoor = new ArrayList<Integer[]>();
		for(int i=0;i<reg0.length;i++){
			tmask1[reg0[i][0]][reg0[i][1]][reg0[i][2]] = false; //tmask1 is the foreground mask(now remove the candidate first), got by thresholding thr0:thr1
			smask[reg0[i][0]][reg0[i][1]][reg0[i][2]] = 1;
			RegCoor.add(new Integer[] {reg0[i][0], reg0[i][1], reg0[i][2]});
			g1[i] = d[reg0[i][0]][reg0[i][1]][reg0[i][2]];
		}

		boolean done = false;
		int ntry = 1;
		int nNeibPrev = 0;
		ArrayList<Integer[]>NeibCoor = new ArrayList<Integer[]>();
		double [] g2 = new double[1];
		while(!done){
			for(int i=0;i<RegCoor.size();i++){
				for(int hs = -1;hs<2;hs++){
					for(int vs = -1;vs<2;vs++){
						for(int zs = -1;zs<2;zs++) {
							/*if(hs==0 && vs==0) //the point itself needs to be tested too
							continue;*/ 
							int tmp_z = zs+RegCoor.get(i)[0];
							int tmp_y = vs+RegCoor.get(i)[1];
							int tmp_x = hs+RegCoor.get(i)[2];
							if(tmp_z>=tmask.length || tmp_y>=tmask[0].length || tmp_x>=tmask[0][0].length || tmp_z<0 || tmp_y<0 || tmp_x<0)
								continue;
							if(smask[tmp_z][tmp_y][tmp_x]==0 && bmask[tmp_z][tmp_y][tmp_x] && !tmask1[tmp_z][tmp_y][tmp_x]){
								NeibCoor.add(new Integer[] {tmp_z, tmp_y, tmp_x});
							}
						}
					}
				}
			}
			ArrayList<Integer[]> uniqueNeiCoor = new ArrayList<Integer[]>();
			if(NeibCoor.size()>0)
				uniqueNeiCoor.add(NeibCoor.get(0));
			for(int i=1;i<NeibCoor.size();i++){
				boolean add_flag = true;
				for(int j=0;j<uniqueNeiCoor.size();j++){
					if(uniqueNeiCoor.get(j)[0]==NeibCoor.get(i)[0] && uniqueNeiCoor.get(j)[1]==NeibCoor.get(i)[1] && uniqueNeiCoor.get(j)[2]==NeibCoor.get(i)[2]){
						add_flag = false;
						break;
					}
				}
				if(add_flag){
					uniqueNeiCoor.add(NeibCoor.get(i));
				}
			}
			int nNeib = uniqueNeiCoor.size();
			if(nNeib>=nPix || nNeibPrev==nNeib || ntry==qNtry){
				g2 = new double [nNeib];
				java.util.Iterator<Integer[]> it = uniqueNeiCoor.iterator();
				int g2_cnt = 0;
				while(it.hasNext()){
					//System.out.println(""+d[it.next()[0]][it.next()[1]]);
					Integer[] tmp_it = it.next();
					g2[g2_cnt++]=d[tmp_it[0]][tmp_it[1]][tmp_it[2]];
					//
				}
				done = true;
			}
			else {
				java.util.Iterator<Integer[]> it = NeibCoor.iterator();
				while(it.hasNext()){
					RegCoor.add(it.next());
				}
			}
			ntry++;
			nNeibPrev = nNeib;
		}
		//test statistics: g2 save the neighboring pixel coordinates, g1 is the pixels inside ROI
		BasicMath mBM = new BasicMath();
		double t0 = mBM.sampleMean(g1)-mBM.sampleMean(g2);
		int M = g1.length;
		int N = g2.length;

		// test the size constrain
		if(M < pMin)
			M = pMin;
		if(N<=pMin/2 || N < (M/10))
			return 0;
		if(M>pMax)
			M = pMax;
		if(N>pMax)
			N = pMax;
		double mu = pMu[M-1][N-1];
		double sigma = pSigma[M-1][N-1];
		mu = mu*Math.sqrt(qVar);
		sigma = sigma*Math.sqrt(qVar);

		double zScore = (t0-mu)/sigma;
		//t0 = t0/Math.sqrt(qVar);

		return zScore;
	}
	//output for double matrix
	public void IMwrite(double[][] imarray, String filepath){
		///////////////////// 	  Output intermediate result image (for inspection)     ///////////////////////

		int[][] showImg = new int[imarray.length][imarray[0].length];
		int numRows = showImg.length;
		int numCols = showImg[0].length;
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				if(imarray[y][x]<0)
					showImg[y][x] = 0;
				else{
					showImg[y][x] = (int)(imarray[y][x]);
				}
			}
		}

		int pixel;
		BufferedImage bufferedImage = new BufferedImage(numCols, numRows, BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster raster = bufferedImage.getRaster();
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				pixel = showImg[y][x];
				//pixel =  0 + (pixel << 8)  + (0 <<16);
				//		            System.out.println("The pixel in Matrix: "+pixel);
				raster.setSample(x, y, 0, pixel);
				//		            System.out.println("The pixel in BufferedImage: "+bufferedImage.getRGB(x, y));
			}
		}
		System.out.println("BufferedImage got.");

		File outputfile = new File(filepath);
		try {
			ImageIO.write(bufferedImage, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("BufferedImage saved.");

		//////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	//output for bool matrix
	public void IMwrite(boolean[][] imarray, String filepath){
		///////////////////// 	  Output intermediate result image (for inspection)     ///////////////////////

		int[][] showImg = new int[imarray.length][imarray[0].length];
		int numRows = showImg.length;
		int numCols = showImg[0].length;
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				if(imarray[y][x])
					showImg[y][x] = 255;
				else
					showImg[y][x] = 0;
			}
		}

		int pixel;
		BufferedImage bufferedImage = new BufferedImage(numCols, numRows, BufferedImage.TYPE_INT_RGB);
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				pixel = showImg[y][x];
				pixel =  pixel + (pixel << 8)  + (pixel <<16);
				//		            System.out.println("The pixel in Matrix: "+pixel);
				bufferedImage.setRGB(x, y, pixel);
				//		            System.out.println("The pixel in BufferedImage: "+bufferedImage.getRGB(x, y));
			}
		}
		System.out.println("BufferedImage got.");

		File outputfile = new File(filepath);
		try {
			ImageIO.write(bufferedImage, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("BufferedImage saved.");

		//////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	//output for int matrix
	public void IMwrite(int[][] imarray, String filepath){
		///////////////////// 	  Output intermediate result image (for inspection)     ///////////////////////

		int[][] showImg = new int[imarray.length][imarray[0].length];
		int numRows = showImg.length;
		int numCols = showImg[0].length;
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				if(imarray[y][x]>255)
					showImg[y][x] = 0;
				else
					showImg[y][x] = imarray[y][x];
			}
		}

		int pixel;
		BufferedImage bufferedImage = new BufferedImage(numCols, numRows, BufferedImage.TYPE_INT_RGB);
		for (int x = 0; x < numCols; x++) {
			for (int y = 0; y < numRows; y++) {
				pixel = showImg[y][x];
				pixel =  pixel + (pixel << 8)  + (pixel <<16);
				//		            System.out.println("The pixel in Matrix: "+pixel);
				bufferedImage.setRGB(x, y, pixel);
				//		            System.out.println("The pixel in BufferedImage: "+bufferedImage.getRGB(x, y));
			}
		}
		System.out.println("BufferedImage got.");

		File outputfile = new File(filepath);
		try {
			ImageIO.write(bufferedImage, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("BufferedImage saved.");

		//////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	//For grow neurite
	boolean[][] imopen(boolean[][] inputIm, int connect){
		boolean[][] outIm = new boolean[inputIm.length][inputIm[0].length];
		if(connect==4)
			for(int y =0;y<inputIm.length;y++){
				for(int x=0;x<inputIm[0].length;x++){
					if(inputIm[y][x]){
						outIm[y][x] = true;
						if(y>0)
							outIm[y-1][x] = true;
						if(x>0)
							outIm[y][x-1] = true;
						if(y<inputIm.length-1)
							outIm[y+1][x] = true;
						if(x<inputIm[0].length-1)
							outIm[y][x+1] = true;
					}
					else
						outIm[y][x] = false;
				}
			}
		if(connect==8)
			for(int y =0;y<inputIm.length;y++){
				for(int x=0;x<inputIm[0].length;x++){
					if(inputIm[y][x]){
						outIm[y][x] = true;
						if(y>0)
							outIm[y-1][x] = true;
						if(x>0)
							outIm[y][x-1] = true;
						if(y<inputIm.length-1)
							outIm[y+1][x] = true;
						if(x<inputIm[0].length-1)
							outIm[y][x+1] = true;

						if(y>0 & x>0)
							outIm[y-1][x-1] = true;
						if(y>0 & x<inputIm[0].length-1)
							outIm[y-1][x+1] = true;
						if(y<inputIm.length-1 & x>0)
							outIm[y+1][x-1] = true;
						if(y<inputIm.length-1 & x<inputIm[0].length-1)
							outIm[y+1][x+1] = true;
					}
					else
						outIm[y][x] = false;
				}
			}
		return outIm;
	}
	boolean[][] Skeleton(boolean[][] inIm){

		ByteProcessor outIP =  new ByteProcessor(inIm[0].length, inIm.length);

		for (int i = 0; i < inIm[0].length; i++) {
			for(int j = 0;j<inIm.length;j++){
				// put channel values in an integer
				if(inIm[j][i])
					outIP.set(i, j, 0);
				else
					outIP.set(i, j, 255);
			}
		}
		BinaryProcessor bi_IP = new BinaryProcessor(outIP);
		bi_IP.skeletonize();
		boolean[][] outIm = new boolean[inIm.length][inIm[0].length]; 
		for (int i = 0; i < inIm[0].length; i++) {
			for(int j = 0;j<inIm.length;j++){
				if(bi_IP.get(i, j)>0)
					outIm[j][i] = false;
				else
					outIm[j][i] = true;
			}
		}
		return outIm;
	}
	ArrayList<Integer[]> branchpts(boolean[][] skel){
		ArrayList<Integer[]> branchPts = new ArrayList<Integer[]>();
		for(int y=0;y<skel.length;y++){
			for(int x=0;x<skel[0].length;x++){

				if(!skel[y][x])
					continue;
				int neiCnt = 0;
				for(int ni=-1; ni<=1; ni++) {
					for(int nj=-1; nj<=1; nj++) {
						if(x+ni<0 || y+nj<0 || x+ni>skel[0].length-1 || y+nj>skel.length-1) 
							continue;
						else{
							if(ni == 0 && nj == 0) continue;
							if(skel[y+nj][x+ni]) neiCnt++;
						}
					}
				}
				if(neiCnt>2)
					branchPts.add(new Integer[] {y,x});

			}
		}
		return branchPts;
	}
	ArrayList<Integer[]> endpts(boolean[][] skel){
		ArrayList<Integer[]> endPts = new ArrayList<Integer[]>();
		for(int y=0;y<skel.length;y++){
			for(int x=0;x<skel[0].length;x++){

				if(!skel[y][x])
					continue;
				int neiCnt = 0;
				for(int ni=-1; ni<=1; ni++) {
					for(int nj=-1; nj<=1; nj++) {
						if(x+ni<0 || y+nj<0 || x+ni>skel[0].length-1 || y+nj>skel.length-1) 
							continue;
						else{
							if(ni == 0 && nj == 0) continue;
							if(skel[y+nj][x+ni]) neiCnt++;
						}
					}
				}
				if(neiCnt==1)
					endPts.add(new Integer[] {y,x});

			}
		}
		return endPts;
	}
	ArrayList<Integer[]> branchendpts(boolean[][] skel){
		ArrayList<Integer[]> endPts = new ArrayList<Integer[]>();
		for(int y=0;y<skel.length;y++){
			for(int x=0;x<skel[0].length;x++){

				if(!skel[y][x])
					continue;
				int neiCnt = 0;
				for(int ni=-1; ni<=1; ni++) {
					for(int nj=-1; nj<=1; nj++) {
						if(x+ni<0 || y+nj<0 || x+ni>skel[0].length-1 || y+nj>skel.length-1) 
							continue;
						else{
							if(ni == 0 && nj == 0) continue;
							if(skel[y+nj][x+ni]) neiCnt++;
						}
					}
				}
				if(neiCnt==1 || neiCnt>2)
					endPts.add(new Integer[] {y,x});

			}
		}
		return endPts;
	}
	//bwlabel with regularization
	public int[][] Regbwlabel(boolean[][] BWIM, int CONNECT, int regularization) {
		int[][] labels = new int[BWIM.length][BWIM[0].length];

		NextLabel = 1;

		//First pass
		if(CONNECT==8){
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]==true & labels[i][j]==0) {
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(roi_pos.size()>regularization || loop_cnt>=roi_pos.size())
								break;

							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							for(int ni=-1; ni<=1; ni++) {
								for(int nj=-1; nj<=1; nj++) {
									if(cur_i+ni<0 || cur_j+nj<0 || cur_i+ni>labels.length-1 || cur_j+nj>labels[0].length-1) {
										continue;
									}
									else {
										if(ni == 0 && nj == 0) continue;
										if(labels[cur_i+ni][cur_j+nj] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+ni][cur_j+nj]){
											roi_pos.add(new Integer[] {cur_i+ni,cur_j+nj});
											labels[cur_i+ni][cur_j+nj] = NextLabel;
										}
									}
								}
							}
						}
						NextLabel++;
					}
				}
			}
		}
		else{
			for(int i=0; i<BWIM.length; i++) {
				for(int j=0; j<BWIM[0].length; j++) {
					if(BWIM[i][j]==true & labels[i][j]==0) {
						//Labels of neighbors
						ArrayList<Integer[]> roi_pos = new ArrayList<Integer[]>();
						roi_pos.add(new Integer[] {i,j});
						labels[i][j]=NextLabel;
						int loop_cnt = 0;
						while(true){
							if(roi_pos.size()>regularization || loop_cnt>=roi_pos.size())
								break;
							//System.out.println(""+loop_cnt);
							//if(loop_cnt==443)
							//	loop_cnt = 443;
							int cur_i = roi_pos.get(loop_cnt)[0];
							int cur_j = roi_pos.get(loop_cnt++)[1];
							if(cur_i+1<BWIM.length){
								if(labels[cur_i+1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i+1][cur_j]){
									roi_pos.add(new Integer[] {cur_i+1,cur_j});
									labels[cur_i+1][cur_j] = NextLabel;
								}
							}
							if(cur_j+1<BWIM[0].length){
								if(labels[cur_i][cur_j+1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j+1]){
									roi_pos.add(new Integer[] {cur_i,cur_j+1});
									labels[cur_i][cur_j+1] = NextLabel;
								}
							}
							if(cur_i-1>0){
								if(labels[cur_i-1][cur_j] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i-1][cur_j]){
									roi_pos.add(new Integer[] {cur_i-1,cur_j});
									labels[cur_i-1][cur_j] = NextLabel;
								}
							}
							if(cur_j-1>0){
								if(labels[cur_i][cur_j-1] == 0 & BWIM[cur_i][cur_j] == BWIM[cur_i][cur_j-1]){
									roi_pos.add(new Integer[] {cur_i,cur_j-1});
									labels[cur_i][cur_j-1] = NextLabel;
								}
							}
							//deeplabel(BWIM, labels, i, j ,NextLabel);
						}
						NextLabel++;
					}
				}
			}
		}
		NextLabel--;
		return labels;
	}
	public void maskDisplay3D(int[][][] particleIdx, String title) {
		int width = particleIdx[0][0].length;
		int height = particleIdx[0].length;
		int zSlice = particleIdx.length;
		ImagePlus outputImp = IJ.createHyperStack(title, width, height, 1, zSlice, 1,8/*bitdepth*/);

		ImageStack curStack = outputImp.getStack();
		for (int i = 0; i < width; i++) {
			for(int j = 0;j<height;j++){
				for(int k = 0;k<zSlice;k++){
					int stackNum = outputImp.getStackIndex(1, k, 1);
					// put channel values in an integer
					if(particleIdx[k][j][i]>0) {
						//System.out.println("Number of Synapse: "+curStack.getPixels(stackNum + 1).getClass());
						((byte[]) curStack.getPixels(stackNum + 1))[i + j * width] =  (byte)(particleIdx[k][j][i]);//255
					}
				}
			}
		}		
		//stack.getBitDepth());
		outputImp.show();
	}
	public void maskDisplay3D(boolean[][][] particleIdx, String title) {
		int width = particleIdx[0][0].length;
		int height = particleIdx[0].length;
		int zSlice = particleIdx.length;
		ImagePlus outputImp = IJ.createHyperStack(title, width, height, 1, zSlice, 1,8/*bitdepth*/);

		ImageStack curStack = outputImp.getStack();
		for (int i = 0; i < width; i++) {
			for(int j = 0;j<height;j++){
				for(int k = 0;k<zSlice;k++){
					int stackNum = outputImp.getStackIndex(1, k, 1);
					// put channel values in an integer
					if(particleIdx[k][j][i]) {
						//System.out.println("Number of Synapse: "+curStack.getPixels(stackNum + 1).getClass());
						((byte[]) curStack.getPixels(stackNum + 1))[i + j * width] =  (byte)(255);//255
					}
				}
			}
		}		
		//stack.getBitDepth());
		outputImp.show();
	}
	public void DisplayROI(int nSyn0,int height,int width,int[][] SynR1Idx, ImagePlus imp,String frametitle){
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
		ImagePlus newimp = NewImage.createByteImage (imtitle, width, height, 1,/*stack=1*/ 
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
		double wandVal = 0.0001;
		for(int i=0;i<nSyn0;i++){
			roi2fiu[roi_cnt] = i;
			roi_cnt++;
			Wand w = new Wand(ip);
			w.autoOutline(roi_pt1[i][1],roi_pt1[i][0],wandVal,Wand.EIGHT_CONNECTED); 
			//System.out.println("roi size:"+w.npoints);
			if (w.npoints>0) { // we have an roi from the wand... 
				Roi roi = new PolygonRoi(w.xpoints, w.ypoints, w.npoints, Roi.TRACED_ROI);
				imp.setRoi(roi);
				manager.addRoi(roi);
			}
		}
		manager.setTitle(frametitle);
		manager.runCommand("show all");
		manager.setSize(300, 400);
		
	}

	/***********
	 * functions for noise estimation and variance stablization
	 * *************/
	/// estimate noise level: the estimated parameters are saved in qq
	//// More details can be found: http://www.cs.tut.fi/~foi/sensornoise.html
	void getNoise(double[][] Org_Im, paraQ3D qq) {
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
		double phi22 = 0.250141019881;// sum(conv2(lp.',hp).^2); //
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
		if(qq.pgEst==null) {
			qq.pgEst = new double [pEst0.length];
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = pEst0[i];
		}
		else {
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = qq.pgEst[i]+pEst0[i];
		}
		System.out.println("alpha: "+pEst0[0]+" beta:"+pEst0[1]);
		//// get the variance after stabilization, this can be done by selecting a point from the curve y = alpha*x+beta
		double alpha = pEst0[0];
		double sigma2 = Math.max(pEst0[1],2e-4); // this is the var of gausssian noise
		double t0 = 2/alpha*Math.sqrt(alpha*0+3/8*alpha*alpha+sigma2);
		double t1 = 2/alpha*Math.sqrt(alpha*1+3/8*alpha*alpha+sigma2);
		double var1 = 1/((t1-t0)*(t1-t0));
		if (Double.isNaN(qq.var))
			qq.var = var1*255*255;
		else
			qq.var = qq.var + var1*255*255;
	}
	// Variance stabilization assuming the parameters are already known (noise estimation already done)
	public double [][] getNoiseVSTwithAlpha(double[][] Org_Im, paraQ3D qq){
		double alpha = qq.pgEst[0];
		double sigma2 = qq.pgEst[1]; // this is the var of gausssian noise
		double t0 = 2/alpha*Math.sqrt(Math.max(0, alpha*0+3/8*alpha*alpha+sigma2));
		double t1 = 2/alpha*Math.sqrt(alpha*1+3/8*alpha*alpha+sigma2);

		double[][] outputGt = new double[Org_Im.length][Org_Im[0].length];
		for(int i=0;i<Org_Im.length;i++){
			for(int j=0;j<Org_Im[0].length;j++){
				outputGt[i][j] = 2/alpha*Math.sqrt(Math.max(0, alpha * Org_Im[i][j]/255+3/8*alpha*alpha+sigma2));
				outputGt[i][j] = 255*(outputGt[i][j]-t0)/(t1-t0);
			}
		}
		double var1 = 1/((t1-t0)*(t1-t0));
		System.out.println("alpha: "+qq.pgEst[0]+" beta:"+qq.pgEst[1]+" Variance:"+var1);

		return outputGt;
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
		System.out.println("alpha: "+pEst0[0]+" beta:"+pEst0[1]+" Variance:"+var1);

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
	//// simplified noise estimation Boyu Lyu's idea:
	void getNoiseRobust(double[][] Org_Im, paraQ3D qq) {
		// step 1: input is 0-255 image, normalize data to to 0-1
		double[][] double_G = new double[Org_Im.length][Org_Im[0].length];
		for (int i = 0; i < Org_Im.length; i++)
			for (int j = 0; j < Org_Im[0].length; j++)
				double_G[i][j] = Org_Im[i][j] / 255;

		// step 2. smoothing
//		double[][] avgfilter = { {1.0/9, 1.0/9,  1.0/9}, {1.0/9, 1.0/9,  1.0/9},{1.0/9, 1.0/9,  1.0/9}};
//		double[][] zsmo = imfilter(double_G, avgfilter);
//
//		double[] diff = new double [(Org_Im.length-1)*(Org_Im[0].length-1)]; //save difference without image Edges
//		int idx;
//		double min_sm_val=1, max_sm_val = 0;
//		for (int i = 1; i < Org_Im.length-1; i++) {
//			for (int j = 1; j < Org_Im[0].length-1; j++) {
//				idx = (i-1) * (Org_Im[0].length-1) + (j-1);
//				diff[idx] = (double_G[i][j] - zsmo[i][j]); // assume average does not have noise; indeed: var(a-1/9*sum(a1..a9) = 10/9 var ==> var = 9/10*x
//				if (zsmo[i][j] > max_sm_val)
//					max_sm_val = zsmo[i][j];
//				if (zsmo[i][j] < min_sm_val)
//					min_sm_val = zsmo[i][j];
//			}
		//		}

		int delta = 2;
		double[][] zsmo = new double[Org_Im.length][Org_Im[0].length];
		for (int i = 2; i < Org_Im.length-2; i++) {
			for (int j = 2; j < Org_Im[0].length-2; j++) {
				// 5*5 -3*3 neighbors
				for (int ii = i-2; ii<= i+2; ii++) {
					for (int jj = j-2; jj<= j+2; jj++) {
						if (ii == i-2 | ii == i+2 | jj==j-2 | jj==j+2) {
							zsmo[i][j] += double_G[ii][jj];
						}
					}
				}
				zsmo[i][j] = zsmo[i][j] / 16;
			}
		}		
		int idx = 0;
		double[][] diff = new double [Org_Im.length][Org_Im[0].length]; //save difference without image Edges
		double min_sm_val=1, max_sm_val = 0;
		for (int i = 2; i < Org_Im.length-2; i++) {
			for (int j = 2; j < Org_Im[0].length-2; j++) {
				// 5*5 -3*3 neighbors
				diff[i][j] = (double_G[i][j] - zsmo[i][j]);
				if (zsmo[i][j] > max_sm_val)
					max_sm_val = zsmo[i][j];
				if (zsmo[i][j] < min_sm_val)
					min_sm_val = zsmo[i][j];
			}
		}	

		// step 3. get the intensities and variances
		int nLevels = 256;
		double delta0 = (max_sm_val-min_sm_val) / (nLevels-1);
		
		double[] xMean = new double[nLevels];
		int[] eleCnt = new int[nLevels];
		double[] xVar = new double[nLevels];
		double[] diffMean = new double[nLevels];
		for (int i=0;i<nLevels;i++) {
			xMean[i] = 0;
			eleCnt[i] = 0;
			xVar[i] = 0;
			diffMean[i] = 0;
		}
		int curLevel;
		for (int i = delta; i < Org_Im.length-delta; i++) {
			for (int j = delta; j < Org_Im[0].length-delta; j++) {
				//idx = (i-1) * (Org_Im[0].length-1) + (j-1);
				curLevel = (int)Math.floor((zsmo[i][j]-min_sm_val)/delta0);
				xMean[curLevel] += zsmo[i][j];
				diffMean[curLevel] += diff[i][j];
				eleCnt[curLevel] ++;
			}
		}
		for (int i=0;i<nLevels;i++) {
			xMean[i] /= eleCnt[i];
			diffMean[i] /= eleCnt[i];
		}
		for (int i = delta; i < Org_Im.length-delta; i++) {
			for (int j = delta; j < Org_Im[0].length-delta; j++) {
				idx = (i-delta) * (Org_Im[0].length-delta) + (j-delta);
				curLevel = (int)Math.floor((zsmo[i][j]-min_sm_val)/delta0);
				xVar[curLevel] += Math.pow(diff[i][j] - diffMean[curLevel], 2);
			}
		}
		for (int i=0;i<nLevels;i++) {
			xVar[i] /= (eleCnt[i]-1);
		}
		
		int nPair = 200; // we only use the first 200 values to do estimation
		// LS Estimation use Apache
		/***way 1***/
//		WeightedObservedPoints obs = new WeightedObservedPoints();
//		int pairCnt = 0;
//		for (int i = 3; i < nLevels; i++) {
//			//System.out.println("mean: "+xMean[i]+" var:"+xVar[i]);
//			if (Double.isNaN(xMean[i]) | xVar[i]<=0)
//				continue;
//			obs.add(xMean[i], xVar[i]);
//			pairCnt++;
//			if (pairCnt >= nPair)
//				break;
//		}
//		PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);
//	    double[] pEst = fitter.fit(obs.toList());
		/***way 2***/
		int pairCnt = 0;
	    double[][] A = new double[nPair][2];
	    double[][] A_Inv = new double[2][nPair];
		double[] b_var = new double[nPair];
		for (int i = 0; i < nPair; i++) {
			if (Double.isNaN(xMean[i]) | Double.isNaN(xVar[i]) | xVar[i]<=0)
				continue;
			//System.out.println("mean: "+xMean[i]+" var:"+xVar[i]);
			A[pairCnt][0] = xMean[i];
			A[pairCnt][1] = 1;
			A_Inv[0][pairCnt] = xMean[i];
			A_Inv[1][pairCnt] = 1;
			b_var[pairCnt] = xVar[i];
			pairCnt++;
			if (pairCnt >= nPair)
				break;
		}
		RealMatrix A_real = MatrixUtils.createRealMatrix(A);
		RealMatrix A_Inv_real = MatrixUtils.createRealMatrix(A_Inv);
		RealMatrix AtA = A_Inv_real.multiply(A_real);
		RealVector b = MatrixUtils.createRealVector(b_var);
		RealVector Atb = A_Inv_real.operate(b);
		DecompositionSolver solver = new LUDecomposition(AtA).getSolver();
		RealVector solution = solver.solve(Atb);
		
	    double[] pEst0 = new double[2];
		pEst0[0] = solution.getEntry(0);//0.0080;//
		pEst0[1] = solution.getEntry(1);//-1.4e-4;//
		//System.out.println("alpha: "+pEst0[0]+" sigma:"+pEst0[1]);
		double alpha = pEst0[0];
		double sigma2 = pEst0[1];//Math.max(pEst0[1],2e-4); // this is the var of gausssian noise
		double t0 = 2/alpha*Math.sqrt(Math.max(0, alpha*0+3/8*alpha*alpha+sigma2));
		double t1 = 2/alpha*Math.sqrt(alpha*1+3/8*alpha*alpha+sigma2);

		double var1 = 1/((t1-t0)*(t1-t0));
		

		if(qq.pgEst==null) {
			qq.pgEst = new double [pEst0.length];
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = pEst0[i];
		}
		else {
			for (int i=0;i<qq.pgEst.length;i++)
				qq.pgEst[i] = qq.pgEst[i]+pEst0[i];
		}
		
		qq.var = var1*255*255;
		/***way 1***/
//		qq.varRatio = 0.9;
//		double cumuEleCnt = 0;
//		double totalEleCnt = (double) (Org_Im[0].length * Org_Im.length - eleCnt[0]);
//		int tgLevel = (int)Math.round(nLevels * 0.05);
//		for (int i=1; i<eleCnt.length;i++) {
//			if ((cumuEleCnt / totalEleCnt) > qq.varRatio & i >= tgLevel) {
//				tgLevel = i;
//				break;
//			}
//			cumuEleCnt += (double)eleCnt[i];
//		}
//		qq.varRatioBased = ((A[tgLevel][0])*alpha+sigma2)*255*255;
		/***way 2***/
		double cumuEleCnt = 0;
		double totalEleCnt = (double) (Org_Im[0].length * Org_Im.length - eleCnt[0]);
		int tgLevel = 0;
		for (int i=1; i<eleCnt.length;i++) {
			if ((cumuEleCnt / totalEleCnt) > (qq.varRatio / 2)) {
				tgLevel = i;
				break;
			}
			cumuEleCnt += (double)eleCnt[i];
		}
		qq.varRatioBased = xVar[tgLevel]*255*255;
//		if (qq.varRatioBased < 1) {
//			qq.varRatioBased = 1;
//		}
		System.out.println("alpha: "+pEst0[0]+" beta:"+pEst0[1]+" Variance:"+qq.var+"Variance 0.1:"+qq.varRatioBased);

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
}
