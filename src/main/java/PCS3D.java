import java.util.ArrayList;
import java.util.HashMap;

public class PCS3D {
	//int[][][] labelMap=segPC3d(input,zmap);
	
	public static void main(String[] args) {

		int[][][] input = new int[][][] {
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,3,1,0,0,3,0,0,0},
			{0,0,0,0,0,1,2,5,3,1,2,4,2,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
			{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,1,1,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,1,1,0,0,1,0,0,0},
			{0,0,7,0,0,0,0,6,1,0,0,6,0,0,0},
			{0,0,0,0,0,1,8,9,8,6,8,9,8,1,0},
			{0,0,0,0,0,0,0,6,1,0,0,6,0,0,0},
			{0,0,0,0,0,0,0,3,1,0,0,1,0,0,0},
			{0,0,7,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,1,2,4,3,1,2,4,2,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		};
		
		int[][][] zscore = new int[][][] {
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,3,1,0,0,3,0,0,0},
			{0,0,0,0,0,1,2,5,3,1,2,4,2,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},{
			{0,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,1,1,0,0,0,0,0,0,0},
			{0,1,1,0,0,0,0,1,1,0,0,1,0,0,0},
			{0,1,7,0,0,0,0,6,1,0,0,6,0,0,0},
			{0,0,0,0,0,1,8,9,8,6,8,9,8,1,0},
			{0,0,0,0,0,0,0,6,1,0,0,6,0,0,0},
			{0,0,1,0,0,0,0,3,1,0,0,1,0,0,0},
			{0,1,7,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
			{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,1,2,4,3,1,2,4,2,0,0},
			{0,0,0,0,0,0,0,2,1,0,0,3,0,0,0},
			{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
			{0,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		{
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		},
		};
		
		
		int H = input[0].length;
		int W = input[0][0].length;
		int T=input.length;
		
		boolean[][][] zmap = new boolean[T][H][W];
		
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				for(int z=0;z<T;z++) {
					if(zscore[z][x][y]>0) {
						zmap[z][x][y]=true;
					}else {
						zmap[z][x][y]=false;
					}
				}
			}
		}
		
		//segmentation by principal curvature
		
		double[][][] img3d = new double[T][H][W];
		
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				for(int z=0;z<T;z++) {
					img3d[z][x][y]=new Double(input[z][x][y]);
				}
			}
		}
		
		int[][][] labelMap=segPC3dMap(img3d,zmap);
		
		System.out.print("label:");
		System.out.println();
		for(int z=0;z<T;z++) {
			System.out.print("z="+z);
			System.out.println();
			for(int x=0;x<H;x++) {
				for(int y=0;y<W;y++) {
					System.out.print(labelMap[x][y][z]+" ");
				}
				System.out.println();
			}
			System.out.println();
		}
	}
	
	public static double[][][] shortArrayConvert(short[][] imArray,int height,int width){
		int zSlice = imArray.length;
		double[][][] G = new double[zSlice][height][width];
		for (int zz = 0; zz<zSlice; zz++) {
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					G[zz][i][j] = (double) imArray[zz][i * width + j] ;
				}
			}
		}
		return G;
	}
	
	public static boolean[][][] zmapConvert(double[][][] zmap, double thres){
		int zSlice =zmap.length;
		int height =zmap[0].length;
		int width = zmap[0][0].length;
		
		boolean[][][] z_positive = new boolean[zSlice][height][width];
		for (int zz = 0; zz<zSlice; zz++) {
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					z_positive[zz][i][j] = zmap[zz][i][j]>=thres;
				}
			}
		}
		return z_positive;
	}
	
	public static int[][][] segPC3dMap(double[][][]input,boolean[][][]zmap_o){
		int H = input[0].length;
		int W = input[0][0].length;
		int T = input.length;
		HashMap<Integer, ArrayList<int[]>> hashall=segPC3d(input,zmap_o);
		int[][][]remap3d=hash2map(hashall,H,W,T);
		return remap3d;
	}
	
	
	public static HashMap<Integer, ArrayList<int[]>> segPC3d(double[][][]input,boolean[][][]zmap_o){
		
		int H = input[0].length;
		int W = input[0][0].length;
		int T = input.length;
		HashMap<Integer, ArrayList<int[]>> hashall=new HashMap<Integer, ArrayList<int[]>>();
		int[][][]remap3d=new int[T][H][W];
		if(T>1) {
			HashMap<Integer, ArrayList<int[]>> hashu=new HashMap<Integer, ArrayList<int[]>>();
			HashMap<Integer, ArrayList<int[]>> hashc=new HashMap<Integer, ArrayList<int[]>>();
			
			HashMap<Integer, HashMap<Integer, ArrayList<int[]>>> chash=new HashMap<Integer,HashMap<Integer, ArrayList<int[]>>>();
			int[][][] nodeTree_all= new int[T][][];
			for(int z=0;z<T;z++) {
				double[][]img=input[z];
				boolean[][]zmap=zmap_o[z];
				
				HashMap<Integer, ArrayList<int[]>> PChash=getPC(img,zmap);
				ArrayList<int[]> seed=getPCseed(PChash);
				HashMap<Integer, ArrayList<int[]>> zhash= twoPassConnect2D(zmap);
				HashMap<Integer, ArrayList<int[]>> rehash=recover(seed,zhash);
				chash.put(z, hash2Dto3D(rehash,z));
				if(z==0) {
					hashu=rehash;
				}else {
					if(z==1) {
						hashc=rehash;
						nodeTree_all[z-1]=nodeTree(hashu,hashu,hashc,z-1,T);
					}else {
						nodeTree_all[z-1]=nodeTree(hashc,hashu,rehash,z-1,T);
						hashu=hashc;
						hashc=rehash;
					}
				}
			}
			nodeTree_all[T-1]=nodeTree(hashc,hashu,hashc,T-1,T);
			
			int k=1;
			int[][] nodeclass=new int[T][];
			for(int z=0;z<T;z++) {
				HashMap<Integer, ArrayList<int[]>> hash_z=chash.get(z);
				nodeclass[z]=new int[hash_z.size()];
				for(int num_c=1;num_c<=hash_z.size();num_c++) {
					if(nodeTree_all[z][num_c-1][0]==-1) {
						hashall.put(k, hash_z.get(num_c));
						nodeclass[z][num_c-1]=k;
						k++;
					}
				}
			}
			
			for(int z=0;z<T;z++) {
				HashMap<Integer, ArrayList<int[]>> hash_z=chash.get(z);
				for(int num_c=1;num_c<=hash_z.size();num_c++) {
					if(nodeTree_all[z][num_c-1][0]!=-1) {
						int[] pnode=findPnode(nodeTree_all,z,num_c-1);
						k=nodeclass[pnode[0]][pnode[1]];
						hashall.get(k).addAll(hash_z.get(num_c));
					}
				}
			}
			
			
		}else {
			HashMap<Integer, ArrayList<int[]>> hash2d=PCS.segPC(input[0],zmap_o[0]);
			hashall=hash2Dto3D(hash2d,0);
		}
		return hashall;
	}
	
	public static int[] findPnode(int[][][] nodeTree_all,int z,int num_c) {
		int[] Pnode= {-1,-1};
		if(nodeTree_all[z][num_c][0]!=-1) {
			Pnode=findPnode(nodeTree_all,nodeTree_all[z][num_c][0],nodeTree_all[z][num_c][1]);
		}else {
			Pnode= new int[] {z,num_c};
		}
		return Pnode;
	}
	
	public static int [][] nodeTree(HashMap<Integer, ArrayList<int[]>> hashc,
			HashMap<Integer, ArrayList<int[]>> hashu,
			HashMap<Integer, ArrayList<int[]>> hashl,int z,int maxz){
		int [][] nodeTree_c=new int[hashc.size()][];
		for(int i=1;i<=hashc.size();i++) {
			ArrayList<int[]> cl=hashc.get(i);
			int[] ucont= {0,-1};
			int[] lcont= {0,-1};
			int usize=0;
			int lsize=0;
			if(z<maxz-1) {
				for(int j=1;j<=hashl.size();j++) {
					ArrayList<int[]> ll=hashl.get(j);
					int loverlap=0;
					
					for(int[] cpix:cl) {
						for(int[] lpix:ll) {
							if(cpix[0]==lpix[0]&&cpix[1]==lpix[1]) {loverlap++;}
						}
					}
					if(loverlap>lcont[0]) {
						lcont[0]=loverlap;
						lcont[1]=j-1;
						lsize=ll.size();
					}
				}
			}
			
			if(z>0) {
				for(int k=1;k<=hashu.size();k++) {
					ArrayList<int[]> ul=hashu.get(k);
					int uoverlap=0;
					for(int[] cpix:cl) {
						for(int[] upix:ul) {
							if(cpix[0]==upix[0]&&cpix[1]==upix[1]) {uoverlap++;}
						}
					}
					if(uoverlap>ucont[0]) {
						ucont[0]=uoverlap;
						ucont[1]=k-1;
						usize=ul.size();
					}
				}
			}
			
			if(lcont[0]>ucont[0]&&cl.size()<lsize) {
				if(lcont[0]*2>cl.size()) {nodeTree_c[i-1]= new int[] {z+1,lcont[1]}; }else {nodeTree_c[i-1]= new int[] {-1,-1};}	
			}else if(lcont[0]<=ucont[0]&&cl.size()<=usize) {
				if(ucont[0]*2>cl.size()) {nodeTree_c[i-1]= new int[] {z-1,ucont[1]}; }else {nodeTree_c[i-1]= new int[] {-1,-1};}
			}else {
				nodeTree_c[i-1]= new int[] {-1,-1};
			}
		}
		
		return nodeTree_c;
	}
	
	public static HashMap<Integer, ArrayList<int[]>> hash2Dto3D(HashMap<Integer, ArrayList<int[]>> hash2d,int z){
		HashMap<Integer, ArrayList<int[]>> hash3d=new HashMap<Integer, ArrayList<int[]>>();
		
		for(int i=1;i<=hash2d.size();i++) {
			ArrayList<int[]> l=hash2d.get(i);
			ArrayList<int[]> ll=new ArrayList<int[]>();
			for(int[] j:l) {
				int[] jj= new int[] {j[0],j[1],z};
				ll.add(jj);
			}
			hash3d.put(i, ll);
		}
		return hash3d;
	}
	
	public static HashMap<Integer, ArrayList<int[]>> combineHash(HashMap<Integer, ArrayList<int[]>> chash,
			HashMap<Integer, ArrayList<int[]>> hash2d,int z){
		
		int startN=chash.size();
		
		for(int i=1;i<=hash2d.size();i++) {
			ArrayList<int[]> l=hash2d.get(i);
			ArrayList<int[]> ll=new ArrayList<int[]>();
			for(int[] j:l) {
				int[] jj= new int[] {j[0],j[1],z};
				ll.add(jj);
			}
			chash.put(startN+i, ll);
		}	
		return chash;
	}
	
	public static HashMap<Integer, ArrayList<int[]>> recover(ArrayList<int[]> seed, HashMap<Integer, ArrayList<int[]>> zhash) {
		HashMap<Integer, ArrayList<int[]>> remap=new HashMap<Integer, ArrayList<int[]>>();
		HashMap<Integer, ArrayList<int[]>> seedmap=new HashMap<Integer, ArrayList<int[]>>();
		
		int k=1;
		for(int i=1;i<=zhash.size();i++) {
			ArrayList<int[]> l =zhash.get(i);
			ArrayList<int[]> sdl=new ArrayList<int[]>();
			
			for(int[] sd:seed) {
				for(int[] pix:l) {
					if(pix[0]==sd[0]&&pix[1]==sd[1]) {
						sdl.add(sd);
					}
				}
			}
			seedmap.put(i, sdl);
			
			if(sdl.size()<2) {
				remap.put(k, l);
				k++;
			} else {
				HashMap<Integer, ArrayList<int[]>> localSeedMap=new HashMap<Integer, ArrayList<int[]>>();
				for(int[] j:l) {
					double dist=0;
					double dist_min=Double.MAX_VALUE;
					int idx=0;
					for(int ii=0;ii<sdl.size();ii++) {
						int[] sd= sdl.get(ii);
						dist=Math.pow(j[0]-sd[0], 2)+Math.pow(j[1]-sd[1], 2);
						if(dist<dist_min) {
							dist_min=dist;
							idx=ii;
						}
					}
					
					if(localSeedMap.containsKey(idx)) {
						localSeedMap.get(idx).add(j);
					}else {
						ArrayList<int[]> lst_t=new ArrayList<int[]>();
						lst_t.add(j);
						localSeedMap.put(idx,lst_t);
					}
				}
				
				for(int jj=0;jj<sdl.size();jj++) {
					remap.put(k, localSeedMap.get(jj));
					k++;
				}
				
			}
	  }
	return remap;
	}
	
	
	public static ArrayList<int[]> getPCseed( HashMap<Integer, ArrayList<int[]>> hash){
		ArrayList<int[]> seed = new ArrayList<int[]>();
		for(int i=1;i<=hash.size();i++) {
			  ArrayList<int[]> l = hash.get(i);
			  double[] center= new double[]{0.0,0.0};
			  for(int[] p:l) {
				  center[0]=center[0]+p[0];
				  center[1]=center[1]+p[1];
			  }
			  int[] idx=new int[2];
			  idx[0]=(int) Math.round(center[0]/l.size());
			  idx[1]=(int) Math.round(center[1]/l.size());
			  seed.add(idx);
		}
		return seed;
	}
	
	
	
	public static HashMap<Integer, ArrayList<int[]>> getPC(double[][] img,boolean[][] zmap) {
		
		int H=img.length;
		int W=img[0].length;
		boolean[][] PCmap=new boolean[H][W];
		
		double[][][] Grad=getGradient(img);
		double[][][] Gradx=getGradient(Grad[0]);
		double[][][] Grady=getGradient(Grad[1]);
		
		for(int i=0;i<H;i++) {
			for(int j=0;j<W;j++) {
				if((Gradx[0][i][j]+Grady[1][i][j])<0 && (Gradx[0][i][j]*Grady[1][i][j]-Gradx[1][i][j]*Grady[0][i][j])>0&&zmap[i][j]==true) {
					PCmap[i][j]=true;
				}else {
					PCmap[i][j]=false;
				}
			}
		}
		HashMap<Integer, ArrayList<int[]>> PChash= twoPassConnect2D(PCmap);
		return PChash;
	}
	
	public static double[][][] getGradient(double[][] im) {
		
		int H=im.length;
		int W=im[0].length;
		double[][] GradX=new double[H][W];
		double[][] GradY=new double[H][W];
		double[][] sobel= {{1,2,1},{0,0,0},{-1,-2,-1}};//gradient filter by sobel
		
		for(int i=0;i<H;i++) {
			for(int j=0;j<W;j++) {
				GradX[i][j]=0;
				GradY[i][j]=0;
			}
		}
		
		for(int i=1;i<H-1;i++) {
			for(int j=1;j<W-1;j++) {
				for(int ii=0;ii<3;ii++) {
					if(ii==1) {continue;}
					for(int jj=0;jj<3;jj++) {
						GradX[i][j]=GradX[i][j]+sobel[ii][jj]*im[i+ii-1][j+jj-1];
						GradY[i][j]=GradY[i][j]+sobel[ii][jj]*im[i+jj-1][j+ii-1];
					}
				}
			}
		}
		
		double[][][] Grad=new double[2][H][W];
		Grad[0]=GradX;
		Grad[1]=GradY;
		
		return Grad;
		
	}
	
	public static int[][][] hash2map(HashMap<Integer, ArrayList<int[]>> map,int H,int W,int T) {
		
		
		  int[][][] labelMap = new int[H][W][T];
		  for(int i=1;i<=map.size();i++) {
			  ArrayList<int[]> l = map.get(i);
			  for(int[] p:l) {
				  labelMap[p[0]][p[1]][p[2]] = i;
			  }
		  }
		  return labelMap;
	}

	public static HashMap<Integer, ArrayList<int[]>> twoPassConnect2D(boolean[][] input) {  
		  int width = input.length;
		  int height = input[0].length;
		  
		  int[][] label = new int[width][height];
		  ArrayList<Integer> list = new ArrayList<Integer>();
		  list.add(0);
		  int curLabel = 1;
		  for(int i=0;i<width;i++) {
		   for(int j=0;j<height;j++) {
		    if(input[i][j]) {
		     int[] labels = new int[4];
		     labels[0] = i>0 && j>0? label[i-1][j-1]:0;//4 or 8 connected
		     labels[1] = i>0? label[i-1][j]:0;
		     labels[2] = i>0 && j<height-1? label[i-1][j+1]:0;
		     labels[3] = j>0? label[i][j-1]:0;
		     
		     ArrayList<Integer> labelList = new ArrayList<Integer>();
		     int min = Integer.MAX_VALUE;
		     for(int ii=0;ii<4;ii++) {
		      if(labels[ii]!=0) {
		       min = Math.min(min, labels[ii]);
		       labelList.add(labels[ii]);
		      }
		     }
		     if(labelList.size()==0) {
		      label[i][j] = curLabel;
		      list.add(0);
		      curLabel++;
		     }else {
		      label[i][j] = min;
		      for(int ii=0;ii<labelList.size();ii++) {
		       union_connect(min,labelList.get(ii),list);
		      }
		     }
		    }
		   }
		  }
		  
		  HashMap<Integer,Integer> rootMap = new HashMap<Integer, Integer>();
		  int cnt = 1;
		  HashMap<Integer, ArrayList<int[]>> map = new HashMap<Integer, ArrayList<int[]>>();
		  for(int i=0;i<width;i++) {
		   for(int j=0;j<height;j++) {
		    if(label[i][j]!=0) {
		     int root = union_find(label[i][j], list);
		     int value;
		     if(rootMap.get(root)!=null) {
		      value = rootMap.get(root);
		     }else {
		      value = cnt;
		      rootMap.put(root, value);
		      cnt++;
		     }
		     
		     label[i][j] = value;
		     ArrayList<int[]> l = map.get(value);
		     if(l==null) {
		      l= new ArrayList<int[]>();
		      map.put(value, l);
		     }
		     l.add(new int[] {i,j});
		    }
		   }
		  }
		  
		  int[][] labelMap = new int[width][height];
		  for(int i=1;i<=map.size();i++) {
			  ArrayList<int[]> l = map.get(i);
			  for(int[] p:l) {
				  labelMap[p[0]][p[1]] = i;
			  }
		  }
		  return map;
	 }
	
	public static int union_find(int label, ArrayList<Integer> list){
		  int i = label;
		  
		  while(list.get(i)!=0) {
		   i = list.get(i);
		  }
		  if(i!=label)
		   list.set(label, i);
		  return i;
		 }
		 
	 public static void union_connect(int label1, int label2, ArrayList<Integer> list) {
		  if(label1==label2)
		   return;
		  int i = union_find(label1,list);
		  int j = union_find(label2,list);
		  if(i!=j)
		   list.set(j, i);
	 }
}
