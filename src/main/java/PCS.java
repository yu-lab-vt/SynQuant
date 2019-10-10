import java.util.ArrayList;
import java.util.HashMap;

public class PCS {
	//int[][]labelMap=segPC(input,zmap);
	
	
	public static void main(String[] args) {

		int[][] input = new int[][] {
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
		};
		
		int[][] zscore = new int[][] {
			{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,0,1,0,0,0,0,0,0,0,0,0,1,1,0},
			{0,1,1,0,0,0,1,1,0,0,0,0,0,1,0},
			{0,0,0,0,0,0,1,1,1,0,1,1,0,0,0},
			{0,0,0,0,1,1,1,1,1,1,1,1,1,0,0},
			{0,0,0,1,1,1,1,1,1,1,1,1,1,1,0},
			{0,0,0,0,1,1,1,1,1,0,1,1,1,0,0},
			{0,0,0,0,0,0,1,1,1,0,1,1,0,0,0},
			{0,1,1,0,0,0,1,1,0,0,0,0,0,0,0},
			{0,1,1,0,0,0,0,0,0,0,0,0,0,1,1},
		};
		
		
		int H = input.length;
		int W = input[0].length;
		
		boolean[][] zmap = new boolean[H][W];
		
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				if(zscore[x][y]>0) {
					zmap[x][y]=true;
				}else {
					zmap[x][y]=false;
				}
			}
		}
		
		
		
		// display input
		System.out.print("Original:");
		System.out.println();
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				System.out.print(input[x][y]+" ");
			}
			System.out.println();
		}
		System.out.println();
		
		// display zmap
		System.out.print("zmap:");
		System.out.println();
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				System.out.print(zscore[x][y]+" ");
			}
			System.out.println();
		}
		System.out.println();
		
		//segmentation by principal curvature
		

		double[][] img = new double[H][W];
		
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				img[x][y]=new Double(input[x][y]);
			}
		}
		
		int[][]labelMap=segPCmap(img,zmap);
		
		System.out.print("zmap:");
		System.out.println();
		
		for(int x=0;x<H;x++) {
			for(int y=0;y<W;y++) {
				System.out.print(labelMap[x][y]+" ");
			}
			System.out.println();
		}
		System.out.println();
		
	}
	
	public static int[][] segPCmap(double[][]input,boolean[][]zmap){
		
		HashMap<Integer, ArrayList<int[]>> rehash=segPC(input,zmap);
		int[][] remap=hash2map(rehash,zmap);
		return remap;
	}
	
	public static HashMap<Integer, ArrayList<int[]>> segPC(double[][]input,boolean[][]zmap){
		
		
		HashMap<Integer, ArrayList<int[]>> PChash=getPC(input,zmap);
		ArrayList<int[]> seed=getPCseed(PChash);
		HashMap<Integer, ArrayList<int[]>> zhash= twoPassConnect2D(zmap);
		HashMap<Integer, ArrayList<int[]>> rehash=recover(seed,zhash);
		
		return rehash;
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
	
	public static boolean[][] getPC_all(double[][] img) {
		
		int H=img.length;
		int W=img[0].length;
		boolean[][] PCmap=new boolean[H][W];
		
		double[][][] Grad=getGradient(img);
		double[][][] Gradx=getGradient(Grad[0]);
		double[][][] Grady=getGradient(Grad[1]);
		
		for(int i=0;i<H;i++) {
			for(int j=0;j<W;j++) {
				if((Gradx[0][i][j]+Grady[1][i][j])<=0 && (Gradx[0][i][j]*Grady[1][i][j]-Gradx[1][i][j]*Grady[0][i][j])>=0) {
					PCmap[i][j]=true;
				}else {
					PCmap[i][j]=false;
				}
			}
		}
		return PCmap;
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
	
	public static int[][] hash2map(HashMap<Integer, ArrayList<int[]>> map,boolean[][] input) {
		
		  int width = input.length;
		  int height = input[0].length;
		
		  int[][] labelMap = new int[width][height];
		  for(int i=1;i<=map.size();i++) {
			  ArrayList<int[]> l = map.get(i);
			  for(int[] p:l) {
				  labelMap[p[0]][p[1]] = i;
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
