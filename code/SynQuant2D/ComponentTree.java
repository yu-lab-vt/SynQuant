import java.util.ArrayList;
/**
 * This class build a component tree from a image array. Also contain another two functions "cpt2map" which transfers component tree to 
 * zscore maps and "extractBestSyn" which extracts highest-score region.
 * 
 * Major idea of tree building and zscore calculation is from Dr.Petter Ranefall(petter.ranefall@it.uu.se).
 * For details of creating component tree, please refer to “Fast Adaptive Local Thresholding Based on Ellipse fit”(ISBI'16) and 
 * “Building the Component Tree in Quasi-Linear Time” (TIP 2006 Vol.15 No.11).
 * @version 1.0
 * @date 2016-07-11
 *
 */


public class ComponentTree{
	protected int[] sortedIndex;//store the elements in a new ordered array
	protected int[] parNode; //Parent nodes
	
	//for ScanAll3 use, initialization
	protected int [][][] kMapx;//regions numbering for each threshold
	protected double [][][] zMapx;//zscore maps for each threshold
	protected int Nthr; //number of thresholds used
	protected int[] thrrg; // thresholds
	
	//for extract highest zscore
	protected ArrayList<Double> Zscore_Vec;//zscore vector
	protected ArrayList<Integer> levels;
	double C;//highest zscore
	int I;//the intensity level that C lies in 
	//for prior knowledge test
	int[][] BxCor;//each area's bounding box coordinates
	//for component tree build
	
	short[] imArray;
    byte[] outputArray;
    double[] diffN; //Neighborhood difference
    double[] zscore;
	
	double[] pixSum; //Pixel sum inside objects
	double[] pixSumN; //Pixel sum in neighbourhood
	int[] areas; //Areas inside objects
	int[] areasN; //Areas of neighborhoods
	byte[] usedN; //Flag if used in neighborhood
	int width;
	int height;
	int nSlices;
	//CONSTANTS
	public static final byte UNDEFINED = (byte)0;
	public static final byte NOT_OBJECT = (byte)1;
	public static final byte OBJECT = (byte)2;
	public static final byte NOT_USED_AS_N = (byte)0;
	public static final byte USED_AS_N_ONCE = (byte)1;
	public static final byte USED_AS_N_MORE = (byte)2;
	/* 1 build the component tree of input image array (8-connected)
	 * 2 label objects by size test
	 * 3 calculate z-score of each object pixel */
	public ComponentTree(short[] imArrayIn,byte[] usedNIn, int w, int h,ParaP p, ParaQ q) 
	{
		imArray = imArrayIn;
		width = w;
		height = h;
		int nPixels=width*height;
		//Create nodes
        parNode=new int[nPixels];
        BxCor = new int[nPixels][4];//ymin,xmin,ymax,xmax.
		int x,y,x0,x2,y0,y2;
		int i,j,k;
		outputArray = new byte[nPixels];
		diffN = new double[nPixels];
		pixSum = new double[nPixels];
		pixSumN = new double[nPixels];
		usedN = new byte[nPixels]; 
		areas = new int[nPixels];
		areasN = new int[nPixels];
		for (i=0; i<nPixels; i++)
		{
			y=i/width;
			x=i-y*width;
			areas[i] = 1;
			areasN[i] = 0;
			pixSum[i] = imArray[i];
			pixSumN[i] = 0;
			diffN[i] = 0;
			//usedN[i] = NOT_USED_AS_N;
			usedN[i] = usedNIn[i];
			BxCor[i][0] = y;
			BxCor[i][1] = x;
			BxCor[i][2] = y;
			BxCor[i][3] = x;
		}
		
		//Sort points
		// create a counting array, counts, with a member for 
		// each possible discrete value in the input.  
		// initialize all counts to 0.
		BasicMath BM = new BasicMath();
		int maxP = BM.vectorMax(imArray);
		int minP = BM.vectorMin(imArray);
		int nLevels = maxP-minP + 1;
		int[] counts = new int[nLevels];
		// for each value in the unsorted array, increment the
		// count in the corresponding element of the count array
		for (i=0; i<nPixels; i++)
		{
			counts[imArray[i]-minP]++;
		}
		// accumulate the counts - the result is that counts will hold
		// the offset into the sorted array for the value associated with that index
		for (i=1; i<nLevels; i++)
		{
			counts[i] += counts[i-1];
		}
		// store the elements in a new ordered array
		sortedIndex = new int[nPixels];
		for (i = nPixels-1; i >= 0; i--)
		{
			// decrementing the counts value ensures duplicate values in A
			// are stored at different indices in sorted.
			sortedIndex[--counts[imArray[i]-minP]] = i;                
		}

		//Init nodes
		for (i=0; i<nPixels; i++)
		{
			parNode[i]=i;
		}
		//Search in decreasing order
		int curNode;
		int adjNode;
		boolean found;
		for (i = nPixels-1; i >= 0; i--)
		{
			j=sortedIndex[i];
			curNode=j;
			//System.out.println("Image Value"+imArray[j]);
			y=j/width;
			x=j-y*width;			

			y0=y-1;
			y2=y+1;
			x0=x-1;
			x2=x+1;

			//Later neigbours x2,y2
			found = false;
			if(y2<height)
			{
				k=x+width*y2;
				if(imArray[k]>=imArray[j])
				{
					adjNode=findNode(k);
					if(curNode!=adjNode)
					{
						curNode=mergeNodes(adjNode,curNode);
						found = true;
					}
				}
				if(x2<width)
				{
					k=x2+width*y2;
					if(imArray[k]>=imArray[j])
					{
						adjNode=findNode(k);
						if(curNode!=adjNode)
						{
							curNode=mergeNodes(adjNode,curNode);
							found = true;
						}
					}
				}
				if(x0>=0)
				{
					k=x0+width*y2;
					if(imArray[k]>=imArray[j])
					{
						adjNode=findNode(k);
						if(curNode!=adjNode)
						{
							curNode=mergeNodes(adjNode,curNode);
							found = true;
						}
						
					}
				}
			}
			if(x2<width)
			{
				k=x2+width*y;
				if(imArray[k]>=imArray[j])
				{
					adjNode=findNode(k);
					if(curNode!=adjNode)
					{
						curNode=mergeNodes(adjNode,curNode);
						found = true;
					}
				}
			}
			//Earlier neighbours x0,y0. No need to check =
			if(x0>=0)
			{
				k=x0+width*y;
				if(imArray[k]>imArray[j])
				{
					adjNode=findNode(k);
					if(curNode!=adjNode)
					{
						curNode=mergeNodes(adjNode,curNode);
						found = true;
					}
					
				}
			}
			if (y0 >= 0) {
				k = x + width * y0;
				if (imArray[k] > imArray[j]) {
					adjNode = findNode(k);
					if (curNode != adjNode) {
						curNode = mergeNodes(adjNode,curNode);
						found = true;
					}
				}
				if(x2<width)
				{
					k=x2+width*y0;
					if(imArray[k]>imArray[j])
					{
						adjNode=findNode(k);
						if(curNode!=adjNode)
						{
							curNode=mergeNodes(adjNode,curNode);
							found = true;
						}
					}
				}
				if(x0>=0)
				{
					k=x0+width*y0;
					if(imArray[k]>imArray[j])
					{
						adjNode=findNode(k);
						if(curNode!=adjNode)
						{
							curNode=mergeNodes(adjNode,curNode);
							found = true;
						}
						
					}
				}
			}
			
			if (!found)
			{
				//Later neigbours x2,y2
				if(y2<height)
				{
					k=x+width*y2;
					if (usedN[k] == NOT_USED_AS_N)
					{
						pixSumN[j] += imArray[k];
						areasN[j]++;
						usedN[k] = USED_AS_N_ONCE;
					}
					if(x2<width)
					{
						k=x2+width*y2;
						if (usedN[k] == NOT_USED_AS_N)
						{
							pixSumN[j] += imArray[k];
							areasN[j]++;
							usedN[k] = USED_AS_N_ONCE;
						}
					}
					if(x0>=0)
					{
						k=x0+width*y2;
						if (usedN[k] == NOT_USED_AS_N)
						{
							pixSumN[j] += imArray[k];
							areasN[j]++;
							usedN[k] = USED_AS_N_ONCE;
						}
					}
				}
				if(x2<width)
				{
					k=x2+width*y;
					if (usedN[k] == NOT_USED_AS_N)
					{
						pixSumN[j] += imArray[k];
						areasN[j]++;
						usedN[k] = USED_AS_N_ONCE;
					}
				}
				//Earlier neighbours x0,y0.
				if(x0>=0)
				{
					k=x0+width*y;
					if (usedN[k] == NOT_USED_AS_N)
					{
						pixSumN[j] += imArray[k];
						areasN[j]++;
						usedN[k] = USED_AS_N_ONCE;
					}
				}
				if(y0>=0)
				{
					k=x+width*y0;
					if (usedN[k] == NOT_USED_AS_N)
					{
						pixSumN[j] += imArray[k];
						areasN[j]++;
						usedN[k] = USED_AS_N_ONCE;
					}
					if(x2<width)
					{
						k=x2+width*y0;
						if (usedN[k] == NOT_USED_AS_N)
						{
							pixSumN[j] += imArray[k];
							areasN[j]++;
							usedN[k] = USED_AS_N_ONCE;
						}
					}
					if(x0>=0)
					{
						k=x0+width*y0;
						if (usedN[k] == NOT_USED_AS_N)
						{
							pixSumN[j] += imArray[k];
							areasN[j]++;
							usedN[k] = USED_AS_N_ONCE;
						}
					}
				}
				usedN[j] = USED_AS_N_ONCE;
				diffN[j] = pixSum[j]/(double)areas[j];
				if (areasN[j] > 0)
					diffN[j] -= pixSumN[j]/(double)areasN[j];
			}
		}
		ObjLabel(p.min_size,p.max_size);
		zscore = new double[nPixels];
		for (i=0; i<nPixels; i++)
		{
			y=i/width;
			x=i-y*width;
			double LH = (double)BxCor[i][3]-BxCor[i][1]+1;
			double LW = (double)BxCor[i][2]-BxCor[i][0]+1;
			double ratio = LH>LW? LH/LW: LW/LH;
			if(areas[i]<p.min_size || ratio>p.maxWHratio || (areas[i]/(double)(LH*LW))<p.minfill){
				zscore[i] = 0;
				continue;
			}
			if(outputArray[i] == OBJECT)
			{
				zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p.mu, p.sigma,q.var);
			}
		}
	}
	/*Label object or not: here we simply  test the size  to decide object or not*/
	public void ObjLabel(int minSize ,int maxSize)
	{
		int i,j;
		int nPixels = sortedIndex.length;
		for (i = nPixels-1; i >= 0; i--)
		{
			j=sortedIndex[i];
			if (areas[j]<=maxSize)// && areas[j]>=minSize)
				outputArray[j] = OBJECT;
			else
				outputArray[j] = NOT_OBJECT;
		}
	}
	public int findNode(int e)
	{
		if(parNode[e]!=e)
		{
			int root = findNode(parNode[e]);
			//parNode[e] = root; //This cannot be used here
			return root;
		}
		else
		{
			return e;
		}
	}	
	public int mergeNodes(int e1,int e2)//, double[][] pMu,double[][] pSigma, double qVar)
	{
		//e1-adjacent; e2-current
		int res;
		int m;
		
		if(imArray[e1]==imArray[e2])
		{
			res=Math.max(e1,e2);
			m=Math.min(e1,e2);
		}
		else
		{
			res=e2;
			m=e1;
		}
		//Compute new neighbours
		int y=e2/width;
		int x=e2-y*width;			
		
		int y0=y-1;
		int y2=y+1;
		int x0=x-1;
		int x2=x+1;
		int k;
		if(y2<height)
		{
			k=x+width*y2;
			if (usedN[k] == NOT_USED_AS_N)
			{
				pixSumN[res] += imArray[k];
				areasN[res]++;
				usedN[k] = USED_AS_N_ONCE;
			}
			if(x2<width)
			{
				k=x2+width*y2;
				if (usedN[k] == NOT_USED_AS_N)
				{
					pixSumN[res] += imArray[k];
					areasN[res]++;
					usedN[k] = USED_AS_N_ONCE;
				}
			}
			if(x0>=0)
			{
				k=x0+width*y2;
				if (usedN[k] == NOT_USED_AS_N)
				{
					pixSumN[res] += imArray[k];
					areasN[res]++;
					usedN[k] = USED_AS_N_ONCE;
				}
			}
		}
		if(x2<width)
		{
			k=x2+width*y;
			if (usedN[k] == NOT_USED_AS_N)
			{
				pixSumN[res] += imArray[k];
				areasN[res]++;
				usedN[k] = USED_AS_N_ONCE;
			}
		}
		//Earlier neighbours x0,y0.
		if(x0>=0)
		{
			k=x0+width*y;
			if (usedN[k] == NOT_USED_AS_N)
			{
				pixSumN[res] += imArray[k];
				areasN[res]++;
				usedN[k] = USED_AS_N_ONCE;
			}
		}
		if(y0>=0)
		{
			k=x+width*y0;
			if (usedN[k] == NOT_USED_AS_N)
			{
				pixSumN[res] += imArray[k];
				areasN[res]++;
				usedN[k] = USED_AS_N_ONCE;
			}
			if(x2<width)
			{
				k=x2+width*y0;
				if (usedN[k] == NOT_USED_AS_N)
				{
					pixSumN[res] += imArray[k];
					areasN[res]++;
					usedN[k] = USED_AS_N_ONCE;
				}
			}
			if(x0>=0)
			{
				k=x0+width*y0;
				if (usedN[k] == NOT_USED_AS_N)
				{
					pixSumN[res] += imArray[k];
					areasN[res]++;
					usedN[k] = USED_AS_N_ONCE;
				}
			}
		}		
		
		areasN[res] += areasN[m];
		pixSumN[res] += pixSumN[m];
		if (usedN[e2] == USED_AS_N_ONCE)
		{
			areasN[res] -= areas[e2];
			pixSumN[res] -= pixSum[e2];
		}
		areas[res] += areas[m];
		pixSum[res] += pixSum[m];
		parNode[m]=res;
		
		usedN[e2] = USED_AS_N_MORE;
		
		diffN[res] = pixSum[res]/(double)areas[res];
		
		if (areasN[res] > 0)
			diffN[res] -= pixSumN[res]/(double)areasN[res];
		//System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
		//System.out.println("BC:" +BxCor[m][3]+" "+BxCor[m][2]+" "+BxCor[m][1]+" "+BxCor[m][0]);
		if (BxCor[res][0] > BxCor[m][0])
			BxCor[res][0] = BxCor[m][0];
		if (BxCor[res][1] > BxCor[m][1])
			BxCor[res][1] = BxCor[m][1];
		if (BxCor[res][2] < BxCor[m][2])
			BxCor[res][2] = BxCor[m][2];
		if (BxCor[res][3] < BxCor[m][3])
			BxCor[res][3] = BxCor[m][3];
		//System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
		return res;
	}
	/*Calculate the z-score of each pixel*/
	public double zscoreCal(double t0, int M/*in*/, int N/*nei*/, double[][] pMu,double[][] pSigma, double qVar){
		if(M < 4)
	        M = 4;
	    if(N<=3 || N < (M/10))
	        return 0;
		if(N<4)
			N = 4;
	    if(M>100)
	        M = 100;
	    if(N>100)
	        N = 99;
	    double mu = pMu[M-1][N-1];
	    double sigma = pSigma[M-1][N-1];
	    mu = mu*Math.sqrt(qVar);
	    sigma = sigma*Math.sqrt(qVar);
	    double zScore = (t0-mu)/sigma;
	    t0 = t0/Math.sqrt(qVar);
		
		return zScore;
	}
	/* This function is to connect the component tree to the following functions like scanAll3Updt in scanAll3 class.
	 * It transfer the component tree into Nthr zscore maps (zMapx) and index maps (kMapx).
	 * Each map corresponding to a threshold, the connected region in one map share same zscore and index.
	 * 
	 * Also this function build two arraylists for the following FDR control: one is Zscore_Vec, which saves all zscores.
	 * The other is levels, which saves the map level(just the threshold level) of each zscore.*/
	public void cpt2map(ParaP p){
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		kMapx = new int[Nthr][height][width];
		boolean[][][] bMapx = new boolean[Nthr][height][width];
		zMapx = new double[Nthr][height][width];
		int nPixels = width*height;
		ImageHandling IH = new ImageHandling();
		
		Zscore_Vec = new ArrayList<Double>();
		levels = new ArrayList<Integer>();
		for (int i = nPixels-1; i >=0; i--)
		{
			int j=sortedIndex[i];
			if(outputArray[j] != OBJECT)
				continue;
			
			int y=j/width;
			int x=j-y*width;	
			int e = j;
			int intensity = imArray[e];
			if(intensity<p.thr0)
				continue;
			int Inl = (intensity-p.thr0)/p.thrg;
			int Inlevel = Inl;
			zMapx[Inlevel][y][x] = zscore[e];
			bMapx[Inlevel][y][x] = true;

			while (parNode[e] != e){
				intensity = imArray[parNode[e]];
				if(intensity<p.thr0 | outputArray[parNode[e]] != OBJECT){
					break;
				}
				int InlevelPar = (intensity-p.thr0)/p.thrg;
				zMapx[InlevelPar][y][x] = zscore[parNode[e]];
				bMapx[InlevelPar][y][x] = true;
				e = parNode[e];
			}
		}
		int lvCnt = 0;
		for(int i=0;i<Nthr;i++){
			kMapx[i] = IH.bwlabel(zMapx[i], 8,Zscore_Vec);
			for(int j=lvCnt;j<Zscore_Vec.size();j++){
				levels.add(i);
			}
			lvCnt = Zscore_Vec.size();
		}
	}
	/* The same function as BestSyn.extractBestSyn_fdr. Extract the region with highest zscore. 
	 * Now we just need to scan the ArrayList Zscore_Vec and levels. Now need to go through all maps.*/
	public void extractBestSyn(ParaP p){
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		int nPixels = width*height;
		boolean[] checkedMark = new boolean[nPixels];
		if(Zscore_Vec==null){//if the vector haven't been build, this will not be used
			Zscore_Vec = new ArrayList<Double>();
			levels = new ArrayList<Integer>();
			for (int i = 0; i <=nPixels; i++)
			{
				int j=sortedIndex[i];
				if(outputArray[j] != OBJECT)
					continue;
				if (checkedMark[j])
					continue;
				int e = j;
				int intensity = imArray[e];
				if(intensity<p.thr0)
					continue;
				int Inlevel = (intensity-p.thr0)/p.thrg;
				checkedMark[e] = true;
				while (parNode[e] != e){
					intensity = imArray[parNode[e]];
					if(intensity<p.thr0)
						break;
					int newInlevel = (int) Math.floor((intensity-p.thr0)/(double)p.thrg);
					checkedMark[e] = true;
					if(newInlevel!=Inlevel){
						Zscore_Vec.add(zscore[e]);
						levels.add(Inlevel);
						Inlevel = newInlevel;
					}
					e = parNode[e];
					
				}
				Zscore_Vec.add(zscore[e]);
				levels.add(Inlevel);
			}
		}
		C = 0;
		I = 0;
		int CI_cnt = 0;
		
		for(int i=0;i<Zscore_Vec.size();i++)
		{
			if(Zscore_Vec.get(i)>C){
				C = Zscore_Vec.get(i);
				I = levels.get(i);
				CI_cnt = i;
			}
		}
		Zscore_Vec.remove(CI_cnt);
		levels.remove(CI_cnt);
		
	}
}