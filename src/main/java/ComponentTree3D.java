import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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


public class ComponentTree3D{
	protected int[] sortedIndex;//store the elements in a new ordered array
	protected int[] parNode; //Parent nodes

	//for ScanAll3 use, initialization
	protected int [][][][] kMapx;//regions numbering for each threshold
	protected double [][][][] zMapx;//zscore maps for each threshold
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

	double[] voxSum; //Pixel sum inside objects
	double[] voxSumN; //Pixel sum in neighbourhood
	int[] areas; //Areas inside objects
	int[] areasN; //Areas of neighborhoods
	byte[] usedN; //Flag if used in neighborhood
	int width;
	int height;
	int zSlice;
	//CONSTANTS
	public static final byte UNDEFINED = (byte)0;
	public static final byte NOT_OBJECT = (byte)1;
	public static final byte OBJECT = (byte)2;
	public static final byte NOT_USED_AS_N = (byte)0;
	public static final byte USED_AS_N_ONCE = (byte)1;
	public static final byte USED_AS_N_MORE = (byte)2;
	//////////////////////Not finished/////////////////////////////////////
	/* 1 build the component tree of input image array (8-connected)
	 * 2 label objects by size test
	 * 3 calculate z-score of each object pixel */
	public ComponentTree3D(short[][] imArrayIn,byte[] usedNIn, paraP3D p, paraQ3D q) 
	{
		//imArray = imArrayIn;
		zSlice = imArrayIn.length;
		width = q.Nx;
		height = q.Ny;
		int nVoxels=width*height*zSlice; // voxels in a 3D image
		int nPixels = width*height; //pixels in a single slice
		imArray = new short[nVoxels];
		int tmpCnt = 0;
		
		for(int i = 0; i < imArrayIn.length; i++) {
			for(int j = 0; j < imArrayIn[0].length; j++)
				imArray[tmpCnt++] = imArrayIn[i][j]; // java.util.Arrays
		}
		//Create nodes
		parNode=new int[nVoxels];
		BxCor = new int[nVoxels][6];//ymin,xmin,zmin,  ymax,xmax,zmax.
		int x,y,z, rmder, x0,x2,y0,y2,z0, z2;
		int i,j,k;
		outputArray = new byte[nVoxels];
		diffN = new double[nVoxels];
		voxSum = new double[nVoxels];
		voxSumN = new double[nVoxels];
		usedN = new byte[nVoxels]; 
		areas = new int[nVoxels];
		areasN = new int[nVoxels];
		for (i=0; i<nVoxels; i++)
		{
			rmder = i % nPixels;
			z = i / nPixels;
			y=rmder/width;
			x=rmder-y*width;
			areas[i] = 1;
			areasN[i] = 0;
			voxSum[i] = imArray[i];
			voxSumN[i] = 0;
			diffN[i] = 0;
			//usedN[i] = NOT_USED_AS_N;
			usedN[i] = usedNIn[i];
			BxCor[i][0] = y;
			BxCor[i][1] = x;
			BxCor[i][2] = z;
			BxCor[i][3] = y;
			BxCor[i][4] = x;
			BxCor[i][5] = z;
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
		for (i=0; i<nVoxels; i++)
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
		sortedIndex = new int[nVoxels];
		for (i = nVoxels-1; i >= 0; i--)
		{
			// decrementing the counts value ensures duplicate values in A
			// are stored at different indices in sorted.
			sortedIndex[--counts[imArray[i]-minP]] = i;                
		}

		//Init nodes
		for (i=0; i<nVoxels; i++)
		{
			parNode[i]=i;
		}
		//Search in decreasing order
		int curNode;
		int adjNode;
		int ii,jj, kk, tmpIdx;
		boolean found;
		for (i = nVoxels-1; i >= 0; i--)
		{
			j=sortedIndex[i];
			curNode=j;
			//System.out.println("Image Value"+imArray[j]);
			rmder = j % nPixels;
			z = j / nPixels;
			y=rmder/width;
			x=rmder-y*width;

			
			found = false;
			/* ROI selection: for the same z-stack, we use 8 neighbors, but for other z-stacks, we only consider two direct neighbors
			* the order of the neighbors are very important
			* we go through the larger neighbors first, then lower ones
			*/
			y0=y-1;
			y2=y+1;
			x0=x-1;
			x2=x+1;
			z0=z-1;
			z2=z+1;
			/**for debug*
			if (imArray[j]>0) {
				System.out.print(j+" " + imArray[j]+"\n");
			}
			if(imArray[j]>=255)//usedN[28+width*28+6*nPixels] != NOT_USED_AS_N)
				System.out.println(" "+x+" "+y+" "+z);
			*/
			//Later neigbours x2,y2
			if(z2<zSlice) {
				k = x+width*y+z2*nPixels;
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
			if(y2<height)
			{
				k=x+width*y2+z*nPixels;
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
					k=x2+width*y2+z*nPixels;
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
					k=x0+width*y2+z*nPixels;
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
				k=x2+width*y+z*nPixels;
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
			if(z0>=0) {
				k = x+width*y+z0*nPixels;
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
				k=x0+width*y+z*nPixels;
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
				k = x + width * y0+z*nPixels;
				if (imArray[k] > imArray[j]) {
					adjNode = findNode(k);
					if (curNode != adjNode) {
						curNode = mergeNodes(adjNode,curNode);
						found = true;
					}
				}
				if(x2<width)
				{
					k=x2+width*y0+z*nPixels;
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
					k=x0+width*y0+z*nPixels;
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
				/*****Debug***
				if(j==13979) {
					System.out.print("neighbor: "+voxSumN[j]+" "+areasN[j]+" self: "+voxSum[j]+" "+areas[j]+" "+"\n");
				}***/
				y0=Math.max(y-1,0);
				y2=Math.min(y+1, height-1);
				x0=Math.max(x-1,0);
				x2=Math.min(x+1,width-1);
				z0=Math.max(z-1,0);
				z2=Math.min(z+1,zSlice-1);
				// for neighboring pixels' value we consider 26 neighbors
				for (ii=z2;ii>=z0;ii--) {
					for (jj=y2;jj>=y0;jj--) {
						for(kk=x2;kk>=x0;kk--) {
							if( ii==z & jj==y & kk==x)
								continue;
								tmpIdx = kk+width*jj+ii*nPixels;
								if (usedN[tmpIdx] == NOT_USED_AS_N)
								{
									voxSumN[j] += imArray[tmpIdx];
									areasN[j]++;
									usedN[tmpIdx] = USED_AS_N_ONCE;
									/*if(j==13979) {
										System.out.print("neighbor val: "+imArray[tmpIdx]+" "+ii +" "+jj+" "+kk+"\n");
									}*/
								}
						}
					}
				}
				usedN[j] = USED_AS_N_MORE;
				diffN[j] = voxSum[j]/(double)areas[j];
				if (areasN[j] > 0)
					diffN[j] -= voxSumN[j]/(double)areasN[j];
			}
			
		}
		ObjLabel(p.min_size,p.max_size);
		zscore = new double[nVoxels];
		for (i=0; i<nVoxels; i++)
		{
			rmder = i % nPixels;
			z = i / nPixels;
			y=rmder/width;
			x=rmder-y*width;
			double LH = (double)BxCor[i][4]-BxCor[i][1]+1;
			double LW = (double)BxCor[i][3]-BxCor[i][0]+1;
			double LZ = (double)BxCor[i][5]-BxCor[i][2]+1;
			double ratio = LH>LW? LH/LW: LW/LH;
			/* for debug
			if(i==13775) {
				System.out.println(""+outputArray[i]+" "+areas[i]+" "+voxSum[i]+" "+areasN[i]+" "+voxSumN[i]+" ");
			}*/
			if(areas[i]<p.min_size || ratio>p.maxWHratio || (areas[i]/(double)(LH*LW*LZ))<p.minfill){
				zscore[i] = 0;
				continue;
			}
	
			if(outputArray[i] == OBJECT)
			{
				zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p,q.var);
			}
		}

		
		/******Debug for the correctness of component tree building****
		for(i=2; i<6;i++) {
			for(j = 21;j<30;j++){
				for(kk = 21; kk<30;kk++) {
					tmpIdx = kk+j*width+i*nPixels;
					adjNode=findNode(tmpIdx);
					System.out.print(adjNode+" "+areas[adjNode]+"; ");
					//System.out.print(imArray[tmpIdx]+"; ");
				}
			}
			System.out.print("\n\n");
		}
		System.out.print("Debug in Component Tree Done!\n");***/
	}
	/*Label object or not: here we simply  test the size  to decide object or not*/
	public void ObjLabel(int minSize ,int maxSize)
	{
		int i,j;
		int nVoxels = sortedIndex.length;
		for (i = nVoxels-1; i >= 0; i--)
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
	public int mergeNodes(int e1,int e2)/*e1 adjacent node, e2 current node*/
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
		int curNeiCnt = areasN[res];
		if (curNeiCnt==53)
			curNeiCnt = 53;
		/*****Debug****
		if(res==13979) {
			System.out.print("neighbor: "+voxSumN[res]+" "+areasN[res]+" self: "+voxSum[res]+" "+areas[res]+" "+"\n");
		}*/
		//Compute new neighbours
		int z = e2 / (width*height);
		int rmder = e2 % (width*height);
		int y=rmder/width;
		int x=rmder-y*width;
		
		int y0=Math.max(y-1,0);
		int y2=Math.min(y+1, height-1);
		int x0=Math.max(x-1,0);
		int x2=Math.min(x+1,width-1);
		int z0=Math.max(z-1,0);
		int z2=Math.min(z+1,zSlice-1);
		// for neighboring pixels' value we consider 26 neighbors
		//System.out.print("Before Merging"+voxSumN[res]+" "+areasN[res]+" "+"\n");
		int ii, jj, kk, tmpIdx;
		for (ii=z2;ii>=z0;ii--) {
			for (jj=y2;jj>=y0;jj--) {
				for(kk=x2;kk>=x0;kk--) {
					if( ii==z & jj==y & kk==x)
						continue;
					tmpIdx = kk+width*jj+ii*(width*height);
					if (usedN[tmpIdx] == NOT_USED_AS_N)
					{
						voxSumN[res] += imArray[tmpIdx];
						areasN[res]++;
						usedN[tmpIdx] = USED_AS_N_ONCE;
					}
				}
			}
		}
		/*****Debug***
		if(res==136)
			System.out.print("e2: "+areasN[e2]+ "e1: "+areasN[e1]+"\n");
		**/
		areasN[res] += areasN[m];
		voxSumN[res] += voxSumN[m];
		if (usedN[e2] == USED_AS_N_ONCE) // e2 ever be used as neighbors, now we need to remove them
		{
			areasN[res] -= areas[e2];
			voxSumN[res] -= voxSum[e2];
		}
		/*****Debug***
		if(res==136)
			System.out.print("Before Merging res: "+curNeiCnt+", but after Merging res: "+areasN[res]+"\n");
		*/
		areas[res] += areas[m];
		voxSum[res] += voxSum[m];
		parNode[m]=res;

		usedN[e2] = USED_AS_N_MORE;

		diffN[res] = voxSum[res]/(double)areas[res];

		if (areasN[res] > 0)
			diffN[res] -= voxSumN[res]/(double)areasN[res];

		//System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
		//System.out.println("BC:" +BxCor[m][3]+" "+BxCor[m][2]+" "+BxCor[m][1]+" "+BxCor[m][0]);
		if (BxCor[res][0] > BxCor[m][0])
			BxCor[res][0] = BxCor[m][0];
		if (BxCor[res][1] > BxCor[m][1])
			BxCor[res][1] = BxCor[m][1];
		if (BxCor[res][2] > BxCor[m][2])
			BxCor[res][2] = BxCor[m][2];
		if (BxCor[res][3] < BxCor[m][3])
			BxCor[res][3] = BxCor[m][3];
		if (BxCor[res][4] < BxCor[m][4])
			BxCor[res][4] = BxCor[m][4];
		if (BxCor[res][5] < BxCor[m][5])
			BxCor[res][5] = BxCor[m][5];
		//System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
		return res;
	}
	/*Calculate the z-score of each pixel*/
	public double zscoreCal(double t0, int M/*in*/, int N/*nei*/, paraP3D p, double qVar){
		double[][] pMu = p.mu;
		double[][] pSigma = p.sigma;
		
		if(M < p.min_size)
			M = p.min_size;
		if(N<=p.min_size || N < (M/10))
			return 0;
		if(N<p.min_size)
			N = p.min_size;
		if(M>p.max_size)
			M = p.max_size;
		if(N>p.max_size)
			N = p.max_size;
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
	public void cpt2map(paraP3D p){
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		kMapx = new int[Nthr][zSlice][height][width];
		boolean[][][][] bMapx = new boolean[Nthr][zSlice][height][width];
		zMapx = new double[Nthr][zSlice][height][width];
		int nVoxels = width*height*zSlice;
		ImageHandling IH = new ImageHandling();

		Zscore_Vec = new ArrayList<Double>();
		levels = new ArrayList<Integer>();
		for (int i = nVoxels-1; i >=0; i--)
		{
			int j=sortedIndex[i];
			if(outputArray[j] != OBJECT)
				continue;

			//Compute new neighbours
			int z = j / (width*height);
			int rmder = j % (width*height);
			int y=rmder/width;
			int x=rmder-y*width;

			int e = j;
			int intensity = imArray[e];
			if(intensity<p.thr0)
				continue;
			int Inl = (intensity-p.thr0)/p.thrg;
			int Inlevel = Inl;
			zMapx[Inlevel][z][y][x] = zscore[e];
			bMapx[Inlevel][z][y][x] = true;

			while (parNode[e] != e){
				intensity = imArray[parNode[e]];
				if(intensity<p.thr0 | outputArray[parNode[e]] != OBJECT){
					break;
				}
				int InlevelPar = (intensity-p.thr0)/p.thrg;
				zMapx[InlevelPar][z][y][x] = zscore[parNode[e]];
				bMapx[InlevelPar][z][y][x] = true;
				e = parNode[e];
			}
		}
		int lvCnt = 0;
		BasicMath bm = new BasicMath();
		for(int i=0;i<Nthr;i++){
			System.out.println("Thresholded region's zscore: "+bm.matrix3DMax(zMapx[i])); 
			kMapx[i] = IH.bwlabel3D(zMapx[i], 26,Zscore_Vec);
			for(int j=lvCnt;j<Zscore_Vec.size();j++){
				levels.add(i);
			}
			lvCnt = Zscore_Vec.size();
		}
	}
	/* The same function as BestSyn.extractBestSyn_fdr. Extract the region with highest zscore. 
	 * Now we just need to scan the ArrayList Zscore_Vec and levels. No need to go through all maps.*/
	public void extractBestSyn(paraP3D p){
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		int nVoxels = width*height*zSlice;
		boolean[] checkedMark = new boolean[nVoxels];
		System.out.println("size of candidate particle set: "+Zscore_Vec.size());
		if(Zscore_Vec==null){//if the vector haven't been build, this will not be used
			Zscore_Vec = new ArrayList<Double>();
			levels = new ArrayList<Integer>();
			for (int i = 0; i <=nVoxels; i++)
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