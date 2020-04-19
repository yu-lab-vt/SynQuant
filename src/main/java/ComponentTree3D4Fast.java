import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
/**
 * This class build a component tree from a image array. Also contain another two functions "cpt2map" which transfers component tree to 
 * zscore maps and "extractBestSyn" which extracts highest-score region.
 * 
 * Major idea of tree building and zscore calculation is from Dr.Petter Ranefall(petter.ranefall@it.uu.se).
 * For details of creating component tree, please refer to "Fast Adaptive Local Thresholding Based on Ellipse fit"(ISBI'16) and 
 * "Building the Component Tree in Quasi-Linear Time" (TIP 2006 Vol.15 No.11).
 * @version 1.0
 * @date 2016-07-11
 * Handle the region with large number of voxels
 * @version 1.1
 * @date 2019-06-08
 * 
 * @contributors Congchao Wang, Boyu Lyu, Petter Ranefall
 * @contact ccwang@vt.edu, boyu93@vt.edu, petter.ranefall@it.uu.se
 */


public class ComponentTree3D4Fast{
	protected int[] sortedIndex;//store the elements in a new ordered array
	protected int[] parNode; //Parent nodes
	
	//for extract highest zscore
	protected ArrayList<Double> Zscore_Vec;//zscore vector
	protected ArrayList<Integer> zScoreIdx;
	protected ArrayList<Integer> levels;

	/***for component tree build***/
	// type 1: saved after initialization
	short[] imArray;
	double[] outputArray;
	double[] zscore;
	int width, height, zSlice;
	int minSize, maxSize;
	int nPixels, nVoxels;
	int cpt_nei; // neighbors considered
	// type 2: deleted after initialization
	double[] diffN; //Neighborhood difference
	double[] voxSum; //Pixel sum inside objects
	double[] voxSumN; //Pixel sum in neighbourhood
	int[] areas; //Areas inside objects
	int[] areasN; //Areas of neighborhoods
	byte[] usedN; //Flag if used in neighborhood
	int[][] BxCor;//each area's bounding box coordinates

	//CONSTANTS
	public static final double UNDEFINED = 0;
	public static final double NOT_OBJECT = -1;
	public static final double OBJECT = 1;
	public static final byte NOT_USED_AS_N = (byte)0;
	public static final byte USED_AS_N_ONCE = (byte)1;
	public static final byte USED_AS_N_MORE = (byte)2;
	//////////////////////Finished/////////////////////////////////////
	/* 1 build the component tree of input image array (8-connected)
	 * 2 label objects by size test
	 * 3 calculate z-score of each object pixel */
	public ComponentTree3D4Fast(short[][] imArrayIn,byte[] usedNIn, paraP3D p, paraQ3D q) 
	{
		zSlice = imArrayIn.length;
		width = q.Nx;
		height = q.Ny;
		minSize = p.min_size;
		maxSize = p.max_size;
		cpt_nei = p.cpt_nei;
		nVoxels=width*height*zSlice; // voxels in a 3D image
		nPixels = width*height; //pixels in a single slice
		imArray = new short[nVoxels];
		int tmpCnt = 0;
		
		for(int i = 0; i < imArrayIn.length; i++) {
			for(int j = 0; j < imArrayIn[0].length; j++)
				imArray[tmpCnt++] = imArrayIn[i][j]; //java.util.Arrays
		}
		//Create nodes
		parNode=new int[nVoxels];
		BxCor = new int[nVoxels][6];//ymin,xmin,zmin,  ymax,xmax,zmax.
		int x,y,z, rmder, x0,x2,y0,y2,z0, z2;
		int i,j,k;
		outputArray = new double[nVoxels];
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
			usedN[i] = NOT_USED_AS_N;
			//usedN[i] = usedNIn[i];
			BxCor[i][0] = y;
			BxCor[i][1] = x;
			BxCor[i][2] = z;
			BxCor[i][3] = y;
			BxCor[i][4] = x;
			BxCor[i][5] = z;
		}

		/***** Sort points:
		* create a counting array, counts, with a member for 
		* each possible discrete value in the input.  
		* initialize all counts to 0.
		* ***/
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
		//endTime1=System.nanoTime();
		//System.out.println("Initialization running time: "+(endTime1-startTime1)/1e9);
		//Search in decreasing order
		//startTime1=System.nanoTime();
		int curNode;
		int adjNode;
		int ii,jj, kk, tmpIdx;
		boolean found;
		for (i = nVoxels-1; i >= 0; i--)
		{
			j=sortedIndex[i];
			//if (j==37)
			//	System.out.println(imArray[j]);
			curNode=j;
			//System.out.println("Idx "+j+" Image Value"+imArray[j]);
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
			/*** for debug*
			if (imArray[j]>0) {
				System.out.print(j+" " + imArray[j]+"\n");
			}
			if(imArray[j]>=255)//usedN[28+width*28+6*nPixels] != NOT_USED_AS_N)
				System.out.println(" "+x+" "+y+" "+z);
			***/
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
				if (cpt_nei == 26) {// if we believe neighbors in z-direction
					z0=Math.max(z-1,0);
					z2=Math.min(z+1,zSlice-1);
				}else { //default: we do not fully believe neighbors in z-direction
					z0=z;
					z2=z;
				}
				// for neighboring pixels' value
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

		zscore = new double[nVoxels];
		Zscore_Vec = new ArrayList<Double>();
		zScoreIdx = new ArrayList<Integer>();
		for (j=nVoxels-1; j>=0; j--)
		{
			i=sortedIndex[j];
//			rmder = i % nPixels;
//			z = i / nPixels;
//			y=rmder/width;
//			x=rmder-y*width;
//			if ((x<=17 & x>=15) & (y<=18 & y>=16) & z==0) {
//				int e = i;
//				System.out.println(e+" summary: "+z+" "+y+" "+x+" "+"score: "+zscore[e]+" Other "+areas[e]+" "+voxSum[e]+" "+areasN[e]+" "+voxSumN[e]+" ");
////				while (imArray[e] > 1 & e != parNode[e]) {
////					rmder = e % nPixels;
////					z = e / nPixels;
////					y=rmder/width;
////					x=rmder-y*width;
////					System.out.println(e+" summary: "+z+" "+y+" "+x+" "+"score: "+zscore[e]+" Other "+areas[e]+" "+voxSum[e]+" "+areasN[e]+" "+voxSumN[e]+" ");
////					e = parNode[e];
////				}
//			}
			double LH = (double)BxCor[i][4]-BxCor[i][1]+1;
			double LW = (double)BxCor[i][3]-BxCor[i][0]+1;
			double LZ = (double)BxCor[i][5]-BxCor[i][2]+1;
			double ratio = LH>LW? LH/LW: LW/LH;
			/* for debug
			if(y==276 && x==73) {
				System.out.println(""+outputArray[i]+" "+areas[i]+" "+voxSum[i]+" "+areasN[i]+" "+voxSumN[i]+" ");
			}*/
			if(areas[i]>=p.max_size) {
				zscore[i] = -1;
				outputArray[i] = NOT_OBJECT;
			}else if(areas[i]<p.min_size || ratio>p.maxWHratio || (areas[i]/(double)(LH*LW*LZ))<p.minfill) {
				zscore[i] = -1;
			}else {
//				rmder = i % nPixels;
//				z = i / nPixels;
//				y=rmder/width;
//				x=rmder-y*width;
//				if (i==6640) {//(x<=26 & x>=25) & (y<=29 & y>=27) & z==0) {
//					int e = i; 
//					System.out.println("summary: "+z+" "+y+" "+x+" "+"score: "+zscore[e]+" Other "+areas[e]+" "+voxSum[e]+" "+areasN[e]+" "+voxSumN[e]+" ");
//				}
				zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p,q.var);
				Zscore_Vec.add(zscore[i]);
				zScoreIdx.add(i);
			}
			//System.out.println("summary: "+z+" "+y+" "+x+" "+"score: "+zscore[i]+" Other "+areas[i]+" "+voxSum[i]+" "+areasN[i]+" "+voxSumN[i]+" ");
			/*if (zscore[i]>255) {
				zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p,q.var);
				System.out.println("zScore: "+zscore[i]);
			}*/
				
		}
//		i = 51485;
//		while (i != parNode[i]) {
//			System.out.println(i+" summary: score: "+zscore[i]+" Other "+areas[i]+" "+voxSum[i]+" "+areasN[i]+" "+voxSumN[i]+" ");
//			i = parNode[i];
//		}
//		System.out.println(" ----------------- ");
		/***
		 * Generate the zscore map
		 * 1. if there is input zscore threshold use the zscore threshold
		 * 2. otherwise use fdr control with 0.05 threshold
		 */
		if (p.zscore_thres>0) {  // we do not use this
			objLabel_zscore(p.zscore_thres);
		}
		else {
			objLabel_bestlevel();
		}
		////Debug for the correctness of component tree building
//		for(i=2; i<6;i++) {
//			for(j = 21;j<30;j++){
//				for(kk = 21; kk<30;kk++) {
//					tmpIdx = kk+j*width+i*nPixels;
//					adjNode=findNode(tmpIdx);
//					System.out.print(adjNode+" "+areas[adjNode]+"; ");
//					//System.out.print(imArray[tmpIdx]+"; ");
//				}
//			}
//			System.out.print("\n\n");
//		}
//		System.out.print("Debug in Component Tree Done!\n");
		
		//// remove the variables
		diffN = null;
		voxSum = null;
		voxSumN = null;
		areas = null;
		areasN = null;
		usedN = null;
		BxCor = null;
	}
	///// label the objects by z-score threshold set by users
	public void objLabel_zscore(double zscore_thres) {
		int i,j;
		for (i=nVoxels-1; i>=0; i--)
		{
			j=sortedIndex[i];
			if (outputArray[j] == UNDEFINED)
			{
				int e = j;
				while(zscore[e] < zscore_thres && outputArray[e] == UNDEFINED) {
					e = parNode[e];
				}
				if (zscore[e] >= zscore_thres) { // this region is good
					double cur_zscore = zscore[e];
					double cur_label = outputArray[e];
					while(imArray[e] == imArray[parNode[e]] & outputArray[e] == UNDEFINED)
						e = parNode[e];
					if (cur_label == UNDEFINED) { // this region haven't been labeled
						int e1 = j;
						while(e1 != e)
						{
							outputArray[e1] = OBJECT;
							zscore[e1] = cur_zscore;
							e1 = parNode[e1];
						}
						if (outputArray[e1] == UNDEFINED) {
							outputArray[e1] = OBJECT;
							zscore[e1] = cur_zscore;
						}
						e1 = parNode[e];
						while(e1 != parNode[e1] & outputArray[e1] == UNDEFINED)
						{
							outputArray[e1] = NOT_OBJECT;
							zscore[e1] = -1;
							e1 = parNode[e1];
						}
						if (outputArray[e1] == UNDEFINED) {
							outputArray[e1] = NOT_OBJECT;
							zscore[e1] = -1;
						}
						if (outputArray[e1] == OBJECT) { // we did wrong thing if the root is OBJECT
							// under such condition, we need to correct previous from j to e1
							int e2 = j;
							while(e2 != e1)
							{
								outputArray[e2] = OBJECT;
								zscore[e2] = zscore[e1];
								e2 = parNode[e2];
							}
						}
					}
					else{ // this region has been labeled (should only have NOT_OBJECT label)
						int e1 = j;
						while(e1 != e) // label j to e as the same
						{
							outputArray[e1] = cur_label;
							zscore[e1] = cur_zscore;
							e1 = parNode[e1];
						}
					}
				}else { // this region is bad
					int e1 = j;
					while(e1 != e) // label j to e as the same
					{
						outputArray[e1] = outputArray[e];
						zscore[e1] = zscore[e];
						e1 = parNode[e1];
					}
				}
			}
		}
	}
	///// label the objects by finding the best level based on z-score
	public void objLabel_bestlevel() {
		int i,j;
		// sort all pixels with their zscore values with descending order
		int[] zScoreSortedIdx = new int[zscore.length];
		for (i=0; i<zscore.length;i++)
			zScoreSortedIdx[i] = i;
		ArrayList<Double> tmpZV = new ArrayList<Double>(Zscore_Vec);
		Collections.sort(tmpZV, Collections.reverseOrder());
		for(i=0;i<tmpZV.size();i++){
			zScoreSortedIdx[i] = zScoreIdx.get(Zscore_Vec.indexOf(tmpZV.get(i)));
			zScoreSortedIdx[zScoreSortedIdx[i]] = i;
		}

		// label the objects by finding the best level based on z-score
		for (i = 0; i < zScoreSortedIdx.length; i++)
		{
			j=zScoreSortedIdx[i];
			if (outputArray[j] == UNDEFINED)
			{
				int e = j;
				//System.out.println("Idx "+j+" Image Value"+imArray[j]+" "+zscore[j]);
				while(imArray[e] == imArray[parNode[e]] && outputArray[e] == UNDEFINED)
					e = parNode[e];
				if (outputArray[e] == UNDEFINED)
				{
					int e1 = parNode[e];				
					//					if (e1!=e && zscore[e]>0) {
					//						outputArray[e] = OBJECT;
					//						continue;
					//					}
					double ObjOrNot = NOT_OBJECT;
					double tmpZscore = -1;
					while(e1 != parNode[e1])
					{
						if (outputArray[e1]!=UNDEFINED) {
							ObjOrNot = outputArray[e1];
							tmpZscore = zscore[e1];
							break;
						}
						e1 = parNode[e1];
					}


					if (ObjOrNot==OBJECT) { //label j to e1 to object
						e1 = j;
						while(e1 != parNode[e1] && outputArray[e1] == UNDEFINED){
							outputArray[e1] = OBJECT;
							//diffN[e1] = diffN[e];
							zscore[e1] = tmpZscore;
							e1 = parNode[e1];
						}
						if (outputArray[e1] == UNDEFINED) {
							outputArray[e1] = OBJECT;
							zscore[e1] = tmpZscore;
						}
						//						if (outputArray[e] == UNDEFINED) {
						//							outputArray[e] = OBJECT;
						//							zscore[e] = tmpZscore;
						//						}
					}
					else { // label parNode[e] to e1 NOT_OBJECT, label j to e object (only e and j is enough)
						e1 = parNode[e];
						while(e1 != parNode[e1] && outputArray[e1] == UNDEFINED){
							outputArray[e1] = NOT_OBJECT;
							//diffN[e1] = diffN[e];
							zscore[e1] = -1;
							e1 = parNode[e1];
						}
						if (outputArray[e1] == UNDEFINED) {
							outputArray[e1] = NOT_OBJECT;
							zscore[e1] = -1;
						}
						outputArray[j] = OBJECT;
						outputArray[e] = OBJECT;
						zscore[e] = zscore[j];
					}

					//findBestLevel(j, -1);
				}
				else
				{
					int e1 = j;
					while(e1 != e)
					{
						outputArray[e1] = outputArray[e];
						//diffN[e1] = diffN[e];
						zscore[e1] = zscore[e];
						e1 = parNode[e1];
					}
				}
			}
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
//		if(e1==6640 | e2==6640) {
//			System.out.print("e1:"+imArray[e1]+" e2:"+imArray[e1] +" "+imArray[e2]+"\n");
//		}
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
		//int curNeiCnt = areasN[res];
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
		int z0,z2;
		if (cpt_nei == 26) {// if we believe neighbors in z-direction
			z0=Math.max(z-1,0);
			z2=Math.min(z+1,zSlice-1);
		}else { //default: we do not fully believe neighbors in z-direction
			z0=z;
			z2=z;
		}

		
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
		//System.out.println("parentVal:"+imArray[res]+" sonVal:"+imArray[m]);
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
	// find the best threshold for each location
	public int findBestLevel(int e, int optE)
	{
		//if (e==82)
		//	System.out.println(imArray[e]);
		double optDiffN;
		if (optE >= 0)
		{
			optDiffN = zscore[optE];
		}
		else
		{
			optDiffN = 0;
		}
		int startE = e;
		while(imArray[e] == imArray[parNode[e]] && parNode[e] != e)
			e = parNode[e];
		if (outputArray[e] > 0)//==OBJECT
		{
			return e;
		}
		if (outputArray[e] == NOT_OBJECT)
		{
			return optE;
		}
		if (parNode[e] == e)
		{
			outputArray[e] = NOT_OBJECT;
			return optE;
		}
		//System.out.println("size "+areas[e]+"zscore "+zscore[e]);
		if (areas[e] > maxSize)
		{
			outputArray[startE] = NOT_OBJECT;
			while(startE != parNode[startE] && outputArray[parNode[startE]] == UNDEFINED)
			{
				startE = parNode[startE];
				outputArray[startE] = NOT_OBJECT;
			}
			return optE;
		}
		if (areas[e] >= minSize && zscore[e] > optDiffN)
		{
			optDiffN = zscore[e];
			optE = e;
		}
		optE = findBestLevel(parNode[e], optE);
		if (optE >= 0)
		{
			optDiffN = zscore[optE];
		}
		else
		{
			optDiffN = 0;
		}
		if (optE >= 0 && imArray[e]  >= imArray[optE]) //Opt found below
		{
			while(startE != e)
			{
				outputArray[startE] = OBJECT;
				zscore[startE] = zscore[optE];
				startE = parNode[startE];
			}
			outputArray[e] = OBJECT;
			zscore[e] = zscore[optE];
			//Zscore_Vec.add(zscore[optE]);
			return optE;
		}
		
		// Below best level
		outputArray[startE] = NOT_OBJECT;
		
		while(startE != e)
		{
			startE = parNode[startE];
			outputArray[startE] = NOT_OBJECT;
		}
		return optE;	
	}
	/*Calculate the z-score of each pixel*/
	public double zscoreCal(double t0, int M/*in*/, int N/*nei*/, paraP3D p, double qVar){
		double[][] pMu = p.mu;
		double[][] pSigma = p.sigma;
		// our look up table start from 10 pixel, so the min_size = max(10, p.min_size);
		if(M < Math.max(p.min_size, 10))
			M = Math.max(p.min_size, 10);
		if(N<=p.min_size || N < (M/10))
			return -1;
		if(N < Math.max(p.min_size, 10))
			N = Math.max(p.min_size, 10);
		
		double mu, sigma;
//		double scalefactor = 1;
//		while(M>=pMu.length || N>=pMu[0].length) {
//			M = M/10;
//			scalefactor = scalefactor*Math.sqrt(10);
//		}
//		mu = pMu[M-1][N-1]*scalefactor;
//		
//		scalefactor = 1;
//		while(N>=pMu[0].length) {
//			
//		}
/**way 1: truncated gaussian for approximation**/
		double sigmaScl = 1;
		if (M>=pMu.length || N>=pMu[0].length) { //approximation of super large particle
		    sigmaScl = Math.sqrt((double)(M+N)/500);
		    M = (int) Math.floor(((double)M)/(M+N)*500);
		    N = 500 - M;//If we use previous function, the approximation is larger than expected
		}       
		mu = pMu[M-1][N-1];
		sigma = pSigma[M-1][N-1];
		mu = mu*Math.sqrt(qVar);
		sigma = sigma*Math.sqrt(qVar)/sigmaScl;
		double zScore = (t0-mu)/sigma;
		return zScore;
/**way 2: use integral to re-estimate parameters based on order-statistics**/
//		if (M>=pMu.length || N>=pMu[0].length) {
//			mu = p.CalMu(M,N,1000);
//			sigma = p.CalSigma(M,N,1000);
//		}else {
//			mu = pMu[M-1][N-1];
//			sigma = pSigma[M-1][N-1];
//		}
//
//		mu = mu*Math.sqrt(qVar);
//		sigma = sigma*Math.sqrt(qVar);
//		return (t0-mu)/sigma;
	}
	
	
	public ArrayList<Double> removeDuplicates(ArrayList<Double> list) {
        // Store unique items in result.
        ArrayList<Double> result = new ArrayList<Double>();

        // Record encountered Strings in HashSet.
        HashSet<Double> set = new HashSet<Double>();

        // Loop over argument list.
        for (Double item : list) {

            // If String is not in set, add it to the list and the set.
            if (!set.contains(item)) {
                result.add(item);
                set.add(item);
            }
        }
        return result;
    }
}
