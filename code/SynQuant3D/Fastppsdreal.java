import java.util.ArrayList;
import java.util.Collections;
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


public class Fastppsdreal{
	protected int[] sortedIndex;//store the elements in a new ordered array
	protected int[] parNode; //Parent nodes
	
	protected boolean[][] SynR1;
	protected int nSyn0;
	//for prior knowledge test
	int[][] BxCor;//each area's bounding box coordinates
	int minSize, maxSize;
	double maxWHratio, minfill;
	
	
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
	public Fastppsdreal(short[] imArrayIn, int w, int h,ParaP p, ParaQ q) 
	{
		imArray = imArrayIn;
		width = w;
		height = h;
		minSize = p.min_size;
		maxSize = p.max_size;
		maxWHratio = p.maxWHratio;
		minfill = p.minfill;
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
			usedN[i] = NOT_USED_AS_N;
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
		//ObjLabel();
		zscore = new double[nPixels];
		for (i=0; i<nPixels; i++)
		{
			/*y=i/width;
			x=i-y*width;
			double LH = (double)BxCor[i][3]-BxCor[i][1]+1;
			double LW = (double)BxCor[i][2]-BxCor[i][0]+1;
			double ratio = LH>LW? LH/LW: LW/LH;
			if(areas[i]<p.min_size || ratio>p.maxWHratio || (areas[i]/(double)(LH*LW))<p.minfill){
				zscore[i] = 0;
				continue;
			}
			else{
				zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p.mu, p.sigma,q.var);
			}*/
			zscore[i] = zscoreCal(diffN[i],areas[i],areasN[i],p.mu, p.sigma,q.var);
		}
		/*scan */
		for (i = nPixels-1; i >= 0; i--)
		{
			j=sortedIndex[i];
			if (outputArray[j] == UNDEFINED)
			{
				int e = j;
				while(imArray[e] == imArray[parNode[e]] && outputArray[e] == UNDEFINED)
					e = parNode[e];
				if (outputArray[e] == UNDEFINED)
				{
					findBestLevel(j, -1);
				}
				else
				{
					int e1 = j;
					while(e1 != e)
					{
						outputArray[e1] = outputArray[e];
						zscore[e1] = zscore[e];
						e1 = parNode[e1];
					}
				}
			}
		}
		double[][] zMap = new double[height][width];
		for (i = 0; i < height; i++)
			for (j = 0; j < width; j++)
				zMap[i][j] = zscore[i * width + j] ;// the data is saved by line
		//ImageHandling IH = new ImageHandling();
		ArrayList<Double> zscore_vec = new ArrayList<Double>();
		//int[][] kMap = IH.bwlabel(zMap, 8, zscore_vec);
		
		double threshold = FDR_Control(zscore_vec,p);
		SynR1 = new boolean [height][width];
		for (i = 0; i < height; i++)
			for (j = 0; j < width; j++)
				if(zMap[i][j]>threshold && outputArray[i * width + j]==OBJECT)
					SynR1[i][j] = true;
		//IH.bwlabel(SynR1,8);
		//nSyn0 = IH.NextLabel;
	}
	/*Label object or not*/
	public int findBestLevel(int e, int optE)
	{
		double optZscore;
		if (optE >= 0)
		{
			optZscore = zscore[optE];
		}
		else
		{
			optZscore = 0;
		}
		int startE = e;
		while(imArray[e] == imArray[parNode[e]] && parNode[e] != e)
			e = parNode[e];
		if (outputArray[e] == OBJECT)
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

		double LH = (double)BxCor[e][3]-BxCor[e][1]+1;
		double LW = (double)BxCor[e][2]-BxCor[e][0]+1;
		double ratio = LH>LW? LH/LW: LW/LH;
		if (areas[e] >= minSize && zscore[e] > optZscore && ratio<=maxWHratio && (areas[e]/(double)(LH*LW))>=minfill)
		{			
			optZscore = zscore[e];
			optE = e;
		}
		optE = findBestLevel(parNode[e], optE);
		if (optE >= 0)
		{
			optZscore = zscore[optE];
		}
		else
		{
			optZscore = 0;
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
	
	public double FDR_Control(ArrayList<Double> z_vec, ParaP p) {
		/*extract the best map*/
		BasicMath mBM = new BasicMath();

		/*FDR Control*/
		ArrayList<Double> pvalue = new ArrayList<Double>();
		for(int i=0;i<z_vec.size();i++){
			pvalue.add(mBM.zTop(z_vec.get(i)));
		}
		int total_len = pvalue.size();
		double[] FDR_Con = new double[total_len];
		for(int i=0;i<total_len;i++)
			FDR_Con[i] = p.fdr*((double)(i+1))/total_len;
		if(p.fdr_den){
			for(int i=0;i<total_len;i++){
				FDR_Con[i] = FDR_Con[i]/Math.log(total_len);
			}
		}
		Collections.sort(pvalue);// Ascend; add Collections.reverseOrder() for descend
		double OK_pvalue = 0;
		for(int i=pvalue.size()-1;i>=0;i--){
			if(pvalue.get(i)<FDR_Con[i]){
				OK_pvalue = pvalue.get(i);
				break;
			}
		}
		return mBM.pToZ(OK_pvalue);
	}
}
