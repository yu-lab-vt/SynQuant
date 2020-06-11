import java.lang.*;


//Parameter initialization, mainly image properties
public class paraQ3D{
	public double maxVal;  // image max intensity
	public int Pix_per_synapse;  // image channel for dendrite length calculation
    //public double [][] Gt; // a stablized version of original image G.
    public double[] pgEst;
    public double var = 0;
    public int Nx,Ny, Nz;
    public int ntry;
    public double[][][][] synZscore = null; // if we have some prior infor for z-score, this is generally from pre-channel
    public int _NumChannel;
    public int _wayCombinePrePost;
    public int NumChannelProcessed = 0; // if we have multi-channel, when NumChannelProcessed==_NumChannel, we can do dfr control
    public int curTps = 0; // indicate which time point we are processing if each channel is a video. start from 0 to t-1
    public double varRatio = 0.05; // use the 0.1 location in y=a*x+b to estimate the noise
    public double varRatioBased = 0;
    
    public double ExtendedDistance=0; // extended distance
    public double zAxisMultiplier=1; // z axis extended distance multiplier
    
    public void init(short [][] pixels, double Lx,int inntry, int height, int width) {
    	double tmpMax=0;
    	for(int i = 0;i<pixels.length;i++)
    		for(int j = 0; j < pixels[0].length; j++)
	    		if(tmpMax<pixels[i][j])
	    			tmpMax = pixels[i][j];

    	maxVal = (double) tmpMax;
		Ny = height;
		Nx = width;
    	Nz = pixels.length;
    	//double Lx = 2.076*1e-7;
    	double Lx2 = (Lx*1e6)*(Lx*1e6);
    	Pix_per_synapse = (int)Math.round(1/Lx2);
    	ntry = inntry;
    }
    public paraQ3D(int NumChannel, int wayCombinePrePost, double varRatioIn){
    	_wayCombinePrePost = wayCombinePrePost;
    	_NumChannel = NumChannel;
    	varRatio = varRatioIn;
    };

}