import java.lang.*;


//Parameter initialization, mainly image properties
public class ParaQ{
	public double maxVal;  // image max intensity
	public int Pix_per_synapse;  // image channel for dendrite length calculation
    //public double [][] Gt; // a stablized version of original image G.
    public double[] pgEst;
    public double var = 0;
    int Nx,Ny;
    int ntry;
    public ParaQ(short [] pixels,double Lx,int inntry) {
    	int length = pixels.length;
    	double tmpMax=0;
    	for(int i = 0;i<length;i++)
    		if(tmpMax<pixels[i])
    			tmpMax = pixels[i];

    	maxVal = (double) tmpMax;
    	
    	//double Lx = 2.076*1e-7;
    	double Lx2 = (Lx*1e6)*(Lx*1e6);
    	Pix_per_synapse = (int)Math.round(1/Lx2);
    	ntry = inntry;
    }

}