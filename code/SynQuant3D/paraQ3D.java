import java.lang.*;


//Parameter initialization, mainly image properties
public class paraQ3D{
	public double maxVal;  // image max intensity
	public int Pix_per_synapse;  // image channel for dendrite length calculation
    //public double [][] Gt; // a stablized version of original image G.
    public double[] pgEst;
    public double var;
    int Nx,Ny, Nz;
    int ntry;
    public paraQ3D(short [][] pixels, double Lx,int inntry, int height, int width) {
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

}