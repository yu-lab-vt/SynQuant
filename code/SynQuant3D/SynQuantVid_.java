import java.awt.Color;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import ij.*;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.*;
import ij.process.*;

/**
 * This plugin 
 * 1. segments synapse/puncta using Probability principled synapse detection
 * algorithm modified from "PPSD: Probability Principled Synapse Detection"
 * 2. extract dendrite from dendrite channel and use linear regression to anaylyse
 * relationships between synapse and 
 * 
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @version 1.0
 * @date 2016-05-31
 *
 *
 */

public class SynQuantVid_ implements PlugIn {
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp; //Current Image
	protected int ntry=1; //decide max number of neighbors considered in zscore calculation
	protected int width; //image width
	protected int height; //image height
	protected int timePts; //image height
	protected int zSlice; //image height
	protected int numChannels; //image height
	protected int presyn_chl;
	protected int postsyn_chl;
	protected int den_chl;
	protected double fdr;
	protected int MinSize,MaxSize;
	protected int NumSynSite=0;
	protected ImagePlus outputImp=null;
	
	boolean fastflag = false;

	public void run(String arg) {
		try {
			if (showDialog())
				synQuant3D_real();
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}
	public boolean showDialog() 
	{
		// Get input parameter
		GenericDialog gd = new GenericDialog("3D Particles - Data and Parameter Setting");
		gd.addNumericField("FDR Control for Particle Detection: ", 0.05, 3);//2.5-3.5 good
		gd.addNumericField("Min Particle Size: ", 20, 0);//2.5-3.5 good
		gd.addNumericField("Max Particle Size: : ", 500, 0);//2.5-3.5 good
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		fdr = gd.getNextNumber();
		MinSize = (int) gd.getNextNumber();
		MaxSize = (int) gd.getNextNumber();
		if(Double.isNaN(fdr) || (fdr<=0) || (fdr>=1)){
			IJ.showMessage("Invalid parameter(s).\n" + "0-1 needed.");
			return false;
		}
		return true;
	}
	public void synQuant3D_real() {
		imp = WindowManager.getCurrentImage();
		stack = imp.getStack();
		
		timePts = imp.getNFrames();
		numChannels = imp.getNChannels();
		zSlice = imp.getNSlices();
		width = imp.getWidth();
		height = imp.getHeight();
		
		if (zSlice==1) {//imp.isHyperStack() == false){
			IJ.log("Current only support 3D data!");
			return;
		}
		double vox_x = imp.getCalibration().pixelWidth; 
		if(vox_x==1)//simulated data
			vox_x = 2.0757e-7;
		else
			vox_x = vox_x*1e-6;//real data
		
		int type = imp.getType();
		
		ppsd3D particle3D_det = null;
		//String Roititle = "3D particle detection";
		outputImp = IJ.createHyperStack("Found Particles", width, height, numChannels, zSlice, timePts,16/*bitdepth*/);
				//stack.getBitDepth());
		outputImp.show();
		String saveFilePath = "C:\\Users\\Congchao\\Desktop\\Probjects_Google_Drive\\SynQuant3D_\\src\\res\\";
		//BasicMath bm = new BasicMath();
		long startTime1=System.nanoTime();
		//particle detection
		for (int i=1; i <= timePts; i++){
			short[][] Arr3D = stack2array(type, stack, i); // #zstack*#pixels in one slice
			particle3D_det = new ppsd3D(Arr3D, width, height,vox_x, fdr,fastflag,MinSize,MaxSize);
			int tmpMax = particle3D_det.nSyn0;//bm.matrix3DMax(particle3D_det.SynR1Idx);
			System.out.println("Number of Synapse: "+tmpMax);
			//display data
			SynapticDisplay(particle3D_det.SynR1Idx, i);
			//write results in txt e.g.z-score
			zScoreWriter(particle3D_det.zscoreList, saveFilePath);
			outputImp.updateAndDraw();
		}
		long endTime1=System.nanoTime();
		System.out.println("Finished. Data size: "+zSlice+" * "+height+" * "+width+" * "+timePts+" Total running time: "+(endTime1-startTime1)/1e9); 
	}
	// save the z-score, corresponding particle index and threshold into txt files 
	public void zScoreWriter(ArrayList<double[]> zscoreList, String fileRootName) {
		String zscoreTxt = fileRootName+"zScore.txt";
		String thresTxt = fileRootName+"threshold.txt";
		String pIdxTxt = fileRootName+"particleIndex.txt";
		
        try{
            // Create index file
            File file = new File(pIdxTxt);
            FileWriter fw;
            // If file doesn't exists, then create it
            if (!file.exists()) {
            	file.createNewFile();    
            	fw = new FileWriter(file.getAbsoluteFile());
            }
            else {
            	fw = new FileWriter(file.getAbsoluteFile(), true);
            }

            BufferedWriter bw = new BufferedWriter(fw);

            // Write index file
            String numberAsString = String.format ("%d", (int)zscoreList.get(0)[0]);
            bw.append(numberAsString);
            for(int i=1;i<zscoreList.size();i++) {
            	numberAsString = String.format ("%d", (int)zscoreList.get(i)[0]);
            	bw.append("," + numberAsString);
            }
            bw.append("\n");
            // Close connection
            bw.close();
            
            // Create zscore file
            file = new File(zscoreTxt);
            // If file doesn't exists, then create it
            if (!file.exists()) {
            	file.createNewFile();    
            	fw = new FileWriter(file.getAbsoluteFile());
            }
            else {
            	fw = new FileWriter(file.getAbsoluteFile(), true);
            }

            bw = new BufferedWriter(fw);

            // Write zscore file
            numberAsString = String.format ("%.3f", zscoreList.get(0)[1]);
            bw.append(numberAsString);
            for(int i=1;i<zscoreList.size();i++) {
            	numberAsString = String.format ("%.3f", zscoreList.get(i)[1]);
            	bw.append("," + numberAsString);
            }
            bw.append("\n");
            // Close connection
            bw.close();
            
            // Create threshold file
            file = new File(thresTxt);
            // If file doesn't exists, then create it
            if (!file.exists()) {
            	file.createNewFile();    
            	fw = new FileWriter(file.getAbsoluteFile());
            }
            else {
            	fw = new FileWriter(file.getAbsoluteFile(), true);
            }
            
            bw = new BufferedWriter(fw);

            // Write threshold file
            numberAsString = String.format ("%d", (int)zscoreList.get(0)[2]);
            bw.append(numberAsString);
            for(int i=1;i<zscoreList.size();i++) {
            	numberAsString = String.format ("%d", (int)zscoreList.get(i)[2]);
            	bw.append("," + numberAsString);
            }
            bw.append("\n");
            // Close connection
            bw.close();
        }
        catch(Exception e){
        	System.out.println(e);
        }
	}
	public void SynapticDisplay(int [][][] kSynR1, int curTimePt){
		//ImagePlus outimp = NewImage.createImage(null, width, height, 1,16, NewImage.FILL_BLACK);
		//ImageProcessor outIP =  outimp.getProcessor();
		ImageStack curStack = outputImp.getStack();
		for (int i = 0; i < width; i++) {
			for(int j = 0;j<height;j++){
				for(int k = 0;k<zSlice;k++){
					int stackNum = outputImp.getStackIndex(numChannels, k, curTimePt);
					// put channel values in an integer
					if(kSynR1[k][j][i]>0) {
						//System.out.println("Number of Synapse: "+curStack.getPixels(stackNum + 1).getClass());
						((short[]) curStack.getPixels(stackNum + 1))[i + j * width] =  (short)kSynR1[k][j][i];
					}
				}
			}
		}
		//return outimp;
	}

	public short[][] stack2array(int type, ImageStack stack, int frameNum){
		int mask=0xff;
		int nPixels=width*height;
		short[][] imArray = new short[zSlice][nPixels];
		if (type == ImagePlus.GRAY16)
		{
			mask=0xffff;
			for(int zz=1;zz<=zSlice;zz++) {
				int curSliceNum = imp.getStackIndex(numChannels, zz, frameNum);
				short[] pixels = (short[])stack.getPixels(curSliceNum);
				int intP = (int)(mask&pixels[0]);
				for (int i=0; i<nPixels; i++)
				{
					intP=(int)(mask&pixels[i]);
					short p = (short)(intP/2);
					imArray[zz][i]=p;
				}
			}
		}  
		else if (type == ImagePlus.GRAY8) 
		{
			mask=0xff;
			for(int zz=1;zz<=zSlice;zz++) {
				int curSliceNum = imp.getStackIndex(numChannels, zz, frameNum);
				byte[] pixels = (byte[])stack.getPixels(curSliceNum);
				for (int i=0; i<nPixels; i++)
				{
					short p=(short)(mask&pixels[i]);
					imArray[zz-1][i]=p;
				}
			}
		}
		else
		{
			IJ.log("Pixel format not supported");
			return null;
		}
		return imArray;
	}
	
	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ,
	 * loads an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args
	 *            unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins
		// menu
		Class<?> clazz = SynQuantVid_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		/*ImagePlus image = IJ
				.openImage("C:\\Users\\Congchao\\Desktop\\test series for process tracking-1.tif");//C2-3-weak-z_Maximum intensity projection.tif");
		image.show();

		 //run the plugin
		IJ.runPlugIn(clazz.getName(), "");*/
	}

	void error() {
		IJ.showMessage("3D particle", "Error");
	}
}
