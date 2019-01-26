import java.awt.Color;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
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

public class SynQuant_ implements PlugIn {
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp; //Current Image
	protected int ntry=1; //decide max number of neighbors considered in zscore calculation
	protected int width; //image width
	protected int height; //image height
	protected int presyn_chl;
	protected int postsyn_chl;
	protected int den_chl;
	protected double fdr;
	protected int MinSize,MaxSize;
	protected int NumSynSite=0;
	
	boolean fastflag = false;

	public void run(String arg) {
		try {
			if (showDialog())
				synQuant_real();
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}
	public boolean showDialog() 
	{
		// Get input parameter
		GenericDialog gd = new GenericDialog("SynQuant - Data and Parameter Setting");
		gd.addNumericField("Channel for Pre-synaptic Puncta (if NA, input 0)#: ", 1, 0);//2.5-3.5 good
		gd.addNumericField("Channel for Post-synaptic Puncta (if NA, input 0)#: ", 2, 0);//2.5-3.5 good
		gd.addNumericField("Channel for Dendrite (if NA, input 0)#: ", 3, 0);//2.5-3.5 good
		gd.addNumericField("FDR Control for Synapse Detection: ", 0.05, 3);//2.5-3.5 good
		gd.addNumericField("Min Synapse Size: ", 4, 0);//2.5-3.5 good
		gd.addNumericField("Max Synapse Size: : ", 100, 0);//2.5-3.5 good
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		presyn_chl = (int) gd.getNextNumber();
		postsyn_chl = (int) gd.getNextNumber();
		den_chl = (int) gd.getNextNumber();
		fdr = gd.getNextNumber();
		MinSize = (int) gd.getNextNumber();
		MaxSize = (int) gd.getNextNumber();
		if(Double.isNaN(fdr) || (fdr<=0) || (fdr>=1)){
			IJ.showMessage("Invalid parameter(s).\n" + "0-1 needed.");
			return false;
		}
		return true;
	}
	public void synQuant_real() {
		imp = WindowManager.getCurrentImage();
		stack = imp.getStack();
		width = imp.getWidth();
		height = imp.getHeight();
		
		double vox_x = imp.getCalibration().pixelWidth; 
		if(vox_x==1)//simulated data
			vox_x = 2.0757e-7;
		else
			vox_x = vox_x*1e-6;//real data
		int type = imp.getType();
		
		boolean displayROI = true;
		ppsd postsyn_det = null;
		ppsd presyn_det = null;
		GrowNeurite den_det = null;
		String Roititle = ""; 
		//postsynapse detection
		if(postsyn_chl>0 & presyn_chl==0 & den_chl==0){
			short[] postsynArr = stack2array(type, stack, postsyn_chl);
			Roititle = "Post-Synaptic Puncta";
			postsyn_det = new ppsd(postsynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
		}
		//presynapse detection
		if(postsyn_chl==0 & presyn_chl>0 & den_chl==0){
			short[] presynArr = stack2array(type, stack, presyn_chl);
			Roititle = "Pre-Synaptic Puncta";
			presyn_det = new ppsd(presynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
		}
		//dendrite detection
		if(postsyn_chl==0 & presyn_chl==0 & den_chl>0){
			short[] denArr = stack2array(type, stack, den_chl);
			Roititle = "Dendrite";
			den_det = new GrowNeurite(denArr, width, height,vox_x,displayROI,Roititle);
			ResultsTable den_table = new ResultsTable();
			den_table.incrementCounter();
			den_table.addLabel("Total Length of Dendrite:");
			den_table.addValue("Value", den_det.totallength);
			den_table.showRowNumbers(false);
			den_table.show("Dendrite Detection Summary Table");
		}
		//synaptic site detection
		if(postsyn_chl>0 & presyn_chl>0 & den_chl==0){
			displayROI = false;
			short[] postsynArr = stack2array(type, stack, postsyn_chl);
			Roititle = "Post-Synaptic Puncta";
			postsyn_det = new ppsd(postsynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
			
			short[] presynArr = stack2array(type, stack, presyn_chl);
			Roititle = "Pre-Synaptic Puncta";
			presyn_det = new ppsd(presynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
			Roititle = "Synaptic Sites";
			SynapticDisplay(presyn_det, postsyn_det,Roititle);
			
			ResultsTable synapseSite_table = new ResultsTable();
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Pre-synaptic Puncta number: ");
			synapseSite_table.addValue("Value", presyn_det.nSyn0);
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Post-synaptic Puncta number: ");
			synapseSite_table.addValue("Value", postsyn_det.nSyn0);
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Synaptic site number: ");
			synapseSite_table.addValue("Value", NumSynSite);
			synapseSite_table.showRowNumbers(false);
			synapseSite_table.show("Synaptic sites Detection Summary Table");
		}
		//synapse quantification with post-synaptic puncta
		if(postsyn_chl>0 & presyn_chl==0 & den_chl>0){
			//synapse detection
			displayROI = true;
			short[] postsynArr = stack2array(type, stack, postsyn_chl);
			Roititle = "Post-Synaptic Puncta";
			postsyn_det = new ppsd(postsynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
			//dendrite extraction
			displayROI = false;
			Roititle = "Dendrite";
			short[] denArr = stack2array(type, stack, den_chl);
			den_det = new GrowNeurite(denArr, width, height,vox_x,displayROI,Roititle);
			
			LinearTest LT = new LinearTest(postsyn_det.kSynR1,den_det);
			//show features
			ResultsTable Ft_table = new ResultsTable();
			int DenCnt = 0;
			for (int i=0;i<den_det.CurNum;i++) {
				Ft_table.incrementCounter();
				DenCnt = DenCnt+1;
				Ft_table.addLabel("Dendrite piece #"+DenCnt);
				Ft_table.addValue("Number of synapse on it",LT.denSynNum[i]);
				Ft_table.addValue("Dendrite length",LT.DenFeatures[i][0]);
				Ft_table.addValue("Dendrite scale",LT.DenFeatures[i][1]);
				Ft_table.addValue("Dendrite intensity",LT.DenFeatures[i][2]);
			}
			Ft_table.showRowNumbers(false);
			Ft_table.show("Synapse Quantification Feature Table");
			
			ResultsTable esm_table = new ResultsTable();

			esm_table.incrementCounter();
			esm_table.addLabel("Coefficents:");
			esm_table.addValue("Length",LT.betas[0]);
			esm_table.addValue("Scale",LT.betas[1]);
			esm_table.addValue("Intensity",LT.betas[2]);
			esm_table.addValue("Intercept",LT.betas[3]);
			esm_table.showRowNumbers(false);
			esm_table.show("Synapse Quantification Coefficents Table");
			
			ResultsTable sm_table = new ResultsTable();

			sm_table.incrementCounter();
			sm_table.addLabel("Total # of synapses:");sm_table.addValue("Value",postsyn_det.nSyn0);
			sm_table.incrementCounter();
			sm_table.addLabel("Total Length of Dendrite:");sm_table.addValue("Value",den_det.totallength);
			sm_table.incrementCounter();
			sm_table.addLabel("Puncta Density per Unit Length:");sm_table.addValue("Value",(double)postsyn_det.nSyn0/den_det.totallength);
			sm_table.showRowNumbers(false);
			sm_table.show("Synapse Quantification Summary Table");
		}
		
		//synapse quantification with pre-synaptic puncta
		if (postsyn_chl == 0 & presyn_chl > 0 & den_chl > 0) {
			// synapse detection
			displayROI = true;
			short[] presynArr = stack2array(type, stack, presyn_chl);
			Roititle = "Pre-Synaptic Puncta";
			presyn_det = new ppsd(presynArr, width, height, vox_x, fdr, fastflag, MinSize, MaxSize, displayROI,
					Roititle);
			// dendrite extraction
			displayROI = false;
			Roititle = "Dendrite";
			short[] denArr = stack2array(type, stack, den_chl);
			den_det = new GrowNeurite(denArr, width, height, vox_x, displayROI, Roititle);

			LinearTest LT = new LinearTest(presyn_det.kSynR1, den_det);
			// show features
			ResultsTable Ft_table = new ResultsTable();
			int DenCnt = 0;
			for (int i = 0; i < den_det.CurNum; i++) {
				Ft_table.incrementCounter();
				DenCnt = DenCnt + 1;
				Ft_table.addLabel("Dendrite piece #" + DenCnt);
				Ft_table.addValue("Number of synapse on it", LT.denSynNum[i]);
				Ft_table.addValue("Dendrite length", LT.DenFeatures[i][0]);
				Ft_table.addValue("Dendrite scale", LT.DenFeatures[i][1]);
				Ft_table.addValue("Dendrite intensity", LT.DenFeatures[i][2]);
			}
			Ft_table.showRowNumbers(false);
			Ft_table.show("Synapse Quantification Feature Table");

			ResultsTable esm_table = new ResultsTable();

			esm_table.incrementCounter();
			esm_table.addLabel("Coefficents:");
			esm_table.addValue("Length", LT.betas[0]);
			esm_table.addValue("Scale", LT.betas[1]);
			esm_table.addValue("Intensity", LT.betas[2]);
			esm_table.addValue("Intercept", LT.betas[3]);
			esm_table.showRowNumbers(false);
			esm_table.show("Synapse Quantification Coefficents Table");

			ResultsTable sm_table = new ResultsTable();

			sm_table.incrementCounter();
			sm_table.addLabel("Total # of synapses:");
			sm_table.addValue("Value", presyn_det.nSyn0);
			sm_table.incrementCounter();
			sm_table.addLabel("Total Length of Dendrite:");
			sm_table.addValue("Value", den_det.totallength);
			sm_table.incrementCounter();
			sm_table.addLabel("Puncta Density per Unit Length:");
			sm_table.addValue("Value", (double) presyn_det.nSyn0 / den_det.totallength);
			sm_table.showRowNumbers(false);
			sm_table.show("Synapse Quantification Summary Table");
		}
		//synaptic site quantification
		if(postsyn_chl>0 & presyn_chl>0 & den_chl>0){
			displayROI = false;
			short[] postsynArr = stack2array(type, stack, postsyn_chl);
			Roititle = "Post-Synaptic Puncta";
			postsyn_det = new ppsd(postsynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
			
			short[] presynArr = stack2array(type, stack, presyn_chl);
			Roititle = "Pre-Synaptic Puncta";
			presyn_det = new ppsd(presynArr, width, height,vox_x, fdr,fastflag,MinSize,MaxSize,displayROI,Roititle);
			Roititle = "Synaptic Sites";
			boolean[][] SynapticSite = SynapticDisplay(presyn_det, postsyn_det,Roititle);
			
			short[] denArr = stack2array(type, stack, den_chl);
			den_det = new GrowNeurite(denArr, width, height,vox_x,displayROI,Roititle);
			LinearTest LT = new LinearTest(SynapticSite,den_det);
			
			
			//show features
			ResultsTable Ft_table = new ResultsTable();
			int DenCnt = 0;
			for (int i=0;i<den_det.CurNum;i++) {
				Ft_table.incrementCounter();
				DenCnt = DenCnt+1;
				Ft_table.addLabel("Dendrite piece #"+DenCnt);
				Ft_table.addValue("Number of synapse on it",LT.denSynNum[i]);
				Ft_table.addValue("Dendrite length",LT.DenFeatures[i][0]);
				Ft_table.addValue("Dendrite scale",LT.DenFeatures[i][1]);
				Ft_table.addValue("Dendrite intensity",LT.DenFeatures[i][2]);
			}
			Ft_table.showRowNumbers(false);
			Ft_table.show("Synapse Quantification Feature Table");
			
			ResultsTable esm_table = new ResultsTable();

			esm_table.incrementCounter();
			esm_table.addLabel("Coefficents:");
			esm_table.addValue("Length",LT.betas[0]);
			esm_table.addValue("Scale",LT.betas[1]);
			esm_table.addValue("Intensity",LT.betas[2]);
			esm_table.addValue("Intercept",LT.betas[3]);
			esm_table.showRowNumbers(false);
			esm_table.show("Synapse Quantification Coefficents Table");			
			
			ImageHandling IH = new ImageHandling();
			IH.bwlabel(SynapticSite, 8);
			ResultsTable ssq_table = new ResultsTable();
			ssq_table.incrementCounter();
			ssq_table.addLabel("Total # of synaptic sites:");ssq_table.addValue("Value",IH.NextLabel);
			ssq_table.incrementCounter();
			ssq_table.addLabel("Total Length of Dendrite:");ssq_table.addValue("Value",den_det.totallength);
			ssq_table.incrementCounter();
			ssq_table.addLabel("Puncta Density per Unit Length:");ssq_table.addValue("Value",(double)IH.NextLabel/den_det.totallength);

			ssq_table.showRowNumbers(false);
			ssq_table.show("Synaptic Site Quantification Summary Table");
			
			ResultsTable synapseSite_table = new ResultsTable();
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Pre-synaptic Puncta number: ");
			synapseSite_table.addValue("Value", presyn_det.nSyn0);
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Post-synaptic Puncta number: ");
			synapseSite_table.addValue("Value", postsyn_det.nSyn0);
			synapseSite_table.incrementCounter();
			synapseSite_table.addLabel("Synaptic site number: ");
			synapseSite_table.addValue("Value", NumSynSite);
			synapseSite_table.showRowNumbers(false);
			synapseSite_table.show("Synaptic sites Detection Summary Table");
		}
	}
	public short[] stack2array(int type, ImageStack stack, int channel){
		int mask=0xff;
		int nPixels=width*height;
		short[] imArray = new short[nPixels];
		if (type == ImagePlus.GRAY16)
		{
			mask=0xffff;
			//Find min, max. Copy to imArray
			short[] pixels = (short[])stack.getPixels(channel);
			int intP = (int)(mask&pixels[0]);
			for (int i=0; i<nPixels; i++)
			{
				intP=(int)(mask&pixels[i]);
				short p = (short)(intP/2);
				imArray[i]=p;
			}
		}  
		else if (type == ImagePlus.GRAY8) 
		{
			mask=0xff;
			//Find min, max. Copy to imArray
			byte[] pixels = (byte[])stack.getPixels(channel);
			for (int i=0; i<nPixels; i++)
			{
				short p=(short)(mask&pixels[i]);
				imArray[i]=p;
			}
		}
		else
		{
			IJ.log("Pixel format not supported");
			return null;
		}
		return imArray;
	}
	public boolean[][] SynapticDisplay(ppsd presyn, ppsd postsyn,String RoiTitile){
		ImagePlus outimp = NewImage.createRGBImage("Synaptic Sites, G:PreSynaptic Puncta, B:PostSynaptic Puncta", width, height, 1,NewImage.FILL_BLACK);
		ImageProcessor outIP =  outimp.getProcessor();
		boolean[][] overlappedSyn = new boolean[height][width]; 
		for (int i = 0; i < width; i++)
			for(int j = 0;j<height;j++)
				if(presyn.kSynR1[j][i] & postsyn.kSynR1[j][i])
					overlappedSyn[j][i] = true;
		for (int i = 0; i < width; i++) {
			for(int j = 0;j<height;j++){
			// put channel values in an integer
			if(overlappedSyn[j][i])
				outIP.set(i, j, (((int)200 & 0xff) << 16)
				          + (((int)presyn.G[j][i] & 0xff) << 8)
				          +  ((int)postsyn.G[j][i] & 0xff));
			else
				outIP.set(i, j, (((int)0 & 0xff) << 16) 
				          + (((int)presyn.G[j][i] & 0xff) << 8)
				          +  ((int)postsyn.G[j][i] & 0xff));

			}
		}
		outimp.show();
		outimp.updateAndDraw();
		ImageHandling IH = new ImageHandling();
		int[][] ovSynLabel = IH.bwlabel(overlappedSyn, 8);
		NumSynSite = IH.NextLabel;
		IH.DisplayROI(IH.NextLabel,height,width,ovSynLabel, outimp,RoiTitile);
		return overlappedSyn;
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
		Class<?> clazz = SynQuant_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		//ImagePlus image = IJ
		//		.openImage("C:\\Users\\ccwang\\Desktop\\Microglia\\3-weak-z_Maximum intensity projection_Crop.tif");//C2-3-weak-z_Maximum intensity projection.tif");
		//image.show();

		// run the plugin
		//IJ.runPlugIn(clazz.getName(), "");
	}

	void error() {
		IJ.showMessage("PPSD", "Error");
	}
}
