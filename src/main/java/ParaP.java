import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;
import java.io.FileNotFoundException;  
import java.lang.*;
import java.util.Scanner;

import javax.swing.JOptionPane;

//parameter initialization, mainly for synapse's prior knowledge
public class ParaP {
    public String trans = "anscombe"; // not used here, for anscombe transform
    public String test = "order"; // order statistics is used here
    public String correct = "bonf"; // bonferroni correction
    public double thrSig = 0.05;  // corrected significance threshold, this is not used in fact
    public int thr0;
    public int thr1;
    public int thrg=10;  // step size
    public double thrZ = NormInv(thrSig,0,1);
    public int min_size = 4;  // !! minimun synapse size in pixels, should be image specific
    public int max_size = 100; // !! maxmum synapse size in pixels, should be image specific
    public int minIntensity = 10; 
    public int max_ratio = 2;  //
    public int conn_direc = 8;  // define neighbors
	public double minfill = 0.5; //for convex shape of synapse
	public double maxWHratio = 2.0; //remove too long or too short regions, spatial-ratio of detected synapse
	public final double [][] mu = InitialMu();//look up table for Mu
	public double [][] sigma=InitialSigma();//look up table for sigma
	public double fdr;
	public boolean fdr_den = true;
	
	public ParaP(double InputFdr, int Lowthreshold, int Highthreshold,int MinSize,int MaxSize){
		fdr = InputFdr;
		thr0 = Math.max(minIntensity,Lowthreshold);
		thr1 = Math.min(250,Highthreshold);
		min_size = MinSize;
		max_size = MaxSize;
	}
	
	/*read mu from lookup table*/
	/***read mu from lookup table***/
	public double[][] InitialMu(){
		double [][] matrix = new double[max_size][max_size];
		matrix = InitialMuSigmaInJar("mu.txt");
        return matrix;
	}
	/***read sigma from lookup table***/
	public double[][] InitialSigma(){
		double [][] matrix = new double[max_size][max_size];
		matrix = InitialMuSigmaInJar("sigma.txt");
        return matrix;                             
	}
	/***Read mu.txt and sigma.txt, which will be included in jar file***/
	public double[][] InitialMuSigmaInJar(String FileName){
		InputStream stream = paraP3D.class.getResourceAsStream(FileName);
		if (stream == null) JOptionPane.showMessageDialog(null, "Resource not located.");
		Scanner input = null;
		try {
		 input = new Scanner (stream);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			JOptionPane.showMessageDialog(null, "Scanner error");
		}
		
		//while (input.hasNextLine()) JOptionPane.showMessageDialog(null, input.nextLine());
		

		double [][] matrix = new double[max_size][max_size];
		int line = 0;
		String str = null;
		// read by line
		while(input.hasNextLine()){
			str = input.nextLine();
			String[] Linemu = str.split(",");
			for (int i = 0; i < Linemu.length; i++) {  
				matrix[line][i] = Double.parseDouble(Linemu[i]);  
				//System.out.println(matrix[line][i]);  
			}
			line++;
		}
		input.close();
        return matrix;                             
		
	}
	/***read sigma and mu from txt file***/
	public double[][] InitialMuSigma(String Filepath){
		ClassLoader classLoader = getClass().getClassLoader();
		File file = new File(classLoader.getResource(Filepath).getFile());
		double [][] matrix = new double[max_size][max_size];
        BufferedReader br = null;
        int line = 0;
		// read
		try {
			br  = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String str = null;
		// read by line
		try {
			while((str=br.readLine())!=null){
			
				String[] Linemu = str.split(",");
			    for (int i = 0; i < Linemu.length; i++) {  
			    	matrix[line][i] = Double.parseDouble(Linemu[i]);  
			        //System.out.println(matrix[line][i]);  
			    }
			    line++;
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return matrix;                             
		
	}
	/*norm inverse, only use here, in the other file use the function zTop and pToz in BasicMath*/
	public double NormInv(double p, double mu, double sigma) {
		/**
		 * Java Implementation of http://support.microsoft.com/kb/827358/en-us implementing 
		 * Apache POIs {@link FreeRefFunction} to be used as a user defined function.
		 * 
		 * Mathematics found here: 
		 * https://gist.github.com/kmpm/1211922/
		 * https://gist.github.com/kmpm/1211922/raw/a11e0dfc9fab493bcdadc669f3213d11f1897ebf/norminv.js
		 * 
		 * Register for a given workbook with:
		 * 
		    <code>
				final UDFFinder udfs = new DefaultUDFFinder(new String[]{ "MS_NormInv" }, new FreeRefFunction[]{ new NormInv() }) ;        
				workbook.addToolPack(new AggregatingUDFFinder(new UDFFinder[] {udfs}));
		    </code>
		 * @author michael.simons, 2013-02-21
		 */
		
		if(p < 0 || p > 1) 
			throw new RuntimeException("The probality p must be bigger than 0 and smaller than 1");		
		if(sigma < 0)
			throw new RuntimeException("The standard deviation sigma must be positive");
		if(p == 0)
			return Double.NEGATIVE_INFINITY;		
		if(p == 1)
			return Double.POSITIVE_INFINITY;		
		if(sigma == 0)
			return mu;		
		double  q, r, val;
 
		q = p - 0.5;
 
		/* 0.075 <= p <= 0.925 */
		if(Math.abs(q) <= .425) {
			r = .180625 - q * q;
			val =
		         q * (((((((r * 2509.0809287301226727 +
		                    33430.575583588128105) * r + 67265.770927008700853) * r +
		                  45921.953931549871457) * r + 13731.693765509461125) * r +
		                1971.5909503065514427) * r + 133.14166789178437745) * r +
		              3.387132872796366608)
		         / (((((((r * 5226.495278852854561 +
		                  28729.085735721942674) * r + 39307.89580009271061) * r +
		                21213.794301586595867) * r + 5394.1960214247511077) * r +
		              687.1870074920579083) * r + 42.313330701600911252) * r + 1);
		}
		/* closer than 0.075 from {0,1} boundary */
		else {
		     /* r = min(p, 1-p) < 0.075 */
			 if (q > 0) {
				 r = 1 - p;
			 } else {
				 r = p;
			 }
 
			 r = Math.sqrt(-Math.log(r));
			 /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
 
			 if (r <= 5) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
				 r += -1.6;
				 val = (((((((r * 7.7454501427834140764e-4 +
		                     .0227238449892691845833) * r + .24178072517745061177) *
		                   r + 1.27045825245236838258) * r +
		                  3.64784832476320460504) * r + 5.7694972214606914055) *
		                r + 4.6303378461565452959) * r +
		               1.42343711074968357734)
		              / (((((((r *
		                       1.05075007164441684324e-9 + 5.475938084995344946e-4) *
		                      r + .0151986665636164571966) * r +
		                     .14810397642748007459) * r + .68976733498510000455) *
		                   r + 1.6763848301838038494) * r +
		                  2.05319162663775882187) * r + 1);
			 } else { /* very close to  0 or 1 */
				 r += -5;
				 val = (((((((r * 2.01033439929228813265e-7 +
		                     2.71155556874348757815e-5) * r +
		                    .0012426609473880784386) * r + .026532189526576123093) *
		                  r + .29656057182850489123) * r +
		                 1.7848265399172913358) * r + 5.4637849111641143699) *
		               r + 6.6579046435011037772)
		              / (((((((r *
		                       2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
		                      r + 1.8463183175100546818e-5) * r +
		                     7.868691311456132591e-4) * r + .0148753612908506148525)
		                   * r + .13692988092273580531) * r +
		                  .59983220655588793769) * r + 1);
		      }
 
		      if (q < 0.0) {
		          val = -val;
		      }
		  }
 
		  return mu + sigma * val;		
	}
	

}