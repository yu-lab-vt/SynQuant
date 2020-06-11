import java.lang.*;

public class ExtendPreChannelResult {
	public paraQ3D q;
	public int ValidCnt=0;
	public double[][][][] ExtendedResult;
	public int [][] KernelIdx;
	
	
    public ExtendPreChannelResult(paraQ3D Inputq){
    	q=Inputq;
    	GetExtendedResult();
    	
    	//TestGetExtendedResult();
    };
    
    public void GetKernel() {
		double ExtendedDistance=q.ExtendedDistance;
		double zAxisMultiplier=q.zAxisMultiplier;
		
		int max_xy=(int) Math.ceil(ExtendedDistance);
		int max_z=(int) Math.ceil(ExtendedDistance/zAxisMultiplier);
		KernelIdx=new int [(2*max_xy+1)*(2*max_xy+1)*(2*max_z+1)][3];
		
		double zAxisMultiplierP2=zAxisMultiplier*zAxisMultiplier;
		double DistanceP2=ExtendedDistance*ExtendedDistance;
		
		
		ValidCnt=0;
		for (int i = -max_xy; i <= max_xy; i++) {
			for (int j = -max_xy; j <= max_xy; j++) {
				for (int k = -max_z; k <= max_z; k++) {
					if(i*i+j*j+k*k*zAxisMultiplierP2<=DistanceP2) {
						KernelIdx[ValidCnt][0]=k;
						KernelIdx[ValidCnt][1]=i;
						KernelIdx[ValidCnt][2]=j;
						ValidCnt=ValidCnt+1;
					}
					
				}
			}
		}
    }
    
    public void GetExtendedResult() {
    	
    	GetKernel();

    	ExtendedResult=new double[q.synZscore.length][q.synZscore[0].length][q.synZscore[0][0].length][q.synZscore[0][0][0].length];
		for(int i=0; i<q.synZscore[0][0].length; i++) {
			for(int j=0; j<q.synZscore[0][0][0].length; j++) {
				for (int k=0; k<q.synZscore[0].length; k++) {
					for(int timeCnt=0;timeCnt<q.synZscore.length;timeCnt++) {
						ExtendedResult[timeCnt][k][i][j]=0;
					}
				}
			}
		}
    	
    	int Ci;
    	int Cj;
    	int Ck;
		for(int i=0; i<ExtendedResult[0][0].length; i++) {
			for(int j=0; j<ExtendedResult[0][0][0].length; j++) {
				for (int k=0; k<ExtendedResult[0].length; k++) {
					if (q.synZscore[0][k][i][j]>0) {
						for (int IdxCnt=0; IdxCnt<ValidCnt;IdxCnt++) {
							Ci=i+KernelIdx[IdxCnt][1];
							if(Ci<0) {continue;}
							if(Ci>=ExtendedResult[0][0].length) {continue;}
							Cj=j+KernelIdx[IdxCnt][2];
							if(Cj<0) {continue;}
							if(Cj>=ExtendedResult[0][0][0].length) {continue;}
							Ck=k+KernelIdx[IdxCnt][0];
							if(Ck<0) {continue;}
							if(Ck>=ExtendedResult[0].length) {continue;}
							
							for(int timeCnt=0;timeCnt<q.synZscore.length;timeCnt++) {
								ExtendedResult[timeCnt][Ck][Ci][Cj]=q.synZscore[timeCnt][k][i][j];
							}
						}
					}
				}
			}
		}
    	
    }

}
