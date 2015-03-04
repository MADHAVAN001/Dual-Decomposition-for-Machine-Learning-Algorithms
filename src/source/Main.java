package source;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import Jama.Matrix;

public class Main {
	public static void main(String[] args) {
		int m =6000;
		int n = 279;
		int cycle = 10;
		final double mu = 10;
		final double sigma = 25;
		String opfile = "dd-Re-2-m.txt";
		Matrix m_op = new Matrix(300,1);
		Matrix time = new Matrix(300,1);
		for(m = 500;m<=6000;m+=500)
		{
		Matrix X = new Matrix(m, n);
		Matrix Y = new Matrix(m, 1);
		/*for (int i = 0; i < X.getRowDimension(); i++)
			for (int j = 0; j < X.getColumnDimension(); j++) {
				double random = Math.random();
				double s = random - (int) random;
				X.set(i, j, mu + sigma * s);
			}
		for (int i = 0; i < Y.getRowDimension(); i++)
			for (int j = 0; j < Y.getColumnDimension(); j++) {
				double random = Math.random();
				double s = random - (int) random;
				Y.set(i, j, mu + sigma * s);
			}
			*/
		//Matrix beta_accurate = ((((X.transpose()).times(X)).inverse()).times(X
			//	.transpose())).times(Y);
		String Path = "blogData_train.csv";
		csvfileinput.input(X, Y, Path);
		for (int y = 0; y < X.getRowDimension(); y++)
			for (int t = 0; t < X.getColumnDimension(); t++) {
				double rand = Math.random();
				X.set(y, t, X.get(y, t) + rand * 0.01);
			}
		String Path1 = "/home/hduser/Desktop/Data/";
		//datafileinput.input(X, Y, Path1);
		long startTime = System.currentTimeMillis();
		Matrix beta = dualdecomposition(X, Y, cycle);
		long stopTime = System.currentTimeMillis();
	    long elapsedTime = stopTime - startTime;
		//Matrix abs_performance = (((X.times(beta_accurate).minus(Y))
			//	.transpose()).times(X.times(beta_accurate).minus(Y)));
		Matrix abs_performance = new Matrix(cycle,1);
		 //Matrix beta_accurate = new Matrix(m,1);
		 for(int tmp = 0;tmp<abs_performance.getRowDimension();tmp++)
		 {
			 abs_performance.set(tmp, 0, 0);
		 }
		Matrix performance = mse_performance.mse_problem(beta, X, Y, cycle);
		System.out.println(abs_performance.get(0, 0));
		
		m_op.set((m-500)/500, 0, performance.get(9, 0));
		time.set((m-500)/500, 0, elapsedTime);
		}
		try {
			FileWriter fw = new FileWriter(opfile);
			for (int j = 0; j < 300; j++) {

				double diff = m_op.get(j, 0)
						/ m;
				fw.append(String.valueOf(diff));
				fw.append(',');

				fw.append(String.valueOf(time.get(j, 0)));
				fw.append('\n');
			}
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	public static Matrix dualdecomposition(Matrix X, Matrix Y, int cycle) {
		
		Matrix performance = new Matrix(cycle,1);
		Matrix beta = new Matrix(X.getColumnDimension(), cycle);
		Matrix gamma = new Matrix(X.getRowDimension(), 1);
		Matrix eps = new Matrix(X.getRowDimension(), 1);
		Matrix beta_i = new Matrix(X.getColumnDimension(),1);
		
		for (int i = 0; i < beta.getRowDimension(); i++)
			for (int j = 0; j < beta.getColumnDimension(); j++)
				beta.set(i, j, 0);
		for (int i = 0; i < gamma.getRowDimension(); i++)
			for (int j = 0; j < gamma.getColumnDimension(); j++)
				gamma.set(i, j, 1);
		
		for(int k = 0;k<eps.getRowDimension();k++)
		{
			eps.set(k, 0, -Y.get(k, 0));
		}
		double max1 = 1;
		double max2 = 1;
		
		for (int i = 1; i < cycle; i++) {
			double alpha = 0.1/max1;
			double alpha1 = 0.1/max2;
			
			max1 = Math.abs(eps.get(0, 0));
			for(int p = 0;p<eps.getRowDimension();p++)
			{
				if(Math.abs(eps.get(p, 0))>max1)
					max2 = Math.abs(eps.get(p, 0));
			}
			System.out.println("Cycle: "+ i + " Max1: "+max1+"Max2: "+max2);
			for (int j = 0; j < eps.getRowDimension(); j++) {
				//eps.set(j,0,eps.get(j, 0)- (alpha1 * (eps.get(j, 0) - gamma.get(j, 0))));

				eps.set(j,0,gamma.get(j,0));
			}
			double sum_eps = 0;
			for(int k = 0;k<eps.getRowDimension();k++)
				sum_eps+=(eps.get(k, 0)*eps.get(k, 0));
			Matrix sub_perf = (gamma.transpose()).times((X.times(beta_i).minus(Y)).minus(eps));
			double sum_i = 0;
			Matrix perf = (X.times(beta_i)).minus(eps).minus(Y);
			for(int p = 0;p<X.getRowDimension();p++)
				sum_i += Math.abs(perf.get(p, 0)); 
			
			System.out.println("sum_eps: "+sum_eps+"   sub_perf: "+ sub_perf.get(0, 0)+ "  perf"+sum_i);
			performance.set(i, 0, 0.5*sum_eps+sub_perf.get(0, 0));
			for (int k = 0; k < X.getColumnDimension(); k++) {
				beta.set(k,i,beta.get(k, i - 1)- Math.pow(10,-11.5) * ((((gamma.transpose()).times(X)).transpose()).get(k, 0)));
			}
			for(int k = 0;k<X.getColumnDimension();k++)
				beta_i.set(k, 0,beta.get(k, i));
			for (int j = 0; j < X.getRowDimension(); j++) {
				gamma.set(j,0,gamma.get(j, 0)+ alpha * ((((X.times(beta_i)).minus(Y)).minus(eps)).get(j,0)));
			}
			Matrix tmp1 = ((X.times(beta_i)).minus(eps)).minus(Y);
			max1 = Math.abs(tmp1.get(0, 0));
			for(int p = 0;p<tmp1.getRowDimension();p++)
			{
				if(Math.abs(tmp1.get(p, 0))>max1)
					max1 = Math.abs(tmp1.get(p, 0));
			}
		}
		return beta;
	}
}