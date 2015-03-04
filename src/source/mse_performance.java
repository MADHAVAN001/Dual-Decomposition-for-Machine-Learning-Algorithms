package source;
import Jama.Matrix;

public class mse_performance{
	public static Matrix mse_problem(Matrix beta_cd, Matrix X, Matrix Y, int num_cycle)
	{
	Matrix performance = new Matrix(num_cycle,1);
	for(int i = 0;i<num_cycle;i++)
	{
		Matrix beta_cd_i =  new Matrix(beta_cd.getRowDimension(),1);
		for(int j = 0;j<beta_cd.getRowDimension();j++)
		{
			beta_cd_i.set(j, 0, beta_cd.get(j, i));
		}
		Matrix perform = ((X.times(beta_cd_i).minus(Y)).transpose())
				.times(X.times(beta_cd_i).minus(Y));
		performance.set(i,0,perform.get(0, 0));
	}
	return performance;
		
}
}