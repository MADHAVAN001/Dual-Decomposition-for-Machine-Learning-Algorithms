package source;

import Jama.Matrix;
import java.io.*;

public class csvfileinput {
	public static void input(Matrix X, Matrix Y, String S) {
		BufferedReader br = null;
		String line = "";
		String cvsSplitBy = ",";
		String minus="-";

		int j = 0;
		try {
			
			br = new BufferedReader(new FileReader(S));
			br.readLine();
			while ((line = br.readLine()) != null) {
				j++;
				String[] tmp = line.split(cvsSplitBy);
				if (j == X.getRowDimension())
					break;
				for (int i = 1; i < X.getColumnDimension()-1; i++) {
					if(tmp[i].charAt(0)=='-')
					{
						String[] tmpif = tmp[i].split(minus);
						X.set(j - 1, i-1, -1*Double.parseDouble(tmpif[1]));
							
					}
					else
						X.set(j - 1, i-1, Double.parseDouble(tmp[i]));
				}
				
				Y.set(j - 1, 0, Double.parseDouble(tmp[280]));
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null)
				try {
					br.close();

				} catch (IOException e) {
					e.printStackTrace();
				}
		}
	}

}