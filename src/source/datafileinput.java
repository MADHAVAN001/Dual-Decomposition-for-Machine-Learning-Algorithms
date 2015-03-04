package source;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import Jama.Matrix;

public class datafileinput{
	public static void input(Matrix X, Matrix Y, String Path)
	{
		
			BufferedReader br = null;
			String line = "";
			String cvsSplitBy = " ";
			String colon=":";

			int j = 0;
			try {
				for(int file = 0;file<=120;file++)
				{
				br = new BufferedReader(new FileReader(Path+file));
				while ((line = br.readLine()) != null) {
					j++;
					if(X.getRowDimension()==j)
						return;
					String[] tmp = line.split(cvsSplitBy);
					Y.set(j - 1, 0, Double.parseDouble(tmp[0]));
					for(int i = 1;i<tmp.length;i++)
					{
						String[] tmpstr = tmp[i].split(colon);
						X.set(j-1, Integer.parseInt(tmpstr[0]), Double.parseDouble(tmpstr[1]));
					}
					
				}}
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