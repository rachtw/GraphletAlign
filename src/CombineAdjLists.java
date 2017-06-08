import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;
import java.util.regex.*;
//java CombineAdjLists [line num] [input: 2nd species.neighbor] [input: 3rd species.neighbor] ... [output: 1st species.neighbor]
public class CombineAdjLists {
  public static void main(String[] args){
    try {
      BufferedReader breader;
      BufferedWriter bwriter=new BufferedWriter(new FileWriter(args[args.length-1]));
      int base=Integer.parseInt(args[0]);
      Pattern p = Pattern.compile("[\\[\\];\\-\t, ]");
      Matcher m;
      StringTokenizer st;
      String str;
      for (int i=1;i<args.length-1;i++){
        breader=new BufferedReader(new FileReader(args[i]));
        int linecount=0;
        while ((str=breader.readLine())!=null) {
          st=new StringTokenizer(str,p.toString(),true);
          while (st.hasMoreTokens()) {
            m = p.matcher(str=st.nextToken());
            if (m.matches())
              bwriter.write(str);
            else
              bwriter.write(Integer.parseInt(str)+base+"");
          }
          bwriter.write("\n");
          linecount++;
        }
        breader.close();
        base+=linecount;
      }
      bwriter.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }    
}