import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.Collections;
//java ParseNetwork [*.net] [output neighbor file] [output name file] [threshold]
public class ParseNetwork {
  public static void main(String[] args){
    try {     
      LinkedHashMap<String, TreeSet<String>> neighbors=new LinkedHashMap<String, TreeSet<String>>();
      HashMap<String,Integer> serialNo=new HashMap<String,Integer>();
      HashMap<String,String> geneName=new HashMap<String,String>();
      StringTokenizer st;
      String str,str1,str2;
      TreeSet<String> neighborSet;
      double threshold = 0;
      if (args.length==4) threshold=Double.parseDouble(args[3]);
      int i=1;
      BufferedReader breader=new BufferedReader(new FileReader(args[0]));
			/*
      if (args[0].indexOf("1.0")==-1) {
        breader.readLine();
        breader.readLine();
        //System.err.println("2.0");
      }
			*/
      while ((str=breader.readLine())!=null) {
        st=new StringTokenizer(str,"\t");
        str1=st.nextToken();
        str2=st.nextToken();
        if (str2.equals("pp")) str2=st.nextToken();
        if (args.length<4 || Double.parseDouble(st.nextToken()) >= threshold) {
          if ((neighborSet=neighbors.get(str1))==null) {
            neighbors.put(str1,neighborSet=new TreeSet<String>());
            serialNo.put(str1,new Integer(i));
						//if (!str1.equals(i+""))
						  //System.err.println("assign "+str1+" to "+i);
            i++;
          }
          neighborSet.add(str2);
          if ((neighborSet=neighbors.get(str2))==null) {
            neighbors.put(str2,neighborSet=new TreeSet<String>());
            serialNo.put(str2,new Integer(i));
						//if (!str2.equals(i+""))
							//System.err.println("assign "+str2+" to "+i);
            i++;
          }
          neighborSet.add(str1);
        }
      }
      breader.close();
      
      BufferedWriter bwriter=new BufferedWriter(new FileWriter(args[1]));
      BufferedWriter bwriter2=new BufferedWriter(new FileWriter(args[2]));
      ArrayList<Integer> serialNoList=new ArrayList<Integer>();
      Iterator<String> it=neighbors.keySet().iterator();
      Iterator<String> it2;
      while (it.hasNext()) {
        bwriter.write(serialNo.get((str=it.next()))+"[");
        bwriter2.write(str);
        bwriter2.newLine();
        
        serialNoList.clear();
        it2=neighbors.get(str).iterator();
        while (it2.hasNext())
          serialNoList.add(serialNo.get(it2.next()));
        Collections.sort(serialNoList);
        for (i=0;i<serialNoList.size();i++)
          bwriter.write(serialNoList.get(i)+";");
        bwriter.write("]");
        bwriter.newLine();
      }
      bwriter.close();
      bwriter2.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }    
}