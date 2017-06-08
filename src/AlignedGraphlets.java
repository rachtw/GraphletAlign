import java.util.ArrayList;
import java.io.Serializable;
import java.util.TreeSet;
import java.util.Iterator;
public class AlignedGraphlets implements Serializable, Comparable<AlignedGraphlets> {
  public ArrayList<TreeSet<Integer>> hitGroupNoList;
  public ArrayList<ArrayList<Integer>> alignedProteinListList;
  public int numEdges=0;
  public boolean distance2=false;
  public AlignedGraphlets(TreeSet<Integer> i, ArrayList<Integer> s) {
    hitGroupNoList=new ArrayList<TreeSet<Integer>>();
    alignedProteinListList=new ArrayList<ArrayList<Integer>>();
    hitGroupNoList.add(i);
    alignedProteinListList.add(s);
  }
  // copy constructor
  AlignedGraphlets(AlignedGraphlets dg) {
    hitGroupNoList=new ArrayList<TreeSet<Integer>>();
    alignedProteinListList=new ArrayList<ArrayList<Integer>>();
    
    Iterator<TreeSet<Integer>> i1=dg.hitGroupNoList.iterator();
    while (i1.hasNext()) {
      TreeSet<Integer> set1=i1.next();
      TreeSet<Integer> set2=new TreeSet<Integer>();
      for(Integer intObj: set1) set2.add(new Integer(intObj)); 
      hitGroupNoList.add(set2);
    }
    Iterator<ArrayList<Integer>> i2=dg.alignedProteinListList.iterator();
    while (i2.hasNext()) {
      ArrayList<Integer> list1=i2.next();
      ArrayList<Integer> list2=new ArrayList<Integer>();
      for(Integer intObj: list1) list2.add(new Integer(intObj));     
      alignedProteinListList.add(list2);
    }
    numEdges=dg.numEdges;
    distance2=dg.distance2;
  }  
  private int compareTreeSet(TreeSet<Integer> a1, TreeSet<Integer> a2) {
    Integer int1=null,int2=null;
    Iterator<Integer> i1=a1.iterator(),i2=a2.iterator();
    while (i1.hasNext() && i2.hasNext() && (int1=i1.next()).equals(int2=i2.next())) {
      ;
    }
    if (!int1.equals(int2)) {
      return int1.compareTo(int2);
    } else {
      if (!i2.hasNext()) {
        if (!i1.hasNext()) {
          return 0;
        }
        return 1;
      }
      if (!i1.hasNext()) {return -1;}
    }
    return 0;
  }
  private int compareArrayList(ArrayList<Integer> a1, ArrayList<Integer> a2) {
    Integer s1,s2;
    for (int i=0;i<a1.size();i++) {
      if ((s1=a1.get(i)) < (s2=a2.get(i))) return -1;
      if (s1 > s2) return 1;
    }
    return 0;
  }
  //compare last group nos & group contents
  public int compareTo(AlignedGraphlets dg) {
    int j;
    int last=alignedProteinListList.size()-1;
    if ((j=compareTreeSet(hitGroupNoList.get(last),dg.hitGroupNoList.get(last)))!=0)
      return j;
    if ((j=compareArrayList(alignedProteinListList.get(last),dg.alignedProteinListList.get(last)))!=0)
      return j;
    return 0;
  }
  //compare by group content
  public boolean equals(AlignedGraphlets dg) { //not used currently
//incorrect-begin
/*
    if (this.hitGroupNoList.size()!=dg.hitGroupNoList.size())
      return false;
    int i;
    if (!this.hitGroupNoList.equals(dg.hitGroupNoList)) return false;
    for (i=0;i<alignedProteinListList.size();i++)
      if (!new TreeSet(alignedProteinListList.get(i)).equals(new TreeSet(dg.alignedProteinListList.get(i))))
        return false;
    return true;
*/
//incorrect-end
    return this.hitGroupNoList.size()==dg.hitGroupNoList.size() &&
           //this.hitGroupNoList.equals(dg.hitGroupNoList) && this.alignedProteinListList.equals(dg.alignedProteinListList);
           this.alignedProteinListList.equals(dg.alignedProteinListList);
  }
  public String toString() {
    String str="";
    Iterator<Integer> it;
    for (int i=0;i<alignedProteinListList.size();i++) {
      it=hitGroupNoList.get(i).iterator();
      while (it.hasNext())
        str+=it.next()+";";
      str=str.substring(0,str.length()-1);
      str+=alignedProteinListList.get(i).toString()+"-";
    }
    return str;
  }
  public void clear() {
    hitGroupNoList.clear();
    alignedProteinListList.clear();
  }
  public String prefix() {
    String str="";
    Iterator<Integer> it;
    int limit=(int)Math.ceil((double)alignedProteinListList.size()/2);
    for (int i=0;i<limit;i++) {
      it=hitGroupNoList.get(i).iterator();
      while (it.hasNext())
        str+=it.next()+";";
      str=str.substring(0,str.length()-1);
      str+=alignedProteinListList.get(i).toString()+"-";
    }
    return str;
  }
  public String proteinSet(int index) {
    String str="";
    for (int i=0;i<alignedProteinListList.size();i++) {
      str+=alignedProteinListList.get(i).get(index)+" ";
    }
    str=str.substring(0,str.length()-1);
    return str;
  }
  public void proteinSet(TreeSet<Integer> s,int index) {
    for (int i=0;i<alignedProteinListList.size();i++) {
      s.add(alignedProteinListList.get(i).get(index));
    }
  }
  public void restProteins(TreeSet<Integer> s,int index) {
    int limit=(int)Math.ceil((double)alignedProteinListList.size()/2);
    for (int i=limit;i<alignedProteinListList.size();i++) {
      s.add(alignedProteinListList.get(i).get(index));
    }
  }
}