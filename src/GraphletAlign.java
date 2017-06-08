import java.util.StringTokenizer;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.ListIterator;
import java.io.*;
import java.math.BigInteger;
import java.util.LinkedList;
import java.util.Stack;
import java.util.Collections;
import java.lang.Math;
import java.util.TreeSet;
import java.util.TreeMap;
import java.util.Set;
import java.text.DecimalFormat;
import java.math.BigDecimal;
//java GraphletAlign
//[*.hits] [*.adjl] [*.blastp.no] [*.adj] [*.auto] [# disjoint vertices] [# accumulated proteins from each species] [output folder]

public class GraphletAlign {
  final static boolean DEBUG=false;
  final static boolean INDUCED_SUBGRAPH=true;
  final static boolean REPEAT_GROUP=true;
  final static boolean OUTPUT_PROTEIN_SETS=false;
  final static boolean ABBREVIATION=false;
  final static boolean DISTANCE2=true;
  final static boolean DISTANCE2_COUNT=false; // not used
  final static int DISTANCE2_NUM_THRESHOLD=2; // not used
  final static boolean BLAST_SCORE=true;
  final static boolean EXPECTATION_OR_BIT_SCORE=true;
  final static boolean OUTPUT_ALIGNMENT=true;
  //final static boolean EXAMINE_SUBGRAPH=false; // given up
  //final static boolean EXACT_SUBSET=false; // given up
  //final static boolean ENABLE_GROUP_THRESHOLD=false; // given up
  //final static int GROUP_THRESHOLD=10; // given up
  
  //have not testd on DISTANCE2=false;
  //need to debug when more subgraph copies than species  

  private String[] args;
  private boolean SINGLE_SPECIES=true;
  private int NUM_DISJOINT_VERTICES=3;
  private int NUM_EDGES;
  private int NUM_VERTICES=0;
  private int SPECIES_GRAPH_SIZE; // from argument
  private ArrayList<Integer> SPECIES_NUM_LINES=new ArrayList<Integer>(); // from arguments
    //private String[] species;
  private int num_species=1;
  private int numSubgraphsInOneSpecies;
  private TreeMap<Integer,Integer> subgraphNo_to_speciesNo;
  private TreeMap<Integer,Integer> speciesNo_to_subgraphNo;
  private ArrayList<ArrayList<Integer>> hitGroupList; // hit protein group list
  private ArrayList<ArrayList<Integer>> proteinToHitGroupNosList;
  private boolean[][] sameGroup; 
  private ArrayList<ArrayList<Integer>> proteinToNeighborArray; // neighbor
  private ArrayList<ArrayList<Integer>> proteinToNeighborArray2; //@ DISTANCE2
  private ArrayList<TreeMap<Integer,BigDecimal>> blastScores; //@ BLAST_SCORE

/*
  // for reducing graph amount
  private int order[];
  private ArrayList<ArrayList<AlignedGraphlets>> arrayArrayDG;
  private int[][][] adjacencyMatrix=new int[60][][];
  private ArrayList<TreeMap<Integer, TreeSet<Integer>>> dataIndex;
  private ArrayList<TreeMap<Integer, TreeSet<Integer>>> vertex_vdgNo;
  private int fileNo;
*/
  class ProteinListComparator implements java.util.Comparator<ArrayList<Integer>> {
    public int compare(ArrayList<Integer> a1, ArrayList<Integer> a2) {
      Integer s1,int2;
      int i;
      for (i=0;i<a1.size();i++) {
        if (i>=a2.size()) return 1;
        if ((s1=a1.get(i)) < (int2=a2.get(i))) return -1;
        if (s1 > int2) return 1;
      }
      if (i<a2.size()) return -1;
      return 0;
    }
    public int compareCombination(ArrayList<Integer> a1, ArrayList<Integer> a2) {
      ArrayList<Integer> a1_=new ArrayList<Integer>(a1);
      Collections.sort(a1_);
      ArrayList<Integer> a2_=new ArrayList<Integer>(a2);
      Collections.sort(a2_);
      return compare(a1_,a2_);
    }
  }

  // compare the last added group, considering only the group content
  class AlignedGraphletsComparator implements java.util.Comparator<AlignedGraphlets> {
    TreeSet<ArrayList<Integer>> otherLevelRedundancy;
    public AlignedGraphletsComparator(TreeSet<ArrayList<Integer>> otherLevelRedundancy) { this.otherLevelRedundancy=otherLevelRedundancy; }
    public int compare(AlignedGraphlets dg1, AlignedGraphlets dg2) {
      int j;
      int last=dg1.alignedProteinListList.size()-1;
      if (ABBREVIATION) {
        if ((j=(new ProteinListComparator()).compareCombination(dg1.alignedProteinListList.get(last),dg2.alignedProteinListList.get(last)))!=0) return j;
        //if (DEBUG) System.err.println(dg1.alignedProteinListList.get(last)+" equals "+dg2.alignedProteinListList.get(last));
        ArrayList<Integer> lastList=new ArrayList<Integer>(dg1.alignedProteinListList.get(last));
        Collections.sort(lastList);
        otherLevelRedundancy.add(lastList); //otherLevelRedundancy stores sorted list;
      } else {
        if ((j=(new ProteinListComparator()).compare(dg1.alignedProteinListList.get(last),dg2.alignedProteinListList.get(last)))!=0) return j;
        ArrayList<Integer> lastList=new ArrayList<Integer>(dg1.alignedProteinListList.get(last));
        otherLevelRedundancy.add(lastList);
      }
      return 0;
    }
  }

  private void removeRedundancy(LinkedList<AlignedGraphlets> queue) {
    //if (DEBUG) System.err.println("starting queue="+queue);
    TreeSet<ArrayList<Integer>> otherLevelRedundancy=new TreeSet<ArrayList<Integer>>(new ProteinListComparator());
    Collections.sort(queue, new AlignedGraphletsComparator(otherLevelRedundancy));
    if (otherLevelRedundancy.size()!=0) {
      //if (DEBUG) System.err.println("sorted queue="+queue);
      ListIterator<AlignedGraphlets> itDG=queue.listIterator();
      if (itDG.hasNext()) {
        AlignedGraphlets dg2=null;
        AlignedGraphlets dg=itDG.next();
        int last=dg.hitGroupNoList.size()-1;
        TreeSet<Integer> noSet;
        while (itDG.hasNext()) {
          ArrayList<Integer> a=new ArrayList<Integer>(dg.alignedProteinListList.get(last));
          if (ABBREVIATION) Collections.sort(a);
          if (otherLevelRedundancy.contains(a)) {
            //if (DEBUG) System.err.println("redundant: "+dg.alignedProteinListList.get(last)+" => "+a);
            noSet=dg.hitGroupNoList.get(last);
            while (itDG.hasNext()) {
              ArrayList<Integer> a2=new ArrayList<Integer>((dg2=itDG.next()).alignedProteinListList.get(last));
              //if (DEBUG) System.err.println("check "+dg2.alignedProteinListList.get(last)+" => "+a2);
              if (ABBREVIATION) Collections.sort(a2);
              if (a2.equals(a)) {
                //if (DEBUG) System.err.println("remove it");
                noSet.addAll(dg2.hitGroupNoList.get(last));
                itDG.remove();
              } else break;
            }
            otherLevelRedundancy.remove(a);
            dg=dg2;
          } else
            dg=itDG.next();
        }
      }
    }
    //use AlignedGraphlets's comparator, comparing both group no's & content of last group
    Collections.sort(queue);
    //if (DEBUG) System.err.println("after removing redundancy, queue="+queue);
  }

  public GraphletAlign(String[] args) {this.args=args;}

  private int getCombination (int n, int r) {
    double c = (double) n/r;
    for (int i=1;i<r;i++) {
      c *= ((double)(n-i))/(r-i);
    }
    return (int) c;
  }

  private BigDecimal calculateBlastScore(AlignedGraphlets dg) {
    int i,j,k;
    ArrayList<Integer> a;
    BigDecimal total_score=new BigDecimal(1);
    BigDecimal score;
    TreeMap<Integer,BigDecimal> m;
    if (EXPECTATION_OR_BIT_SCORE) {
      for (i=0;i<dg.alignedProteinListList.size();i++) {
        a=dg.alignedProteinListList.get(i);
        for (j=0;j<a.size();j++) {
          for (k=0;k<a.size();k++) {
            if (j==k) continue;
            m=blastScores.get(a.get(j));
            if (m==null) {
              if (DEBUG) System.err.println(dg);
              if (DEBUG) System.err.println(a.get(j)+","+a.get(k));
              score=new BigDecimal(1);
            } else {
              score=m.get(a.get(k));
              if (score==null) {
                if (DEBUG) System.err.println(dg);
                if (DEBUG) System.err.println(a.get(j)+","+a.get(k));
                score=new BigDecimal(1);
              }
            }
            if (score.equals(BigDecimal.valueOf(0))) {
              score=new BigDecimal("1E-180");
              if (DEBUG) System.err.println("set 1E-180");
            }
            if (DEBUG) System.err.print(a.get(j)+","+a.get(k)+"=");
            if (DEBUG) System.err.println(score+"");
            total_score=total_score.multiply(score);
            if (DEBUG) System.err.println("=>"+total_score);
          }
        }
      }
    } else {
      for (i=0;i<dg.alignedProteinListList.size();i++) {
        a=dg.alignedProteinListList.get(i);
        for (j=0;j<a.size();j++) {
          for (k=0;k<a.size();k++) {
            if (j==k) continue;
            m=blastScores.get(a.get(j));
            if (m==null) {
              //if (DEBUG) System.err.println(dg);
              //if (DEBUG) System.err.println(a.get(j)+","+a.get(k));
              score=new BigDecimal(0);
            } else {
              score=m.get(a.get(k));
              if (score==null) {
                //if (DEBUG) System.err.println(dg);
                //if (DEBUG) System.err.println(a.get(j)+","+a.get(k));
                score=new BigDecimal(0);
              }
            }
            //if (DEBUG) System.err.print(a.get(j)+","+a.get(k)+"=");
            //if (DEBUG) System.err.println(score+"");
            total_score=total_score.add(score);
            //if (DEBUG) System.err.println("=>"+total_score);
          }
        }
      }
    }
    return total_score;
  }

  private void calculateHitGroupIntersection(TreeMap<ArrayList<Integer>,TreeSet<Integer>> hitGroupIntersection) {
    int i,j;
    TreeSet<Integer> t;
    for (i=0;i<hitGroupList.size();i++) {
      for (j=i+1;j<hitGroupList.size();j++) {
        ArrayList<Integer> a=new ArrayList<Integer>(hitGroupList.get(i));
        a.retainAll(hitGroupList.get(j));
        if (a.size()!=0) {
          if ((t=hitGroupIntersection.get(a))==null)
            hitGroupIntersection.put(a,t=new TreeSet<Integer>());
          t.add(i); t.add(j);
        }
      }
    }
  }

  class RedundantHitGroupInfo {
    boolean occurs;
    TreeSet<Integer> hitGroupNos;
    public RedundantHitGroupInfo(boolean occurs_,TreeSet<Integer> hitGroupNos_)
      {occurs=occurs_; hitGroupNos=hitGroupNos_;}
    public String toString() {
      return hitGroupNos.toString();
    }
  }

  //output of enumerateAlignedProteinSets()
  private TreeMap<ArrayList<Integer>,RedundantHitGroupInfo> firstLevelRedundancy
    =new TreeMap<ArrayList<Integer>,RedundantHitGroupInfo>(new ProteinListComparator()); //for first level in dg&dg tree
  private LinkedList<AlignedGraphlets> queue=new LinkedList<AlignedGraphlets>(); //queue of dg&dg tree

  private int getInitialCombinationArray(int[] proteinNoCombinationArray,
                                  int[] boundary,
                                  int speciesNo,
                                  ArrayList<Integer> group) {
    int i;
    if (SINGLE_SPECIES) {
      for (i=0;i<NUM_DISJOINT_VERTICES;i++)
        proteinNoCombinationArray[i]=i;
      //if (DEBUG) System.err.println("proteinNoCombinationArray="+Arrays.toString(proteinNoCombinationArray));

      return getCombination(group.size(),NUM_DISJOINT_VERTICES);
    } else {
      //find the boundary
      boundary[0]=0;
      int pre_i=0;
      int j;
      for (j=1;j<boundary.length-1;j++) {
        for (i=pre_i;i<group.size();i++) {
          if (group.get(i)>SPECIES_NUM_LINES.get(j-1)) {
            boundary[j]=i;
            break;
          }
        }
        if (i==group.size()) break;
        pre_i=i;
      }
      for (i=j;i<boundary.length;i++)
        boundary[i]=group.size();
/* // errorneous
      int j=0;
      boundary[0]=0;
      for (i=0;i<group.size();i++)
        if (group.get(i)>SPECIES_NUM_LINES[j]) {
          boundary[j+1]=i;
          //if (DEBUG) System.err.println("boundry["+(j+1)+"]="+i);
          j++;
          if (j+1==boundary.length-1) break;
        }
      for (i=j+1;i<boundary.length;i++)
        boundary[i]=group.size();
*/

      if (DEBUG) {
        System.err.print("boundary:");
        for (i=0;i<boundary.length;i++)
          System.err.print(boundary[i]+",");
        System.err.println("");
      }

      if (boundary[speciesNo+1]-boundary[speciesNo]<numSubgraphsInOneSpecies)
        return 0;

      //calculate # combinations
      int numCombinations=1;
      for (i=0;i<boundary.length-1;i++)
        if (i!=speciesNo)
          numCombinations*=(boundary[i+1]-boundary[i]);
      numCombinations*=getCombination(boundary[speciesNo+1]-boundary[speciesNo],numSubgraphsInOneSpecies);

      for (i=0;i<num_species;i++)
        if (i!=speciesNo)
          proteinNoCombinationArray[speciesNo_to_subgraphNo.get(i)]=boundary[i];
      for (i=0;i<numSubgraphsInOneSpecies;i++)
        proteinNoCombinationArray[i+speciesNo]=boundary[speciesNo]+i;
      if (DEBUG) System.err.println("initial proteinNoCombinationArray="+Arrays.toString(proteinNoCombinationArray));

      return numCombinations;
    }
  }

  private boolean enumerateAlignedProteinSets(TreeMap<ArrayList<Integer>,TreeSet<Integer>> in,
                               int queryDegree,
                               boolean removeRedundancy,
                               int[] proteinNoCombinationArray,
                               int[] boundary,
                               int speciesNo) {

    int numCombinations;
    ArrayList<Integer> tmpProteinList;
    int i,j,k,l;
    int int1;
    ArrayList<Integer> group;
    boolean lastCombination=false;

    Iterator<ArrayList<Integer>> itArrayInteger=in.keySet().iterator();
    while (itArrayInteger.hasNext()) {
      group=itArrayInteger.next();
      //if (DEBUG && group.size()!=0) System.err.println(group);

      if (group.size()<NUM_DISJOINT_VERTICES) {
        lastCombination=true;
        continue;
      }
      //======== fill the first group of proteins ========
      if (in.size()==1)
        numCombinations=1;
      else
        numCombinations=getInitialCombinationArray(proteinNoCombinationArray,boundary,speciesNo,group);

      if (DEBUG) System.err.println("group="+group.toString());
      if (DEBUG && in.size()!=1) System.err.println("numCombinations="+numCombinations);

      if (SINGLE_SPECIES) {

        for (i=0;i<numCombinations;i++) {
          //if (DEBUG && (i+1)%100000==0) System.err.println("Combination "+i);

          tmpProteinList=new ArrayList<Integer>();
          for (j=0;j<NUM_DISJOINT_VERTICES;j++)
            tmpProteinList.add(group.get(proteinNoCombinationArray[j]));

          //check degree
          if (DISTANCE2) {
            for (j=0;j<tmpProteinList.size();j++)
              if (proteinToNeighborArray.get(int1=tmpProteinList.get(j)).size()
                  + proteinToNeighborArray2.get(int1).size() < queryDegree)
                break;
          } else {
            for (j=0;j<tmpProteinList.size();j++)
              if (proteinToNeighborArray.get(tmpProteinList.get(j)).size() < queryDegree) {
/*
                if (DEBUG) {
                  System.err.println(tmpProteinList.get(j)+"'s degree < queryDegree "+queryDegree);
                  System.err.println("doesn't add combination "+tmpProteinList);
                }
*/
                break;
              }
          }
          if (j==tmpProteinList.size()) {
            RedundantHitGroupInfo rgi;
            //generate a set of redundant groups
            if (!removeRedundancy) {
              if ((rgi=firstLevelRedundancy.get(tmpProteinList))==null)
                firstLevelRedundancy.put(tmpProteinList,new RedundantHitGroupInfo(false,in.get(group)));
              else
                rgi.hitGroupNos.addAll(in.get(group));
            //generate first level groups
            } else if ((rgi=firstLevelRedundancy.get(tmpProteinList))==null) {
              queue.add(new AlignedGraphlets(new TreeSet<Integer>(in.get(group)),tmpProteinList));
              //if (DEBUG) System.err.println("add combination "+tmpProteinList);
            } else if (rgi.occurs==false) {
              queue.add(new AlignedGraphlets(new TreeSet<Integer>(rgi.hitGroupNos),tmpProteinList));
              //if (DEBUG) System.err.println("add combination "+tmpProteinList);
              rgi.occurs=true;
            } //else if (DEBUG) System.err.println("doesn't add combination "+tmpProteinList);

          }

          // start from last index, find the element that does not have the largest possible value
          for (k=NUM_DISJOINT_VERTICES-1; k>=0 && proteinNoCombinationArray[k]==group.size()-(NUM_DISJOINT_VERTICES-k); k--)
            ;
          if (k==-1) { // all elements have their largest possible values
            //if (DEBUG) System.err.println("lastCombination=true");
            lastCombination=true;
            break;
          }

          proteinNoCombinationArray[k]++; // increment the element just found
          for (j=k+1;j<NUM_DISJOINT_VERTICES;j++)
            proteinNoCombinationArray[j]=proteinNoCombinationArray[k]+(j-k); // the element is followed by an increasing sequence

        } // for each combination i

      } else { //@ !SINGLE_SPECIES

        for (i=0;i<numCombinations;i++) {
          //if ((i+1)%100000==0) System.err.println("Combination "+i);
          // allow multiple proteins in a species
          /*
          k=0;
          for (j=0;j<NUM_DISJOINT_VERTICES;j++)
            if (proteinNoCombinationArray[j]>=boundary[k]) {
              if (proteinNoCombinationArray[j]<boundary[k+1]) {
                k++;
                if (k==boundary.length-1) break;
              } else
                break;
            }

          if (k==boundary.length-1) {
          */

          tmpProteinList=new ArrayList<Integer>();
          if (DEBUG) System.err.print("add ");
          for (j=0;j<NUM_DISJOINT_VERTICES;j++) {
            tmpProteinList.add(group.get(proteinNoCombinationArray[j]));
            if (DEBUG) System.err.print(proteinNoCombinationArray[j]+":"+group.get(proteinNoCombinationArray[j])+",");
          }
          if (DEBUG) System.err.println("");
          if (DEBUG) System.err.println("add combination "+tmpProteinList);

          //check degree
          if (DISTANCE2) {
            for (j=0;j<tmpProteinList.size();j++)
              if (proteinToNeighborArray.get(int1=tmpProteinList.get(j)).size()
                  + proteinToNeighborArray2.get(int1).size() < queryDegree)
                break;
          } else {
            for (j=0;j<tmpProteinList.size();j++)
              if (proteinToNeighborArray.get(tmpProteinList.get(j)).size() < queryDegree)
                break;
          }
          if (j==tmpProteinList.size()) {
            RedundantHitGroupInfo rgi;
            //generate a set of redundant groups
            if (!removeRedundancy) {
              if ((rgi=firstLevelRedundancy.get(tmpProteinList))==null)
                firstLevelRedundancy.put(tmpProteinList,new RedundantHitGroupInfo(false,in.get(group)));
              else
                rgi.hitGroupNos.addAll(in.get(group));
            //generate first level groups
            } else if ((rgi=firstLevelRedundancy.get(tmpProteinList))==null) {
              queue.add(new AlignedGraphlets(new TreeSet<Integer>(in.get(group)),tmpProteinList));
              if (DEBUG) System.err.println("add combination "+tmpProteinList);
            } else if (rgi.occurs==false) {
              queue.add(new AlignedGraphlets(new TreeSet<Integer>(rgi.hitGroupNos),tmpProteinList));
              if (DEBUG) System.err.println("add combination "+tmpProteinList);
              rgi.occurs=true;
            } else if (DEBUG) System.err.println("doesn't add combination "+tmpProteinList);

          }

          //increment indices
          for (j=num_species-1;j>=0;j--) {
            if (j==speciesNo) {
              // start from last index, find the element that does not have the largest possible value
              for (k=speciesNo+numSubgraphsInOneSpecies-1; k>=speciesNo
                && proteinNoCombinationArray[k]==boundary[subgraphNo_to_speciesNo.get(k)+1]+(k-speciesNo-numSubgraphsInOneSpecies); k--)
                ;
              if (DEBUG) System.err.println("k="+k);
              if (k==speciesNo-1) { // all elements have their largest possible values
                //reset index
                for (k=0;k<numSubgraphsInOneSpecies;k++)
                  proteinNoCombinationArray[k+speciesNo]=boundary[speciesNo]+k;
                if (DEBUG) System.err.println("lastCombination=true");
                if (speciesNo==0) lastCombination=true;
              } else {
                proteinNoCombinationArray[k]++; // increment the element just found
                for (l=k+1;l<speciesNo+numSubgraphsInOneSpecies;l++)
                  proteinNoCombinationArray[l]=proteinNoCombinationArray[k]+(l-k); // the element is followed by an increasing sequence
                break;
              }
            } else {
              proteinNoCombinationArray[speciesNo_to_subgraphNo.get(j)]++;
              if (proteinNoCombinationArray[speciesNo_to_subgraphNo.get(j)]==boundary[j+1]) {
                proteinNoCombinationArray[speciesNo_to_subgraphNo.get(j)]=boundary[j];
                if (j==0) {
                  if (DEBUG) System.err.println("lastCombination=true");
                  lastCombination=true;
                }
              } else {
                break;
              }
            }
          }

        } // for each combination i
      } //if (SINGLE_SPECIES)
      if (DEBUG && in.size()!=1) System.err.println("leaving function: proteinNoCombinationArray="+Arrays.toString(proteinNoCombinationArray));
    } // while
    return lastCombination;
  }

  public void branchNBound() {
    try {
      BufferedReader breader;
      String str,str2;
      int i,j,k, groupNo;
      int int1,int2;
      StringTokenizer st,st2;

      // ========== read arguments ==========
/*
      System.err.println("args.length="+args.length);
      for (i=0;i<args.length;i++)
        System.err.println(args[i]);
*/

      File f=new File(args[1]);
      
      /*
      // find species names (no use)
      st=new StringTokenizer(str,"_");
      species=new String[st.countTokens()];
      i=0;
      System.err.print("species: ");
      while (st.hasMoreTokens()) {
        species[i++]=st.nextToken();
        System.err.print(species[i-1]+" ");
      }
      System.err.println("");
      */

      // need change
      NUM_DISJOINT_VERTICES=Integer.parseInt(args[5]);
      System.err.println("Find "+NUM_DISJOINT_VERTICES+" aligned graphlets");
      
      breader=new BufferedReader(new FileReader(args[6]));
      i=0;
      while ((str=breader.readLine())!=null) {
        SPECIES_NUM_LINES.add(Integer.parseInt(str));
        System.err.println("Network "+i+": protein "+(i==0?1:SPECIES_NUM_LINES.get(i-1)+1)+" to "+SPECIES_NUM_LINES.get(i));
        i++;
      }
      breader.close();
      
      if (SPECIES_NUM_LINES.size()>1) {
        SINGLE_SPECIES=false;
        num_species=SPECIES_NUM_LINES.size();
        System.err.println(num_species+" species");
      } else
        System.err.println("single species");
        
      SPECIES_GRAPH_SIZE=SPECIES_NUM_LINES.get(SPECIES_NUM_LINES.size()-1);

      if (!SINGLE_SPECIES) {
        numSubgraphsInOneSpecies=NUM_DISJOINT_VERTICES-num_species+1;
        //System.err.println("numSubgraphsInOneSpecies="+numSubgraphsInOneSpecies);
      } else
        numSubgraphsInOneSpecies=-1;
        
      // ========== read data ==========

      // read group.list
      hitGroupList=new ArrayList<ArrayList<Integer>>(SPECIES_GRAPH_SIZE); // data structure

      //hit group no's
      proteinToHitGroupNosList=new ArrayList<ArrayList<Integer>>(SPECIES_GRAPH_SIZE+1);
      for (i=0;i<SPECIES_GRAPH_SIZE+1;i++)
        proteinToHitGroupNosList.add(new ArrayList<Integer>());
      sameGroup=new boolean[SPECIES_GRAPH_SIZE+1][SPECIES_GRAPH_SIZE+1];
      for (i=0;i<sameGroup.length;i++)
        for (j=0;j<sameGroup[0].length;j++)
          sameGroup[i][j]=false;
      ArrayList<Integer> hitGroupNoList; //tmp variable

      ArrayList<Integer> tmpHitGroup;
      breader=new BufferedReader(new FileReader(args[0]));
      groupNo = 0;
      while ((str=breader.readLine())!=null) {
        st=new StringTokenizer(str,"\\[\\], ");
        tmpHitGroup=new ArrayList<Integer>();
        if (st.countTokens() >= NUM_DISJOINT_VERTICES) {
          while (st.hasMoreTokens()) {
            tmpHitGroup.add(int1=Integer.parseInt(st.nextToken()));
            proteinToHitGroupNosList.get(int1).add(new Integer(groupNo));
          }
        } /* else {
          while (st.hasMoreTokens()) {
            int1=Integer.parseInt(st.nextToken());
            proteinToHitGroupNosList[int1]=-1;
          }
        } */
        hitGroupList.add(tmpHitGroup);
        groupNo++;
        for (i=0;i<tmpHitGroup.size();i++) {
          for (j=0;j<tmpHitGroup.size();j++) {
            sameGroup[tmpHitGroup.get(i)][tmpHitGroup.get(j)]=true;
          }
        }
      }
      breader.close();
/*
      if (DEBUG) {
        for (i=0;i<proteinToHitGroupNosList.size();i++)
          System.err.println(proteinToHitGroupNosList.get(i).toString());
      }
*/
      //neighbor
      breader=new BufferedReader(new FileReader(args[1]));
      proteinToNeighborArray=new ArrayList<ArrayList<Integer>>(SPECIES_GRAPH_SIZE+1); // data structure

      //@ DISTANCE2
      BufferedReader breader2=null;
      if (DISTANCE2) {
        breader2=new BufferedReader(new FileReader(args[1]+".2")); // distance-2 neighbors
        proteinToNeighborArray2=new ArrayList<ArrayList<Integer>>(SPECIES_GRAPH_SIZE+1); // data structure
      }

      Integer intObject;
      Integer groupNoObject;

      i=1;
      proteinToNeighborArray.add(new ArrayList<Integer>());
      if (DISTANCE2) proteinToNeighborArray2.add(new ArrayList<Integer>());
      while ((str=breader.readLine())!=null) {
        proteinToNeighborArray.add(new ArrayList<Integer>());        
        if (DISTANCE2) {
          str2=breader2.readLine();
          proteinToNeighborArray2.add(new ArrayList<Integer>());
        }

        //if (proteinToHitGroupNosList[i]==-1) {
        if (proteinToHitGroupNosList.get(i).size()==0) {
          i++;
          continue;
        }

        // read neighbors

        st=new StringTokenizer(str,"\\[\\]");
        str=st.nextToken();
        int2=Integer.valueOf(str);
        if (st.hasMoreTokens()) {
          st2=new StringTokenizer(st.nextToken(),";, ");
          while (st2.hasMoreTokens()) {
            int1=Integer.valueOf(st2.nextToken());
            //if (proteinToHitGroupNosList[int1]!=-1)
            if (proteinToHitGroupNosList.get(i).size()!=0 && int1!=int2)
              proteinToNeighborArray.get(i).add(int1);
          }
        }

        //@ DISTANCE2
        if (DISTANCE2) {
          st=new StringTokenizer(str2,"\\[\\]");
          st.nextToken();
          if (st.hasMoreTokens()) {
            st2=new StringTokenizer(st.nextToken(),";, ");
            while (st2.hasMoreTokens()) {
              int1=Integer.valueOf(st2.nextToken());
              //if (proteinToHitGroupNosList[int1]!=-1)
              if (proteinToHitGroupNosList.get(i).size()!=0)
                proteinToNeighborArray2.get(i).add(int1);
            }
          }
        }

        i++;
      } //while
      breader.close();
      if (DISTANCE2) breader2.close(); //@ DISTANCE2

      //@ BLAST_SCORE
      BufferedReader breader3=null;
      if (BLAST_SCORE) {
        //[*.hits] [*.adjl] [*.blastp.no] [*.adj] [*.auto] [# disjoint vertices] [# accumulated proteins from each species]
        breader3=new BufferedReader(new FileReader(args[2])); // blastp.pairs
        blastScores=new ArrayList<TreeMap<Integer,BigDecimal>>(SPECIES_GRAPH_SIZE+1);
        for (i=0;i<SPECIES_GRAPH_SIZE+1;i++)
          blastScores.add(new TreeMap<Integer,BigDecimal>());
        while ((str=breader3.readLine())!=null) {
          st=new StringTokenizer(str,"\t ");
          TreeMap<Integer,BigDecimal> blastScoreMap=blastScores.get(i=Integer.parseInt(st.nextToken()));
          if (blastScoreMap==null)
            blastScoreMap=new TreeMap<Integer,BigDecimal>();
          blastScoreMap.put(Integer.parseInt(st.nextToken()),new BigDecimal(st.nextToken()));
        }
        breader3.close(); //@ BLAST_SCORE
      }

      //AlignedGraphlets & bound
      //data structures
      Stack<LinkedList<AlignedGraphlets>> stack=new Stack<LinkedList<AlignedGraphlets>>();

      //persistent variables
      LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>> neighboringGroupNoToProteinsMap= new LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>>();
      //@ DISTANCE2
      LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>> neighboringGroupNoToProteinsMap2=null;
      if (DISTANCE2) neighboringGroupNoToProteinsMap2 = new LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>>();

      AlignedGraphlets dg,dg2;
      int edgeIndex;
      boolean addNewProtein;

      //temporary variables
      TreeSet<Integer> proteinSet,proteinSet2,proteinSet3;
      Set<Integer> groupNoSet;
      ArrayList<TreeSet<Integer>> proteinSetList=new ArrayList<TreeSet<Integer>>(NUM_DISJOINT_VERTICES);
      ArrayList<Integer> neighboringGroupProteinsSubset;
      ArrayList<ArrayList<Integer>> neighboringGroupProteins,neighboringGroupProteins2;
      TreeSet<Integer> previousProteinSet=new TreeSet<Integer>();
      Iterator<Integer> itInteger;
      Iterator<Integer> itInt;
      ListIterator<Integer> litInt;
      ListIterator<Integer> litInteger;
      ListIterator<AlignedGraphlets> litDG;
      RedundantHitGroupInfo rgi;
      TreeSet<Integer> groupNos;
      Integer previousProtein;
      int queryDegree,previousDegreeReduced,previousGroupIndex;
      ArrayList<Integer> group, group2;
      boolean hasNeighbors;
      boolean toSkip;

      TreeMap<ArrayList<Integer>,TreeSet<Integer>> in=new TreeMap<ArrayList<Integer>,TreeSet<Integer>>(new ProteinListComparator());
      boolean lastCombination;
      int numCombinations;
      int totalNumPairs,numPairs;

      // ====== read topology ======

      ArrayList<ArrayList<Integer>> topology=new ArrayList<ArrayList<Integer>>();
      TreeMap<Integer,ArrayList<Integer>> rules=new TreeMap<Integer,ArrayList<Integer>>();
      String strAdj,strRule,fileSuffix;
      subgraphNo_to_speciesNo=new TreeMap<Integer,Integer>();
      speciesNo_to_subgraphNo=new TreeMap<Integer,Integer>();

      for (int speciesNo=0;speciesNo<num_species;speciesNo++) {
      if (speciesNo==1 && (SINGLE_SPECIES || NUM_DISJOINT_VERTICES==num_species)) break;
      //if (DEBUG) System.err.println("|||||||||||||||||||||speciesNo="+speciesNo+"||||||||||||||||||||||");

      subgraphNo_to_speciesNo.clear();
      speciesNo_to_subgraphNo.clear();
      for (i=0;i<=speciesNo;i++) {
        subgraphNo_to_speciesNo.put(i,i);
        speciesNo_to_subgraphNo.put(i,i);
      }
      for (i=speciesNo+1;i<speciesNo+numSubgraphsInOneSpecies;i++)
        subgraphNo_to_speciesNo.put(i,speciesNo);
      for (i=speciesNo+numSubgraphsInOneSpecies;i<NUM_DISJOINT_VERTICES;i++) {
        subgraphNo_to_speciesNo.put(i,i-NUM_DISJOINT_VERTICES+num_species);
        speciesNo_to_subgraphNo.put(i-NUM_DISJOINT_VERTICES+num_species,i);
      }
      //if (DEBUG) System.err.println("subgraphNo_to_speciesNo="+subgraphNo_to_speciesNo.toString());
      //if (DEBUG) System.err.println("speciesNo_to_subgraphNo="+speciesNo_to_subgraphNo.toString());

      // ====== calculate intersection ======

      //for first level in dg&dg tree
      TreeMap<ArrayList<Integer>,TreeSet<Integer>> hitGroupIntersection=new TreeMap<ArrayList<Integer>,TreeSet<Integer>>(new ProteinListComparator()); // group list
      int[] proteinNoCombinationArray=new int[NUM_DISJOINT_VERTICES];
      int[] boundary=null; //@MULTI_SPECIES
      if (!SINGLE_SPECIES)
        boundary=new int[SPECIES_NUM_LINES.size()+1];

      calculateHitGroupIntersection(hitGroupIntersection);
      //if (DEBUG) System.err.println("hitGroupIntersection="+hitGroupIntersection.toString());
      if (DEBUG) System.err.println("hitGroupIntersection's size="+hitGroupIntersection.size());
      enumerateAlignedProteinSets(hitGroupIntersection,(int)0,false,proteinNoCombinationArray,boundary,speciesNo); //generates firstLevelRedundancy

      //if (DEBUG) System.err.println("firstLevelRedundancy = "+firstLevelRedundancy);
      if (DEBUG) System.err.println("# redundant sets of NUM_DISJOINT_VERTICES vertices = "+firstLevelRedundancy.size());

      // ====== read topology and rules ======
      //[*.hits] [*.adjl] [*.blastp.no] [*.adj] [*.auto] [# disjoint vertices] [# accumulated proteins from each species]
      if (DEBUG) System.err.println("Open topology file "+args[3]);
      BufferedReader breaderAdj=new BufferedReader(new FileReader(args[3])); //*.adj file
      if (DEBUG) System.err.println("Open automorphism file "+args[4]);
      BufferedReader breaderAuto=new BufferedReader(new FileReader(args[4])); //*.auto
      
	    System.err.println("Reading topology...");
      while ((strAdj=breaderAdj.readLine())!=null) {

      if (DEBUG) System.err.println(strAdj);
      st=new StringTokenizer(strAdj,"-");
      NUM_VERTICES=Integer.parseInt(st.nextToken());
      System.err.println(NUM_VERTICES+" vertices");
      NUM_EDGES=Integer.parseInt(st.nextToken());
      System.err.println(NUM_EDGES+" edges");

      fileSuffix=strAdj;

      topology.clear();
      strAdj=breaderAdj.readLine();
      st=new StringTokenizer(strAdj," ");
      int1=Integer.parseInt(st.nextToken());
      int2=Integer.parseInt(st.nextToken());
      topology.add(new ArrayList<Integer>());
      topology.add(new ArrayList<Integer>());
      topology.get(int1).add(int2);
      topology.get(int2).add(int1);
      while (!(strAdj=breaderAdj.readLine()).equals("")) {
        st=new StringTokenizer(strAdj," ");
        int1=Integer.parseInt(st.nextToken());
        int2=Integer.parseInt(st.nextToken());
        //if (DEBUG) System.err.println(int1+" "+int2);
        if (int1>=topology.size()) {
          System.err.println("traversal bug");
          return;
        }
        if (int2==topology.size()) {
          topology.add(new ArrayList<Integer>());
        } else if (int2>topology.size()) {
          System.err.println("traversal bug");
          return;
        }
        topology.get(int1).add(int2);
        topology.get(int2).add(int1);
      }

      if (DEBUG) {
        System.err.print("degree:");
        for (i=0;i<NUM_VERTICES;i++)
          System.err.print(i+"("+topology.get(i).size()+"):"+topology.get(i).toString());
        System.err.println("");
      }
      
      for (i=0;i<NUM_VERTICES;i++) {
        System.err.print("Node "+i+": ");
        for (j=0;j<topology.get(i).size();j++) {
          System.err.print(topology.get(i).get(j)+" ");
        }
        System.err.println("");
      }
      
      Integer intObj,intObj2;

      rules.clear();
      for (int1=1;int1<NUM_VERTICES;int1++)
        rules.put(new Integer(int1), new ArrayList<Integer>());
      System.err.println("Reading rules...");
      while (!(strRule=breaderAuto.readLine()).equals("")) {
        st=new StringTokenizer(strRule," ");
        intObj=new Integer(st.nextToken());
        intObj2=new Integer(st.nextToken());
        rules.get(intObj2).add(intObj);
      }
      Iterator<Integer> it=rules.keySet().iterator();
      while (it.hasNext()) {
        int1=it.next();
        System.err.println(int1+">"+rules.get(int1).toString());
      }

      //if (fileSuffix.equals("4-3-1") || fileSuffix.equals("4-3-2")) continue;
      //File file = new File(args[args.length-1]+"/"+(str=f.getName().substring(0,f.getName().length()-5))+"."+NUM_DISJOINT_VERTICES+"."+speciesNo+"."+fileSuffix);
      File file = new File(args[args.length-1]+"/"+(str=f.getName().substring(0,f.getName().length()-5))+"."+fileSuffix);
      BufferedWriter bwriter=new BufferedWriter(new FileWriter(file));
      System.err.println("Write "+file.getName());

      //BufferedWriter[] bwriterProteinSet=new BufferedWriter[NUM_DISJOINT_VERTICES];
      String pre_prefix,prefix; //@ OUTPUT_PROTEIN_SETS
      //ArrayList<TreeSet<Integer>> proteinSets; //@ OUTPUT_PROTEIN_SETS
      ArrayList<Integer> tmpVertexList;

/*
      if (OUTPUT_PROTEIN_SETS) {
        pre_prefix="";
        proteinSets=new ArrayList<TreeSet<Integer>>(NUM_DISJOINT_VERTICES);
        for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
          bwriterProteinSet[i]=new BufferedWriter(new FileWriter(args[args.length-1]+"/"+(str=f.getName().substring(0,f.getName().length()-5))+"."+NUM_DISJOINT_VERTICES+".set."+i,true)); // protein set;
          proteinSets.add(new TreeSet<Integer>());
        }
      }
*/

      totalNumPairs=0;

      // ====== go over every protein group ======
/*
      queryDegree=degree[edges[0][0]];
      for (int g=0;g<hitGroupList.size();g++) {
        enumerateAlignedProteinSets(g,queryDegree,queue);
      } // for g

      if (!queue.isEmpty()) {
        Collections.sort(queue, new AlignedGraphletsComparator(firstLevelRedundancy));
        //if (DEBUG) System.err.println("unique queue="+queue);
        //stack.push(queue);
        queue.clear();
      }
*/

      //for (int g=201;g<=201;g++) {
      System.err.println("Total "+hitGroupList.size()+" blastp hit groups");
      for (int g=0;g<hitGroupList.size();g++) {
        System.err.println("Blastp hit group "+g);
        //if (DEBUG) System.err.println(hitGroupList.get(g).toString());

        in.clear();
        TreeSet<Integer> t=new TreeSet<Integer>();
        t.add(g);
        in.put(group=hitGroupList.get(g),t);
        queue.clear();

        //======== fill the first group of proteins ========

        numCombinations=getInitialCombinationArray(proteinNoCombinationArray,boundary,speciesNo,group);
        System.err.println("  Total "+numCombinations+" sets of aligned nodes for the first node of the graphlet");

        //calculate # pairs
        if (SINGLE_SPECIES) {
          numPairs = getCombination(group.size(),2);
          if (DEBUG) System.err.println("numPairs="+numPairs);
          totalNumPairs += numPairs;
        } else {
          numPairs=0;
          for (i=1;i<boundary.length;i++)
            for (j=i+1;j<boundary.length;j++)
              numPairs+=(boundary[i]-boundary[i-1])*(boundary[j]-boundary[j-1]);
          totalNumPairs+=numPairs;
          if (DEBUG) System.err.println("numPairs="+numPairs);

        }

        if (numCombinations>0) {
          lastCombination=false;
          //continue;
        }
        else lastCombination=true;

        numCombinations=0;
        while (!lastCombination) {
          numCombinations++;
          if (numCombinations%5000==0) System.err.println(numCombinations);
          queryDegree=(int)topology.get(0).size();
          lastCombination=enumerateAlignedProteinSets(in,queryDegree,true,proteinNoCombinationArray,boundary,speciesNo); //considers firstLevelRedundancy
          
          //if (g==11 && numCombinations<5000) continue;
          if (DEBUG) System.err.println(queue.toString());
          if (!queue.isEmpty())
            //if (queue.peek().hitGroupList.get(0).get(0).equals(Integer.valueOf((int)1167))
                //&& queue.peek().hitGroupList.get(0).get(1).equals(Integer.valueOf((int)7116)))
              stack.push(queue);
            //else queue.clear();

          //====== traversing the stack of queues of AlignedGraphlets's ======
          while (!stack.empty()) {

            dg=stack.peek().poll();

            //find from which protein to extend
            tmpVertexList=topology.get(dg.hitGroupNoList.size());
            previousGroupIndex=NUM_VERTICES;
            for (int1=0;int1<tmpVertexList.size();int1++) {
              if ((int2=tmpVertexList.get(int1))<previousGroupIndex)
                previousGroupIndex=int2;
            }

            if (DEBUG) System.err.println("=== edge "+previousGroupIndex+"-"+dg.hitGroupNoList.size()+" === "+dg);

            // new queue in the next level
            queue=new LinkedList<AlignedGraphlets>();

            // determine if the protein to extend to is already in the list
/*
            addNewProtein=true;
            if (edges[edgeIndex][1]<dg.hitGroupNoList.size())
              addNewProtein=false;
            else if (edges[edgeIndex][1]!=dg.hitGroupNoList.size()) {
              System.err.println("BUG in input pattern!!");
              System.exit(0);
            }
*/
            //====== collect the neighbors of the group =====
            //====== and the group which each neighbor belongs to =====
            //====== store them in neighboringGroupNoToProteinsMap, =====
            //====== neighboringGroupNoToProteinsMap2, =====
            //====== neighboringProteinsMap =====

            previousDegreeReduced=(int)topology.get(int1=previousGroupIndex).size();
            if (DEBUG) System.err.println("previousDegreeReduced="+previousDegreeReduced);
            tmpVertexList=topology.get(previousGroupIndex);
            for (int2=0;int2<dg.hitGroupNoList.size();int2++) {
              if (tmpVertexList.contains(int2)) //connected
                previousDegreeReduced--;
            }
            if (DEBUG) System.err.println("previousDegreeReduced="+previousDegreeReduced);

            group=dg.alignedProteinListList.get(previousGroupIndex);
            queryDegree=(int)topology.get(dg.hitGroupNoList.size()).size();

            // a hash map that maps a group no to a set of neighbors
            neighboringGroupNoToProteinsMap.clear();
            if (DISTANCE2) neighboringGroupNoToProteinsMap2.clear();

            //if (addNewProtein) {

              previousProteinSet.clear();
              for (i=0;i<dg.alignedProteinListList.size();i++)
                previousProteinSet.addAll(dg.alignedProteinListList.get(i));

              if (!DISTANCE2) {

                //====== process neighboringGroupNoToProteinsMap ======
                for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                  previousProtein=group.get(i);

                  proteinSet=new TreeSet<Integer>(proteinToNeighborArray.get(previousProtein));
                  //if (DEBUG) System.err.println("previous protein "+previousProtein+"'s neighbors"+proteinSet);

                  //cannot appear in current graph
                  proteinSet.removeAll(previousProteinSet);

                  //remove by degree
                  itInteger=proteinSet.iterator();
                  while (itInteger.hasNext())
                    if (proteinToNeighborArray.get(itInteger.next()).size() < queryDegree)
                      itInteger.remove();

                  //if (DEBUG) System.err.println("=>"+proteinSet);

                  hasNeighbors=false;
                  if (proteinSet.size()!=0) {
                    // check which group each neighbor belongs to
                    itInteger=proteinSet.iterator();
                    while (itInteger.hasNext()) {
                      intObject=itInteger.next();
                      //groupNo=proteinToHitGroupNosList[intObject];
                      hitGroupNoList=proteinToHitGroupNosList.get(intObject);
                      litInt=hitGroupNoList.listIterator();
                      while (litInt.hasNext()) {
                        groupNoObject=litInt.next();
                        //if (DEBUG) System.err.println("Insert "+intObject+" to group "+groupNoObject);

                        // save the neighbor's group no
                        if ((neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(groupNoObject))==null) {
                          if (i==0) { //first graph copy
                            neighboringGroupNoToProteinsMap.put(groupNoObject, neighboringGroupProteins=new ArrayList<ArrayList<Integer>>(NUM_DISJOINT_VERTICES));
                          } else
                            continue;
                        }
    /*
                        if (DEBUG) {
                          System.err.print("*distance 1 map* ");
                          System.err.print(groupNoObject+"=>");
                          for (j=0;j<neighboringGroupProteins.size();j++)
                            System.err.print(j+":"+neighboringGroupProteins.get(j)+",");
                          System.err.println("");
                        }
    */

                        // insert to the i-th tree set
                        if (neighboringGroupProteins.size()==i) {
                          neighboringGroupProteins.add(neighboringGroupProteinsSubset=new ArrayList<Integer>());
                          neighboringGroupProteinsSubset.add(intObject);
                          hasNeighbors=true;
                        } else if (neighboringGroupProteins.size()==i+1) {
                          neighboringGroupProteins.get(i).add(intObject);
                          hasNeighbors=true;
                        } else if (neighboringGroupProteins.size()>i+1)
                          System.err.println("bugbugbugbugbugbug neighboringGroupProteins.size()>i+1");
    /*
                        if (DEBUG) {
                          System.err.print("*distance 1 map* ");
                          System.err.print(groupNoObject+"=>");
                          for (j=0;j<neighboringGroupProteins.size();j++)
                            System.err.print(j+":"+neighboringGroupProteins.get(j)+",");
                          System.err.println("");
                        }
    */
                      }
                    }
                  } else {
                    neighboringGroupNoToProteinsMap.clear();
                    //f (DEBUG) System.err.println("one of the previous proteins has no neighbors");
                    break;
                  }

                  if (!hasNeighbors) {
                    neighboringGroupNoToProteinsMap.clear();
                    //if (DEBUG) System.err.println("one of the previous proteins has no usable neighbors");
                    break;
                  }
                } // for (i=0;i<NUM_DISJOINT_VERTICES;i++)

              } else { // if @ DISTANCE2

                //====== process neighboringGroupNoToProteinsMap ======
                proteinSetList.clear();
                for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                  previousProtein=group.get(i);

                  proteinSetList.add(proteinSet=new TreeSet<Integer>(proteinToNeighborArray.get(previousProtein)));
                  if (DEBUG) System.err.println("previous protein "+previousProtein+"'s neighbors"+proteinSet);

                  //cannot appear in current graph
                  proteinSet.removeAll(previousProteinSet);

                  //remove by degree
                  itInteger=proteinSet.iterator();
                  while (itInteger.hasNext())
                    if (proteinToNeighborArray.get(int1=itInteger.next()).size()
                        + proteinToNeighborArray2.get(int1).size() < queryDegree) {
                      itInteger.remove();
                    }

                  if (DEBUG) System.err.println("filtering by vertex degree=>"+proteinSet);

                  if (proteinSet.size()!=0) {

                    //check which group each neighbor belongs to
                    itInteger=proteinSet.iterator();
                    while (itInteger.hasNext()) {
                      intObject=itInteger.next();
                      //groupNo=proteinToHitGroupNosList[intObject];
                      hitGroupNoList=proteinToHitGroupNosList.get(intObject);

                      litInt=hitGroupNoList.listIterator();
                      while (litInt.hasNext()) {
                        groupNoObject=litInt.next();
                        if (DEBUG) System.err.println("Insert "+intObject+" to group "+groupNoObject);

                        // save the neighbor's group no
                        if ((neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(Integer.valueOf(groupNoObject)))==null) {
                          neighboringGroupNoToProteinsMap.put(new Integer(groupNoObject), neighboringGroupProteins=new ArrayList<ArrayList<Integer>>(NUM_DISJOINT_VERTICES));
                          neighboringGroupNoToProteinsMap2.put(new Integer(groupNoObject), new ArrayList<ArrayList<Integer>>(NUM_DISJOINT_VERTICES));
                        }
    
                        if (DEBUG) {
                          System.err.print("*distance 1 map* ");
                          System.err.print(groupNoObject+"=>");
                          for (j=0;j<neighboringGroupProteins.size();j++)
                            System.err.print(j+":"+neighboringGroupProteins.get(j)+",");
                          System.err.println("");
                        }
    
                        // insert to the i-th tree set; don't skip any proteins
                        if (neighboringGroupProteins.size()>i+1)
                          System.err.println("bugbugbugbugbugbug neighboringGroupProteins.size()>i+1");
                        for (j=neighboringGroupProteins.size();j<=i;j++)
                          neighboringGroupProteins.add(new ArrayList<Integer>());
                        neighboringGroupProteins.get(i).add(intObject);

                        if (DEBUG) {
                          System.err.print("*distance 1 map* ");
                          System.err.print(groupNoObject+"=>");
                          for (j=0;j<neighboringGroupProteins.size();j++)
                            System.err.print(j+":"+neighboringGroupProteins.get(j)+",");
                          System.err.println("");
                        }

                      } //while
                    } //while
                  } //if
                } // for (i=0;i<NUM_DISJOINT_VERTICES;i++)

                groupNoSet=neighboringGroupNoToProteinsMap.keySet();
                if (DEBUG) System.err.println("group no=>neighbor proteins:"+neighboringGroupNoToProteinsMap);

                //====== process neighboringGroupNoToProteinsMap2 ======

                for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                  previousProtein=group.get(i);

                  proteinSet2=new TreeSet<Integer>(proteinToNeighborArray2.get(previousProtein));
                  if (DEBUG) System.err.println("previous protein "+previousProtein+"'s neighbors in distance 2"+proteinSet2);

                  //cannot appear in distance 1
                  proteinSet2.removeAll(proteinSetList.get(i));
                  if (DEBUG) System.err.println("removing distance 1 neighbors=>"+proteinSet2);
                  //cannot appear in current graph
                  proteinSet2.removeAll(previousProteinSet);
                  if (DEBUG) System.err.println("removing neighbors already involved=>"+proteinSet2);

                  // check which group each neighbor belongs to
                  itInteger = proteinSet2.iterator();
                  while (itInteger.hasNext()) {
                    intObject=itInteger.next();
                    //groupNo=proteinToHitGroupNosList[intObject];
                    hitGroupNoList=proteinToHitGroupNosList.get(intObject);

                    litInt=hitGroupNoList.listIterator();
                    while (litInt.hasNext()) {
                      groupNoObject=litInt.next();
                      if (DEBUG) System.err.println("Insert "+intObject+" to group "+groupNoObject);

                      // must be in the neighboringGroupNoToProteinsMap already
                      if (!groupNoSet.contains(groupNoObject)) continue;

                      neighboringGroupProteins2=neighboringGroupNoToProteinsMap2.get(groupNoObject);
                      neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(groupNoObject);
    
                      if (DEBUG) {
                        System.err.print("*distance 2 map* ");
                        System.err.print(groupNoObject+"=>");
                        for (j=0;j<neighboringGroupProteins2.size();j++)
                          System.err.print(j+":"+neighboringGroupProteins2.get(j)+",");
                        System.err.println("");
                      }
    
                      // insert to the i-th tree set
                      if (neighboringGroupProteins2.size()>i+1) {
                        System.err.println("bugbugbugbugbugbug neighboringGroupProteins2.size()>i+1");
                        return;
                      }
                      // fill the space before i
                      if (neighboringGroupProteins2.size()<i && neighboringGroupProteins.size()>=i) {
                        for (j=neighboringGroupProteins2.size();j<i;j++) {
                          if (neighboringGroupProteins.get(j).size()>0)
                            neighboringGroupProteins2.add(new ArrayList<Integer>());
                          else
                            break;
                        }
                      }
                      if (neighboringGroupProteins2.size()==i) {
                        neighboringGroupProteins2.add(neighboringGroupProteinsSubset=new ArrayList<Integer>());
                        neighboringGroupProteinsSubset.add(intObject);
                      } else if (neighboringGroupProteins2.size()==i+1) {
                        neighboringGroupProteins2.get(i).add(intObject);
                      } else {
                        neighboringGroupNoToProteinsMap.remove(groupNoObject);
                        neighboringGroupNoToProteinsMap2.remove(groupNoObject);
                      }
    
                      if (DEBUG) {
                        System.err.print("*distance 2 map* ");
                        System.err.print(groupNoObject+"=>");
                        for (j=0;j<neighboringGroupProteins2.size();j++)
                          System.err.print(j+":"+neighboringGroupProteins2.get(j)+",");
                        System.err.println("");
                      }
    
                    } //while
                  } //while
                } // for (i=0;i<NUM_DISJOINT_VERTICES;i++)
              } // if @ DISTANCE2

              if (DEBUG) System.err.println("group no=>neighbor proteins:"+neighboringGroupNoToProteinsMap);
              if (DEBUG) if (DISTANCE2) System.err.println("group no=>neighbor proteins 2:"+neighboringGroupNoToProteinsMap2);

              //====== eleminate neighboring groups that have occurred =====
              if (!REPEAT_GROUP) {
                for (i=0;i<dg.hitGroupNoList.size();i++) {
                  itInt=dg.hitGroupNoList.get(i).iterator();
                  while (itInt.hasNext()) {
                    neighboringGroupNoToProteinsMap.remove(j=itInt.next().intValue());
                    if (DISTANCE2) neighboringGroupNoToProteinsMap2.remove(j);
                  }
                }
              }

              //====== eleminate neighboring groups =====
              //====== which are not adjacent to every protein ======

              if (!DISTANCE2) {
                itInt=neighboringGroupNoToProteinsMap.keySet().iterator();
                while (itInt.hasNext()) {
                  if ((neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(itInt.next())).size() < NUM_DISJOINT_VERTICES)
                    itInt.remove();
                }
              } else {
                itInt=neighboringGroupNoToProteinsMap.keySet().iterator();
                while (itInt.hasNext()) {
                  groupNoObject=itInt.next();

                  neighboringGroupProteins2=neighboringGroupNoToProteinsMap2.get(groupNoObject);
                  for (j=neighboringGroupProteins2.size();j<NUM_DISJOINT_VERTICES;j++)
                    neighboringGroupProteins2.add(new ArrayList<Integer>());

                  neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(groupNoObject);
                  for (j=neighboringGroupProteins.size();j<NUM_DISJOINT_VERTICES;j++)
                    neighboringGroupProteins.add(new ArrayList<Integer>());

                  for (j=0;j<NUM_DISJOINT_VERTICES;j++) {
                    if (neighboringGroupProteins.get(j).size()+neighboringGroupProteins2.get(j).size()==0)
                      break;
                  }
                  if (j!=NUM_DISJOINT_VERTICES) {
                    itInt.remove();
                    neighboringGroupNoToProteinsMap2.remove(groupNoObject);
                  }
                }
              }

              if (DEBUG) System.err.println("After removing groups, group no=>neighbor proteins="+neighboringGroupNoToProteinsMap);
              if (DEBUG) if (DISTANCE2) System.err.println("After removing groups, group no=>neigbhor proteins 2="+neighboringGroupNoToProteinsMap2);

              //====== eleminate certain neighboring groups or proteins with smaller no=====
              if (neighboringGroupNoToProteinsMap.size()!=0) {
                int nIndex=(int)(dg.hitGroupNoList.size()),nSym;
                ArrayList<Integer> arrayList=rules.get(nIndex); // arrayList[i] < nIndex
                for (i=0;i<arrayList.size();i++) {
                  nSym=arrayList.get(i);
                  //if (DEBUG) System.err.println(nSym+"<"+nIndex);
                  //remove vertices with smaller sequence no. in certain neighboring group
                  int1=Collections.min(dg.alignedProteinListList.get(nSym)).shortValue(); //min element
                  //if (DEBUG) System.err.println("remove vertices (< "+int1+") ");
                  itInt=neighboringGroupNoToProteinsMap.keySet().iterator();
                  while (itInt.hasNext()) {
                    neighboringGroupProteins=neighboringGroupNoToProteinsMap.get(intObj=itInt.next());
                    for (j=0;j<neighboringGroupProteins.size();j++) {
                      litInteger=neighboringGroupProteins.get(j).listIterator();
                      while (litInteger.hasNext()) {
                        intObject=litInteger.next();
                        if (intObject.shortValue()<int1) {
                          litInteger.remove();
                          //if (DEBUG) System.err.print(intObject+" ");
                        }
                      }
                      //if (DEBUG) System.err.println("from group "+dg.hitGroupNoList.get(nSym)+":"+j+":"+neighboringGroupProteins.get(j));
                      if (!DISTANCE2 && neighboringGroupProteins.get(j).size()==0) break;
                    }
                    if (!DISTANCE2 && j!=neighboringGroupProteins.size()) { // neighbors of some protein are empty
                      itInt.remove();
                      //if (DEBUG) System.err.println("remove group "+intObj);
                      continue;
                    }
                    if (DISTANCE2) {
                      neighboringGroupProteins=neighboringGroupNoToProteinsMap2.get(intObj);
                      for (j=0;j<neighboringGroupProteins.size();j++) {
                        litInteger=neighboringGroupProteins.get(j).listIterator();
                        while (litInteger.hasNext()) {
                          intObject=litInteger.next();
                          if (intObject.shortValue()<int1) {
                            litInteger.remove();
                            // if (DEBUG) System.err.print(intObject+",");
                          }
                        }
                        // if (DEBUG) System.err.println(" from group "+dg.hitGroupNoList.get(nSym)+"-2:"+j+":"+neighboringGroupProteins.get(j));
                        if (neighboringGroupNoToProteinsMap.get(intObj).get(j).size()
                            +neighboringGroupProteins.get(j).size()==0) break;
                      }
                      if (j!=neighboringGroupProteins.size()) { // neighbors of some protein are empty
                        itInt.remove();
                        //if (DEBUG) System.err.println("remove group "+intObj);
                        continue;
                      }
                    } // if (DISTANCE2)
                  } // while (itInteger.hasNext())
                } // for (i
              } // if (neighboringGroupNoToProteinsMap.size()!=0


              if (DEBUG) System.err.println("After eliminating smaller proteins, group no=>neighbor proteins="+neighboringGroupNoToProteinsMap);
              if (DEBUG) if (DISTANCE2) System.err.println("After eliminating smaller proteins, group no=>neigbhor proteins 2="+neighboringGroupNoToProteinsMap2);


              if (neighboringGroupNoToProteinsMap.size()!=0)
                goOverNeighboringGroup
                    (queue,
                    neighboringGroupNoToProteinsMap,
                    neighboringGroupNoToProteinsMap2,
                    dg,
                    group,
                    previousDegreeReduced,
                    topology);

              if (!queue.isEmpty()) removeRedundancy(queue);

              //====== eleminate neighboring groups that have occurred =====
              if (!REPEAT_GROUP) {
                if (DEBUG) System.err.println("###remove alignments with overlapping group nos###");

                litDG=queue.listIterator();
                while (litDG.hasNext()) {
                  dg2=litDG.next();
                  if (DEBUG) System.err.println(dg2);
                  if ((rgi=firstLevelRedundancy.get(dg2.alignedProteinListList.get(dg.hitGroupNoList.size())))!=null) {
                    dg2.hitGroupNoList.get(dg2.hitGroupNoList.size()-1).addAll(rgi.hitGroupNos);
                    if (DEBUG) System.err.println("new group no's="+rgi.hitGroupNos);
                    for (i=0;i<dg2.hitGroupNoList.size()-1;i++) {
                      //groupNos=(TreeSet<Integer>)DeepCopy.copy(dg2.hitGroupNoList.get(i));
                      groupNos.clear();
                      TreeSet<Integer> set1=dg2.hitGroupNoList.get(i);
                      for(Integer intObj3: set1) groupNos.add(new Integer(intObj));                      
                      
                      if (DEBUG) System.err.println("group no's at "+i +"="+groupNos);
                      groupNos.retainAll(rgi.hitGroupNos);
                      if (!groupNos.isEmpty()) {
                        if (DEBUG) System.err.println("...overlap");
                        break;
                      }
                    }
                    if (i!=dg.hitGroupNoList.size()) litDG.remove();
                  }
                }
              }

              if (dg.hitGroupNoList.size()==NUM_VERTICES-1) {
              //if (edgeIndex==NUM_EDGES-1) {
                //if (!queue.isEmpty() && EXAMINE_SUBGRAPH) queue=examineSubgraph(queue);

                double queueSize=0;
                while (!queue.isEmpty()) {
                  dg2=queue.poll();

                  if (DISTANCE2) {
                    //if (dg2.distance2) {
                      if (OUTPUT_ALIGNMENT) {
                        bwriter.write(dg2.toString());
                        if (BLAST_SCORE)  bwriter.write("\t");
                      }
                      if (BLAST_SCORE) bwriter.write(calculateBlastScore(dg2).toString());
                      bwriter.newLine();
  /*
                      if (OUTPUT_PROTEIN_SETS) {
                        if ((prefix=dg2.prefix()).equals(pre_prefix)) {
                          for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                            dg2.restProteins(proteinSets.get(i),i); // protein set;
                          }
                        } else {
                          for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                            str=proteinSets.get(i).toString();
                            bwriterProteinSet[i].write(str.substring(1,str.length()-1)); // protein set;
                            bwriterProteinSet[i].newLine();
                            proteinSets.get(i).clear();
                            dg2.proteinSet(proteinSets.get(i),i);
                          }
                          pre_prefix=prefix;
                        }
                      }
                      queueSize++;
  */
                    //}
                  } else {
                    if (OUTPUT_ALIGNMENT) {
                      bwriter.write(dg2.toString());
                      if (BLAST_SCORE)  bwriter.write("\t");
                    }
                    if (BLAST_SCORE) bwriter.write(calculateBlastScore(dg2).toString());
                    bwriter.newLine();
  /*
                    if (OUTPUT_PROTEIN_SETS) {
                      if ((prefix=dg2.prefix()).equals(pre_prefix)) {
                        for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                          dg2.restProteins(proteinSets.get(i),i); // protein set;
                        }
                      } else {
                        for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                          str=proteinSets.get(i).toString();
                          bwriterProteinSet[i].write(str.substring(1,str.length()-1)); // protein set;
                          bwriterProteinSet[i].newLine();
                          proteinSets.get(i).clear();
                          dg2.proteinSet(proteinSets.get(i),i);
                        }
                        pre_prefix=prefix;
                      }
                    }

                    queueSize++;
  */
                  }
                }
  /*
                if (queueSize!=0) {
                  DecimalFormat df = new DecimalFormat("#");
                  System.out.println(df.format(queueSize));
                }
  */
              } else
                stack.push(queue);

            while (!stack.empty() && stack.peek().isEmpty())
              stack.pop();

          } //while (!stack.empty())
        } // while (!lastCombination)
      }// for g
      System.err.println("totalNumPairs="+totalNumPairs);
      bwriter.close();
      if (file.length()==0) {
        System.err.println(file.getName()+" is empty");
        file.delete();
      }
/*
      if (OUTPUT_PROTEIN_SETS)
        for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
          bwriterProteinSet[i].newLine();
          bwriterProteinSet[i].close(); // protein set;
        }
*/

      Iterator<ArrayList<Integer>> itArrayInteger=firstLevelRedundancy.keySet().iterator();
      while (itArrayInteger.hasNext())
        firstLevelRedundancy.get(itArrayInteger.next()).occurs=false;

      } //adj
      breaderAdj.close();
      breaderAuto.close();

      } //for speciesNo
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
/*
  private void generateCombinations(int[][] proteinNoCombinationArray, //will change values
                                    int m,int n) {
    int i,j,k;
    for (i=0;i<n;i++)
      proteinNoCombinationArray[0][i]=i;

    for (i=1;i<proteinNoCombinationArray.length;i++) {
      //copy values in the previous combination
      for (j=0;j<n;j++)
        proteinNoCombinationArray[i][j]=proteinNoCombinationArray[i-1][j];

      k=n-1;
      while (k>=0 && proteinNoCombinationArray[i-1][k]==m-n+k)
        k--;
      proteinNoCombinationArray[i][k]=proteinNoCombinationArray[i][k]+1;
      for (j=k+1;j<n;j++)
        proteinNoCombinationArray[i][j]=proteinNoCombinationArray[i][k]+j-k;
    }
  }
*/
  private void goOverNeighboringGroup
      (LinkedList<AlignedGraphlets> queue, //will change values
       LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>> neighboringGroupNoToProteinsMap,
       LinkedHashMap<Integer,ArrayList<ArrayList<Integer>>> neighboringGroupNoToProteinsMap2,
       AlignedGraphlets dg,
       ArrayList<Integer> group,
       int previousDegreeReduced, ArrayList<ArrayList<Integer>> topology)
  {

    Iterator<Integer> itInteger,itInteger2,itInteger3;
    ListIterator<Integer> litInt;

    ArrayList<Integer> neighboringGroupNos;
    int groupNo;
    Integer intObject;
    Integer groupNoObject;
    ArrayList<ArrayList<Integer>> neighboringGroupProteins=new ArrayList<ArrayList<Integer>>(), neighboringGroupProteins2=new ArrayList<ArrayList<Integer>>();
    ListIterator<ArrayList<Integer>> itList;

    TreeSet<Integer> groupNoTreeSet;
    TreeSet<Integer> proteinSet;
    ArrayList<Integer> newAlignedProteins=new ArrayList<Integer>(NUM_DISJOINT_VERTICES);
    Integer[][] neighboringProteinsMatrix;
    int i,j,k;
    int[] indices=new int[NUM_DISJOINT_VERTICES];
    boolean overlap,someHasNoNeighbor,allDistance2;
    AlignedGraphlets dg2;

    int numPermutations,size,degree=0;
    int distance2count;
    int advance_position,pre_advance_position;

    ArrayList<Integer> oldAlignedProteins,topologyNeighbors;
    int int1,int2;
    boolean addNewProtein;
    
    //go over the neighboring groups
    //NOTE: dg, group remain constant in this loop
    neighboringGroupNos=new ArrayList<Integer>(neighboringGroupNoToProteinsMap.keySet());
    Collections.sort(neighboringGroupNos);
   
    litInt=neighboringGroupNos.listIterator();
    //if (DEBUG) System.err.println("neighboring group="+neighboringGroupNos.toString());
    while (litInt.hasNext()) {
      groupNo=(groupNoObject=litInt.next()).intValue();
      if (DEBUG) System.err.println("** neighboring group="+groupNo);
      //neighboringGroupProteins=(ArrayList<ArrayList<Integer>>)DeepCopy.copy(neighboringGroupNoToProteinsMap.get(groupNoObject));
      neighboringGroupProteins.clear();
      itList=neighboringGroupNoToProteinsMap.get(groupNoObject).listIterator();
      while (itList.hasNext()) {
        ArrayList<Integer> list1=itList.next();
        ArrayList<Integer> list2=new ArrayList<Integer>(list1.size());
        for(Integer intObj: list1) list2.add(new Integer(intObj));
        neighboringGroupProteins.add(list2);
      }
        
      if (DEBUG) System.err.println("neighboringGroupProteins="+neighboringGroupProteins);
      if (DISTANCE2) {
        //neighboringGroupProteins2=(ArrayList<ArrayList<Integer>>)DeepCopy.copy(neighboringGroupNoToProteinsMap2.get(groupNoObject));
        neighboringGroupProteins2.clear();
        itList=neighboringGroupNoToProteinsMap2.get(groupNoObject).listIterator();
        while (itList.hasNext()) {
          ArrayList<Integer> list1=itList.next();
          ArrayList<Integer> list2=new ArrayList<Integer>(list1.size());
          for(Integer intObj: list1) list2.add(new Integer(intObj));        
          neighboringGroupProteins2.add(list2);      
        }
        if (DEBUG) System.err.println("neighboringGroupProteins2="+neighboringGroupProteins2);
      }

      numPermutations=1;
      if (!DISTANCE2)
        for (i=0;i<NUM_DISJOINT_VERTICES;i++)
          numPermutations*=neighboringGroupProteins.get(i).size();
      else
        for (i=0;i<NUM_DISJOINT_VERTICES;i++)
          numPermutations*=neighboringGroupProteins.get(i).size()+neighboringGroupProteins2.get(i).size();
      //if (DEBUG) System.err.println("# permutations="+numPermutations);

      //reset the indices
      for (i=0;i<NUM_DISJOINT_VERTICES;i++)
        indices[i]=0;
/*
      if (DEBUG) {
        System.err.print("index ");
        for (i=0;i<NUM_DISJOINT_VERTICES;i++)
          System.err.print(indices[i]+",");
        System.err.println("");
      }
*/
      //go over the permuation of neighbors of the proteins in the neighboring group
      //if a permutation of neighbors satisfy the vertex disjoint contraint,
      //add the combination of the proteins in the neighboring group to the queue
      for (k=0;k<numPermutations;k++) {

        //@ DISTANCE2
        newAlignedProteins.clear();
        if (DISTANCE2) {
          if (DISTANCE2_COUNT) {
            distance2count=0;
            for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
              if (indices[i]>=(size=neighboringGroupProteins.get(i).size())) {
                newAlignedProteins.add(neighboringGroupProteins2.get(i).get(indices[i]-size));
                distance2count++;
              } else
                newAlignedProteins.add(neighboringGroupProteins.get(i).get(indices[i]));
            }
            if (distance2count > DISTANCE2_NUM_THRESHOLD) {
              //if (DEBUG) System.err.println("too many distance 2");
              continue;
            }
          } else {
            //concatenates distance 1 & distance 2 neighbor lists
            //indices[i] indicates the border between the two
            allDistance2=true;
            for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
              if (indices[i]>=(size=neighboringGroupProteins.get(i).size())) {
                newAlignedProteins.add(neighboringGroupProteins2.get(i).get(indices[i]-size));
              } else {
                newAlignedProteins.add(neighboringGroupProteins.get(i).get(indices[i]));
                allDistance2=false;
              }
            }
            if (allDistance2) {
              //if (DEBUG) System.err.println("All distance 2");
              break; //all the following permutations are all distance-2 combinations
            }
          }
        } else {
          for (i=0;i<NUM_DISJOINT_VERTICES;i++)
            newAlignedProteins.add(neighboringGroupProteins.get(i).get(indices[i]));
        }
/*
        if (DEBUG) {
          System.err.print("permutation "+k+":");
          for (i=0;i<NUM_DISJOINT_VERTICES;i++)
            System.err.print(newAlignedProteins.get(i)+",");
          System.err.println("");
        }
*/
        //check if this set of neighbors overlap
        overlap=false;
        for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
          for (j=0;j<i;j++) {
            if (newAlignedProteins.get(i).equals(newAlignedProteins.get(j))) {
              overlap=true;
              break;
            }
          }
          if (overlap) break;
        }

        if (!overlap) {

            // check induced connections
            addNewProtein=true;
            topologyNeighbors=topology.get(int1=(int)dg.hitGroupNoList.size());
            if (!DISTANCE2) {
              for (int2=0;int2<int1;int2++) {
                oldAlignedProteins=dg.alignedProteinListList.get(int2);
                if (topologyNeighbors.contains(int2)) { //connected

                  for (i=0;i<NUM_DISJOINT_VERTICES;i++)
                    if (!new TreeSet<Integer>(proteinToNeighborArray.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i)))
                      break;
                  if (i!=NUM_DISJOINT_VERTICES) {
                    break;
                  }
                } else if (INDUCED_SUBGRAPH) { //not connected
                  for (i=0;i<NUM_DISJOINT_VERTICES;i++)
                    if (!new TreeSet<Integer>(proteinToNeighborArray.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i)))
                      break;
                  if (i==NUM_DISJOINT_VERTICES) {
                    break;
                  }
                } 
              } //for (int2=0;int2<int1-1;int2++)
              if (int2!=int1) addNewProtein=false;

            } else { // if DISTANCE 2

              for (int2=0;int2<int1;int2++) {
                oldAlignedProteins=dg.alignedProteinListList.get(int2);
                //if (DEBUG) System.err.println(dg.hitGroupNoList.get(int2)+":"+oldAlignedProteins);
                if (topologyNeighbors.contains(int2)) { //connected
                  //if (DEBUG) System.err.println("should be connected");
                  j=0;
                  for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                    if (!new TreeSet<Integer>(proteinToNeighborArray.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i))) {
                      if (!new TreeSet<Integer>(proteinToNeighborArray2.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i))) {
                        //if (DEBUG) System.err.println(newAlignedProteins.get(i)+" and "+oldAlignedProteins.get(i)+" are not neighbors");
                        break;
                      } else {
                        j++;
                        //if (DEBUG) System.err.println(newAlignedProteins.get(i)+" has indirect neighbor "+oldAlignedProteins.get(i));
                      }
                    } //else if (DEBUG) System.err.println(newAlignedProteins.get(i)+" has direct neighbor "+oldAlignedProteins.get(i));
                  }
                  if (i!=NUM_DISJOINT_VERTICES || j==NUM_DISJOINT_VERTICES) {
                    break;
                  }
                } else if (INDUCED_SUBGRAPH) { //not connected
                  //if (DEBUG) System.err.println("should not be connected");
                  j=0;
                  for (i=0;i<NUM_DISJOINT_VERTICES;i++) {
                    if (!new TreeSet<Integer>(proteinToNeighborArray.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i))) {
                      if (!new TreeSet<Integer>(proteinToNeighborArray2.get(newAlignedProteins.get(i))).contains(oldAlignedProteins.get(i))) {
                        //if (DEBUG) System.err.println(newAlignedProteins.get(i)+" and "+oldAlignedProteins.get(i)+" are not neighbors");
                        break;
                      } else {
                        j++;
                        //if (DEBUG) System.err.println(newAlignedProteins.get(i)+" has indirect neighbor "+oldAlignedProteins.get(i));
                      }
                    } //else if (DEBUG) System.err.println(newAlignedProteins.get(i)+" has direct neighbor "+oldAlignedProteins.get(i));
                  }
                  if (i==NUM_DISJOINT_VERTICES && j<NUM_DISJOINT_VERTICES) {
                    break;
                  }
                } 
              } //for
              if (int2!=int1) {
                //if (DEBUG) System.err.println("give up this branch");
                addNewProtein=false;
              } 
            } // // if DISTANCE 2

          if (addNewProtein) {

            //modify AlignedGraphlets
            //dg2 = (AlignedGraphlets) DeepCopy.copy(dg);
            dg2 = new AlignedGraphlets(dg);
            
            //if (DISTANCE2) dg2.distance2=true;
            //dg2.numEdges++;
            
            groupNoTreeSet=new TreeSet<Integer>();
            groupNoTreeSet.add(groupNo);
            dg2.hitGroupNoList.add(groupNoTreeSet);
            dg2.alignedProteinListList.add(new ArrayList<Integer>(newAlignedProteins));
            overlap = false;
/*
//incorrect-begin
            ListIterator<AlignedGraphlets> itDG=queue.listIterator();
            while (itDG.hasNext()) {
              if (dg2.equals(itDG.next())) {
                overlap = true;
                break;
              }
            }
            if (!overlap) {
//incorrect-end
*/
              //if (DEBUG) System.err.println("add "+dg2);
              queue.add(dg2);
/*
//incorrect-begin
            } else {
              //if (DEBUG) System.err.println("has identical graph, not adding "+dg2);
            }
//incorrect-end
*/
          }
/*
          } else {
            if (DEBUG) System.err.println("low e value, drop "+newAlignedProteins);
          }
*/
        } else { //if overlap
          //if (DEBUG) System.err.println("...overlap");
        }

        //increment indices
        advance_position=0;
        for (i=NUM_DISJOINT_VERTICES-1;i>=0;i--) {
          indices[i]++;
          if (DISTANCE2) {
            pre_advance_position=advance_position;
            if (indices[i]==(neighboringGroupProteins.get(i).size()+neighboringGroupProteins2.get(i).size()))
              advance_position=2;
            else if (indices[i]==neighboringGroupProteins.get(i).size())
              advance_position=1;
            else
              break;
            if (indices[i]==neighboringGroupProteins.get(i).size()) {
              if (pre_advance_position==2) {
                indices[i+1]=0;
                if (indices[i]==(neighboringGroupProteins.get(i).size()+neighboringGroupProteins2.get(i).size()))
                  indices[i]=0;
                else
                  break;
              } else if (pre_advance_position==1) {
                indices[i]=0;
                indices[i+1]=neighboringGroupProteins.get(i+1).size();
                break;
              } else {
                indices[i]=0;
              }
            } else if (indices[i]==(neighboringGroupProteins.get(i).size()+neighboringGroupProteins2.get(i).size())) {
              if (pre_advance_position==2) {
                indices[i]=0;
                indices[i+1]=0;
              } else if (pre_advance_position==1) {
                indices[i]=neighboringGroupProteins.get(i).size();
                indices[i+1]=neighboringGroupProteins.get(i+1).size();
                break;
              } else {
                indices[i]=neighboringGroupProteins.get(i).size();
              }
            }
          } else if (indices[i]==neighboringGroupProteins.get(i).size())
            indices[i]=0;
          else
            break;
        }
/*
        if (DEBUG) {
          System.err.print("index ");
          for (i=0;i<NUM_DISJOINT_VERTICES;i++)
            System.err.print(indices[i]+",");
          System.err.println("");
        }
*/
      } //for k - go over the permutation of neighbors of the size-combination of proteins in the last group
    } //itInteger - go over the neighboring groups
  }

  public static void main(String[] args) {
    GraphletAlign ga = new GraphletAlign(args);
    ga.branchNBound();
  }
}

  // NUM_LINES is defined when EXAMINE_SUBGRAPH = true
  // for allocating space in advance

  /*
  //biogrid Dmela.2
  final static int[] NUM_LINES = {30, 1, 5, 5, 1,
                 4, 1, 0, 40, 40,
                 5, 5, 10, 4, 1,
                 6, 3, 1, 3, 1,
                 1, 0, 0, 1, 0,
                0, 0, 0, 0};
                */
/*
  //biogrid Scere.2.9
  final static int[] NUM_LINES = {4225, 24, 48869, 259, 5736,
                 15892, 7146, 1496, 1000, 81117,
                 0, 10350, 35310, 12861, 215,
                 2658, 39248, 1977, 7917, 114,
                 1500, 2981, 10113, 4788, 7670,
                 6830, 921, 26616, 2200};
*/
/*
  //dip Scere.2
  final static int[] NUM_LINES = {86, 4, 76, 71, 4,
    13, 1, 0, 68, 180,
    0, 14, 1, 14, 0,
    17, 4, 0, 2, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0};
*/
    /*
  //dip Dmela.2
  final static int[] NUM_LINES = {45, 3, 42, 0, 2,
    2, 2, 0, 49, 49,
    9, 6, 12, 8, 1,
    6, 5, 1, 3, 1,
    1, 0, 0, 1, 0,
    0, 0, 0, 0};
  */