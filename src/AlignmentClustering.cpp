#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include "Utility.h"
using namespace std;
//AlignmentClustering [protein_names] [output directory] [topology dir] [SINGLE_SPECIES] [threshold] [max proteins] [alignment files]
//
//Clustering on protein set combination of all species
//**with preclustering**
//output protein no
bool SINGLE_SPECIES=true;
double THRESHOLD=0.50;
int MAX_PROTEIN=100;
//const int MIN_PROTEIN=10;
const bool AVERAGE_SCORE=true;
const int NUMBER_CLUSTERS_PER_PROTEIN=1;
const int STARTING_INDEX_ALIGNMENT_FILE=7;
//1: for each species, overlapping part of larger module < THRESHOLD, module size < MAX_PROTEIN
//2: for each species, overlapping part of average module size < THRESHOLD, module size < MAX_PROTEIN
//3: for each species, overlapping part of smaller module < THRESHOLD, module size < MAX_PROTEIN
//4: for all species together, overlapping part of average module size < THRESHOLD, module size < MAX_PROTEIN
int num_copies;

vector<struct alignment> alignments; //data structure
vector< vector<string> > protein_names; //data structure

//sort by e-values
struct alignment {
  long e;
  vector<vector<int> > proteins; //sorted by protein number when inputting
  int num_edges;
  alignment() {num_edges=0;}  
  void merge(const struct alignment& a) {
    vector<vector<int> > tmp_proteins = proteins;
    
    for (int i=0;i<num_copies;i++) {
      vector<int> v;
      set_union(proteins[i].begin(),proteins[i].end(),a.proteins[i].begin(),a.proteins[i].end(),insert_iterator< vector<int> >(v,v.begin()));
      proteins[i]=v;
    }

    if (a.e < e) e=a.e; //stores best e
      
    for (int i=0;i<proteins.size();i++)
      if (proteins[i].size() > MAX_PROTEIN) {
        proteins = tmp_proteins;
        return;
      }
    
  }
  // return false if num_copies!=2 or reaching the max size
  bool reverse_merge(const struct alignment& a) {
    if (num_copies!=2) {
      cerr << "reverse_merge:num_copies!=2" << endl;
      return false;
    }
    
    vector<vector<int> > tmp_proteins = proteins;
    
    for (int i=0;i<num_copies;i++) {
      vector<int> v;
      set_union(proteins[i].begin(),proteins[i].end(),a.proteins[1-i].begin(),a.proteins[1-i].end(),insert_iterator< vector<int> >(v,v.begin()));
      proteins[i]=v;
    }

    if (a.e < e) e=a.e; //stores best e
    for (int i=0;i<proteins.size();i++)
      if (proteins[i].size() > MAX_PROTEIN) {
        proteins = tmp_proteins;
        return false;
      }
    return true;
  }    
  bool empty() const {
    return proteins.empty();
  }
  void clear() {
    proteins.clear();
  }
  bool operator==(const struct alignment& a) const {
    int i;
    for (i=0;i<num_copies;i++) {
      if (proteins[i]!=a.proteins[i]) break;
    }
    if (i!=num_copies)
      return false;
    return true;
  }
  
  // protein #, edge #, e-value
  // protein # (fewer first), edge # (fewer first), e-value (less first)
  // protein # (fewer first), edge # (fewer first), e-value (more first)
  bool operator<(const struct alignment& a) const
  {
    //cerr << "compare to " << a << endl;
    int i,size=0,size2=0;
    for (i=0;i<num_copies;i++) {
      size+=proteins[i].size();
      size2+=a.proteins[i].size();
    }
    if (size==size2) {
      if (num_edges == a.num_edges) {
        if (e == a.e) {
          for (i=0;i<num_copies;i++) {
            if (proteins[i]!=a.proteins[i])
              return proteins[i]<a.proteins[i];
          }
          return false;
        }
        return e < a.e;
      }
      return num_edges > a.num_edges;
    }
    return size > size2;
  }
  friend std::ostream& operator<<(std::ostream& out, const struct alignment& a)
  {
    int i;
    for (i=0;i<num_copies;i++) {
      vector<int>::const_iterator it=a.proteins[i].begin();
      out << protein_names[*it][0];
      for (++it;it!=a.proteins[i].end();++it)
        out << "," << protein_names[*it][0];
      //copy(a.proteins[i].begin(),a.proteins[i].end(),ostream_iterator<int,char>(out,","));
      out << "\t";
    }
    /*
    for (i=0;i<num_copies;i++) {
      out << a.proteins[i].size() << " ";
    }
    */
    out << a.num_edges;
    out << "\t" << a.e;
/*
    for (i=0;i<a.files.size();i++)
      out << a.files[i] << " " << a.lines[i] << " " << a.es[i] << "|";
*/
    return out;
  }
};
//sort by module sizes
struct lt_cluster
{
  bool operator()(const alignment& a1, const alignment& a2) const
  {
    int size1=0,size2=0;
    for (int i=0;i<num_copies;i++) {
      size1+=a1.proteins[i].size();
      size2+=a2.proteins[i].size();
    }
    if (size1==size2)
      return a1 < a2;
    return size1 > size2;
  }
};
struct protein_to_alignments {
  int p;
  vector<int> as;
  protein_to_alignments() { p = -1; }
  bool operator<(const struct protein_to_alignments& pa) const {
    return as.size() > pa.as.size();
  }
  //overload << for printing VertexDisjointGraph object
  friend std::ostream& operator<<(std::ostream& out, const struct protein_to_alignments& pa)
  {
    out << pa.p << ":";
    copy(pa.as.begin(),pa.as.end(),ostream_iterator<int,char>(out, ", "));
    return out;
  }
};
bool all_overlap;
bool overlap2(const struct alignment& a1, const struct alignment& a2) {
  //go over set pairs of each species
  int i;
  double average_size;
  double p1,p2;
  all_overlap=true;
  for (i=0;i<num_copies;i++) {
    set<int> s;
    average_size=(a1.proteins[i].size()+a2.proteins[i].size())/(double)2;
    set_intersection(a1.proteins[i].begin(),a1.proteins[i].end(),a2.proteins[i].begin(),a2.proteins[i].end(),insert_iterator< set<int> >(s,s.begin()));
    if (s.size()/(double)average_size < THRESHOLD) break;
    if (s.size()/(double)a1.proteins[i].size()!=1 && s.size()/(double)a2.proteins[i].size()!=1) all_overlap=false;
  }
  if (i!=num_copies) return false;
  /*
  for (i=0;i<num_copies;i++)
    if (a1.proteins[i].size()>=MAX_PROTEIN || a2.proteins[i].size()>=MAX_PROTEIN) return false;
    */
  return true;
}
bool overlap2_reverse(const struct alignment& a1, const struct alignment& a2) {
  if (num_copies!=2) {
    cerr << "overlap2_reverse:num_copies!=2" << endl;
    return false;
  }
  //go over set pairs of each copy
  int i;
  double average_size;
  double p1,p2;
  all_overlap=true;
  for (i=0;i<num_copies;i++) {
    set<int> s;
    average_size=(a1.proteins[i].size()+a2.proteins[1-i].size())/(double)2;
    set_intersection(a1.proteins[i].begin(),a1.proteins[i].end(),a2.proteins[1-i].begin(),a2.proteins[1-i].end(),insert_iterator< set<int> >(s,s.begin()));
    if (s.size()/(double)average_size < THRESHOLD) break;
    if (s.size()/(double)a1.proteins[i].size()!=1 && s.size()/(double)a2.proteins[1-i].size()!=1) all_overlap=false;
  }
  if (i!=num_copies) return false;
  /*
  for (i=0;i<num_copies;i++)
    if (a1.proteins[i].size()>=MAX_PROTEIN || a2.proteins[i].size()>=MAX_PROTEIN) return false;
    */
  return true;    
}
int main(int argc, char * argv[]) {
  // determine num_copies

  ifstream fin;
  string str;
  fin.open(argv[STARTING_INDEX_ALIGNMENT_FILE]);
  if (!fin) {
    cerr << "files does not exist " << argv[STARTING_INDEX_ALIGNMENT_FILE] << endl;
    fin.clear();
    return 1;
  }
  getline(fin, str);
  stringtokenizer st(str,"[]-\t");
  stringtokenizer::const_iterator stit = st.begin();
  stit++; //num_vertices++;
  num_copies=numtok(*stit," ,");

  fin.close();
  fin.clear();
  //cerr << "num_copies=" << num_copies << endl;
  //cerr << "num_vertices=" << num_vertices << endl;

  
  SINGLE_SPECIES=(string(argv[4])=="true"?true:false);
  THRESHOLD=atof(argv[5]);
  MAX_PROTEIN=atoi(argv[6]);
  
  // ====== read species name ======

  ofstream fout;
  string file(argv[1]);
  int i;
  while (file.at(file.length()-1)=='/')
    file=file.substr(0,file.length()-1);
  file=file.substr(i=file.find_last_of("/")+1,file.length()-i);
  string species_name[num_copies];
  int j,k=-1;
  for (i=0;i<num_copies;i++) {
    j=k+1;
    k=file.find_first_of("_",j);
    species_name[i]=file.substr(j,k-j);
  }

  // ====== read protein names ======

  fin.open(argv[1]);
  if (!fin) {
    cerr << "file does not exist " << argv[1] << endl;
    fin.clear();
    return -1;
  }
  cerr << "open " << argv[1] << endl;

  protein_names.push_back(vector<string>());
  while (getline(fin, str)) {
    vector<string> v;
    if (!str.empty()) {
      stringtok<string>(v,str,";");
      if (v[v.size()-1]=="")
        v.erase(v.end()-1);
    }
    protein_names.push_back(v);
  }
  fin.close();
  fin.clear();

  vector<struct protein_to_alignments> proteins_to_alignments(protein_names.size());
  cerr << "read " << proteins_to_alignments.size() << " proteins" << endl << endl;
  
  // read alignments
  int line;
  int l;
  string pre_prefix,prefix;
  vector< set<int> > protein_sets(num_copies);
  long best_e=180,tmp_e;
  long double d;
  string e;
  int v1, v2;
  vector< set<int> > last_proteins(num_copies);
  int s;

  // read alignments
  for (i=STARTING_INDEX_ALIGNMENT_FILE;i<argc;i++) {

    best_e=180;

    if (byte_size_of_file(argv[i]) == 0) continue;

    file = string(argv[i]);
    k=file.find_last_of(".");
    file=file.substr(k+1,file.length()-k-1);
    
    int num_vertices;
    stot<int>(num_vertices,file.substr(0,k=file.find_first_of("-")));
    int num_edges;
    stot<int>(num_edges,file.substr(k+1,file.find_first_of("-",k+1)-k-1));

    //cerr << "num_vertices=" << num_vertices << endl;
    //cerr << "num_edges=" << num_edges << endl;

    std::stringstream ss;
    ss << argv[3] << "/" << num_vertices << ".adj";
    fin.open(ss.str().c_str());
    if (!fin) {
      cerr << "files does not exist " << ss.str() << endl;
      fin.clear();
      return -1;
    }

    //cerr << "open " << ss.str() << endl;

    int num_edges_last_vertex = 0;
    while (getline(fin, str)) {
      if (str==file) {
        //cerr << str << endl;
        for (j=0;j<num_edges;j++) {
          getline(fin,str);
          //cerr << str << endl;
          stot<int>(v1,str.substr(0,k=str.find_first_of(" ")));
          stot<int>(v2,str.substr(k+1,str.find_first_of(" ",k+1)-k-1));
            //cerr << v1 << " " << v2 << endl;
          if (v1==num_vertices-1 || v2==num_vertices-1)
            num_edges_last_vertex++;
        }
        break;
      }
    }

    fin.close();
    fin.clear();
    //cerr << "num_edges_last_vertex=" << num_edges_last_vertex << endl;
    
    fin.open(argv[i]);
    if (!fin) {
      cerr << "files does not exist " << argv[i] << endl;
      fin.clear();
      continue;
    }
    
    cerr << "open " << argv[i] << endl;

    getline(fin, str);
    str=str.substr(0,str.find_last_of("\t"));
    pre_prefix=str.substr(0,str.find_last_of("-",str.length()-2));
    //cerr << "pre_prefix=" << pre_prefix << endl;

    /*
    int num_vertices2=0;

    stringtokenizer st(str,"[]-\t");
    stringtokenizer::const_iterator stit = st.begin();
    //cerr << *stit << endl;
    stit++; num_vertices2++;
    for (; stit != st.end(); stit++) {
      //cerr << *stit << endl;
      num_vertices2++;
    }
    num_vertices2/=2;

    cerr << "num_vertices2=" << num_vertices2 << endl;

    if (num_vertices2!=num_vertices) {
      cerr << "#vertices wrong" << endl;
      return -1;
    }
    */
    fin.close();
    fin.clear();
    
    
    fin.open(argv[i]);

    //line=1;
    while (getline(fin, str)) {
      l=str.find_last_of("\t");
      e=str.substr(l+1,str.length()-l-1);
      //cerr << "e=" << e << endl;
      str=str.substr(0,l);
      prefix=str.substr(0,str.find_last_of("-",str.length()-2));
      //cerr << "prefix=" << prefix << endl;

      if (prefix!=pre_prefix) {
        //cerr << prefix << endl;
  
        pre_prefix=prefix;
        struct alignment a;
        if (AVERAGE_SCORE) a.e=best_e/num_vertices;
        else a.e=best_e;
        //cerr << "a.e=" << a.e << endl;
        a.num_edges+=num_edges * num_copies;
        //cerr << "precluster" << endl;
        for (k=0;k<num_copies;k++) {
          a.proteins.push_back(vector<int>(protein_sets[k].begin(),protein_sets[k].end()));
          //copy(protein_sets[k].begin(),protein_sets[k].end(),ostream_iterator<int,char>(cerr," "));
          //cerr << endl;
    
          //build index
          vector<int>& vs=a.proteins[k];
          //copy(vs.begin(),vs.end(),ostream_iterator<int,char>(cerr," "));
          //cerr << endl;    
          for (j=0;j<vs.size();j++) {
            proteins_to_alignments[vs[j]].p=vs[j];
            proteins_to_alignments[vs[j]].as.push_back(alignments.size());
            //cerr << proteins_to_alignments[vs[j]] << endl;
          }
          protein_sets[k].clear();

          //copy(last_proteins[k].begin(),last_proteins[k].end(),ostream_iterator<int,char>(cerr," "));
          //cerr << endl;
    
          a.num_edges+=(last_proteins[k].size()-1) * num_edges_last_vertex;
          last_proteins[k].clear();
        }
        //cerr << a << endl;
        //alignments.insert(a);
        alignments.push_back(a);
        //cerr << alignments[alignments.size()-1] << endl;
        best_e=180;
      }

      stringtokenizer st(str,"[]-");
      stringtokenizer::const_iterator stit = st.begin();
      for (j=0;j<num_vertices;j++,stit++) {
        //cerr << (*stit) << endl;
        stit++;
        //cerr << (*stit) << endl;
        stringtokenizer st2(*stit,", ");
        stringtokenizer::const_iterator stit2 = st2.begin();
        for (k=0;k<num_copies;k++,stit2++) {
          protein_sets[k].insert(s=(int)atoi((*stit2).c_str()));
          if (j==num_vertices-1)
            last_proteins[k].insert(s);
        }
      }
/*
      cerr << "line" << endl;
      for (k=0;k<num_copies;k++) {
        copy(protein_sets[k].begin(),protein_sets[k].end(),ostream_iterator<int,char>(cerr," "));
        cerr << endl;
      }
*/

      if ((l=e.find_last_of("E"))!=-1)
        stot<long>(tmp_e,e.substr(l+1,e.length()-l-1));
      else {
        stot<long double>(d,e);
        if (d == 0) tmp_e = -180;
        else tmp_e=(long)log10(d);
      }
      if (tmp_e<best_e) best_e=tmp_e;
      //line++;
    }
    struct alignment a;
    if (AVERAGE_SCORE) a.e=best_e/num_vertices;
    else a.e=best_e;
    //cerr << "a.e=" << a.e << endl;
    a.num_edges+=num_edges * num_copies;
    //cerr << "precluster" << endl;
    for (k=0;k<num_copies;k++) {
      a.proteins.push_back(vector<int>(protein_sets[k].begin(),protein_sets[k].end()));
      //copy(protein_sets[k].begin(),protein_sets[k].end(),ostream_iterator<int,char>(cerr," "));
      //cerr << endl;
      
      //build index
      vector<int>& vs=a.proteins[k];
      //copy(vs.begin(),vs.end(),ostream_iterator<int,char>(cerr," "));
      //cerr << endl;
      for (j=0;j<vs.size();j++) {
        proteins_to_alignments[vs[j]].p=vs[j];
        proteins_to_alignments[vs[j]].as.push_back(alignments.size());
        //cerr << proteins_to_alignments[vs[j]] << endl;
      }
      protein_sets[k].clear();

      //copy(last_proteins[k].begin(),last_proteins[k].end(),ostream_iterator<int,char>(cerr," "));
      //cerr << endl;
      
      a.num_edges+=(last_proteins[k].size()-1) * num_edges_last_vertex;
      last_proteins[k].clear();
    }
    //cerr << a << endl;
    //alignments.insert(a);
    alignments.push_back(a);
    best_e=180;

    fin.close();
    fin.clear();
  }
  /*
  for (vector<struct alignment>::iterator ait=alignments.begin();
       ait!=alignments.end();ait++) {
    cerr << *ait << endl;
  }
  */
  cerr << "read " << alignments.size() << " alignments" << endl << endl;
  //cerr << "size of index=" << proteins_to_alignments.size() << endl;
  /*
  i = 0;
  for (vector<protein_to_alignments >::iterator vit=proteins_to_alignments.begin();
       vit!=proteins_to_alignments.end();vit++) {
    cerr << i << ":" << *vit << endl;
    i++;
  }
  */
  sort(proteins_to_alignments.begin(),proteins_to_alignments.end());
    
  map<int, int> proteins_to_index;
  i = 0;
  for (vector<protein_to_alignments >::iterator vit=proteins_to_alignments.begin();
       vit!=proteins_to_alignments.end();vit++) {
    //cerr << i << ":" << endl;//<< *vit << endl;
    proteins_to_index[vit->p] = i;
    /*
    for (vector<int>::iterator vit2=vit->as.begin(); vit2!=vit->as.end(); vit2++) {
      cerr << alignments[*vit2] << endl;
    }
    */
    i++;
  }
  
  //============clustering===============

  set<struct alignment,lt_cluster> clusters;
  set<struct alignment> tmp_alignments;
  int number_clusters_per_protein;
  bool protein_removed;
  for (vector<protein_to_alignments >::iterator vit=proteins_to_alignments.begin();
       vit!=proteins_to_alignments.end();vit++) {
    if (vit->as.size() > 0 && vit->p != -1) {
      //cerr << "protein " << vit->p << endl;

      // create tmp list
      for (vector<int>::iterator vit2=vit->as.begin();vit2!=vit->as.end();vit2++) {
        struct alignment& a = alignments[*vit2];
        tmp_alignments.insert(a);
        //cerr << alignments[*vit2] << endl;
      }
      // clear the protein
      vit->as.clear();
      /*
      int size=tmp_alignments.size();
      float percentage=0.9;
      */
      //cerr << tmp_alignments.size() << " tmp_alignments" << endl;
      number_clusters_per_protein=0;
      for (set<struct alignment>::iterator sit=tmp_alignments.begin();sit!=tmp_alignments.end();) {
  
        struct alignment a;
        a=*sit;
        /*
        if (tmp_alignments.size()/(long double)size < percentage) {
          cerr << "compare to " << tmp_alignments.size() << " tmp_alignments" << endl;
          percentage-=0.1;
        }
        */
        //cerr << "prepare to merge to " << a << endl;

        set<struct alignment>::iterator sit2=sit; sit2++;
        for (;sit2!=tmp_alignments.end();) {
          //cerr << "compare to " << *sit2;
          if (overlap2(a,*sit2)) {
            //cerr << "...overlaps => ";
            a.merge(*sit2);
            //cerr << a;
            //saved_alignments.insert(*sit2);    
            //sit2++;
            tmp_alignments.erase(sit2++);
          }  else if (SINGLE_SPECIES && overlap2_reverse(a,*sit2)) {
            a.reverse_merge(*sit2);
            tmp_alignments.erase(sit2++);
          } else
            ++sit2;
          //cerr << endl;
        }

        // if the module is too small, don't add it to the cluster list
        /*
        for (i=0;i<num_copies;i++) {
          if (a.proteins[i].size() < MIN_PROTEIN)
            break;
        }
        if (i==num_copies) {
          for (sit2 = saved_alignments.begin();
               sit2 != saved_alignments.end(); sit2++) {
            tmp_alignments.erase(tmp_alignments.find(*sit2));
          }
          */    
        tmp_alignments.erase(sit++);

        //cerr << a << endl;
        clusters.insert(a);
        // clear proteins
        i=0;
        for (vector<vector<int> >::iterator pit = a.proteins.begin(); pit != a.proteins.end(); pit++) {
          for (vector<int>::iterator pit2 = pit->begin(); pit2 != pit->end(); pit2++) {
            //if (i==1 && *pit2 < 9000)
              //cerr << a << endl;
            //cerr << "clear protein " << *pit2 << endl;
            proteins_to_alignments[proteins_to_index[*pit2]].as.clear();
          }
          i++;
        }
        //} else
        //  sit++;

        //variation 1: give up the alignments that were not merged with the first one
        //break;
        number_clusters_per_protein++;
        if (number_clusters_per_protein == NUMBER_CLUSTERS_PER_PROTEIN) break;
      }
      
      tmp_alignments.clear();

      
    } else {
      //cerr << "protein " << vit->p << " has been assigned a cluster" << endl;
    }
  } // end for
  cerr << clusters.size() << " clusters generated" << endl;

  //==============output==============

  string filepath;
  set<int> protein_set;
  ofstream fout2;
  
  ostringstream strm3;
  strm3 << argv[2] << "/clusters";  
  ofstream fout3;
  fout3.open(strm3.str().c_str());
  
  for (i=0;i<num_copies;i++) {
    protein_set.clear();

    /*
    file=species_name[i];
    ostringstream strm4;
    strm4 << "rm -rf " << argv[2] << "/species_" << i;
    system(strm4.str().c_str());
    
    ostringstream strm;
    strm << "mkdir " << argv[2] << "/species_" << i;
    system(strm.str().c_str());

    ostringstream strm2;
    strm2 << "mkdir " << argv[2] << "/species_" << i << "/1-1000";
    system(strm2.str().c_str());
    */

    j=1;
    for (set<struct alignment,lt_cluster>::iterator sit=clusters.begin();
         sit!=clusters.end();sit++,j++) {
      if (i==0) fout3 << *sit << endl;

      const vector<int>& proteins=sit->proteins[i];
      protein_set.insert<vector<int>::const_iterator>(proteins.begin(),proteins.end());

      /*
      ostringstream strm2;
      strm2 << argv[2] << "/species_" << i << "/1-1000/" << j;
      filepath=strm2.str();

      fout.open(filepath.c_str());
      if (!fout) {
        cerr << "file does not exist " << filepath << endl;
        fout.clear();
        return -1;
      }

      strm2 << ".no";
      fout2.open(strm2.str().c_str());
      if (!fout2) {
        cerr << "file does not exist " << filepath << endl;
        fout2.clear();
        return -1;
      }
      */
      
      //for (vector<int>::const_iterator it=proteins.begin();it!=proteins.end();it++) {
        //fout << *it << endl;
      //}

      /*
      if (file=="Vchol") {
        //cerr << "Vchol" << endl;
        for (vector<int>::const_iterator it=proteins.begin();it!=proteins.end();it++) {
          vector<string>& v=protein_names[(int)*it];
          for (vector<string>::iterator it2=v.begin();it2!=v.end();it2++)
            fout << "VC_" << it2->substr(2,it2->length()-2) << endl;
          fout2 << *it << endl;
        }
      } else if (file=="Cjeju") {
        //cerr << "Cjeju" << endl;
        for (vector<int>::const_iterator it=proteins.begin();it!=proteins.end();it++) {
          vector<string>& v=protein_names[(int)*it];
          for (vector<string>::iterator it2=v.begin();it2!=v.end();it2++)
            if ((*it2)[it2->length()-1]=='c') {
              fout << "CJE_" << it2->substr(2,it2->length()-3) << endl;
            } else {
              fout << "CJE_" << it2->substr(2,it2->length()-2) << endl;
            }
          fout2 << *it << endl;
        }
      } else if (file=="Hpylo") {
        //cerr << "Hpylo" << endl;
        for (vector<int>::const_iterator it=proteins.begin();it!=proteins.end();it++) {
          vector<string>& v=protein_names[(int)*it];
          for (vector<string>::iterator it2=v.begin();it2!=v.end();it2++)
            fout << "HP_" << it2->substr(2,it2->length()-2) << endl;
          fout2 << *it << endl;
        }
      } else */
      /*
      for (vector<int>::const_iterator it=proteins.begin();it!=proteins.end();it++) {
        vector<string>& v=protein_names[(int)*it];
        fout << v[0] << endl;
        //fout2 << *it << endl;
      }
      
      fout.close();
      fout.clear();
      */
      //fout2.close();
      //fout2.clear();
    }
    cerr << "species " << i << ": total " << protein_set.size() << " proteins in the cluster" << endl;

  }
  fout3.close();
  fout3.clear();

}
