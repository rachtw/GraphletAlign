#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include "Utility.h"
#include <map>
#include <cstdlib>
#include <cassert>
#include <sstream>
using namespace std;
//BlastpTopHits *name... 
//argc-6: [*blastp] 
//argc-5: [OUTPUT blastp with protein no] 
//argc-4: [OUTPUT hit groups] 
//argc-3: [E-VALUE THRESHOLD] 
//argc-2: [# TOP_HITS] 
//argc-1: [BIDIRECTION]

//write *.groups
//Assume input blastp file contains sorted hits

long double THRESHOLD=1e-6; // must < THRESHOLD, ex. 1e-6 means allowing Xe-7
int TOP_HITS=INT_MAX; //Top TOP_HITS Bidirectional Hits "among all species" 
const bool EXPECTATION_VALUE_OR_BIT_SCORE=true;
bool BI_DIRECTION=true;
const bool WITHIN_SPECIES=true; // if WITHIN_SPECIES is set, the top hits can be within the same species

struct protein {
  string name;
  string no;
  bool operator==(const struct protein& p2) const {
    return name==p2.name && no==p2.no;
  }
  bool operator<(const struct protein& p2) const {
    if (name.size()==p2.name.size())
      return name < p2.name;
    return name.size()<p2.name.size();
  }
  //overload << for printing VertexDisjointGraph object
  friend std::ostream& operator<<(std::ostream& out, const protein& p)
  {
    //out << "{" << mp.p << "," << mp.e << "}";
    out << p.name << "(" << p.no << ")";
    return out;
  }
};

struct matching_protein {
  struct protein p;
  long double e;
  bool operator==(const matching_protein& m2) const
  {
    return p==m2.p && e==m2.e;
  }
  bool operator<(const matching_protein& m2) const
  {
    if (EXPECTATION_VALUE_OR_BIT_SCORE) {
      if (e==m2.e) {
        return p < m2.p;
      }
      return e < m2.e;
    } else {
      if (e==m2.e) {
        return p < m2.p;
      }
      return e > m2.e;
    }
  }
  //overload << for printing VertexDisjointGraph object
  friend std::ostream& operator<<(std::ostream& out, const matching_protein& mp)
  {
    out << "{" << mp.p << "," << mp.e << "}";
    return out;
  }
};

struct lt_integer
{
  bool operator()(const string& s1, const string& s2) const
  {
    if (s1.size()==s2.size())
      return s1 < s2;
    return s1.size()<s2.size();
  }
};

bool lt_integer_func(const string& s1, const string& s2)
{
  if (s1.size()==s2.size())
    return s1 < s2;
  return s1.size()<s2.size();
}

vector<string> num_proteins_each_species;
int find_species(string protein_no) {
  for (int i=0;i<num_proteins_each_species.size();i++) {
    if (lt_integer_func(protein_no,num_proteins_each_species[i])) {
      //cerr << protein_no << "@" << i << endl;
      return i;
    }
  }
  return -1;
}

map<string,string> acc_to_no;
string get_no(string acc) {
  map<string, string>::const_iterator mit;
  mit=acc_to_no.find(acc);
  if (mit!=acc_to_no.end()) {
    //cerr << first_protein_no;
    return mit->second;
  }
  return "";
}

struct lt_no
{
  bool operator()(const string& s1, const string& s2) const
  {
    return lt_integer_func(get_no(s1),get_no(s2));
  }
};
int main(int argc, char * argv[]) {
  //map<string, string>::const_iterator mit1,mit2;
  set<string> accs;

  set<matching_protein> sorted_proteins;
  map<string, set<matching_protein>,lt_integer > no_to_top_matches;
  map<string, map<string, long double, lt_integer>, lt_integer > blast;

//argc-6: [*blastp] 
//argc-5: [OUTPUT blastp with protein no] 
//argc-4: [OUTPUT hit groups] 
//argc-3: [E-VALUE THRESHOLD] 
//argc-2: [# TOP_HITS] 
//argc-1: [BIDIRECTION]
  if (argv[argc-1][0]=='1')
    BI_DIRECTION=true;
  else
    BI_DIRECTION=false;
  
  THRESHOLD=atof(argv[argc-3]);
  cout << "BLASTP expectation value threshold=" << THRESHOLD << endl;
  
  if (strcmp(argv[argc-2],"INT_MAX")!=0)
    TOP_HITS = atoi(argv[argc-2]);
  cout << "Find top " << TOP_HITS << " hits" << endl;
  
  string str,p1,p2;
  ifstream fin;
  ofstream fout;

  //open *.acc
  int line=1;
  for (int i=1;i<argc-6;i++) {
    fin.open(argv[i]);
    if (!fin) {
      cerr << "file does not exist " << argv[i] << endl;
      fin.clear();
      return -1;
    }
    cout << "Open " << argv[i] << endl;
    int i;
    while (getline(fin, str)) {
      stringtokenizer st(str,"; ");
      for (stringtokenizer::const_iterator it = st.begin();
           it != st.end(); it++) {
        str=*it;
      /*
        i=str.find_first_of("-");
        if (i!=string::npos)
          str=str.substr(0,i);
        */
        acc_to_no[str]=to_string<int>(line);
        accs.insert(str);
      }
      line++;
    }
    num_proteins_each_species.push_back(to_string<int>(line-1));
    fin.close();
    fin.clear();
  }
  
  cout << "Number of accumulated proteins in species: ";
  copy(num_proteins_each_species.begin(),num_proteins_each_species.end(),ostream_iterator<string,char>(cout," "));
  cout << endl; 
  
  //read blast results
  
//argc-6: [*blastp] 
//argc-5: [OUTPUT blastp with protein no] 
//argc-4: [OUTPUT hit groups] 
//argc-3: [E-VALUE THRESHOLD] 
//argc-2: [# TOP_HITS] 
//argc-1: [BIDIRECTION]
  
  string first_protein="",second_protein,pre_first_protein,pre_second_protein;
  string first_protein_no,second_protein_no,pre_first_protein_no;
  int first_species;
  long double e;
  fin.open(argv[argc-6]);
  if (!fin) {
    cerr << "file does not exist " << argv[argc-6] << endl;
    fin.clear();
    return -1;
  }
	cout << "open " << argv[argc-6] << endl;
  pre_first_protein=""; pre_second_protein="";
  
  stringstream ss;
  ss << argv[argc-5];
  fout.open(ss.str().c_str());
  if (!fout) {
    cerr << "file does not exist " << ss.str() << endl;
    fout.clear();
    return -1;
  }
  cout << "open " << ss.str() << endl;
  
  while (getline(fin, str)) {
    //cerr << str << endl;
    stringtokenizer st(str," \t");
    stringtokenizer::const_iterator it = st.begin();
    first_protein=*it;
    first_protein_no=get_no(first_protein);
    if (first_protein_no=="") {
      //cerr << "cannot find " << first_protein << endl;
      continue;
    }
    if (first_protein!=pre_first_protein) {
  /* copy(sorted_proteins.begin(),sorted_proteins.end(),ostream_iterator<matching_protein,char>(cerr,"\n"));
      cerr << endl;
  */
      
      pre_first_protein_no=get_no(pre_first_protein);
      if (pre_first_protein!="" && pre_first_protein_no=="") {
        cerr << "bug: cannot find " << pre_first_protein << endl;
        return -1;
      } else {
        //cerr << pre_first_protein << "(" << pre_first_protein_no;
        first_species=find_species(pre_first_protein_no);
        //cerr << "@" << first_species << "): ";
        set<matching_protein>& map_first_protein_no_to=no_to_top_matches[pre_first_protein_no];

        //line=1;
        for (set<matching_protein>::const_iterator sit=sorted_proteins.begin();
            sit!=sorted_proteins.end();sit++) {

          second_protein=sit->p.name;
          second_protein_no=sit->p.no;
          if (second_protein_no=="") {
            cerr << "bug: cannot find " << second_protein << endl;
            return -1;
          }
          //cerr << second_protein << "(" << second_protein_no;
          //cerr << "@" << find_species(second_protein_no) << ")";
          //same species
          if (!WITHIN_SPECIES && num_proteins_each_species.size()>1 && find_species(second_protein_no)==first_species) {
            //cerr  << endl;
            continue;
          }
          //same protein
          if (pre_first_protein_no==second_protein_no) {
            cerr  << "same protein" << pre_first_protein_no << endl;
            continue;
          }
          if (EXPECTATION_VALUE_OR_BIT_SCORE && sit->e>THRESHOLD) {
            //cerr  << endl;
            continue;
          }
          //if (line<=TOP_HITS) {
            map_first_protein_no_to.insert((matching_protein){(protein){second_protein,second_protein_no},sit->e});
            //map_first_protein_to.insert((matching_protein){(protein){second_protein,second_protein_no},sit->e});

            blast[pre_first_protein_no][second_protein_no]=sit->e;
          //}
          //line++;
          //cerr  << endl;
        }
        //cerr << endl;
      }
      sorted_proteins.clear();
    }
    it++; second_protein=*it;
    it++; e=atof((*it).c_str());
    //cerr << "examine " << first_protein << " " << second_protein << " " << e << endl;
    second_protein_no=get_no(second_protein);
    if (second_protein_no=="") {
      //cerr << "cannot find " << second_protein << endl;
      continue;
    }
    if (first_protein==pre_first_protein && second_protein==pre_second_protein) {
      //cerr << "Previous pair:" << pre_first_protein << " " << pre_second_protein << endl;
      //cerr << "skip " << first_protein << " " << second_protein << " " << e << endl;
      continue;
    }
    // output blastp pairs
    fout << first_protein_no << "\t" << second_protein_no << "\t" << *it << endl;

    sorted_proteins.insert((matching_protein){(protein){second_protein,second_protein_no},e});

    pre_first_protein=first_protein;
    pre_second_protein=second_protein;
  } //while
  fin.close();
  fin.clear();
  fout.close();
  fout.clear();

  //copy(sorted_proteins.begin(),sorted_proteins.end(),ostream_iterator<matching_protein,char>(cerr,"\n"));
  //cerr << endl;
  
  pre_first_protein_no=get_no(pre_first_protein);
  if (pre_first_protein_no=="") {
    cerr << "cannot find " << pre_first_protein << endl;
  } else {

    //cerr << pre_first_protein << "(" << pre_first_protein_no;
    first_species=find_species(pre_first_protein_no);
    //cerr << "@" << first_species << "): " << endl;
    set<matching_protein>& map_first_protein_no_to=no_to_top_matches[pre_first_protein_no];

    //line=1;
    for (set<matching_protein>::const_iterator sit=sorted_proteins.begin();
    sit!=sorted_proteins.end();sit++) {

      second_protein=sit->p.name;
      second_protein_no=sit->p.no;
      if (second_protein_no=="") {
        cerr << "cannot find " << second_protein << endl;
        continue;
      }
      //cerr << second_protein << "(" << second_protein_no;
      //cerr << "@" << find_species(second_protein_no) << ")";
      //same species
      if (!WITHIN_SPECIES && num_proteins_each_species.size()>1 && find_species(second_protein_no)==first_species) {
          //cerr  << endl;
          continue;
      }
      //same protein
      if (pre_first_protein_no==second_protein_no) {
          cerr  << "same protein" << pre_first_protein_no << endl;
          continue;
      }
      if (EXPECTATION_VALUE_OR_BIT_SCORE && sit->e>THRESHOLD) {
          //cerr  << endl;
          continue;
      }
      //if (line<=TOP_HITS) {
        map_first_protein_no_to.insert((matching_protein){(protein){second_protein,second_protein_no},sit->e});
        //map_first_protein_to.insert((matching_protein){(protein){second_protein,second_protein_no},sit->e});
      //}
      //cerr  << endl;
      //line++;
    }
    //cerr << endl;
  }
  sorted_proteins.clear();
  pre_first_protein=first_protein;

//argc-6: [*blastp] 
//argc-5: [OUTPUT blastp with protein no] 
//argc-4: [OUTPUT hit groups] 
//argc-3: [E-VALUE THRESHOLD] 
//argc-2: [# TOP_HITS] 
//argc-1: [BIDIRECTION]  
  fout.open(argv[argc-4]);
  if (!fout) {
    cerr << "file does not exist " << argv[argc-4] << endl;
    fout.clear();
    return -1;
  }
  cout << "open " << argv[argc-4] << endl;
      
  //erase those that are not bidirectional
  bool not_found;
  map<string, set<matching_protein>,lt_no >::iterator mit2=no_to_top_matches.begin();
  for (mit2=no_to_top_matches.begin();mit2!=no_to_top_matches.end();) {

    /*
    cerr << mit->first << "(" << get_no(mit->first) << "):";
    copy(mit->second.begin(),mit->second.end(),ostream_iterator<matching_protein,char>(cerr,", "));
    cerr << endl;
	*/
    /*
    cerr << mit2->first << ":" ;
    copy(mit2->second.begin(),mit2->second.end(),ostream_iterator<matching_protein,char>(cerr,", "));
    cerr << endl;
    */
    set<matching_protein> map_first_protein_no_to=mit2->second;
    set<matching_protein>::iterator sit2;

    if (BI_DIRECTION) {
      sit2=map_first_protein_no_to.begin();
      for (set<matching_protein>::iterator sit2=map_first_protein_no_to.begin();
        sit2!=map_first_protein_no_to.end();) {
        //cerr << "look for " << *sit2 << endl;
        map<string, set<matching_protein>,lt_integer>::iterator _mit2;
        set<matching_protein>::const_iterator _sit;

        //second of the pair found
        //if ((_mit=name_to_top_matches.find(sit->p.name))!=name_to_top_matches.end()) {
        if ((_mit2=no_to_top_matches.find(sit2->p.no))!=no_to_top_matches.end()) {

          // try to find the first
          //set<matching_protein>& _map_first_protein_to=_mit->second;
          set<matching_protein>& _map_first_protein_no_to=_mit2->second;
          not_found = true;
          //for (set<matching_protein>::iterator __sit=_map_first_protein_to.begin();
            //__sit!=_map_first_protein_to.end();__sit++) {
          for (set<matching_protein>::iterator __sit2=_map_first_protein_no_to.begin();
            __sit2!=_map_first_protein_no_to.end();__sit2++) {
		  
            //if (__sit->p.name==mit->first) { //found
            if (__sit2->p.no==mit2->first) { //found
                //cerr << "found" << endl;
                not_found = false;
                break;
            }
          }
          if (not_found) { 
            //cerr << "erase " << *sit2 << " from " << mit2->first << " => ";
            //map_first_protein_to.erase(sit++);
            map_first_protein_no_to.erase(sit2++);
          } else {
            //++sit;
            ++sit2;
          }
        } else {
          //map_first_protein_to.erase(sit++);
          map_first_protein_no_to.erase(sit2++);
        }
      }
    }
    /*
    cerr << "results: " << endl;
    cerr << mit->first << "(" << get_no(mit->first) << "):";
    copy(mit->second.begin(),mit->second.end(),ostream_iterator<matching_protein,char>(cerr,", "));
    cerr << endl;
    cerr << mit2->first << ":" ;
    copy(mit2->second.begin(),mit2->second.end(),ostream_iterator<matching_protein,char>(cerr,", "));
    cerr << endl;
    */

    //set emptied
    if (map_first_protein_no_to.size()==0) {
    //if (map_first_protein_to.size()==0) {
      //name_to_top_matches.erase(mit++);
      no_to_top_matches.erase(mit2++);
      //cerr << "set emptied" << endl;
    } else {
      
      if (TOP_HITS!=INT_MAX) {
        //remove proteins besides top hits

        int line=1;
        int size=map_first_protein_no_to.size();
        for (sit2=map_first_protein_no_to.begin();sit2!=map_first_protein_no_to.end();sit2++) {
          if (line>TOP_HITS) break;
          line++;
        }
        for (;sit2!=map_first_protein_no_to.end();) {
          map_first_protein_no_to.erase(sit2++);
        }
        if (map_first_protein_no_to.size()!=TOP_HITS && size!=map_first_protein_no_to.size()) {
          cerr << "top " << map_first_protein_no_to.size() << " hits" << endl;
          return -1;
        }
      }
     
      set<string,lt_integer> s;
      s.insert(mit2->first);
      for (set<matching_protein>::const_iterator sit=map_first_protein_no_to.begin();sit!=map_first_protein_no_to.end();sit++) {
        s.insert(sit->p.no);
      }
      set<string,lt_integer>::const_iterator strit=s.end();
      strit--;
      fout << "[";
      copy(s.begin(),strit,ostream_iterator<string,char>(fout,", "));
      fout << (*strit);
      fout << "]" << endl;
      //++mit;
      ++mit2;
      
    }
  }
  fout.close();
  fout.clear();  
}
