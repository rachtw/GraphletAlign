#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include "Utility.h"
#include <map>
#include <algorithm>
using namespace std;
//SortBlastp *.blastp.sorted
//The input is sorted by protein names
//The output is sorted by names & expect values
//write *.blastp.sorted.sorted
struct homology {
  string input;
  string p1,p2;
  long double e;
  bool operator<(const homology& h) const
  {
    if (p1==h.p1) {
      if (e==h.e) {
        return p2 < h.p2;
      }
      return e < h.e;
    }
    return p1 < h.p1;
  }
  //overload << for printing homology object
  friend std::ostream& operator<<(std::ostream& out, const homology& h)
  {
    out << h.input;
    return out;
  }  
};
int main(int argc, char * args[]) {
  ifstream fin;  
  string str;
  int i,j;
  set<struct homology> s;
  
  string p1,p2,pre_p1="",pre_p2="",input,smallest_input;
  long double e,smallest_e=0;
  
  //read first line
  fin.open(args[1]);
  if (!fin) {
    cerr << "file does not exist " << args[1] << endl;
    fin.clear();
    return -1;
  }
  getline(fin, str);
  smallest_input=str;
  stringtokenizer st(str," \t");
  stringtokenizer::const_iterator it = st.begin();
  pre_p1=*it;
  it++;
  pre_p2=*it;
  it++;
  smallest_e=atof((*it).c_str());  
  fin.close();
  fin.clear();
  
  fin.open(args[1]);
  if (!fin) {
    cerr << "file does not exist " << args[1] << endl;
    fin.clear();
    return -1;
  }
  getline(fin, str);
  while (getline(fin, str)) {
    //cerr << str << endl;
    input=str;
    stringtokenizer st(str," \t");
    stringtokenizer::const_iterator it = st.begin();
    p1=*it;
    it++;
    p2=*it;
    it++;
    e=atof((*it).c_str());
    
    if (p1==pre_p1 && p2==pre_p2) {
      if (e<smallest_e) {
        smallest_input=input;
        smallest_e=e;
      }
    } else {
      struct homology h;
      h.input=smallest_input;
      h.e=smallest_e;
      h.p1=pre_p1;
      h.p2=pre_p2;
      s.insert(h);
    
      pre_p1=p1;
      pre_p2=p2;
      smallest_input=input;
      smallest_e=e;
    }    
    
  } //while
  fin.close();
  fin.clear();
  struct homology h;
  h.input=smallest_input;
  h.e=smallest_e;
  h.p1=pre_p1;
  h.p2=pre_p2;
  s.insert(h);

  ofstream fout;
  stringstream ss;
  ss << args[1] << ".sorted";
  fout.open(ss.str().c_str());
  for (set<struct homology>::const_iterator it=s.begin();it!=s.end();it++) {
    fout << *it << endl;
  }
  fout.close();
  fout.clear();
}
