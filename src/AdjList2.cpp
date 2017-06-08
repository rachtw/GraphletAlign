#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include "Utility.h"
using namespace std;
//./AdjList2 [.adjl] [output: .adjl2]
int main(int argc, char * argv[]) {
  // ====== read neighbors ======
  ifstream fin;
  fin.open(argv[1]);
  if (!fin) {
    cout << "file does not exist " << argv[1] << endl;
    fin.clear();
    return -1;
  }

  cout << "open " << argv[1] << endl;

  int vertex;
  string str;
  int line=1;
  vector<set<int> > vect_neighbor(1);
  while (getline(fin, str)) {
    if (str.empty()) break;
    stringtokenizer st(str,"[];");
    set< int > neighbors;
    stringtokenizer::const_iterator it = st.begin();
    for (it++; it != st.end(); it++) {
      stot<int>(vertex,*it);
      neighbors.insert(vertex);
    }
    neighbors.erase(line); //does not allow self loop
/*
    cerr << line << "[";
    copy(neighbors.begin(),neighbors.end(),ostream_iterator<int,char>(cerr,";"));
    cerr << "]" << endl;
*/
    vect_neighbor.push_back(neighbors);
    line++;
  } // while getline

  fin.close();
  fin.clear();

  cout << "#lines=" << line << endl;

  //calculate indirect neighbors

  ofstream fout;
  fout.open(argv[2]);

  set<int> neighbors;
  int j;
  for (int i=1;i<line;i++) {
    set<int>& s=vect_neighbor[i];
    neighbors.clear();
    for (set<int>::const_iterator sit=s.begin();sit!=s.end();sit++) {
      set<int>& s2=vect_neighbor[*sit];
      neighbors.insert<set<int>::const_iterator>(s2.begin(),s2.end());
    }
    neighbors.erase(i);
    for (set<int>::const_iterator sit=s.begin();sit!=s.end();sit++)
      neighbors.erase(*sit);

    if (neighbors.size()>0) {
      fout << i << "[";
      for (set<int>::const_iterator sit=neighbors.begin();sit!=neighbors.end();sit++)
        fout << *sit << ";";
      fout << "]" << endl;
    } else
      fout << i << "[]" << endl;
  }
  fout.close();
  fout.clear();
}
