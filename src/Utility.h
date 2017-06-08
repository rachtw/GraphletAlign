#ifndef UTILITY_H_
#define UTILITY_H_
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <set>
#include <cstring>
#include <bitset>
#include <queue>
#include <climits>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define dimension_of(X) (sizeof(X)/sizeof(decay_array_to_subtype(X)))

struct node_degree {
  unsigned int n;
  unsigned int d;
  bool operator<(const struct node_degree& nd) const
  {
    if (this->d==nd.d)
      return this->n<nd.n;
    return this->d>nd.d;
  }
  friend std::ostream& operator<<(std::ostream& out, const node_degree& nd)
  {
    out << "{" << nd.n << "," << nd.d << "}";
    return out;
  }
};

void bfs(std::vector<std::vector<unsigned int> >& al,
         std::vector<unsigned int*>& tour) {
  std::set<struct node_degree> unvisited;
  std::set<struct node_degree>::iterator sit;
  std::queue<unsigned int> q;
  std::vector<std::set<unsigned int> > parent(al.size(),std::set<unsigned int>());
  int i,source;
  //std::cerr << "adjacency list:" << std::endl;
  for (i=0;i<al.size();i++) {
    //std::cerr << i << ":";
    unvisited.insert((node_degree){i,al[i].size()});
    //copy(al[i].begin(),al[i].end(),std::ostream_iterator<unsigned int,char>(std::cerr," "));
    //std::cerr << std::endl;
  }
  while (unvisited.size()!=0) {
    //std::cerr << *unvisited.begin() << std::endl;
    q.push(unvisited.begin()->n);
    unvisited.erase(unvisited.begin());
    while (q.size()!=0) {
      source=q.front(); q.pop();
      //std::cerr << "deque " << source << std::endl;
      std::set<struct node_degree> neighbors;
      //std::cerr << "neighbors:" << std::endl;
      for (std::vector<unsigned int>::const_iterator it=al[source].begin();
          it!=al[source].end();it++) {
        //std::cerr << *it << std::endl;
        node_degree nd=(node_degree){*it,al[*it].size()};
        if ((sit=unvisited.find(nd))!=unvisited.end()) {
          //std::cerr << "insert " << nd << std::endl;
          neighbors.insert(nd);
          unvisited.erase(sit);
        } else if (parent[source].find(*it)==parent[source].end()) {
          unsigned int* e=new unsigned int[2];
          e[0]=source;
          e[1]=*it;
          parent[*it].insert(source);
          tour.push_back(e);
          //std::cerr << "store " << e[0] << " " << e[1] << std::endl;
        }
      }
      for (std::set<struct node_degree>::const_iterator it=neighbors.begin();
           it!=neighbors.end();it++) {
        unsigned int* e=new unsigned int[2];
        e[0]=source;
        e[1]=it->n;
        parent[it->n].insert(source);
        tour.push_back(e);
        //std::cerr << "store " << e[0] << " " << e[1] << std::endl;
        q.push(it->n);
        //std::cerr << "enque " << it->n << std::endl;
      }
    }
  }
}

void
graph6toal(std::vector<std::vector<unsigned int> >& al, std::string g6)
{
  int i,j,N_byte_length,k,l,n;
  for (i=0;i<2;i++)
    if ((int)g6.at(i)!=126) break;
  if (i==0) {
    N_byte_length=1; const int N_bit_length=6;
    std::bitset<N_bit_length> bN;
    l=N_bit_length-1;
    for (j=i;j<N_byte_length;j++) {
      std::bitset<6> b((unsigned long)g6.at(j)-63);
      for (k=5;k>=0;k--) {
        bN.set(l,b[k]);
        l--;
      }
    }
    n=(int)bN.to_ulong();
  }
  else if (i==1) {
    N_byte_length=4; const int N_bit_length=18;
    std::bitset<N_bit_length> bN;
    l=N_bit_length-1;
    for (j=i;j<N_byte_length;j++) {
      std::bitset<6> b((unsigned long)g6.at(j)-63);
      for (k=5;k>=0;k--) {
        bN.set(l,b[k]);
        l--;
      }
    }
    n=(int)bN.to_ulong();
  }
  else if (i==2) {
    N_byte_length=8; const int N_bit_length=36;
    std::bitset<N_bit_length> bN;
    l=N_bit_length-1;
    for (j=i;j<N_byte_length;j++) {
      std::bitset<6> b((unsigned long)g6.at(j)-63);
      for (k=5;k>=0;k--) {
        bN.set(l,b[k]);
        l--;
      }
    }
    n=(int)bN.to_ulong();
  }

  //std::cerr << n << std::endl;
  al.assign(n,std::vector<unsigned int>());
  k=1;
  l=0;
  for (i=N_byte_length;i<g6.length();i++) {
    std::bitset<6> b((unsigned long)g6.at(i)-63);
    //std::cerr << b << std::endl;
    for (j=5;j>=0;j--) {
      //std::cerr << l << " " << k << ":" << b.test(j) << std::endl;
      if (b.test(j)) {
        //std::cout << l << " " << k << std::endl;
        al[l].push_back(k);
        al[k].push_back(l);
      }
      l++;
      if (l>=k) {
        l=0; k++;
        if (k>=n) break;
      }
    }
  }
  //std::cout << std::endl;
}

int choose (int n, int r) {
  double c = (double) n/r;
  for (int i=1;i<r;i++) {
    c *= ((double)(n-i))/(r-i);
  }
  return (int) c;
}

void choose(std::vector< std::vector<int> >& indices, int n, int r)
{
  int j;
  indices.clear();
  //std::cerr << "n = " << n << std::endl;
  //std::cerr << "r = " << r << std::endl;

  std::vector<int>* v_=new std::vector<int>();
  indices.push_back(*v_);
  for (j=0;j<r;j++) // increasing sequence
    indices[0].push_back(j);

  //std::cerr << "index:" << std::endl;
  //copy(indices[0].begin(),indices[0].end(),std::ostream_iterator<int>(std::cerr,"\n"));

  int k, i = 1;
  while (true) {
    // start from last index, find the element that does not have the largest possible value
    for (k=r-1; k>=0 && indices[i-1][k]==n-(r-k); k--)
      ;
    if (k==-1) break; // all elements have their largest possible values

    // copy values in the previous combination
    std::vector<int>* v_=new std::vector<int>(indices[i-1]);
    indices.push_back(*v_);
    indices[i][k]++; // increment the element just found
    for (j=k+1;j<r;j++)
      indices[i][j]=indices[i][k]+(j-k); // the element is followed by an increasing sequence

    //std::cerr << "index:" << std::endl;
    //copy(indices[i].begin(),indices[i].end(),std::ostream_iterator<int>(std::cerr,"\n"));

    i++;
  }

  std::cerr << "indices.size()=" << indices.size() << std::endl;

}

void permutation(std::vector< std::vector<short> >& indices, short n, std::vector<short>& solution, std::vector<bool>& used, short k) {
    if (k == n) { // it's a solution
      //copy(solution.begin(),solution.end(),std::ostream_iterator<short,char>(std::cerr,";"));
      //std::cerr << std::endl;
      indices.push_back(solution);
    } else {
      for (int i=0; i<n; i++)  // 試著將第 k 格填入各種數字
        if (!used[i]) {
            used[i] = true;     // 紀錄用過的數字

            solution[k] = i;    // 將第 k 格填入數字 k
            permutation(indices,n,solution,used,k+1);    // iterate next position

            used[i] = false;    // 回收用完的數字
        }
    }
}

void permutation(std::vector< std::vector<short> >& indices, short n) {
  std::vector<short> solution(n);
  std::vector<bool> used(n);
  permutation(indices,n,solution,used,0);
  std::cerr << "permutations:" << std::endl;
  for (int i=0; i<indices.size(); i++) {
    copy(indices[i].begin(),indices[i].end(),std::ostream_iterator<short,char>(std::cerr,","));
    std::cerr << std::endl;
  }
}

template<class T, size_t N> T decay_array_to_subtype(T (&a)[N]);

template<typename T>
bool
stot( T& t, const std::string& s )
{
    std::istringstream iss( s );
    return !(iss >> t).fail();
}

template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

class stringtokenizer {
private:
  std::string _s; // the string to parse
  const char * const _delim; // delimiters = " \t\n" in default
  bool _include_delim; // returns the delimiters as tokens

public:

  class const_iterator { // iterator for stringtokenizer
  private:
    int _i; // the index of the beginning of the token
    int _j; // the index of the end of the token + 1
    stringtokenizer* _st; // a pointer to the stringtokenizer

    void next() { // find the next token

      if (_j == _st->_s.length()) { // there is no more token

        _i = _st->_s.length(); // set both _i & _j to the string length

      } else {

        if (!_st->_include_delim) { // if delimiters will not be returned

          // find the beginning of the token
          _i = _st->_s.find_first_not_of(_st->_delim, _j);
          if (_i==std::string::npos)
            _i=_st->_s.length();

          // find the end of the token
          _j = _st->_s.find_first_of(_st->_delim, _i);
          if (_j==std::string::npos)
            _j=_st->_s.length();

        } else { // if delimiters will be returned as tokens

          // set the beginning of the token to the previous end
          _i = _j;

          if (strchr(_st->_delim,_st->_s.at(_i))) // if the current character is a delimiter

            _j++; // find the end of the token

          else {

            // find the end of the token
            _j = _st->_s.find_first_of(_st->_delim, _i);
            if (_j==std::string::npos)
              _j=_st->_s.length();

          }
        }
      }
    }
  public:
    const_iterator(stringtokenizer* st=NULL, int i=0, int j=0)
      : _st(st), _i(i), _j(j) {
      if (st!=NULL && i==0 && j==0)
        next();
    }
    std::string operator*() const {
      if (_i==_st->_s.length())
        return std::string();
      else
        return _st->_s.substr(_i,_j-_i);
    }
    const_iterator& operator++() {
      if (_i!=_st->_s.length())
        next();
      return *this;
    }
    const_iterator operator++(int) {
      const_iterator it = *this;
      if (_i!=_st->_s.length())
        next();
      return it;
    }
    bool operator==(const_iterator it) const
    { return _st==it._st && _i==it._i && _j==it._j; }
    bool operator!=(const_iterator it) const
    { return _st!=it._st || _i!=it._i || _j!=it._j; }
    friend class stringtokenizer;
  }; // end of const_iterator

  stringtokenizer(const std::string s="", const char * const delim = " \t\n",
                  bool include_delim = false) // default constructor
    : _s(s), _delim(delim), _include_delim(include_delim) {}

  std::string str() { return _s; } // returns the string to parse

  int size() { // returns number of tokens, not including the delimiters
    bool include_delim = _include_delim;
    _include_delim = false;
    int count=0;
    for (const_iterator it=begin(); it!=end(); it++) {
      count++;
    }
    _include_delim = include_delim;
    return count;
  }

  const_iterator begin() { return const_iterator(this); }

  const_iterator end() { return const_iterator(this, _s.length(), _s.length()); }

  friend class stringtokenizer::const_iterator;
};

template <typename T>
void
stringtok (std::vector<T> &v, std::string const &in,
           const char * const delimiters = " \t\n")
{
    const std::string::size_type len = in.length();
          std::string::size_type i = 0;
    while ( i < len )
    {
        // eat leading whitespace
        i = in.find_first_not_of (delimiters, i);
        if (i == std::string::npos)
            return;   // nothing left but white space

        // find the end of the token
        std::string::size_type j = in.find_first_of (delimiters, i);

        T t;
        // push token
        if (j == std::string::npos) {
            stot<T>(t,in.substr(i));
            v.push_back(t);
            return;
        } else {
            stot<T>(t,in.substr(i, j-i));
            v.push_back(t);
        }

        // set up for next loop
        i = j + 1;
    }
}

template <typename T>
void
stringtok (std::set<T> &v, std::string const &in,
           const char * const delimiters = " \t\n")
{
    const std::string::size_type len = in.length();
          std::string::size_type i = 0;
    while ( i < len )
    {
        // eat leading whitespace
        i = in.find_first_not_of (delimiters, i);
        if (i == std::string::npos)
            return;   // nothing left but white space

        // find the end of the token
        std::string::size_type j = in.find_first_of (delimiters, i);

        T t;
        // push token
        if (j == std::string::npos) {
            stot<T>(t,in.substr(i));
            v.insert(t);
            return;
        } else {
            stot<T>(t,in.substr(i, j-i));
            v.insert(t);
        }

        // set up for next loop
        i = j + 1;
    }
}


int
numtok (std::string const &in,
        const char * const delimiters = " \t\n")
{
    const std::string::size_type len = in.size();
          std::string::size_type i = 0;
          std::string::size_type j = 0;
    int num = 0;
    while ( i < len )
    {
        // eat leading whitespace
        i = in.find_first_not_of (delimiters, i);
        if (i == std::string::npos)
            break;   // nothing left but white space

        // find the end of the token
        j = in.find_first_of (delimiters, i);
        num++;

        // push token
        if (j == std::string::npos) {
            break;
        }
        // set up for next loop
        i = j + 1;
    }
    return num;
}

bool endsWith(std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

inline void keep_window_open()
{
	std::cin.clear();
	std::cout << "Please enter a character to exit\n";
	char ch;
	std::cin >> ch;
	return;
}

off_t byte_size_of_file(char* c_str) {
  int file_desc( open(c_str, O_RDWR) );
  struct stat buffer;
  int status( fstat(file_desc, &buffer) );
  close( file_desc );
  return buffer.st_size;
}
#endif
