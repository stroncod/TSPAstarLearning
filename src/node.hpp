#include <boost/heap/fibonacci_heap.hpp>
#include <boost/any.hpp>
using namespace std;
using namespace boost::heap;

struct Node_h {
    unsigned int city;  // city or state
    double g;
    double h;
    double f; // f = g+h
    //int deph;

    Node_h* father;
    boost::any handler_open;

    //boost::container::vector<unsigned short> sol;
    Node_h(int i,double j, double h_, double k, Node_h* l) : city(i),g(j),h(h_),f(k),father(l){}
   
    ~Node_h(){}
};

struct compare_states {
    bool operator()(const Node_h* s1, const Node_h* s2) const {
        //return n1.id > n2.id;
        return s1->f > s2->f ;
    }
};

struct Edge {
    double cost;
    int to_city;
    Edge(double c, int i) : cost(c),to_city(i){};
};

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


typedef fibonacci_heap<Node_h*,compare<compare_states> >::handle_type open_handle;
fibonacci_heap<Node_h*, boost::heap::compare<compare_states> > open;
boost::unordered_map<Node_h*, open_handle> open_map;
boost::container::vector<Node_h*> closed;