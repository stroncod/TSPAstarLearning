#include <iostream>
#include <cstdint> //or <stdint.h>
#include <limits>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/any.hpp>
#include <boost/range/adaptor/map.hpp>
#include <random> 
#include <boost/container/vector.hpp>
//#include <boost/container/stable_vector.hpp>
//#include <boost/container/static_vector.hpp>
//#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <stdexcept>
#include <chrono>
#include <functional>

#include "node.hpp"


#define max_problem_size 1000
#define LARGE 100000000
#define MAX_CITIES 1002

using namespace std;
using namespace boost::heap;

double cities_coor[max_problem_size][2];
double distance_matrix[MAX_CITIES][MAX_CITIES];
double paso_matrix[max_problem_size][max_problem_size];
int succ_matrix[max_problem_size][max_problem_size];
vector<vector<Edge>> sorted_edges;
int ncities = 0;
double w;

int initial_city = 0;
int generated_nodes;
int expanded_nodes;
int it_driver;

void read_problem(const char *filename) {

    FILE *f;
    int y,x,num;
    f = fopen(filename,"r");
    if (f==NULL){printf( "No se puede abrir el fichero.\n" );}
    rewind (f);
    for (y = 0; y < ncities; y++) {
      for (x = 0; x < 3; x++) {
		fscanf(f, "%d", &num);
		if (x > 0) cities_coor[y][x-1] = num;
      }
    }
    fclose (f);
}


void distance_matrix_caculation() {
  int y,y1;

  for (y = 0; y < ncities; y++) {
    for (y1 = 0; y1 < ncities; y1++) {
      if (y != y1) distance_matrix[y][y1] = sqrt(pow((cities_coor[y1][1]-cities_coor[y][1]),2) + pow((cities_coor[y1][0]-cities_coor[y][0]),2));
      
      else distance_matrix[y][y1] =0;


      //Se inicializa matrices para ordenar sucesores
      paso_matrix[y][y1] = distance_matrix[y][y1];
      succ_matrix[y][y1] = y1;
     
    }
  }

}

void succ_matrix_caculation() {
  int y,x,j,aux1,aux2;
  for (y = 0; y < ncities; y++) {
    for (x = 0; x < ncities; x++) {
      for (j = 0; j < ncities-1; j++) {
        if (paso_matrix[y][j] > paso_matrix[y][j+1]) {
          aux1 = paso_matrix[y][j];
          aux2 = succ_matrix[y][j];
          paso_matrix[y][j] = paso_matrix[y][j+1];
          succ_matrix[y][j] = succ_matrix[y][j+1];
          paso_matrix[y][j+1] = aux1;
          succ_matrix[y][j+1] = aux2;
        }
      }
    }
  }
}


vector<int> create_solution(Node_h* candidate) {

	vector<int> solution;
	solution.push_back(candidate->city);
	Node_h* parent = candidate->father;
	while(parent != NULL) {
		solution.push_back(parent->city);
		parent = parent->father;
	}

	return solution;
}

vector<short> fill_visited_cities(Node_h* candidate) {

	vector<short> cities(ncities);
	cities[candidate->city] = 1;
	cities[initial_city] = 1;
	Node_h* parent = candidate->father;
	while(parent != NULL) {

		cities[parent->city]= 1;
		parent = parent->father;
	}
	return cities;
}


void print_node(Node_h* node, vector<int> solution, vector<short> cities ){
	cout<<"Node id: "<< node->city<<endl;
	cout<<"g: "<<node->g<<endl;
	cout<<"h: "<<node->h<<endl;
	cout<<"f: "<<node->f<<endl;
	cout<<"Subtour: ";
	for(unsigned i = 0; i < solution.size(); ++i) {
		cout<<solution[i]<<" ";
	}
	cout<<endl;
	cout<<"Visited Cities: ";
	for(unsigned i = 0; i < cities.size(); ++i) {
		if(cities[i]) {
			cout<<i<<" ";
		}
	}
	cout<<endl;
}

void sort_edges() {

	vector<double> v;
	vector<Edge> edges;

	for(unsigned i = 0; i < ncities; ++i) {
	 	for(unsigned j = 0; j < ncities; ++j) {
	 		v.push_back(distance_matrix[i][j]);

	 	}
	 	for( auto k: sort_indexes(v) ){
	 		Edge edge = Edge(v[k],k); 		
	 		edges.push_back(edge);
	 	}
	 	sorted_edges.push_back(edges);
	 	v.clear();
	 	edges.clear();
	 } 

 }

double in_out(short city, vector<short> visited ) {

   double val = 0.0; 
   int count_a=0, count_b=0;

   for(int i=0; i<ncities; i++){ // for the missing cities
   		
        if(i != city && !visited[i] ){ //if is not the same city and is not in the subtour
            int count_a = count_b = 0;
            while(count_a < 2){
            	//only for edges who are not connected to an interior node; initial node is never interior
                if(!visited[sorted_edges[i][count_b].to_city]  || sorted_edges[i][count_b].to_city==initial_city ){
                	if(sorted_edges[i][count_b].to_city != i) { //eliminate edges with the same nodes in both sides
                    	val += sorted_edges[i][count_b].cost;
                    	count_a++;
                	}
                }
                count_b++;
                //cin.get();
            }
        }   
    }
   
    val += sorted_edges[city][1].cost + sorted_edges[initial_city][1].cost;

    return val*0.5;

 }


bool in_closed(Node_h * node ) {

	return find(closed.begin(),closed.end(),node) != closed.end();
} 

bool in_open(Node_h * node ){
	auto it_open = open_map.find(node);
			
	if(it_open != open_map.end())  return true;

	return false;
}

void print_open(){
	cout<<"open: ";
	for (auto &i : open )
	{
		cout<<"c: "<<i->city<<" h:"<<i->h<<" ";
	}
	cout<<endl;
}

void print_worst(){
	if(worst_open.size() < 0) {
		cout<<"Worst is empty"<<endl;
	} else {

		cout<<"worst: ";
		for (auto &i : worst_open )
			cout<<i->city<<" ";
		cout<<endl;
	}
}

void backup(Node_h* current) {
	double min_f = 0.0;
	Node_h * parent = current->father;
	if(parent != NULL) {
		for(auto &past_succ: parent->succs){
			if(min_f < past_succ->f) {
				min_f = past_succ->f;
			}
		}

		if(parent->h < min_f) {
			parent->h = min_f;
		}
		backup(parent);
	}

}

void get_successors(Node_h* current, vector<short> cities_visited){

	//al generar estados tengo que revisar que no esten en closed 
	vector<Node_h*> empty = vector<Node_h*>();
	vector<Node_h*> v = vector<Node_h*>() ;
	if(current->father == NULL) {
		for(short i=0; i<ncities; i++){ 
			if(i != initial_city /*&& !in_closed(i)*/ ) {
				float h,g,f;
				//when all cities are visited
				if(current->depth >= ncities-1){
					g = current->g + distance_matrix[current->city][i] + distance_matrix[i][initial_city];
					h = 0;
					f = g;

				} else {
					h = h = in_out(i,cities_visited);
					g = current->g + distance_matrix[current->city][i]; 
					f = (g+h*w)*1000000+h;
				}
				Node_h* succ = new Node_h(i,g,h,f,current->depth+1,current);
				//cout<<i<<" ";
				//print_node(succ,current_solution,cities_visited);
				generated_nodes++;

				if(!in_closed(succ)) {
					auto prt_open = open.push(succ);
					open_map.emplace(succ,prt_open);
					v.push_back(succ);
				}
			}

			current->succs.insert(current->succs.begin(),v.begin(),v.end());
			v.clear();

		}
	} else {

		//cout<<"parent: "<<current->father->city<<endl;
		for(auto &past_succ: current->father->succs){
			
			if(current->city != past_succ->city && past_succ->city != initial_city /*&& !in_closed(past_succ->city)*/){
				float h,g,f;

				//cout<<"Succ: "<<past_succ->city<<endl;
				//when all cities are visited
				if(current->depth >= ncities-1){
					g = current->g + distance_matrix[current->city][past_succ->city] + distance_matrix[past_succ->city][initial_city];
					h = 0;
					f = g;

				} else {
					double 	h0 = in_out(past_succ->city,cities_visited);

					h = max((current->h-distance_matrix[current->city][past_succ->city]),h0);
					//pathmax
					g = current->g + distance_matrix[current->city][past_succ->city]; 
					f = (g+h*w)*1000000+h;

				}
				
				Node_h* succ = new Node_h(past_succ->city,g,h,f,current->depth+1,current);
				generated_nodes++;
				auto prt_open = open.push(succ);				
				open_map.emplace(succ,prt_open);
				v.push_back(succ);

			}
			current->succs.insert(current->succs.begin(),v.begin(),v.end());
			v.clear();		
		}
	}


	auto prt_worst = worst_open.push(current);
	worst_map.emplace(current,prt_worst);

	auto it_worst = worst_map.find(current->father);
	
	if(it_worst != worst_map.end()) {

		worst_open.erase(it_worst->second);
		
		worst_map.erase(it_worst->first);
	}
	
	
}

void init_first_node(){

	vector<Node_h*> v;
	Node_h* initial_node = new Node_h(initial_city,0,0,0,1,NULL); 
	auto prt_open = open.push(initial_node);
	//auto prt_worst = worst_open.push(initial_node);
	//worst_map.emplace(initial_node,prt_worst);
	open_map.emplace(initial_node,prt_open);
}

Node_h* aStar(int init_city, double w, int lookahead, int backsteps) {
	
	int missing_cities = ncities;
	vector<short> cities_visited;
	vector<int> current_solution;
	Node_h* current;
	int iter = 0;


	if(it_driver > 0) lookahead = backsteps;

	while(!open.empty() && iter < lookahead) {
		current = open.top();
		cout<<"c: "<<current->city<<endl;
		cout<<"h: "<<current->h<<endl;
		if(current->father != NULL)
			cout<<"p: "<<current->father->city<<endl;


		cities_visited = fill_visited_cities(current);

		if(current->depth >= ncities) return current;

		closed.push_back(current);
		open.pop();
		open_map.erase (current);

		expanded_nodes++;
		
		get_successors(current,cities_visited);
		//print_open();
		//print_worst();
		cin.get();
		iter++;	
	}
	return current;
}



void undo_Astar(int backsteps) {

	double min_h = LARGE;

	if(worst_open.size() < backsteps) backsteps = worst_open.size(); 

	for(unsigned i = 0; i < backsteps; ++i) {
	
		Node_h* worst = worst_open.top();

		cout<<"w: "<<worst->city<<endl;

		for (auto &node : worst->succs){
			
			cout<<"succ: "<<node->city<<endl;

			if(in_open(node)){
				open.erase(open_map[node]);
				open_map.erase(node);
			}

			if((node->h + distance_matrix[worst->city][node->city])< min_h ) 
				min_h = node->h + distance_matrix[worst->city][node->city];

		}

		worst->h =  max(worst->h,min_h);
		worst->f =  (worst->g+worst->h*w)*1000000+worst->h;
		auto prt_open = open.push(worst);
		open_map.emplace(worst,prt_open);
		closed.erase(remove(closed.begin(), closed.end(), worst), closed.end());


		worst_open.pop();
		int i_succ = 0;
		min_h = LARGE;
		for (auto &succ : worst->father->succs){
			if(in_open(succ)) {
				//cout<<"entro"<<endl;
				i_succ++;
			}

			if((succ->h + distance_matrix[worst->father->city][succ->city])< min_h ) 
				min_h = succ->h + distance_matrix[worst->father->city][succ->city];
		}

		if(i_succ == worst->father->succs.size() ){
			//cout<<"agregando a worst"<<endl;
			cout<<"wp: "<<worst->father->city<<endl;
			worst->father->h =  max(worst->h,min_h);
			cout<<"new wp h: "<<worst->father->h<<endl;
			cout<<"new w h: "<<worst->h<<endl;
			worst->father->f =  (worst->g+worst->h*w)*1000000+worst->h;
			auto prt_worst = worst_open.push(worst->father);
			worst_map.emplace(worst->father,prt_worst);
		}

	}

	//print_worst();
	//print_open();
	cin.get();

}


int rm_driver (double w, int lookahead, int backsteps) {

	init_first_node();
	while(1) {
		cout<<"Doing A* search"<<endl;
		
		Node_h* result = aStar(initial_city,w,lookahead, backsteps);


		if(result->depth >= ncities) {
			//generate solution report
			vector<int> final_solution = create_solution(result);
			cout<<"Best solution find (cost): " <<result->g<<endl;
			cout<<"Best solution find (path): ";
			for (vector<int>::reverse_iterator i =final_solution.rbegin(); i != final_solution.rend(); ++i )
				cout<<*i<<" "; 
			cout<<final_solution[ncities-1]<<" ";
			cout<<endl;

			return 1;
		}
		cout<<"Undoing A*"<<endl;	
		undo_Astar(backsteps);

		it_driver++;
	}
	return -1;
}


void search_driver(int lookahead, double w, int backsteps) {
	
	//add to change search algorithm

		if(lookahead == 0){

			lookahead = LARGE;
		}
		it_driver = 0;

		chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
		//int sol_response = aStar(initial_city,w,lookahead);
		int sol_response = rm_driver(w,lookahead,backsteps);
		chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();

		if(sol_response != -1){
			chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(end-start);
			cout<<"Resolution Time: "<<time_span.count()<<endl;
			cout<<"Expanded Nodes: "<<expanded_nodes<<endl;
			cout<<"Generated Nodes: "<<generated_nodes<<endl;
		} else {
			cout<<"No solution found!"<<endl;
		}

		open.clear();
		closed.clear();
		worst_map.clear();
		expanded_nodes=0;
		generated_nodes=0;
	
	
}

int main(int argc, char const *argv[])
{
	w = stod(argv[1]);
	int lookahead = stoi(argv[2]);
	cout<<"l:" <<lookahead<<endl;
	int backsteps = stoi(argv[3]);
	ncities = 4;
	read_problem("../problems/AdaptedFormat/14.mtsp");
	distance_matrix_caculation();
	
	int values[16] = {0,4,5,15,4,0,12,2,5,12,0,3,15,2,3,0};
	int k = 0;
	for(unsigned i = 0; i < ncities; ++i) {
		for(unsigned j = 0; j < ncities; ++j) {

			distance_matrix[i][j] = values[k];
			k++;
		}
	}
	
	succ_matrix_caculation();
	sort_edges();

	expanded_nodes = 0;
	generated_nodes = 0;
	search_driver(lookahead,w,backsteps);

	return 0;
}