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
#include <boost/container/stable_vector.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/array.hpp>
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


int initial_city = 0;
int generated_nodes;
int expanded_nodes;

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

vector<short> fill_visited_cities(vector<int> solution) {

	vector<short> cities;
	for(unsigned i = 0; i < ncities; ++i) {
		cities.push_back(0);
	}
	cities[initial_city] = 1;

	if(!solution.empty()) {
		for(unsigned i = 0; i < solution.size(); ++i) {
			cities[solution[i]]= 1;
		}
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

 double in_out(short city, vector<int> subtour, vector<short> visited ) {

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
                    //cout << " - Sumando " << cf[i] << " nodo " << min_edge[cf[i]][count_b].nodo1 << " - Cost: " << min_edge[cf[i]][count_b].cost  << endl;
                }
                count_b++;
                //cin.get();
            }
        }   
    }
   
    val += sorted_edges[city][1].cost + sorted_edges[initial_city][1].cost;

    return val*0.5;

 }



double heuristic1(int city,vector<short> cities_visited) 
{
  int i;
  double value=0;
  for (i=0; i<ncities;i++) 
    if (i!=city && cities_visited[i]==0)
      value+=distance_matrix[i][succ_matrix[i][1]]+distance_matrix[i][succ_matrix[i][2]];

  
  for (i=0; i<ncities;i++)
    if (i!=initial_city && !cities_visited[succ_matrix[initial_city][i]]){
      value+=distance_matrix[initial_city][succ_matrix[initial_city][i]];
      break;
    }
  
  return value/2;
}

int aStar(int init_city, double w, int lookahead) {

	
	int missing_cities = ncities;
	vector<short> cities_visited;
	vector<int> current_solution;
	int iter = 0;

	cout<<"Doing A* search"<<endl;


	
	cities_visited = fill_visited_cities(current_solution);
	
	Node_h* initial_node = new Node_h(init_city,0,0,0,NULL);
	// initial_node->h = in_out(initial_city,cities_visited);

	auto prt_open = open.push(initial_node);
	open_map.emplace(initial_node,prt_open);

	while(!open.empty() && iter < lookahead) {

		
		//cout<<"open size: "<<open.size()<<endl;
		Node_h* current = open.top();
		

		current_solution = create_solution(current);
		cities_visited = fill_visited_cities(current_solution);
		//cout<<"*****Current node****"<<endl;
		//print_node(current,current_solution,cities_visited);
		
		if(current_solution.size() >= ncities) {
			//generate solution report
			//
			
			cout<<"Best solution find (cost): " <<current->g<<endl;
			cout<<"Best solution find (path): ";
			for (vector<int>::reverse_iterator i =current_solution.rbegin(); i != current_solution.rend(); ++i ) { 
				cout<<*i<<" ";
			} 
			cout<<current_solution[ncities-1]<<" ";
			cout<<endl;

			
			return 1;
		}	

		closed.push_back(current);
		//cout<<"size before: "<<open.size()<<endl;
		open.pop();
		//cout<<"size after: "<<open.size()<<endl;

		/*auto node_f = open_map.find(current);
		open_map.erase(node_f->first);
		*/
		missing_cities = ncities - current_solution.size();


		expanded_nodes++;
		
		//getSuccessors
		//cout<<"*****SUCCESSORS****"<<endl;
		for(short i=0; i<ncities; i++){ 
			if(!cities_visited[i] && i != initial_city) {
				int succ_id;
				float h,g,f;
				//when all cities are visited
				if(current_solution.size() >= ncities-1){
					g = current->g + distance_matrix[current->city][i] + distance_matrix[i][initial_city];
					h = 0;
					f = g;

				} else {
					h = in_out(i,current_solution,cities_visited);
					//h = heuristic1(i,cities_visited);

					g = current->g + distance_matrix[current->city][i]; 
					f = (g+h*w)*1000000+h;

				}
				Node_h* succ = new Node_h(i,g,h,f,current);
				//cout<<i<<" ";
				//print_node(succ,current_solution,cities_visited);
				generated_nodes++;
				open.push(succ);
			}

		}
		//cout<<endl; 
		iter++;
		//cin.get();

		
	}
	return -1;
}


void search_driver(int lookahead, double w) {
	
	//add to change search algorithm

		if(lookahead == 0){

			lookahead = LARGE;
		}
		chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
		int sol_response = aStar(initial_city,w,lookahead);
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
		open_map.clear();
		expanded_nodes=0;
		generated_nodes=0;
	
	
}

int main(int argc, char const *argv[])
{
	double w = 1.0;
	int lookahead = 0;
	ncities = 35;
	read_problem("../problems/AdaptedFormat/51.mtsp");
	distance_matrix_caculation();
	/*
	int values[16] = {0,4,5,15,4,0,12,2,5,12,0,3,15,2,3,0};
	int k = 0;
	for(unsigned i = 0; i < ncities; ++i) {
		for(unsigned j = 0; j < ncities; ++j) {

			distance_matrix[i][j] = values[k];
			k++;
		}
	}
	*/
	succ_matrix_caculation();
	sort_edges();

	expanded_nodes = 0;
	generated_nodes = 0;
	search_driver(lookahead,w);

	return 0;
}