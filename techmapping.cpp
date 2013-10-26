/****************************************************************
Author       : Abhishek Bhowmick
Date         : 
Organization : Electrical Engg Dept, IIT Bombay
******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <new>
#include <algorithm>
#include <utility>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

using namespace boost; 
using namespace std;


template <class Graph> struct exercise_vertex {					// this functor does operations on a vertex
    exercise_vertex(Graph& g_) : g(g_) {}
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

    void operator()(const Vertex& v) const					// input argument needed is vertex descriptor
    {
      typedef graph_traits<Graph> GraphTraits;
      typename property_map<Graph, vertex_index_t>::type 
      index = get(vertex_index, g);


// this block will list the out_edges of a vertex
      /*
      std::cout << "out-edges: ";						
      typename GraphTraits::out_edge_iterator out_i, out_end;
      typename GraphTraits::edge_descriptor e;
      for (tie(out_i, out_end) = out_edges(v, g); 
           out_i != out_end; ++out_i) {
        e = *out_i;
        Vertex src = source(e, g), targ = target(e, g);
        std::cout << "(" << index[src] << "," 
                  << index[targ] << ") ";
      */

      std::cout << "in-vertices of " << index[v] << " are ";
      typename GraphTraits::in_edge_iterator in_i, in_end;
      typename GraphTraits::edge_descriptor e;
      for (tie(in_i, in_end) = in_edges(v,g); 			// tie function unpacks the elements in a ...
           in_i != in_end; ++in_i) {				// ... tuple(returned by in_edges) into the variables in_i, in_end
        e = *in_i;
        Vertex src = source(e, g), targ = target(e, g);
        std::cout << index[src] << ", " ;
        }
      std::cout << std::endl;
      }
    Graph& g;
  };

struct list_node {
	int nodeid, nodecost;
	std::string nodeattr;
	struct list_node *next;
	} ;

struct list_node* lists[100];						// memory for 10 graph patterns with 10 linked lists each
int listcount;
int list_id=10;								// this is used to assign 'list_id's to the gate patterns

struct NodeProperty {							// this is a property type for reading graphml attributes
        std::string nattr;
	int type;
        };

template <class Graph> struct extract_lists {				// extract_lists will list out all lists of a graph
    extract_lists(Graph& g_) : g(g_) {}
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    
    void operator()(const Vertex& v) const				// input argument needed is vertex descriptor
    {
      typedef graph_traits<Graph> GraphTraits;
      typename property_map<Graph, vertex_index_t>::type 
      index = get(vertex_index, g);
      typename GraphTraits::out_edge_iterator out_i, out_end;
      typename GraphTraits::edge_descriptor e;
      typename GraphTraits::in_edge_iterator in_i, in_end;      
        
	tie(out_i, out_end) = out_edges(v,g);
	if(out_i==out_end)						// we have an input node, create a linked list for this
	{
		Vertex src;
		
		list_node *temp = new list_node;
		temp -> nodeid = index[v];
		temp -> nodeattr = get(get(&NodeProperty::nattr, g), v);	// second argument is vertex descriptor
		temp -> next = lists[listcount];
		lists[listcount] = temp ;

		tie(in_i, in_end) = in_edges(v,g);

		while(in_i!=in_end)
		{
		e = *in_i;
		src = source(e, g);

		list_node *temp = new list_node;
		temp -> nodeid = index[src];
		temp -> nodeattr = get(get(&NodeProperty::nattr, g), src);
		temp -> next = lists[listcount];
		lists[listcount] = temp;
		
		tie(in_i, in_end) = in_edges(src,g);		
		}
		listcount++;
	}	
      }
    Graph& g;
  };


typedef boost::adjacency_list<vecS, vecS, bidirectionalS, NodeProperty> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
char line[100];


struct lib_gate
{ 
  // int nbr_inputs, bool branches[nbr_inputs] (list the branches for all inputs, two inputs
  // have the same branch, if they are input to a NAND) 
  int nbrpatterns,cost;
  char name[10];  
  char graphfiles[5][100];					// a gate can have atmost 5 patterns
  Graph gr[5];							// 5 graphs for 5 possible patterns
  int gr_list_id[5];						// 5 'list_id's for the 5 possible patterns

} ;

std::vector<lib_gate> gates;
int nbr_library_gates=0;

struct matching_t								// Every node will have one matching
{
	char gate_type[10];							// To store the type of gate
	std::vector<int> coveredVertices;					// To store the vertices covered
	std::vector<int> inputVertices;						// To store the input vertices
	int cost;								// To store the cost of this matching
} ;

matching_t matches[100];						// 100 is the maximum number of nodes
									// the indices of matches will be corresponding to 'nodeid's

std::vector<matching_t> tempMatches;					// The list of temporary matches from which best match will be chosen

void findMatch(Graph& g, const Vertex& v)				// This function finds all matches at node v and pushes them to tempMatches
{
typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);		

	std::vector<int> coveredVertices;				// just store the 'nodeid's of the vertices covered
	std::vector<int> inputVertices;					// just store the 'nodeid's of the input vertices
	std::vector<int> tempVertices;					// create a new temporary list of vertices;

std::string input = "i";						// For the node which is an input
matching_t currentMatch;						// This will be the current match

for(int i=0;i<gates.size();i++)						//for all the libgates
{	
	lib_gate* currentgate;
	currentgate = &gates.at(i);	
	int matchflag = 0;

	for(int j=0;j<gates.at(i).nbrpatterns;j++)			//for all the patterns of the libgate
	{
		int patternmatchflag = 1;

		//create temporary list of strings of subject graph;
		std::vector<list_node*> subjectDAGlists;				// it has to be a vector	 
		int k=0;							// list_id for the subject DAG is 0
		
		while(lists[k]!=NULL)
		{
		subjectDAGlists.push_back(lists[k]);
		k++;
		}
		
		int nbrDAGlists = subjectDAGlists.size();
		
		int m = gates.at(i).gr_list_id[j];
		while(lists[m]!=NULL)						//for all the strings in the pattern graph
		{	
			int stringmatchflag = 0;
			for(int n=0;n<nbrDAGlists;n++)				//for all the strings in the subject graph
			{
				//we have two lists now : subjectDAGlists[n], lists[m] : We have to match them

				//match string of pattern graph to subject graph lists (add the vertices to the temp list of vertices);
				
				int flag = 1;
				list_node *travelnode1, *travelnode2;		// travelnode1 is for DAG, travelnode2 is for pattern graph
				travelnode1 = subjectDAGlists.at(n);
				 
				// start from the point in the DAG list which equals the given node
				while(((travelnode1 -> nodeid)!=index[v])&&((travelnode1 -> next)!=NULL)) // segfault	
				{
				travelnode1 = travelnode1 -> next;
				}
				 
				if((travelnode1 -> next)==NULL)
					{
					continue;					// this DAG string doesn't contain vertex v
					}						

				travelnode2 = lists[m];
				
				while((travelnode2 -> nodeattr)!=input)			// traverse till input of the string of pattern graph
				{							//  string comparison
					if((travelnode2 -> nodeattr)==(travelnode1 -> nodeattr))
						{
							tempVertices.push_back((travelnode1 -> nodeid));//add the vertices to the temp list of vertices
							travelnode1 = travelnode1 -> next;
							travelnode2 = travelnode2 -> next;
						}

					else
						{
							flag = 0;
							break;
						}

				}

			
				if(flag == 1)				// match found equivalent to flag = 1
					{	
					//add the temp list to the coveredVertices list;
					for(int p=0;p < tempVertices.size();p++)	coveredVertices.push_back(tempVertices.at(p));
					
					//add the end of the string to the inputVertices list;
					inputVertices.push_back((travelnode1 -> nodeid));

					//remove string from temporary list of strings of subject graph;
					subjectDAGlists.erase(subjectDAGlists.begin()+n);

					stringmatchflag = 1;
					break;
					}

				else	tempVertices.clear();				//clear tempVertices
					

			}

			
			if(stringmatchflag == 0)	
					{
					patternmatchflag = 0;
					break;	
					}
		m++;												
		}	

		if(patternmatchflag == 1)	
					{
					matchflag = 1;
					break;
					}

	}

	if(matchflag == 1)	
		{
		//add this match (gate type, coveredVertices, inputVertices, cost) to the temporary list of matches;
		for(int r=0;r<10;r++)	currentMatch.gate_type[r] = gates.at(i).name[r];
		for(int s=0;s<coveredVertices.size();s++)	currentMatch.coveredVertices.push_back(coveredVertices.at(s));
		for(int t=0;t<inputVertices.size();t++)	currentMatch.inputVertices.push_back(inputVertices.at(t));
		currentMatch.cost = gates.at(i).cost;
		tempMatches.push_back(currentMatch);

		
		std::cout << " One match for the node " << index[v] << " is " << currentMatch.gate_type << "\n" ;
		std::cout << " The inputs for this match are nodes " ;
		for(int q=0;q<currentMatch.inputVertices.size();q++)	std::cout << currentMatch.inputVertices.at(q) << " " ;
		std::cout << "\n" ;
		
		}

		//clear the coveredVertices list;
		coveredVertices.clear();		
		//clear the inputVertices list;
		inputVertices.clear();		
			
	}

}

void optimalCovering(Graph& g, const Vertex& v) 			// This is the dynamic algorithm that calculates minimum cost cover
{
	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);	
	typedef graph_traits<Graph> GraphTraits;
	
	graph_traits<Graph>::out_edge_iterator out_i, out_end;
	graph_traits<Graph>::edge_descriptor e;     
	tie(out_i, out_end) = out_edges(v,g);			

	while(out_i!=out_end) 							// for each input of node
		{
		Vertex targ;
		e = *out_i;
		targ = target(e, g);
		optimalCovering(g, targ);					// satisfies recurs. assumption : principle of optimal substructure
		out_i++;
		}

	if(out_i == out_end)	
		{

		if(get(get(&NodeProperty::nattr, g), v) == "i")				// this is an input node
			{
			
			char* pin = "input_pin"; 
			for(int i=0;i<10;i++)	matches[index[v]].gate_type[i] = pin[i] ;
			
			matches[index[v]].cost = 0 ;
			matches[index[v]].coveredVertices.push_back(index[v]);
			matches[index[v]].inputVertices.push_back(index[v]);
			}

		else {									// this is a normal node 

			tempMatches.clear();			
			findMatch(g, v);						// this puts all matches in tempMatches
			// Now find the best matching
			int minindex = 0;
			int mincost = 100;						// more than maximum possible cost of a matching
			int tempcost;
			matching_t currentMatch;
			for(int j=0;j<tempMatches.size();j++)		
					{
					tempcost = tempMatches.at(j).cost;
					// add the costs of all inputs
					for(int k=0;k<tempMatches.at(j).inputVertices.size();k++)
						{
						tempcost += matches[tempMatches.at(j).inputVertices.at(k)].cost ;
						}
					if(tempcost < mincost)	
						{
						minindex = j;
						mincost = tempcost;
						}
					}
			// write the gate type
			for(int x=0;x<10;x++)	matches[index[v]].gate_type[x] = tempMatches.at(minindex).gate_type[x];

			// Now we have the best match, write this to the matches array
			matches[index[v]].cost = tempcost ;
			// add the covered vertices
			for(int m=0;m<tempMatches.at(minindex).coveredVertices.size();m++)
					matches[index[v]].coveredVertices.push_back(tempMatches.at(minindex).coveredVertices.at(m));

			// add the input vertices
			for(int n=0;n<tempMatches.at(minindex).inputVertices.size();n++)
					matches[index[v]].inputVertices.push_back(tempMatches.at(minindex).inputVertices.at(n));
			
			}

std::cout << " The match for node " << index[v] << " is " << matches[index[v]].gate_type << " and its cost is " << matches[index[v]].cost << std::endl; 
std::cout << std::endl;
		}

// all the matches have been written 

}



void readGraphMLFile ( Graph& designG, string fileName )
{ 
  boost::dynamic_properties dp;							// dp reads attributes, which are the logic functions at nodes
  dp.property("name", get(&NodeProperty::nattr, designG));			// nattr is the property name, get() gives the property map

  std::ifstream gmlStream;
  gmlStream.open(fileName.c_str(), std::ifstream::in);
 
  boost::read_graphml(gmlStream, designG, dp);
  gmlStream.close();
}

int main(int argc, char* argv[])
{
	
//////////////////// Read graphml file of subject DAG, assumed to be already optimized ////////////////////
	std::ifstream inFile;
	Graph g;
	//std::vector<lib_gate> gates;						// a sequence container of gates
	lib_gate tempgate;
	
	readGraphMLFile(g,"simple.graphml");					// making the graph g from graphml file, which is a tree
	
	typedef property_map<Graph, vertex_index_t>::type IndexMap;		// this is the property map type for the property vertex_index_t
	IndexMap index = get(vertex_index, g);					// get function to get a property map

	//std::cout << "vertices(g) = ";
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;		// this is the vertex_iterator type
	std::pair<vertex_iter, vertex_iter> vp;					// for pair of iterators that are returned by vertices(), avlbl in STL
	for (vp = vertices(g); vp.first != vp.second; ++vp.first)		// vertices(g) returns the vertices
	std::cout << index[*vp.first] <<  " " << 
	get(get(&NodeProperty::nattr, g), *vp.first) << " \n";			// print out indices and properties
	std::cout << std::endl;							// always do this stage as here properties are getting read

/////////////////// Read graphml files of all graph patterns for each gate in the library /////////////////////
// manually write the graphml files for the implementation
	
	inFile.open("library.txt", std::ifstream::in);
	while(inFile.good()) {
	inFile.getline(line,100,'\t');

		if(!strncmp(line,"\nNbr_gates",10))
			{ 
			inFile.getline(line,100,'\t');
			//std::cout << line << std::endl;
			nbr_library_gates = atoi(line);
			//std::cout << "Number of gates is " << nbr_library_gates << std::endl;		
			}

		if(!strncmp(line,"\n\nGATE",6))
			{
			inFile.getline(line,100,'\t');
			inFile.getline(line,100,'\t');
			for(int i=0;i<sizeof(line);i++)	tempgate.name[i] = line[i];
			inFile.getline(line,100,'\t');
			inFile.getline(line,100,'\t');
			tempgate.nbrpatterns = atoi(line);
			
			inFile.getline(line,100,'\t');				// read "Pattern_files"
			for(int j=0;j<(tempgate.nbrpatterns);j++)
				{
			inFile.getline(line,100,'\t');				// read filename
			
			for(int i=0;i<sizeof(line);i++)		tempgate.graphfiles[j][i] = line[i] ;		
			tempgate.gr_list_id[j] = list_id;
			list_id += 10;						// for the next graph pattern
			//std::cout << tempgate.graphfiles[j] << std::endl;	
			}

			inFile.getline(line,100,'\t');				// read "Gate_cost"
			inFile.getline(line,100,'\t');				
			tempgate.cost = atoi(line);
			gates.push_back(tempgate);
			
			}
	}									// file reading complete
	inFile.close();
	
for(int k=0;k<nbr_library_gates;k++){
	
	for(int j=0;j<gates.at(k).nbrpatterns;j++)
				{		
				readGraphMLFile(gates.at(k).gr[j], gates.at(k).graphfiles[j]);		// read the graphml file for the gate patterns
				//index = get(vertex_index, gates.at(k).gr[j]);
				//std::cout << "read " << gates.at(k).graphfiles[j] << std::endl;	
				//std::cout << "vertices(g) of this gate are = \n";
	/*		
	for (vp = vertices(gates[k].gr[j]); vp.first != vp.second; ++vp.first)	
	std::cout << index[*vp.first] <<  " " << get(get(&NodeProperty::nattr, gates[k].gr[j]), *vp.first) << " \n";	
	// printed out indices and properties
	std::cout << std::endl;*/
				}
	}

/*        CODE DISCUSSION
 Find optimal Area Covering
// the final output is to be given as an optimum 'covering' of the output node, which is a 'list' of gates and the total cost
// a gate consists of a type, the input nodes, output node, covered nodes and cost
// the covering is to be gradually built

// ** each node will have its own optimum covering

// for each node, find all 'match'es (same as type 'covering' but with one gate only)
// for each match, the total cost is sum of its own cost and the costs of all the optimum coverings at its inputs
// choose the match with minimum total cost and create the optimum covering for the node

// EDIT : covering also will have only one gate with its inputs, where each input has its own covering, thus avoid redundancy, and its own cost

// the problem is finding the match for a node
// each library gate has a 'complete' list of linked lists (strings), if all the strings lie in the subject graph starting at a node, then we have a 
// match

// create the lists. starting with leaves
// each list will have a matched indicator (bool)
// store the list of linked lists for the subject graph and also the library gates

// Use an array of linked lists (can later convert into a doubly linked list)
// Assumption : each graph pattern has atmost 10 linked lists, so a gate having n patterns can have atmost 10n linked lists
*/


///////////  The next block will extract the linked lists of all the graph patterns /////////////
for(int i=0;i<100;i++)	lists[i] = NULL;

// every graph will have a list_id, the linked lists in the 'lists' array with indices (list_id) to (list_id+9) 
// belong to that particular graph
// ** list_id for the subject graph will be 0 always, so it owns indices 0-9

listcount=0;									
std::for_each(vertices(g).first, vertices(g).second, extract_lists<Graph>(g));		// extract all the linked lists of subject graph

for(int i=0;i<nbr_library_gates;i++)
{
	for(int j=0;j<gates.at(i).nbrpatterns;j++)
	{
		listcount = gates.at(i).gr_list_id[j];
		std::for_each(vertices(gates.at(i).gr[j]).first, vertices(gates.at(i).gr[j]).second, extract_lists<Graph>(gates.at(i).gr[j]));
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////



///////////// this block will just print out all the linked lists (strings)

list_node *travelnode;
for(int i=0;i<listcount;i++)
{
	travelnode = lists[i];
	std::cout << "string " << i << " is " ;
	while(travelnode!=NULL)
	{
		std::cout << travelnode -> nodeattr ;
		travelnode = travelnode -> next;
	}
	std::cout << std::endl;	
}

//////////////////////////////////////////////////////////////////////////

optimalCovering(g, *vertices(g).first);		//find optimal covering of the graph g starting at first node, i.e. output
return 0;
}
