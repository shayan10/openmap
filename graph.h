// graph.h <Starter Code>
// Mohammad Shayan Khan
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_map>
#include <list>

using namespace std;


template <typename VertexT, typename WeightT>
class graph
{
private:
	typedef pair<WeightT, VertexT> edge;
	typedef typename list<edge>::iterator edgeItr;
	typedef typename list<edge>::const_iterator constEdgeItr;
	// Store adjancey 
	unordered_map<VertexT, list<edge>> adj;
	size_t numVerticies;
	size_t numEdges;

	/*
		Purpose: Traverses the edge list for the 'from' vertex
			and returns an iterator to the matching entry in the
			edge list.

		Parameters:
			from - starting vertex
			to - ending vertex
		
		Returns:
			iterator to corresponding entry in the edge list
			If edge not present, returns the end() iterator
			Throws exception if either vertex is not in the graph
	*/
	edgeItr lookupVertex(VertexT from, VertexT to) {
		// If the vertex does not exist, throw exception
		if (adj.count(from) == 0 || adj.count(to) == 0) {
			throw invalid_argument("Vertex does not exist");
		}

		// Search from the start of the edge list
		list<edge>& edgeList = adj.at(from); 
		for (edgeItr itr = edgeList.begin(); itr != edgeList.end(); itr++) {
			// If second parameter of pair matches, return the iterator
			if (itr->second == to) {
				return itr;
			}
		}
		// Otherwise, return the end position
		return edgeList.end();
	}

	/*
		Purpose: Outputs all the edges for the given graph in the
			format (from, weight, from) in a comma-separated format
		Parameters:
			output - reference to the destination output stream
			v - vertex of whose edge list to print
		Returns:
			N/A
	*/
	void printEdgeList(ostream& output, const VertexT v) const {
		// If vertex does not exist
		if (adj.count(v) == 0) {
			return;
		}
		// Otherwise, print any edges
		const list<edge>& edgeList = adj.at(v);
		
		// Print vertex
		output << v << ": ";
		// Print edges
		for (constEdgeItr itr = edgeList.begin(); itr != edgeList.end(); itr++) {
			// Output the edge
			output << "(" << v << ", " << itr->second << ", " << itr->first << ")";
			// If not the last one, output comma
			if (next(itr) != edgeList.end()) {
				output << ", ";
			}
		}
		// End with new line
		output << endl;
	}


public:
	// Default Constructor
	// Initializes adjacency list and set numVerticies to 0
	graph(): adj(), numVerticies(0), numEdges(0) {}

	//
	// NumVertices
	//
	// Returns the # of vertices currently in the graph.
	//
	int NumVertices() const
	{
		return static_cast<int>(this->numVerticies);
	}

	//
	// NumEdges
	//
	// Returns the # of edges currently in the graph.
	//
	int NumEdges() const
	{
		return static_cast<int>(this->numEdges);
	}

	//
	// addVertex
	//
	// Adds the vertex v to the graph if there's room, and if so
	// returns true.  If the graph is full, or the vertex already
	// exists in the graph, then false is returned.
	//
	bool addVertex(VertexT v)
	{
		if (adj.count(v) != 0) {
			return false;
		}
		// Set edge list for vertex v to empty list
		adj[v] = {};
		// Update number of vertices
		numVerticies++;
		// New vertex added
		return true;
	}

	//
	// addEdge
	//
	// Adds the edge (from, to, weight) to the graph, and returns
	// true.  If the vertices do not exist or for some reason the
	// graph is full, false is returned.
	//
	// NOTE: if the edge already exists, the existing edge weight
	// is overwritten with the new edge weight.
	//
	bool addEdge(VertexT from, VertexT to, WeightT weight)
	{
		// If either the from or to vertex is not in the adj list,
		// return false
		if (adj.count(from) == 0 || adj.count(to) == 0) {
			return false;
		}

		edgeItr existingEdge = lookupVertex(from, to);

		// If edge does not exist
		if (existingEdge == adj.at(from).end()) {
			// Add edges between from and to vertex in both directions
			adj[from].push_back({weight, to});			
			// Update number of edges
			numEdges += 1;
			return true;
		}

		// Otherwise, overwrite weight for existing edge
		existingEdge->first = weight;
		// // Also update edge in other direction
		// edgeItr reverseEdge = lookupVertex(to, from);
		// reverseEdge->first = weight;
		return true;

	}

	//
	// getWeight
	//
	// Returns the weight associated with a given edge.  If
	// the edge exists, the weight is returned via the reference
	// parameter and true is returned.  If the edge does not
	// exist, the weight parameter is unchanged and false is
	// returned.
	//
	bool getWeight(VertexT from, VertexT to, WeightT &weight) const {
		
		// If either vertex not in adj, return false
		if (adj.count(from) == 0 || adj.count(to) == 0) {
			return false;
		}

		typename std::list<edge>::const_iterator itr;
		for (itr = adj.at(from).begin(); itr != adj.at(from).end(); itr++) {
			// If to vertex found
			if (itr->second == to) {
				// Return weight in reference parameter
				weight = itr->first; 
				return true;
			}
		}
		// Otherwise, vertex not found, return false
		return false;
	}

	//
	// neighbors
	//
	// Returns a set containing the neighbors of v, i.e. all
	// vertices that can be reached from v along one edge.
	// Since a set is returned, the neighbors are returned in
	// sorted order; use foreach to iterate through the set.
	//
	set<VertexT> neighbors(VertexT v) const
	{
		// If vertex does not exist, return empty set
		if (adj.count(v) == 0) {
			return set<VertexT>();
		}

		set<VertexT> S;

		for (const auto& edge : adj.at(v)) {
			S.insert(edge.second);
		}
		return S;
	}

	//
	// getVertices
	//
	// Returns a vector containing all the vertices currently in
	// the graph.
	//
	vector<VertexT> getVertices() const
	{
		vector<VertexT> verticies;
		for (const auto& v: adj) {
			verticies.push_back(v.first);
		}
		return verticies;
	}

	//
	// dump
	//
	// Dumps the internal state of the graph for debugging purposes.
	//
	// Example:
	//    graph<string,int>  G(26);
	//    ...
	//    G.dump(cout);  // dump to console
	//
	void dump(ostream &output) const
	{
		output << "***************************************************" << endl;
		output << "********************* GRAPH ***********************" << endl;

		output << "**Num vertices: " << this->NumVertices() << endl;
		output << "**Num edges: " << this->NumEdges() << endl;

		output << endl;
		output << "**Vertices:" << endl;
		for (const auto & entry: adj) {
			cout << " " << entry.first;
		}
		output << endl;
		output << "**Edges:" << endl;
		for (const auto & entry: adj) {
			// Print the edge list for every vertex
			printEdgeList(output, entry.first);
		}
		
		output << "**************************************************" << endl;
	}
};
