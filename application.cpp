// application.cpp <Starter Code>
// Mohammad Shayan Khan (UIN: 667707825)
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
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip> /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <functional>

#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;

/*
	Purpose: Search through the buildings vector based on either the full
		name or the abbreviation.
	Parameters:
		buildings - reference to the buildings vector
		searchTerm - abbreviation or name parameter
		searchByAbbreviation - flag to indicate if search should consider the abbreviation
			or the full name
	Returns:
		-1 if building not found
		Otherwise, returns the index of the building in the vector
*/
int searchBuildings(const vector<BuildingInfo> &buildings, string searchTerm,
	bool searchByAbbreviation = true)
{

	for (int i = 0; i < static_cast<int>(buildings.size()); i++)
	{
		// If not searching by abbreviation and partial match found, return index
		if (!searchByAbbreviation)
		{
			if (buildings[i].Fullname.find(searchTerm) != string::npos)
			{
				return i;
			}
		}

		// Otherwise, if abbrevation found, return index
		if (buildings[i].Abbrev == searchTerm)
		{
			return i;
		}
	}

	// If not found, return -1
	return -1;
}

/*
	Purpose: Wrapper function which searches the buildings vector based
		on either the abbreviation or the full name.
	Parameters:
		searchTerm - user-provided string to search the buildings vector with
		buildings - constant reference to the buildings vector
	Returns:
		-1 if not found using the abbreviation or full name
		Otherwise, returns the index of the building in the vector
*/
int findBuilding(const vector<BuildingInfo> &buildings, string searchTerm)
{
	int searchByAbbrev = searchBuildings(buildings, searchTerm);
	// If abbreviation not found, return result of partial match,
	// otherwise return the original result
	return searchByAbbrev == -1 ? searchBuildings(buildings, searchTerm, false) : searchByAbbrev;
}

/*
	Purpose: Helper function to print a building instance in the
		following format to cout:
			name
			(lat, long)
	Parameters:
		building - constant reference to building object to be
			printed
	Returns:
		N/A (outputs to console)
*/
void printBuilding(const BuildingInfo& building) {
	// " (" << building.Abbrev << ")" << 
	cout << " " << building.Fullname << endl;
	cout << " (" << building.Coords.Lat << ", " << building.Coords.Lon << ")" << endl; 
}

/*
	Purpose: Helper function to print a node in the
		following format to cout:
			nodeId
			(lat, long)
	Parameters:
		Nodes - constant reference to the map containing all nodes to be
			considered
		nodeId - integer representing the unique ID of the node
	Returns:
		N/A (outputs to console)
*/
void printNode(const map<long long, Coordinates> & Nodes, long long& nodeId) {
	cout << " " << nodeId << endl;
	cout << " " << "(" << Nodes.at(nodeId).Lat << ", " << Nodes.at(nodeId).Lon << ")" << endl; 
}

/*
	Purpose: Returns a vector containing the
		predecessors for the path from 'start'
		to 'end' in the order 'end', ..., 'start'
	Parameters:
		pred - constant reference to map which maps each vertex
			to its predecessor node in the path (if any)
		start - starting vertex of the path
		end - ending vertex of the path
	Returns:
		vector containing the predeccesors from end -> start
*/
vector<long long> orderPredecessors(const map<long long, long long> & pred, const long long start, 
	const long long end) {

	vector<long long> result;
	long long currNode = end;
	// Until start vertex not reached
	while (currNode != -1) {
		// Add current node to result
		result.push_back(currNode);
		// Process predeccesor of current node
		currNode = pred.at(currNode);
	}
	// Return result
	return result;
}

/*
	Purpose: Outputs a given path from the start
		vertex to the end vertex in the following
		format, given a map of predeccesors for
		each vertex:
			start -> ... -> intermediate node -> ... -> end
	Parameters:
		pred - constant reference to map which maps each vertex
			to its predecessor node in the path (if any)
		start - starting vertex of the path
		end - ending vertex of the path
	Returns:
		N/A (outputs to console)
*/
void printPredecessors(const map<long long, long long> & pred, const long long start, 
	const long long end) {
	// Order the predecessors in reverse order
	vector<long long> orderedPredecessors = orderPredecessors(pred, start, end);
	int n = orderedPredecessors.size();
	// Traverse array in reverse order
	for (int i = n - 1; i >= 0; i--) {
		// Output current element
		cout << orderedPredecessors[i];
		// If not the last index, output arrow
		if (i != 0) {
			cout << "->";
		}
	}
}

/*
	Purpose: Helper function to return the closest node to a building
		amongst all possible footway nodes.
	Parameters:
		Footways - vector containing all footways for the given map
		Nodes - vector containing all possible nodes for the given map
		building - reference to struct object for the given building
	Returns:
		id of the node which is closest to the building
*/
long long findClosestNode(const vector<FootwayInfo>& Footways, const map<long long, Coordinates> & Nodes,
	const BuildingInfo& building) {
	// Initialize min distance to infinity
	double minDistance = numeric_limits<double>::max(), currDistance;
	const Coordinates* currNodeCoordinates = nullptr;
	long long result;
	// For every footway
	for (const auto & footway: Footways) {
		for (const long long& nodeId: footway.Nodes) {
			// Retrieve current node
			currNodeCoordinates = &(Nodes.at(nodeId));
			currDistance = distBetween2Points(building.Coords.Lat, building.Coords.Lon,
				currNodeCoordinates->Lat, currNodeCoordinates->Lon);

			// If smaller distance found
			if (currDistance < minDistance) {
				// Update result id
				result = nodeId;
				// Return 
				minDistance = currDistance;
			}
		}
	}
	return result;
}

/*
	Purpose: Determines the shortest path between the start and 
		destination nodes.
	Parameters:
		G - instance of adjancency-list based Graph object for the given
			map
		start, dest - id's for the start and destination nodes
		distances - reference to map which will contain the distance of all
			nodes from the starting node
		predecessors - reference to map which will contain the predecessors
			of all verticies in the shortest path
	Returns:
		boolean indicating if there was a possible path between the start
			and the destination nodes.
*/
bool shortestPath(const graph<long long, double> & G, long long start, long long dest,
	map<long long, double>& distances, map<long long, long long>& predecessors) {
	
	typedef pair<double, long long> EdgeWeightPair;
	const double INF = numeric_limits<double>::max();

	// Retrieve all verticies
	vector<long long> vertices = G.getVertices();

	// Initialize all distances to infinity
	for (const auto & v: vertices) {
		distances[v] = INF;
		predecessors[v] = -1;
	}

	// Initialize starting distance to 0
	distances[start] = 0;
	// Define priority queue to return vertex with the shortest weight
	auto edgeComp = [](const EdgeWeightPair & a, const EdgeWeightPair &b)
	{
		return a.first > b.first;
	};

	// Defining function to be used as the comparator
	function<bool(const EdgeWeightPair &, const EdgeWeightPair & b)> comparator = edgeComp;
	// Initialize priority queue
	priority_queue<EdgeWeightPair, vector<EdgeWeightPair>, decltype(comparator)> pq(comparator);
	set<long long> neighbors;

	// Add the starting vertex to the priority queue
	pq.emplace(make_pair(distances[start], start));
	double currDist, newDist, weight;
	long long currVertex;
	
	while (!pq.empty()) {
		// Retrieve the weight and vertex of edge
		// at the top of the priority queue
		currVertex = pq.top().second;
		currDist = pq.top().first;
		// Pop topmost entry
		pq.pop();


		// If distance remains infinity, terminate algorithm
		// since shorter paths not possible
		if (distances[currVertex] == INF) {
			break;
		}

		// Otherwise, traverse all neighbors and
		// update priority queue and distances
		neighbors = G.neighbors(currVertex);

		// For all neighbors
		for (const auto & v: neighbors) {
			// Calculate new distance
			G.getWeight(currVertex, v, weight);
			newDist = currDist + weight;
			// If smaller than the current distance
			if (newDist < distances[v]) {
				// Update distance for neighbor
				distances[v] = newDist;
				// Push vertex onto the priority queue
				pq.emplace(make_pair(distances[v], v));
				predecessors[v] = currVertex;
			}
		}
	}

	// Return if path possible from start to end
	return distances[dest] != INF;
}

class DistComp {
	private:
		Coordinates* center;
	public:
		DistComp(): center(nullptr) {}
		DistComp(Coordinates* newCenter): center(newCenter) {}
		bool operator()(const BuildingInfo& a, const BuildingInfo & b) {
			// Return if distance between building A and the center is GREATER than the
			// distance between building B and the center
			return distBetween2Points(a.Coords.Lat, a.Coords.Lon, center->Lat, center->Lon) >
				distBetween2Points(b.Coords.Lat, b.Coords.Lon, center->Lat, center->Lon);
		}
};

//
// Implement your standard application here
//
void application(
	 map<long long, Coordinates> &Nodes, vector<FootwayInfo> &Footways,
	 vector<BuildingInfo> &Buildings, graph<long long, double> &G)
{
	string person1Building, person2Building;
	int person1BuildingIndex, person2BuildingIndex;
	Coordinates building1Coords, building2Coords;			
	Coordinates center;
	BuildingInfo minDistBuilding, origDestBuilding;
	long long p1, p2, dest;											// Placeholder to store start, dest node IDs

	DistComp distComp(&center);

	// Map to store the distances for the path p1 -> dest
	map<long long, double> distFromP1;
	// Map to store the distances for the path p2 -> dest
	map<long long, double> distFromP2;

	// Map to store the predecessors for the path p1 -> dest
	map<long long, long long> predP1;
	// Map to store the predecessors for the path p2 -> dest
	map<long long, long long> predP2;
	
	bool done = false, pathPossibleFromP1 = false, pathPossibleFromP2 = false;

	// Until user enters the termination character
	do
	{
		// Declaring priority queue to order buildings by distance to center
		// Note: Needs to be re-initialized on every iteration since the center will change
		// based on user input
		std::priority_queue<BuildingInfo, std::vector<BuildingInfo>, DistComp> minCenterDistPq(distComp);
		
		cout << endl;
		cout << "Enter person 1's building (partial name or abbreviation), or #> ";
		getline(cin, person1Building);

		// If termination character, break out of the loop
		if (person1Building == "#")
		{
			break;
		}

		cout << "Enter person 2's building (partial name or abbreviation)> ";
		getline(cin, person2Building);
		cout << endl;

		// Search for person1 and person2 buildings in the Buildings vector
		person1BuildingIndex = findBuilding(Buildings, person1Building);	
		person2BuildingIndex = findBuilding(Buildings, person2Building);
		
		// If person 1 building not found, get another pair of inputs
		if (person1BuildingIndex == -1)
		{
			cout << "Person 1's building not found" << endl;
			continue;
		}
		
		// If person 2 building not found, get another pair of inputs
		if (person2BuildingIndex == -1)
		{
			cout << "Person 2's building not found" << endl;
			continue;
		}
		
		// Finding the center between the two buildings
		building1Coords = Buildings[person1BuildingIndex].Coords;
		building2Coords = Buildings[person2BuildingIndex].Coords;
		center = centerBetween2Points(building2Coords.Lat, building2Coords.Lon, building1Coords.Lat,
			building1Coords.Lon);

		// Output person 1 & 2 buildings
		cout << "Person 1's point:" << endl;
		printBuilding(Buildings[person1BuildingIndex]);
		cout << "Person 2's point:" << endl;
		printBuilding(Buildings[person2BuildingIndex]);

		// Insert all buildings into priority queue, arranging them by distance to center
		for (const auto& building: Buildings) {
			minCenterDistPq.push(building);
		}

		origDestBuilding = minCenterDistPq.top();

		// Process nodes
		done = false;

		// While possible path not found or queue is not empty
		while (!done && !minCenterDistPq.empty()) {
			// Retrieve current closest building
			minDistBuilding = minCenterDistPq.top();

			// Determine closest node to the closest and person1 & person2 buildings
			p1 = findClosestNode(Footways, Nodes, Buildings[person1BuildingIndex]);
			p2 = findClosestNode(Footways, Nodes, Buildings[person2BuildingIndex]);
			dest = findClosestNode(Footways, Nodes, minDistBuilding);

			// Determine if paths are possible from p1 -> dest and p2 -> dest
			pathPossibleFromP1 = shortestPath(G, p1, dest, distFromP1, predP1);
			pathPossibleFromP2 = shortestPath(G, p2, dest, distFromP2, predP2);
			
			// Pop minimum element off the stack
			minCenterDistPq.pop();

			// If both possible from both verticies, or impossible from both, exit loop
			if ((pathPossibleFromP1 && pathPossibleFromP2) || (!pathPossibleFromP1 && !pathPossibleFromP2))
			{
				done = true;
				break;
			}
		}

		// If both paths possible, output destination building
		// If no possible, output the original building
		if (pathPossibleFromP1 && pathPossibleFromP2) {
			cout << "Destination Building:" << endl;
			printBuilding(minDistBuilding);
		} else {
			cout << "Destination Building:" << endl;
			printBuilding(origDestBuilding);
		}

		// Output nearest nodes
		cout << endl;
		cout << "Nearest P1 node:" << endl;
		printNode(Nodes, p1);
		cout << "Nearest P2 node:" << endl;
		printNode(Nodes, p2);
		cout << "Nearest destination node:" << endl;
		printNode(Nodes, dest);
		cout << endl;

		// If there was a possible path, output predecessors
		// Otherwise, output that the destination is unreachable
		if (pathPossibleFromP1 && pathPossibleFromP2) {
			cout << "Person 1's distance to dest: " << distFromP1[dest] << " miles" << endl;
			cout << "Path: ";
			printPredecessors(predP1, p1, dest);
			cout << endl << endl;
			cout << "Person 2's distance to dest: " << distFromP2[dest] << " miles" << endl;
			cout << "Path: ";
			printPredecessors(predP2, p2, dest);
			cout << endl << endl;
		} else {
			cout << "Sorry, destination unreachable" << endl;
		}
	} while (person1Building != "#");
}

int main()
{
	graph<long long, double> G;

	// maps a Node ID to it's coordinates (lat, lon)
	map<long long, Coordinates> Nodes;
	// info about each footway, in no particular order
	vector<FootwayInfo> Footways;
	// info about each building, in no particular order
	vector<BuildingInfo> Buildings;
	XMLDocument xmldoc;

	cout << "** Navigating UIC open street map **" << endl;
	cout << endl;
	cout << std::setprecision(8);

	string def_filename = "map.osm";
	string filename;

	cout << "Enter map filename> ";
	getline(cin, filename);

	if (filename == "")
	{
		filename = def_filename;
	}

	//
	// Load XML-based map file
	//
	if (!LoadOpenStreetMap(filename, xmldoc))
	{
		cout << "**Error: unable to load open street map." << endl;
		cout << endl;
		return 0;
	}

	//
	// Read the nodes, which are the various known positions on the map:
	//
	int nodeCount = ReadMapNodes(xmldoc, Nodes);

	//
	// Read the footways, which are the walking paths:
	//
	int footwayCount = ReadFootways(xmldoc, Footways);

	//
	// Read the university buildings:
	//
	int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

	//
	// Stats
	//
	assert(nodeCount == (int)Nodes.size());
	assert(footwayCount == (int)Footways.size());
	assert(buildingCount == (int)Buildings.size());

	// Output the number of nodes, footways, and buildings
	cout << endl;
	cout << "# of nodes: " << Nodes.size() << endl;
	cout << "# of footways: " << Footways.size() << endl;
	cout << "# of buildings: " << Buildings.size() << endl;

	// Adding all the vertices
	for (auto &entry : Nodes)
	{
		G.addVertex(entry.first);
	}

	Coordinates coord1, coord2;
	double distance;
	// For every footway
	for (FootwayInfo &footway : Footways)
	{
		// Consider each pair of nodes
		for (size_t l = 0, r = 1; r < footway.Nodes.size(); l++, r++)
		{
			// Retrieve co-ordinates of both nodes
			coord1 = Nodes[footway.Nodes[l]];
			coord2 = Nodes[footway.Nodes[r]];
			// Calculate distance between the two nodes
			distance = distBetween2Points(coord1.Lat, coord1.Lon, coord2.Lat,
													coord2.Lon);
			// Add edge between nodes in both directions
			G.addEdge(footway.Nodes[l], footway.Nodes[r], distance);
			G.addEdge(footway.Nodes[r], footway.Nodes[l], distance);
		}
	}

	// Output number of verticies and edges
	cout << "# of vertices: " << G.NumVertices() << endl;
	cout << "# of edges: " << G.NumEdges() << endl;
	cout << endl;

	// Execute Application
	application(Nodes, Footways, Buildings, G);

	//
	// done:
	//
	cout << "** Done **" << endl;
	return 0;
}
