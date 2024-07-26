# Open Street Maps Navigation Project

## Overview
This project implements a console-based navigation system using OpenStreetMap data, focusing on the UIC East campus. It demonstrates the practical application of graph theory, data structures, and algorithms in a real-world scenario.

## Key Features
- Custom Graph implementation using adjacency list
- Parsing of OpenStreetMap XML data
- Building navigation between campus locations
- Implementation of Dijkstra's algorithm for shortest path finding
- Efficient nearest node and building search algorithms

## Skills Demonstrated
- Advanced C++ programming
- Graph theory and implementation
- XML parsing
- Geospatial data handling
- Algorithm design and optimization
- Memory management in C++

## Getting Started
1. Clone the repository
2. Compile the project: `make build`
3. Run the program: `make run`

## Usage
The program prompts for two building names or abbreviations on the UIC East campus. It then:
1. Locates the buildings on the map
2. Finds an optimal meeting point
3. Calculates and displays the shortest paths from both starting points to the meeting point

## Project Structure
- `graph.h`: Custom Graph class implementation
- `application.cpp`: Main application logic
- `osm.cpp` & `osm.h`: OpenStreetMap data parsing
- `dist.cpp`: Distance calculation utilities

## Future Improvements
- Implement a graphical user interface
- Extend to cover a larger geographical area
- Optimize for larger datasets
- Add multi-modal transportation options

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments
- OpenStreetMap for providing the map data
- TinyXML2 library for XML parsing
