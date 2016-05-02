# osm-router [![Build Status](https://travis-ci.org/osm-router/osm-router.svg?branch=master)](https://travis-ci.org/osm-router/osm-router)

Probabilistic OpenStreetMap router based on C++ boost::graph. Uses the following hierarchically related header files:

1. 'xml-parser' Returns a query from the [overpass API](http://overpass-api.de) in OSM XML format, dumps to a specified file, and extracts all
   ways and associated nodes from file (stored in class `Xml`). Alternatively, file name may be entered to avoid repeated downloads, in which
   case 'xml-parser' merely extracts ways and nodes.
2. 'Graph' Derived from class `Xml`. Stores the ways and nodes from 'xml-parser' as a `boost::graph`. Can store either a full graph
   representing all data, or the 'compressed' graph containing only junction nodes.
3. 'Router-test' Derived from class `Graph`. Performs routing query from a given origin point and returns distances to all destination points
   within the graph (values <1 for non-reachable nodes).

### build:
1. cd ./build  
2. cmake ..  
3. make
