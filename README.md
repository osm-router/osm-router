# osm-router [![Build Status](https://travis-ci.org/osm-router/osm-router.svg?branch=master)](https://travis-ci.org/osm-router/osm-router)

Probabilistic OpenStreetMap router based on C++ boost::graph. Uses the following hierarchically related header files:

1. 'xml-parser' Returns a query from the [overpass API](http://overpass-api.de) in OSM XML format, dumps to a specified file, and extracts all
   ways and associated nodes from file (stored in class `Xml`). Alternatively, file name may be entered to avoid repeated downloads, in which
   case 'xml-parser' merely extracts ways and nodes.
2. 'Graph' Derived from class `Xml`. Stores the ways and nodes from 'xml-parser' as a `boost::graph`. Can store either a full graph
   representing all data, or the 'compressed' graph containing only junction nodes.
3. 'Router-test' Derived from class `Graph`. Performs routing query from a given origin point and returns distances to all destination points
   within the graph (values <1 for non-reachable nodes).

('.c++' files are included in each case for stand-alone compilation.)

### build:
1. cd ./build  
2. cmake ..  
3. make

------

# To do ...

1. ~~Finish [`osmdatar`](https://github.com/osmdatar/osmdatar/issues/3) (**M**)~~

2. ~~Convert `osm-router` to R package (using current 
    [Router::dijkstra](https://github.com/osm-router/osm-router/blob/master/src/Router-test.h))~~

    a. ~~Use [`osmdatar`](https://github.com/osmdatar/osmdatar) to obtain all data as `R` `sp` structures (**M**)~~

    b. ~~Rewrite the `Graph::makeCompactGraph ()` function in 
    [`src/Graph.h`](https://github.com/osm-router/osm-router/blob/master/src/Graph.h) to read and return `SpatialLinesDataFrame` objects
    directly from `R` (**A**)~~

    c. ~~Build the routing files to pass a `SpatialLinesDataFrame` from `R` back to `Rcpp` (**A**)~~

3. Extend the `R` package to include a vignette and some unit tests (**A**)

4. Incorporate probabilistic routing functionality

    a. Modify [`stochastic-sp.h`](https://github.com/osm-router/osm-router/blob/master/src/stochastic-sp.h) to the first algorithm of 
    [Saerens et al (2009)](http://www.mitpressjournals.org/doi/abs/10.1162/neco.2009.11-07-643) 
    ([download here](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.62.6428&rep=rep1&type=pdf)) (**M**)

    b. Potentially include [this algorithm](http://theory.stanford.edu/~tim/papers/sssr.pdf) too (**M**)

    c. Integrate the different routing algorithms within the `R` package (**A**)

5. Actually do something useful, productive, and interesting withing the `R` package ... (all **A**)

    a. Visualisation

    b. Integrate over all routes to get expected, rather than shortest, distance

6. Apply to case studies to determine such things as ... (all **A**)

    a. Relationship between entropy and actual distances

    b. Performance? (In comparison with standard dijkstra?)

    c. ... ?
