# osm-router [![Build Status](https://travis-ci.org/mpadge/osm-router.svg?branch=master)](https://travis-ci.org/mpadge/osm-router)

C++ implementation of OSM router using boost::graph. Designed to work in a designated area, and so reads data from a planet.osm file. Hard-coded
at present to read data for greater London and greater NYC, and to route between points given in ./data/routing-points-(city).txt.

### build:
1. cd ./build  
2. cmake ..  
3. make
