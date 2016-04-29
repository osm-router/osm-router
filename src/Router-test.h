/***************************************************************************
 *  Project:    osm-router
 *  File:       Router.h
 *  Language:   C++
 *
 *  osm-router is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  osm-router is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  osm-router.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Author:     Mark Padgham / Andreas Petutschnig
 *  E-Mail:     mark.padgham@email.com / andras@petutschnig.de
 *
 *  Description:    C++ implementation of OSM router using boost::graph.
 *                  Designed to work in a designated area, and so reads data
 *                  from a planet.osm file. Hard-coded at present to read data
 *                  for greater London and greater NYC, and to route between
 *                  points given in ./data/routing-points-(city).txt using the
 *                  profile given in profile.cfg
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/


/*
 * Stores a planet-osm file as a boost::graph. The latter rely on property maps,
 * which are statically-typed entities, and so data are initially read and
 * stored as unordered_maps to be subsequently read into the graph. This initial
 * storage as unordered_maps also speeds up mapping node IDs in ways onto
 * lat-lon coordinates.
 *
 * The boost::graphs are constructed in two phases:
 * (1) A graph "gFull" is constructed from all OSM nodes and highways, and used
 * to identify both the largest connected component and all terminal or junction
 * nodes (simply as those that appear in multiple different ways). The routing
 * points are then mapped onto the nearest vertices within this largest
 * component, and these vertices are also added to terminalNodes.  (2) A second
 * reading of the data is used to make "gCompact" which only has terminalNodes
 * as vertices, and edge distances as traced along all intermediate nodes.
 */

#include "Graph.h"

#include <boost/unordered_map.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

struct RoutingPoint
{
    float lon, lat;
    long long node;
    int nodeIndex;
};

// edge and vertex bundles for the boost::graph
struct bundled_edge_type
{ 
    Weight weight; // weighted distance
    float dist;
};

struct bundled_vertex_type
{
    long long id;
    std::string name;
    float lat, lon;
    int component;
};


class Router: public Graph
{
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::directedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        std::string profileName;
    protected:
        float latmin, lonmin, latmax, lonmax;
        std::vector <ProfilePair> profile;

    public:
        std::vector <RoutingPoint> RoutingPointsList;

    Router (std::string xml_file, std::string profile_file,
            float lonmin, float latmin, float lonmax, float latmax)
        : Graph (xml_file, profile_file, lonmin, latmin, lonmax, latmax)
    {
        /*
        err = getBBox ();
        err = readRoutingPoints ();

        err = remapRoutingPoints ();

        std::cout << "Getting distances between routing points";
        std::cout.flush ();
        count = 0;
        for (std::vector<RoutingPoint>::iterator itr=RoutingPointsList.begin();
                itr != RoutingPointsList.end(); itr++)
        {
            err = dijkstra (itr->nodeIndex);
            assert (dists.size () == RoutingPointsList.size ());
            std::cout << "\rGetting distances between routing points " <<
                count << "/" << RoutingPointsList.size () << " ";
            std::cout.flush ();
            count++;
        }
        std::cout << "\rGetting distances between routing points " <<
            RoutingPointsList.size () << "/" << RoutingPointsList.size
            () << std::endl;
        std::cout << "done." << std::endl;
        */
    }
    ~Router ()
    {
    }

    int dijkstra (long long fromNode);
    long long nearestNode (float lon0, float lat0);
};
