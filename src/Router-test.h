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

const float FLOAT_MAX = std::numeric_limits <float>::max ();

struct Bbox
{
    float lonmin, latmin, lonmax, latmax;
};

Bbox get_bbox (float xfrom, float yfrom, float xto, float yto, float expand=0.1);

class Router: public Graph
{
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::directedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        bool _compactGraph;
    protected:
        float _xfrom, _yfrom, _xto, _yto;
    public:
        int err;
        int from_node, to_node;
        std::vector <std::pair <int, float>> dists; // returned from dijkstra

    Router (std::string xml_file, std::string profile_file,
            float xfrom, float yfrom, float xto, float yto,
            float lonmin, float latmin, float lonmax, float latmax,
            bool compact)
        : _xfrom (xfrom), _yfrom (yfrom), _xto (xto), _yto (yto),
            Graph (compact, xml_file, profile_file, 
                lonmin, latmin, lonmax, latmax)
    {
        // neaestNode returns the long long OSM ID, and nodeIndx holds the int
        // index of each ID into the graph.
        from_node = nodeIndx [nearestNode (_xfrom, _yfrom)];
        to_node = nodeIndx [nearestNode (_xto, _yto)];
        err = dijkstra (from_node);
        // err = get_dists (); // demo only
    }
    ~Router ()
    {
    }

    long long nearestNode (float lon0, float lat0);
    int dijkstra (int fromNode);
    int get_dists ();
};

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         ROUTER::NEARESTNODE                        **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


long long Router::nearestNode (float lon0, float lat0)
{
    int id;
    long long node = -INT_MAX;
    float d, lon, lat, xd, yd, dmin = FLOAT_MAX;

    /*
     * This is the only place in the whole program that requires explicit
     * looping over vertices, and takes longer than all other bits. Calculating
     * distances to nearest vertices as sums of the 2 abs distances in lats &
     * lons saves only about 2.5s.
     *
     * nodeNames lists the highway vertices, which must be matched to the full
     * list of allNodes to find corresponding lats&lons.
     */

    boost::property_map< Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gr);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lat = boost::get(&bundled_vertex_type::lat, gr);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lon = boost::get(&bundled_vertex_type::lon, gr);
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gr);

    auto vs = boost::vertices (gr);
    for (auto vit = vs.first; vit != vs.second; ++vit)
    {
        // See Graph::get_connected for why next line == 0
        if (vertex_component [*vit] == 0)
        {
            lat = vertex_lat [*vit];
            lon = vertex_lon [*vit];

            //d = std::abs (lon - lon0) + std::abs (lat - lat0);

            xd = (lon - lon0) * PI / 180.0;
            yd = (lat - lat0) * PI / 180.0;
            d = sin (yd / 2.0) * sin (yd / 2.0) + cos (lat * PI / 180.0) *
                cos (lat0 * PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
            d = 2.0 * atan2 (sqrt (d), sqrt (1.0 - d));
            d = d * 6371.0;

            if (d < dmin)
            {
                dmin = d;
                node = vertex_id [*vit];
            }
        }
    }

    return node;
} // end function calcDist



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::DIJKSTRA                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::dijkstra (int fromNode)
{
    std::vector <Vertex> predecessors (boost::num_vertices (gr)); 
    std::vector <Weight> distances (boost::num_vertices (gr)); 

    boost::property_map < Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gr);
    boost::property_map < Graph_t, float bundled_vertex_type::* >::type 
        vertex_lat = boost::get(&bundled_vertex_type::lat, gr);
    boost::property_map < Graph_t, float bundled_vertex_type::* >::type 
        vertex_lon = boost::get(&bundled_vertex_type::lon, gr);

    auto p_map = boost::make_iterator_property_map
        (&predecessors[0], boost::get(boost::vertex_index, gr));
    auto d_map = boost::make_iterator_property_map
        (&distances[0], boost::get(boost::vertex_index, gr));
    auto w_map = boost::get(&bundled_edge_type::weight, gr); 

    int start = vertex (fromNode, gr);
    boost::dijkstra_shortest_paths (gr, start,
            weight_map (w_map). 
            predecessor_map (p_map).
            distance_map (d_map));

    // Check that all routing points have been reached
    int nvalid = 0;
    for (std::vector <Weight>::iterator itr = distances.begin();
            itr != distances.end (); itr++)
        if ((*itr) < FLOAT_MAX)
            nvalid++;
    assert (nvalid == graph_size);

    // Trace back from each routing point 
    Vertex v;
    dists.resize (0);
    float dist;
    boost::property_map < Graph_t, int bundled_vertex_type::* >::type
        vertex_component = boost::get (&bundled_vertex_type::component, gr);

    auto vs = boost::vertices (gr);
    for (auto vit = vs.first; vit != vs.second; ++vit)
    {
        if (vertex_component [*vit] == 0)
        {
            v = (*vit);
            dist = 0.0;
            for (Vertex u = p_map[v]; u != v; v = u, u = p_map[v]) 
            {
                std::pair<Graph_t::edge_descriptor, bool> edgePair = 
                    boost::edge(u, v, gr);
                Graph_t::edge_descriptor edge = edgePair.first;
                dist += boost::get (&bundled_edge_type::dist, gr, edge);
            }
            dists.push_back (std::make_pair ((*vit), dist));
        }
        else
            dists.push_back (std::make_pair ((*vit), -1.0));
    }
    assert (dists.size () == num_vertices (gr));

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::GET_DISTS                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::get_dists ()
{
    assert (dists.size () > 0); // Asserts that dijkstra has been run.

    auto vs = dists;
    for (auto it = dists.begin (); it != dists.end (); ++it)
    {
        std::cout << "(" << (*it).first << ", " << (*it).second <<
            ")" << std::endl;
    }

    return 0;
}
