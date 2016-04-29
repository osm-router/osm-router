/***************************************************************************
 *  Project:    osm-router
 *  File:       Graph-test.h
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
 *  Author:     Mark Padgham 
 *  E-Mail:     mark.padgham@email.com 
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
 * The boost::graphs are constructed in two phases:
 * (1) A graph "gFull" is constructed from all OSM nodes and highways, and used
 * to identify both the largest connected component and all terminal or junction
 * nodes (simply as those that appear in multiple different ways). The routing
 * points are then mapped onto the nearest vertices within this largest
 * component, and these vertices are also added to terminalNodes.  (2) A second
 * reading of the data is used to make "gCompact" which only has terminalNodes
 * as vertices, and edge distances as traced along all intermediate nodes.
 */

#include "xml-parser.h"

#include <boost/unordered_set.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#define PI 3.1415926535897932384626433832795

typedef boost::unordered_map <long long, int> umapInt; // for node int IDs
typedef boost::unordered_map <long long, int>::iterator umapInt_Itr;

typedef std::vector <Way>::iterator Ways_Itr;
typedef std::pair <std::string, float> ProfilePair;
typedef float Weight;

// See http://theboostcpplibraries.com/boost.unordered
std::size_t hash_value(const int &i)
{
    std::size_t seed = 0;
    boost::hash_combine(seed, i);
    return seed;
}

// edge and vertex bundles for the boost::graph
struct bundled_edge_type
{ 
    long long id;
    std::string name, type;
    // key-val pairs not implemented yet
    //std::vector <std::pair <std::string, std::string>> key_val;
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



class Graph: public Xml
{
    // http://stackoverflow.com/questions/18791319/calculate-number-of-in-and-out-edges-in-a-boostgraph-vertex
    // say only "out_degree" is available for directedS, while bidiretionalS
    // enables also "in_degree"
    //using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
    //      boost::directedS, bundled_vertex_type, bundled_edge_type >;
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::bidirectionalS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        Graph_t gFull, gCompact;

    protected:

    public:
        umapInt nodeNames;
        umapPair_Itr node_itr;
        boost::unordered_set <long long> terminalNodes;

    Graph (std::string file, float lonmin, float latmin, float lonmax, float latmax)
        : Xml (file, lonmin, latmin, lonmax, latmax)
    {
        // Construction fills vectors of "ways" and "nodes" from Xml class
        // These then need to be used to fill the boost::graph
        std::cout << "Parsed " << nodes.size () << " nodes and " << 
            ways.size () << " ways" << std::endl;

        nodeNames.clear ();
        makeFullGraph ();
        dumpGraph ();
        makeCompactGraph ();
    }
    ~Graph ()
    {
        nodeNames.clear ();
        terminalNodes.clear ();
    }

    void dumpGraph ();
    void makeFullGraph ();
    void makeCompactGraph ();
    float calcDist (std::vector <float> x, std::vector <float> y);
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::DUMPGRAPH                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::dumpGraph ()
{
    // Actually dumps the raw ways from which the graph is made. Useful for
    // visually inspecting connected components.
    long long ni;
    std::vector <float> lats, lons;
    umapPair_Itr umapitr;
    typedef std::vector <long long>::iterator ll_Itr;
    std::ofstream out_file;

    out_file.open ("junk.txt", std::ofstream::out);

    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        // Set up first origin node
        ni = (*wi).nodes.front ();
        assert ((umapitr = nodes.find (ni)) != nodes.end ());
        lats.resize (0);
        lons.resize (0);
        lons.push_back ((*umapitr).second.first);
        lats.push_back ((*umapitr).second.second);
        assert (nodeNames.find (ni) != nodeNames.end ());

        // Then iterate over the remaining nodes of that way
        for (ll_Itr it = std::next ((*wi).nodes.begin ());
                it != (*wi).nodes.end (); it++)
        {
            assert ((umapitr = nodes.find (*it)) != nodes.end ());
            lons.push_back ((*umapitr).second.first);
            lats.push_back ((*umapitr).second.second);
            assert (nodeNames.find (*it) != nodeNames.end ());

            assert (lons.size () == lats.size ()); // can't ever fail
            if (lons.size () == 2) // then add edge
            {
                out_file << lons [0] << "," << lats [0] << "," <<
                    lons [1] << "," << lats [1] << std::endl;
                lons.erase (lons.begin ());
                lats.erase (lats.begin ());
            }
        }
    }
    out_file.close ();

} // end function Graph::dumpGraph


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                      FUNCTION::MAKEFULLGRAPH                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::makeFullGraph ()
{
    // The full graph has to be constructed so that distances can be accurately
    // calculated for the compact graph. This routine also fills terminalNodes
    // so that compactGraph can subsequently be made
    int tempi [2];
    float lon1, lat1, lon2, lat2;
    std::vector <float> lats, lons;
    long long ni;
    boost::unordered_set <long long> nodeList;
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;
    typedef std::vector <long long>::iterator ll_Itr;

    // Fill vertices. These are automatically numbered consecutively, so
    // nodeNames is also indexed simultaneously to provide references for edge
    // insertion.
    tempi [0] = 0;
    for (umapPair_Itr it = nodes.begin (); it != nodes.end (); ++it)
    {
        oneVert.id = (*it).first;
        oneVert.lon = (*it).second.first;
        oneVert.lat = (*it).second.second;
        boost::add_vertex (oneVert, gFull);
        nodeNames [oneVert.id] = tempi [0]++;
    }

    nodeList.clear ();
    terminalNodes.clear ();
    // First put all first and last nodes of each way into terminalNodes
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        if (terminalNodes.find ((*wi).nodes.front ()) == terminalNodes.end())
            terminalNodes.insert ((*wi).nodes.front ());
        if (terminalNodes.find ((*wi).nodes.back ()) == terminalNodes.end())
            terminalNodes.insert ((*wi).nodes.back ());
    }
    
    // Then fill edges
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        oneEdge.id = (*wi).id;
        oneEdge.name = (*wi).name;
        oneEdge.type = (*wi).type;

        // Set up first origin node
        ni = (*wi).nodes.front ();
        assert ((umapitr = nodes.find (ni)) != nodes.end ());
        lats.resize (0);
        lons.resize (0);
        lons.push_back ((*umapitr).second.first);
        lats.push_back ((*umapitr).second.second);
        assert (nodeNames.find (ni) != nodeNames.end ());
        tempi [0] = (*nodeNames.find (ni)).second; // int index of node ID

        // Then iterate over the remaining nodes of that way
        for (ll_Itr it = std::next ((*wi).nodes.begin ());
                it != (*wi).nodes.end (); it++)
        {
            assert ((umapitr = nodes.find (*it)) != nodes.end ());
            lons.push_back ((*umapitr).second.first);
            lats.push_back ((*umapitr).second.second);
            assert (nodeNames.find (*it) != nodeNames.end ());
            tempi [1] = (*nodeNames.find (*it)).second;

            if (terminalNodes.find (*it) == terminalNodes.end () &&
                    nodeList.find (*it) == nodeList.end ())
            {
                if (nodeList.find (*it) == nodeList.end ())
                    nodeList.insert (*it);
                else
                    terminalNodes.insert (*it);
            }

            assert (lons.size () == lats.size ()); // can't ever fail
            if (lons.size () == 2) // then add edge
            {
                oneEdge.dist = calcDist (lons, lats);
                boost::add_edge (tempi [0], tempi [1], oneEdge, gFull);
                tempi [0] = tempi [1];
                lons.erase (lons.begin ());
                lats.erase (lats.begin ());
            }
        }
    }

    tempi [0] = terminalNodes.size () + nodeList.size ();
    std::cout << "There are " << terminalNodes.size () <<
        " terminal nodes and " << nodeList.size () << " non-terminal" << 
        "; total = " << tempi [0] << std::endl;
    // This total should be less than the number of vertices in gFull because
    // it only includes nodes which actually arise within ways, while gFull is
    // created from all nodes returned from the overpass query regardless of
    // whether they are part of ways or not.

    std::vector <int> compvec (num_vertices (gFull));
    tempi [0] = boost::connected_components (gFull, &compvec [0]);
    std::cout << "Graph has " << num_vertices (gFull) << " vertices and " <<
        tempi [0] << " connected components" << std::endl;

    nodeList.clear ();
} // end function Graph::makeFullGraph


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     FUNCTION::MAKECOMPACTGRAPH                     **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::makeCompactGraph ()
{
    int tempi [2];
    float dist;
    std::vector <float> lons, lats;
    long long ni;
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;

    // Vertices first
    for (boost::unordered_set <long long>::iterator it = terminalNodes.begin ();
            it != terminalNodes.end (); ++it)
    {
        assert ((umapitr = nodes.find (*it)) != nodes.end ());
        oneVert.id = (*umapitr).first;
        oneVert.lon = ((*umapitr).second).first;
        oneVert.lat = ((*umapitr).second).second;
        boost::add_vertex (oneVert, gCompact);
    }

    std::cout << "Compact Graph has " << num_vertices (gCompact) << 
        " vertices" << std::endl;

    // Then edges
    // from http://theboostcpplibraries.com/boost.graph-vertices-and-edges,
    // but the remainder doesn't work
    typedef boost::graph_traits <Graph_t>::vertex_iterator viter;
    Graph_t::out_edge_iterator eit, eend;
    for (std::pair <viter, viter> vp = vertices (gFull);
            vp.first != vp.second; ++vp.first)
    {
        std::tie (eit, eend) = boost::out_edges (*vp.first, gFull);
    }

    // from
    // https://valelab.ucsf.edu/svn/3rdpartypublic/boost/libs/graph/test/grid_graph_test.cpp
    typedef boost::graph_traits <Graph_t>::vertices_size_type vertices_size_type;
    typedef boost::graph_traits <Graph_t>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits <Graph_t>::edges_size_type edges_size_type;
    typedef boost::graph_traits <Graph_t>::edge_descriptor edge_descriptor;
    const int nc = 7;
    int nedges, counts [nc] = {0, 0, 0, 0, 0, 0};
    BOOST_FOREACH (vertex_descriptor current_vertex, vertices (gFull))
    {
        edges_size_type out_edges_count = 0; 
        edges_size_type in_edges_count = 0;
        std::cout << "[O]";
        BOOST_FOREACH (edge_descriptor out_edge, 
                boost::out_edges (current_vertex, gFull))
        {
            out_edges_count++;
            std::cout << get (boost::vertex_index, gFull, 
                    source (out_edge, gFull)) << "-";
        }
        std::cout << "---[I]";
        BOOST_FOREACH (edge_descriptor in_edge, 
                boost::in_edges (current_vertex, gFull))
        {
            std::cout << get (boost::vertex_index, gFull, 
                    source (in_edge, gFull)) << "-";
        }
        std::cout << std::endl;
        // These include double counts:
        out_edges_count = boost::out_degree (current_vertex, gFull);
        in_edges_count = boost::in_degree (current_vertex, gFull);
        nedges = (int) out_edges_count + (int) in_edges_count;
        if (nedges < nc)
            counts [nedges]++;
    }
    std::cout << "#vertices with (";
    for (int i=0; i<(nc - 1); i++) 
        std::cout << i << ", ";
    std::cout << nc - 1 << ") edges = (";
    for (int i=0; i<(nc - 1); i++)
        std::cout << counts [i] << ", ";
    std::cout << counts [nc - 1] << ")" << std::endl;

    // TODO: Why are there vertices with multiple versions of same out-edges?

} // end function Graph::makeCompactGraph

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         FUNCTION::CALCDIST                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


float Graph::calcDist (std::vector <float> x, std::vector <float> y)
{
    float d, dsum = 0.0, x0, y0, x1 = x[0], y1 = y[0], xd, yd;
    assert (x.size () == y.size ());

    for (int i=1; i<x.size (); i++)
    {
        x0 = x1;
        y0 = y1;
        x1 = x[i];
        y1 = y[i];
        
        xd = (x1 - x0) * PI / 180.0;
        yd = (y1 - y0) * PI / 180.0;

        d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y1 * PI / 180.0) *
            cos (y0 * PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
        d = 2.0 * atan2 (sqrt (d), sqrt (1.0 - d));

        dsum += d * 6371.0;
    }

    return dsum; // in kilometres!
} // end function Graph::calcDist

