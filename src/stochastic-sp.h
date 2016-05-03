/***************************************************************************
 *  Project:    osm-router
 *  File:       stochastic-sp.h
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
 *  Description:    Derived from class 'Graph' of 'Graph.h'. Takes the
 *                  boost::graph from that class and dijkstra-routes all paths
 *                  from a specified origin. Additional function 'make_pmat'
 *                  calculates distances along all possible paths to all nodes
 *                  given an origin-destination pair. 
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

#include <cstdlib> // for std::abs
#include <complex> 

#include <boost/unordered_map.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <Eigen/Eigenvalues>

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
        float _theta;
    protected:
        float _xfrom, _yfrom, _xto, _yto;
        std::vector <int> node_order;
    public:
        int err;
        int from_node, to_node;
        float tempf;
        std::vector <std::pair <int, float>> dists; // returned from dijkstra
        std::vector <std::pair <int, float>> pvec; // returned from make_pmat
        boost::numeric::ublas::matrix <float> cost_mat;
        //Eigen::MatrixXf wmat; // Eigen matrix so that eigenvals can be calculated
        boost::numeric::ublas::matrix <float> wmat, pmat;
        std::vector <float> zn;

    Router (std::string xml_file, std::string profile_file,
            float xfrom, float yfrom, float xto, float yto,
            float lonmin, float latmin, float lonmax, float latmax,
            bool compact, float theta)
        : _xfrom (xfrom), _yfrom (yfrom), _xto (xto), _yto (yto),
            Graph (compact, xml_file, profile_file, 
                lonmin, latmin, lonmax, latmax), _theta (theta)
    {
        from_node = nodeIndx [nearestNode (_xfrom, _yfrom)];
        to_node = nodeIndx [nearestNode (_xto, _yto)];
        tempf = make_cost_mat (from_node, to_node);
        err = calc_zn ();
        err = calc_pmat ();
    }
    ~Router ()
    {
        node_order.resize (0);
        pvec.resize (0);
        cost_mat.resize (0, 0);
        zn.resize (0);
        pmat.resize (0, 0);
    }

    long long nearestNode (float lon0, float lat0);
    int dijkstra (int fromNode);
    int get_dists ();
    int make_pvec (int fromNode, int toNode);
    float make_cost_mat (int fromNode, int toNode);
    int calc_zn ();
    int calc_pmat ();
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



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::MAKE_PMAT                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::make_pvec (int fromNode, int toNode)
{
    /*
     * One SP produces a vector of distances from an origin to all intermediate
     * nodes.  Calculating from the destination produces the complementary
     * vector of distances from each of intermediate nodes to the destination.
     * Thus the total path length through each point can be calculated by adding
     * the two.
     *
     * The probability of travelling from the origin to each intermediate point
     * can then be directly obtained from this vector of total path lengths
     * which pass through those points.
     */
    int err = dijkstra (fromNode);

    pvec.resize (num_vertices (gr));
    for (std::vector <std::pair <int, float>>::iterator it=pvec.begin (); 
            it != pvec.end (); ++it) 
        (*it).second = 0.0;

    assert (dists.size () > 0); // Asserts that dijkstra has been run.
    pvec.resize (0);

    auto vs = dists;
    for (auto it = dists.begin (); it != dists.end (); ++it)
        pvec.push_back (std::make_pair ((*it).first, (*it).second));

    // Then repeat for toNode and add result to pvec
    err = dijkstra (toNode);
    // TODO: boost::zip this with an iterator
    assert (dists.size () == pvec.size ());
    for (int i=0; i<dists.size (); i++)
    {
        assert (dists [i].first == pvec [i].first);
        pvec [i].second += dists [i].second;
    }

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        ROUTER::MAKE_COST_MAT                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

float Router::make_cost_mat (int fromNode, int toNode)
{
    /* Fills both cost matrix (C) and W (Eq. 31) of Saeren et al (2009), which
     * contains finite entries only for direct connections between neighbouring
     * vertices. All other values are set here to FLOAT_MAX. Returns proportion
     * of matrix with finite values 
     *
     * fromNode and toNode are only used to renumber the matrix so the 1st (row,
     * column) is the start node and the last is the destination node. This
     * routine also fills the node_order vector.
     * */

    const int nv = num_vertices (gr);
    int n1, n2, count = 1;

    // fill node_order
    node_order.resize (nv);
    node_order [fromNode] = 0;
    node_order [toNode] = nv - 1;
    for (int i=0; i<nv; i++)
        if (i != fromNode && i != toNode)
            node_order [i] = count++;
    assert (count == (nv - 1));

    //cost_mat.resize (nv, nv);
    cost_mat = boost::numeric::ublas::scalar_matrix <float> (nv, nv, FLOAT_MAX);
    wmat = boost::numeric::ublas::scalar_matrix <float> (nv, nv, 0.0);

    /*
    // iterator over matrix:
    typedef boost::numeric::ublas::matrix <float> matrix;
    for (matrix::iterator1 itr1 = cost_mat.begin1 (); 
            itr1 != cost_mat.end1 (); ++itr1)
    {
        for (matrix::iterator2 itr2 = itr1.begin (); itr2 != itr1.end (); ++itr2)
        {
            std::cout << (*itr2) << ", ";
        }
        std::cout << std::endl;
    }
    */

    typedef boost::graph_traits <Graph_t>::vertex_descriptor vertex_t;
    typedef boost::graph_traits <Graph_t>::edge_descriptor edge_t;
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gr);
    boost::graph_traits <Graph_t>::out_edge_iterator ei, ei_end;

    /* only for Eigen::MatrixXf type
    wmat.resize (nv, nv);
    for (int i=0; i<nv; i++)
        for (int j=0; j<nv; j++)
            wmat (i, j) = 0.0;
    */

    count = 0;
    auto vs = boost::vertices (gr);
    for (auto it = vs.first; it != vs.second; ++it)
    {
        if (vertex_component [*it] == 0)
        {
            boost::tie (ei, ei_end) = out_edges (*it, gr);
            for ( ; ei != ei_end; ++ei)
            {
                edge_t e  = *boost::out_edges (*it, gr).first;
                vertex_t v1 = boost::source (*ei, gr);
                vertex_t v2 = boost::target (*ei, gr);
                n1 = node_order [v1];
                n2 = node_order [v2];
                cost_mat (n1, n2) = gr [e].weight;
                wmat (n1, n2) = exp (-_theta * cost_mat (n1, n2));
                if (!gr [e].oneway)
                    wmat (n2, n1) = wmat (n1, n2);
                count++;
            }
        }
    }

    // Check spectral radius of w
    /*
    Eigen::EigenSolver <Eigen::MatrixXf> es2 (wmat);
    float rho = 0.0, eabs;
    auto evs = es2.eigenvalues ();
    for (int i=0; i<evs.size (); i++)
    {
        eabs = std::abs (evs [i]);
        if (eabs > rho)
            rho = eabs;
    }
    if (rho >= 1.0)
        std::cout << "Spectral radius = " << rho << std::endl;
    assert (rho < 1.0);
    */

    // Equivalent, according to Saerens et al (although numerical comparisons
    // suggest otherwise):
    float rsum, rsum_max = 0.0;
    for (int i=0; i<nv; i++)
    {
        rsum = 0.0;
        for (int j=0; j<nv; j++)
            rsum += wmat (i, j);
        if (rsum > rsum_max)
            rsum_max = rsum;
    }
    assert (rsum_max < 1.0);

    return (float) count / ((float) nv * (float) nv);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           ROUTER::CALC_ZN                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_zn ()
{
    const int max_loops = 1e6;
    const float tol = 1.0e-6;
    int nloops = 0;
    float diff = FLOAT_MAX;

    const int nv = num_vertices (gr) - 1;
    std::vector <float> en, zn_old;

    for (int i=0; i<nv; i++)
    {
        en.push_back (0.0);
        zn_old.push_back (1.0 / (float) nv);
    }

    zn.resize (nv);

    while (diff > tol)
    {
        diff = 0.0;
        for (int i=0; i<nv; i++)
        {
            zn [i] = 0.0;
            // next loop works because wmat (nv-1) is destination point
            for (int j=0; j<nv; j++)
                zn [i] += wmat (i, j) * zn_old [j] + en [i];
            
            diff += std::abs (zn [i] - zn_old [i]);
        }
        for (int i=0; i<nv; i++)
            zn_old [i] = zn [i];
        nloops++;
        if (nloops > max_loops)
            break;
    }

    return nloops;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::CALC_PMAT                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_pmat ()
{
    const int nv = num_vertices (gr);

    pmat.resize (nv, nv);

    // TODO: Delete pmax
    float pmax = 0.0;
    for (int i=0; i<(nv - 1); i++)
        for (int j=0; j<(nv - 1); j++)
        {
            if (zn [i] > 0.0 && zn [j] > 0.0 && wmat (i, j) > 0.0)
                pmat (i, j) = zn [j] * wmat (i, j) / zn [i];
            else
                pmat (i, j) = 0.0;
            if (pmat (i, j) > pmax)
                pmax = pmat (i, j);
        }
    assert (pmax < 1.0);
    // pmat is then the matrix Q of Saeren et al

    for (int i=0; i<nv; i++)
    {
        pmat (nv - 1, i) = 0.0;
        if (wmat (i, nv - 1) > 0.0)
            pmat (i, nv - 1) = wmat (i, 0);
        else
            pmat (i, nv - 1) = 0.0;
    }
    pmat (nv - 1, nv - 1) = 1.0;

    return 0;
}
