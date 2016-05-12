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
#include <boost/numeric/ublas/matrix_proxy.hpp>

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
    public:
        int err;
        int from_node, to_node;
        float tempf;
        std::vector <std::pair <int, float>> dists; // returned from dijkstra
        std::vector <std::pair <int, float>> pvec; // returned from make_pmat
        boost::numeric::ublas::matrix <float> cost_mat;
        boost::numeric::ublas::matrix <float> wmat, pmat, nmat;
        std::vector <float> z1, zn;

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
        err = calc_zn (to_node);
        err = calc_pmat (to_node);
        err = calc_z1 (from_node);
        tempf = calc_energy (from_node, to_node);
        std::cout << "-----energy = " << tempf << "-----" << std::endl;
        err = calc_link_density (from_node, to_node);
        err = dump_routes (from_node, to_node);
    }
    ~Router ()
    {
        pvec.resize (0);
        cost_mat.resize (0, 0);
        z1.resize (0);
        zn.resize (0);
        pmat.resize (0, 0);
        nmat.resize (0, 0);
    }

    long long nearestNode (float lon0, float lat0);
    int dijkstra (int fromNode);
    int make_pvec (int fromNode, int toNode);
    float make_cost_mat (int fromNode, int toNode);
    int calc_z1 (int fromNode);
    int calc_zn (int toNode);
    int calc_pmat (int toNode);
    float calc_energy (int fromNode, int toNode);
    int calc_link_density (int fromNode, int toNode);
    int dump_routes (int fromNode, int toNode);
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
} // end function nearestNode



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
} // end function Router::Dijkstra



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::MAKE_PVEC                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::make_pvec (int fromNode, int toNode)
{
    /*
     * TODO: Delete this!
     *
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
} // end function Router::make_pvec


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
     * column) is the start node and the last is the destination node. */

    const int nv = num_vertices (gr);
    int count = 1;

    //cost_mat.resize (nv, nv);
    cost_mat = boost::numeric::ublas::scalar_matrix <float> (nv, nv, FLOAT_MAX);
    wmat = boost::numeric::ublas::scalar_matrix <float> (nv, nv, 0.0);

    typedef boost::graph_traits <Graph_t>::vertex_descriptor vertex_t;
    typedef boost::graph_traits <Graph_t>::edge_descriptor edge_t;
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gr);
    boost::graph_traits <Graph_t>::out_edge_iterator ei, ei_end;

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
                cost_mat (v1, v2) = gr [e].weight;
                wmat (v1, v2) = exp (-_theta * cost_mat (v1, v2));
                if (!gr [e].oneway)
                    wmat (v2, v1) = wmat (v1, v2);
                count++;
            }
        }
    }

    // Then set cost_mat for start and end nodes
    typedef boost::numeric::ublas::matrix <float> Matf;
    boost::numeric::ublas::matrix_row <Matf> mr (cost_mat, toNode);
    for (unsigned i=0; i<mr.size (); ++i) // mr is a pointer
        mr (i) = FLOAT_MAX;
    /*
    boost::numeric::ublas::matrix_column <Matf> mc (cost_mat, toNode);
    for (unsigned i=0; i<mc.size (); ++i)
        if (mc (i) != FLOAT_MAX)
            mc (i) = 0.0;
    */
    // could use matrix_range, but direct index seems easier
    cost_mat (toNode, toNode) = FLOAT_MAX;

    // wmat needs only one row changed:
    boost::numeric::ublas::matrix_row <Matf> mrw (wmat, toNode);
    for (unsigned i=0; i<mrw.size (); ++i) 
        mrw (i) = 0.0;

    /* Method suggested by Saerens et al to check spectral radius without
     * needing to calculate eigenvalues. */
    float rsum, rsum_max = 0.0;
    for (Matf::iterator2 itr2 = wmat.begin2 (); itr2 != wmat.end2 (); ++itr2)
    {
        rsum = 0.0;
        for (Matf::iterator1 itr1 = itr2.begin (); itr1 != itr2.end (); ++itr1)
            rsum += (*itr1);
        if (rsum > rsum_max)
            rsum_max = rsum;
    }
    // Saeren et al: "the spectral radius ... is always smaller or equal than
    // its maximum absolute row sum norm" - this row sum can therefore be larger
    // than one.  TODO: Devise a better heuristic for spectral radius
    //assert (rsum_max < 1.0);

    return (float) count / ((float) nv * (float) nv);
} // end Router::make_cost_mat


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           ROUTER::CALC_Z1                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_z1 (int fromNode)
{
    const int max_loops = 1e6;
    const float tol = 1.0e-6;
    int nloops = 0;
    float diff = FLOAT_MAX;

    const int nv = num_vertices (gr);
    std::vector <float> z1_old;

    for (int i=0; i<nv; i++)
        z1_old.push_back (1.0 / (float) nv);

    z1.resize (nv);

    while (diff > tol)
    {
        diff = 0.0;
        for (int i=0; i<nv; i++)
        {
            z1 [i] = 0.0;
            for (int j=0; j<nv; j++)
                if (j == fromNode)
                    z1 [i] += wmat (j, i);
                else
                    z1 [i] += wmat (j, i) * z1_old [j];
            
            diff += std::abs (z1 [i] - z1_old [i]);
        }
        for (int i=0; i<nv; i++)
            z1_old [i] = z1 [i];
        nloops++;
        if (nloops > max_loops)
            break;
    }
    z1 [fromNode] = 1.0;

    return nloops;
} // end Router::calc_z1


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           ROUTER::CALC_ZN                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_zn (int toNode)
{
    const int max_loops = 1e6;
    const float tol = 1.0e-6;
    int nloops = 0;
    float diff = FLOAT_MAX;

    const int nv = num_vertices (gr);
    std::vector <float> zn_old;

    for (int i=0; i<nv; i++)
        zn_old.push_back (1.0 / (float) nv);

    zn.resize (nv);

    while (diff > tol)
    {
        diff = 0.0;
        for (int i=0; i<nv; i++)
        {
            zn [i] = 0.0;
            for (int j=0; j<nv; j++)
                if (j == toNode)
                    zn [i] += wmat (i, j);
                else
                    zn [i] += wmat (i, j) * zn_old [j];
            
            diff += std::abs (zn [i] - zn_old [i]);
        }
        for (int i=0; i<nv; i++)
            zn_old [i] = zn [i];
        nloops++;
        if (nloops > max_loops)
            break;
    }
    zn [toNode] = 1.0;

    return nloops;
} // end Router::calc_zn


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          ROUTER::CALC_PMAT                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_pmat (int toNode)
{
    const int nv = num_vertices (gr);

    pmat.resize (nv, nv);

    for (int i=0; i<nv; i++)
        for (int j=0; j<nv; j++)
        {
            if (zn [i] > 0.0 && zn [j] > 0.0 && wmat (i, j) > 0.0)
                pmat (i, j) = zn [j] * wmat (i, j) / zn [i];
            else
                pmat (i, j) = 0.0;
        }
    // pmat is then the matrix Q of Saeren et al

    for (int i=0; i<nv; i++)
    {
        pmat (toNode, i) = 0.0;
        if (wmat (i, toNode) > 0.0)
            pmat (i, toNode) = wmat (i, toNode);
        else
            pmat (i, toNode) = 0.0;
    }
    pmat (toNode, toNode) = 1.0;

    return 0;
} // end Router::calc_pmat


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         ROUTER::CALC_ENERGY                        **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

float Router::calc_energy (int fromNode, int toNode)
{
    float tempf, energy = 0.0;
    const int nv = num_vertices (gr);

    boost::numeric::ublas::matrix <float> zmat, wmat2;
    zmat.resize (nv, nv);
    wmat2.resize (nv, nv);

    for (int i=0; i<nv; i++)
        for (int j=0; j<nv; j++)
        {
            zmat (i, j) = z1 [i] * zn [j];
            wmat2 (i, j) = -wmat (i, j) * cost_mat (i, j);
        }
    // This is Saerens et al Eq. (36) (with extra detail in Eq. 69).
    // The *wrong* way is to exchange (from, to)Node.
    for (int i=0; i<nv; i++)
    {
        tempf = 0.0;
        for (int j=0; j<nv; j++)
            tempf += zmat (fromNode, j) * wmat (j, i);
        energy += zmat (i, toNode) * tempf;
    }
    energy = energy / zn [fromNode];

    zmat.resize (0, 0);
    wmat2.resize (0, 0);

    return energy;
} // end Router::calc_energy


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                      ROUTER::CALC_LINK_DENSITY                     **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::calc_link_density (int fromNode, int toNode)
{
    float tempf;
    const int nv = num_vertices (gr);

    nmat.resize (nv, nv);

    for (int i=0; i<nv; i++)
        for (int j=0; j<nv; j++)
            nmat (i, j) = z1 [i] * zn [j] * wmat (i, j) / z1 [toNode];
            // nmat (i, j) = z1 [i] * zn [j] * wmat (i, j) / zn [fromNode];

    return 0;
} // end Router::calc_link_density


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         ROUTER::DUMP_ROUTES                        **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::dump_routes (int fromNode, int toNode)
{
    // Dumps the link densities
    const int nv = num_vertices (gr);
    std::ofstream  out_file;

    boost::property_map< Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gr);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lat = boost::get(&bundled_vertex_type::lat, gr);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lon = boost::get(&bundled_vertex_type::lon, gr);
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gr);

    std::cout << "from = (" << vertex_lon [fromNode] << ", " <<
        vertex_lat [fromNode] << "); to = (" << vertex_lon [toNode] << 
        ", " << vertex_lat [toNode] << ")" << std::endl;

    out_file.open ("junk.txt", std::ofstream::out);
    
    // TODO: Delete these and write the whole thing better
    assert (pmat.size1() == nmat.size1());
    assert (pmat.size2() == nmat.size2());

    for (int i=0; i<pmat.size1 (); ++i)
        for (int j=0; j<pmat.size2 (); ++j)
            if (pmat (i, j) > 0.0 & nmat (i, j) > nmat (j, i))
                out_file << vertex_id [i] << ", " << vertex_id [j] << ", " <<
                    vertex_lon [i] << "," << vertex_lat [i] << "," <<
                    vertex_lon [j] << "," << vertex_lat [j] << "," <<
                    nmat (i, j) << ", " << nmat (j, i) << std::endl;

    out_file.close ();

    /* R script
    plotgraph <- function ()
    {
        from <- c (-0.120048, 51.5151)
        to <- c (-0.116075, 51.5199)
        fname <- "./build/junk.txt"
        dat <- read.csv (fname, sep=",", header=FALSE)
        xlims <- range (c (dat [,1], dat [,3]))
        ylims <- range (c (dat [,2], dat [,4]))
        plot.new ()
        par (mar=rep (0, 4))
        plot (NULL, NULL, xlim=xlims, ylim=ylims,
              xaxt="n", yaxt="n", xlab="", ylab="")
        x <- log10 (dat [,5]) # passage density
        x [!is.finite (x)] <- NA
        dat [,5] <- 5 * (x - min (x, na.rm=TRUE)) / diff (range (x, na.rm=TRUE))

        junk <- apply (dat, 1, function (i)
                       lines (c (i [1], i [3]), c (i [2], i [4]), lwd=3*i[5]))
        points (from [1], from [2], pch=1, cex=2, col="green")
        points (to [1], to [2], pch=1, cex=2, col="red")
    }
    */
    return 0;
}
