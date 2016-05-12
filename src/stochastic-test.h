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
 * ./test -t 4 starts producing wrong results:
 * [1->2] = 0.2718
 * [2->(3,6)] = 0.065 + 0.196 = 0.2615
 * [(12,15)->16] = 0.1416 + 0.1542 = 0.2958
 * These values should be equal yet aren't.
 *
 * Reverting all edges to one way greatly improves the situation, but there are
 * still slight divergences. There is obviously something causing the entire
 * algorithm to fail.
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


class Test
{
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::directedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        float _theta;
    protected:
        int from_node, to_node;
        Graph_t gr;
    public:
        int err;
        float tempf;
        std::vector <std::pair <int, float>> dists; // returned from dijkstra
        std::vector <std::pair <int, float>> pvec; // returned from make_pmat
        boost::numeric::ublas::matrix <float> cost_mat;
        boost::numeric::ublas::matrix <float> wmat, pmat, nmat;
        std::vector <float> z1, zn;

    Test (float theta)
        : _theta (theta)
    {
        err = make_graph2 ();
        std::cout << "Graph has " << num_vertices (gr) << " vertices" << std::endl;
        from_node = 0;
        to_node = 15;
        tempf = make_cost_mat (from_node, to_node);
        err = calc_zn (to_node);
        err = calc_pmat (to_node);

        // Then expected energy
        err = calc_z1 (from_node);
        tempf = calc_energy (from_node, to_node);
        std::cout << "-----energy = " << tempf << "-----" << std::endl;
        err = calc_link_density (from_node, to_node);
        err = dump_routes (from_node, to_node);
    }
    ~Test ()
    {
        pvec.resize (0);
        cost_mat.resize (0, 0);
        z1.resize (0);
        zn.resize (0);
        pmat.resize (0, 0);
        nmat.resize (0, 0);
    }

    int make_graph2 ();
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
 **                         TEST::MAKE_GRAPH2                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::make_graph2 ()
{
    std::string type = "highway";
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;

    // First network
    /*
    oneVert.component = 0;
    oneVert.id = 1;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 2;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 3;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 4;
    boost::add_vertex (oneVert, gr);

    oneEdge.id = 12;
    oneEdge.dist = oneEdge.weight = 1.0;
    boost::add_edge (0, 1, oneEdge, gr);
    oneEdge.id = 23;
    boost::add_edge (1, 2, oneEdge, gr);
    oneEdge.id = 24;
    boost::add_edge (1, 3, oneEdge, gr);
    oneEdge.id = 32;
    boost::add_edge (2, 1, oneEdge, gr);
    */

    // Second network
    /*
    oneVert.component = 0;
    oneVert.id = 1;
    oneVert.lon = 0;
    oneVert.lat = 3;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 2;
    oneVert.lon = 1;
    oneVert.lat = 3;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 3;
    oneVert.lon = 3.1;
    oneVert.lat = 4;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 4;
    oneVert.lon = 3;
    oneVert.lat = 1;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 5;
    oneVert.lon = 4;
    oneVert.lat = 2;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 6;
    oneVert.lon = 4;
    oneVert.lat = 0;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 7;
    oneVert.lon = 5;
    oneVert.lat = 2;
    boost::add_vertex (oneVert, gr);
    oneVert.id = 8;
    oneVert.lon = 6;
    oneVert.lat = 3;
    boost::add_vertex (oneVert, gr);

    oneEdge.id = 12;
    oneEdge.dist = oneEdge.weight = 0.0;
    boost::add_edge (0, 1, oneEdge, gr);

    oneEdge.dist = oneEdge.weight = 1.0;
    oneEdge.id = 23;
    boost::add_edge (1, 2, oneEdge, gr);
    oneEdge.id = 24;
    boost::add_edge (1, 3, oneEdge, gr);
    oneEdge.id = 34;
    boost::add_edge (2, 3, oneEdge, gr);
    oneEdge.id = 45;
    boost::add_edge (3, 4, oneEdge, gr);
    oneEdge.id = 46;
    boost::add_edge (3, 5, oneEdge, gr);
    oneEdge.id = 67;
    boost::add_edge (5, 6, oneEdge, gr);
    oneEdge.id = 57;
    boost::add_edge (4, 6, oneEdge, gr);
    oneEdge.id = 73;
    boost::add_edge (6, 2, oneEdge, gr);
    oneEdge.id = 78;
    boost::add_edge (6, 7, oneEdge, gr);

    oneEdge.dist = oneEdge.weight = 2.0;
    oneEdge.id = 28;
    boost::add_edge (1, 7, oneEdge, gr);
    oneEdge.id = 38;
    boost::add_edge (2, 7, oneEdge, gr);
    */

    // A more complex network with undirected edges
    oneVert.component = 0;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
        {
            oneVert.id = i * 4 + j + 1;
            oneVert.lon = j;
            oneVert.lat = i;
            boost::add_vertex (oneVert, gr);
        }
    /*
     * 13---14--15--16
     *    /   / |   |
     * 9    10  11--12
     * |  /     |   |
     * 5----6---7   8
     *      |       |
     * 1----2---3---4
     */

    oneEdge.dist = oneEdge.weight = 1.0;
    oneEdge.id = 12;
    boost::add_edge (0, 1, oneEdge, gr);
    oneEdge.id = 23;
    boost::add_edge (1, 2, oneEdge, gr);
    //boost::add_edge (2, 1, oneEdge, gr);
    oneEdge.id = 26;
    boost::add_edge (1, 5, oneEdge, gr);
    //boost::add_edge (5, 1, oneEdge, gr);
    oneEdge.id = 34;
    boost::add_edge (2, 3, oneEdge, gr);
    //boost::add_edge (3, 2, oneEdge, gr);
    oneEdge.id = 48;
    boost::add_edge (3, 7, oneEdge, gr);
    //boost::add_edge (7, 3, oneEdge, gr);
    oneEdge.id = 59;
    boost::add_edge (4, 8, oneEdge, gr);
    //boost::add_edge (8, 4, oneEdge, gr);
    oneEdge.id = 510;
    boost::add_edge (4, 9, oneEdge, gr);
    //boost::add_edge (9, 4, oneEdge, gr);
    oneEdge.id = 65;
    //boost::add_edge (4, 5, oneEdge, gr);
    boost::add_edge (5, 4, oneEdge, gr);
    oneEdge.id = 67;
    boost::add_edge (5, 6, oneEdge, gr);
    //boost::add_edge (6, 5, oneEdge, gr);
    oneEdge.id = 711;
    boost::add_edge (6, 10, oneEdge, gr);
    //boost::add_edge (10, 6, oneEdge, gr);
    oneEdge.id = 812;
    boost::add_edge (7, 11, oneEdge, gr);
    //boost::add_edge (11, 7, oneEdge, gr);
    oneEdge.id = 914;
    boost::add_edge (8, 13, oneEdge, gr);
    //boost::add_edge (13, 8, oneEdge, gr);
    oneEdge.id = 1015;
    boost::add_edge (9, 14, oneEdge, gr);
    //boost::add_edge (14, 9, oneEdge, gr);
    oneEdge.id = 1112;
    boost::add_edge (10, 11, oneEdge, gr);
    //boost::add_edge (11, 10, oneEdge, gr);
    oneEdge.id = 1115;
    boost::add_edge (10, 14, oneEdge, gr);
    //boost::add_edge (14, 10, oneEdge, gr);
    oneEdge.id = 1216;
    boost::add_edge (11, 15, oneEdge, gr);
    oneEdge.id = 1314;
    boost::add_edge (12, 13, oneEdge, gr);
    //boost::add_edge (13, 12, oneEdge, gr);
    oneEdge.id = 1415;
    boost::add_edge (13, 14, oneEdge, gr);
    //boost::add_edge (14, 13, oneEdge, gr);
    oneEdge.id = 1516;
    boost::add_edge (14, 15, oneEdge, gr);

    return 0;
} // end Test::make_graph2



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        ROUTER::MAKE_COST_MAT                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

float Test::make_cost_mat (int fromNode, int toNode)
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
                //if (!gr [e].oneway)
                //    wmat (v2, v1) = wmat (v1, v2);
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
    assert (rsum_max < 1.0);
    //std::cout << "***rmax = " << rsum_max << "***" << std::endl;

    return (float) count / ((float) nv * (float) nv);
} // end Test::make_costmat 


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           ROUTER::CALC_Z1                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::calc_z1 (int fromNode)
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
} // end Test::calc_z1


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           ROUTER::CALC_ZN                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::calc_zn (int toNode)
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
} // end Test::calc_zn


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          TEST::CALC_PMAT                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::calc_pmat (int toNode)
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
} // end Test::calc_pmat


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          TEST::CALC_ENERGY                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

float Test::calc_energy (int fromNode, int toNode)
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
} // end Test::calc_energy


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                      TEST::CALC_LINK_DENSITY                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::calc_link_density (int fromNode, int toNode)
{
    float tempf;
    const int nv = num_vertices (gr);

    nmat.resize (nv, nv);

    for (int i=0; i<nv; i++)
        for (int j=0; j<nv; j++)
            nmat (i, j) = z1 [i] * zn [j] * wmat (i, j) / z1 [toNode];
            // nmat (i, j) = z1 [i] * zn [j] * wmat (i, j) / zn [fromNode];

    return 0;
} // end Test::calc_link_density


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         ROUTER::DUMP_ROUTES                        **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Test::dump_routes (int fromNode, int toNode)
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
        junk <- apply (dat, 1, function (i)
                       lines (c (i [1], i [3]), c (i [2], i [4]), lwd=3*i[5]))
        points (from [1], from [2], pch=1, cex=2, col="green")
        points (to [1], to [2], pch=1, cex=2, col="red")

        x <- c (dat [,1], dat [,3])
        y <- c (dat [,2], dat [,4])
        dat2 <- cbind (x, y)
        dat3 <- dat2 [which (!duplicated (dat2)),]
        points (dat2, pch=1)
    }
    */
    return 0;
}
