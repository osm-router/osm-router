/***************************************************************************
 *  Project:    osm-router
 *  File:       stochastic-sp.c++
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
 *                  given an origin-destination pair. All functionality is
 *                  contained in the header file; this .c++ file just provides a
 *                  stand-along compilable wrapper.
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/

#include "stochastic-sp.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    bool compact;
    float xfrom, yfrom, xto, yto, theta;
    Bbox bbox;
    std::string xml_file, profile_file;

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("graph compact, g", boost::program_options::value <bool>
                (&compact)->default_value (true),
                "use compact graph")
            ("xml_file,f", boost::program_options::value <std::string> 
                (&xml_file)->default_value ("xmldat.xml"), 
                "xml_file name (.xml will be appended)")
            ("profile_file,p", boost::program_options::value <std::string> 
                (&profile_file)->default_value ("../profile.cfg"), 
                "profile file name")
            ("xfrom,a", boost::program_options::value <float> 
                (&xfrom)->default_value (-0.12), "xfrom")
            ("yfrom,b", boost::program_options::value <float> 
                (&yfrom)->default_value (51.515), "yfrom")
            ("xto,c", boost::program_options::value <float> 
                (&xto)->default_value (-0.115), "xto")
            ("yto,d", boost::program_options::value <float> 
                (&yto)->default_value (51.52), "yto")
            ("theta,t", boost::program_options::value <float> 
                (&theta)->default_value (1000.0), "theta")
            ;
        // theta is the spectral radius which is scaled in inverse proportion to
        // the scale of distances. The latter are in km, so theta has to be
        // large.

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        boost::program_options::options_description visible("Allowed options");
        visible.add(generic).add(config);

        boost::program_options::variables_map vm;
        store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).run(), vm);

        notify(vm);

        if (vm.count("help")) {
            std::cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "osm-router, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    

    bbox = get_bbox (xfrom, yfrom, xto, yto);

    Router router (xml_file, profile_file, 
            xfrom, yfrom, xto, yto,
            bbox.lonmin, bbox.latmin, bbox.lonmax, bbox.latmax, 
            compact, theta); 
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                               GETBBOX                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Bbox get_bbox (float xfrom, float yfrom, float xto, float yto, float expand)
{
    // expand is in both directions, so actually twice that value
    float xmin, ymin, xmax, ymax, xrange, yrange;
    Bbox bbox;

    if (xfrom < xto)
    {
        xmin = xfrom;
        xmax = xto;
    } else
    {
        xmin = xto;
        xmax = xfrom;
    }
    if (yfrom < yto)
    {
        ymin = yfrom;
        ymax = yto;
    } else
    {
        ymin = yto;
        ymax = yfrom;
    }

    xrange = xmax - xmin;
    yrange = ymax - ymin;

    bbox.lonmin = xmin - xrange * expand;
    bbox.lonmax = xmax + xrange * expand;
    bbox.latmin = ymin - yrange * expand;
    bbox.latmax = ymax + yrange * expand;

    return bbox;
}
