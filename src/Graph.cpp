/***************************************************************************
 *  Project:    osm-router
 *  File:       Graph.cpp
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
 *  Description:    Derived from class 'Xml' of 'xml-parser.h'. Takes OSM ways
 *                  and nodes from that class and stores them as a boost::graph
 *                  structures.  all highways and associated nodes. All
 *                  functionality is contained in the header file; this .c++
 *                  file just provides a stand-alone compilable wrapped.
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/

#include "Graph.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) 
{
    bool compact;
    float lonmin, latmin, lonmax, latmax;
    std::string xml_file, profile_file;

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("graph compact,g", boost::program_options::value <bool>
                (&compact)->default_value (true),
                "use compact graph")
            ("xml_file,f", boost::program_options::value <std::string> 
                (&xml_file)->default_value ("xmldat.xml"), 
                "xml_file name (.xml will be appended)")
            ("profile_file,p", boost::program_options::value <std::string> 
                (&profile_file)->default_value ("../profile.cfg"), 
                "profile file name")
            ("lonmin,a", boost::program_options::value <float> 
                (&lonmin)->default_value (-0.12), 
                "lonmin (ignored if xml_file exists)")
            ("latmin,b", boost::program_options::value <float> 
                (&latmin)->default_value (51.515), "latmin (ditto)")
            ("lonmax,c", boost::program_options::value <float> 
                (&lonmax)->default_value (-0.115), "lonmax (ditto)")
            ("latmax,d", boost::program_options::value <float> 
                (&latmax)->default_value (51.52), "latmax (ditto)")
            ;

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

    if ((int) xml_file.find (".") < 0)
        xml_file = xml_file + ".xml";

    Graph graph (compact, xml_file, profile_file, 
            lonmin, latmin, lonmax, latmax); 
};
