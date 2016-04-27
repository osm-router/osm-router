/***************************************************************************
 *  Project:    osm-router
 *  File:       xml-test.c++
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
 *  Description:    test of boost xml parser
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/

#include "xml-test.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    std::string city;

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("city,c", boost::program_options::value <std::string> 
                (&city)->default_value ("london"), "city")
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

    std::transform (city.begin(), city.end(), city.begin(), ::tolower);
    if (city.substr (0, 2) == "lo")
        city = "london";
    else if (city.substr (0, 2) == "ny")
        city = "nyc";

    Xml xml (city);

    return 0;
};

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                               GETBBOX                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Xml::getBBox ()
{
    // TODO: Copy old getBBox from Graph.c++
    lonmin = -0.12;
    lonmax = -0.115;
    latmin = 51.515;
    latmax = 51.52;

    return 0;
}; // end function getBBox


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            READOVERPASS                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

std::string Xml::readOverpass ()
{
    const std::string key = "['highway']",
                url_base = "http://overpass-api.de/api/interpreter?data=";
    std::stringstream bbox, query, url;

    bbox << "";
    bbox << "(" << latmin << "," << lonmin << "," << 
        latmax << "," << lonmax << ")";

    query << "";
    query << "(node" << key << bbox.str() << ";way" << key << bbox.str() <<
                    ";rel" << key << bbox.str() << ";";
    url << "";
    url << url_base << query.str() << ");(._;>;);out;";

    CURLplusplus client;
    std::string x = client.Get (url.str().c_str ());

    return x;
}; // end function readOverpass


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              READNODES                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Ways Xml::readNodes ( std::string & is )
{
    // populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;
    std::stringstream istream (is, std::stringstream::in);
    read_xml (istream, pt);

    Ways ways;

    // std::cout << is << std::endl; // The overpass XML data
    traverse (pt);

    return ways;
}

void Xml::traverse (const boost::property_tree::ptree& pt)
{
    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        //std::cout << it->first << ": " << 
        //    it->second.get_value<std::string>() << std::endl;
        if (it->first == "node")
        {
            std::cout << "-----";
            traverseNode (it->second);
        }
        if (it->first == "way")
        {
            std::cout << "-----";
            traverseWay (it->second);
        } else
            traverse (it->second);
    }
}

void Xml::traverseWay (const boost::property_tree::ptree& pt)
{
    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        if (it->first == "k")
            std::cout << it->second.get_value<std::string>();
        else if (it->first == "v")
            std::cout << ": " << it->second.get_value<std::string>() << std::endl;
        else if (it->first == "id" | it->first == "ref")
            std::cout << it->first << ": " << 
                it->second.get_value<std::string>() << std::endl;
        traverseWay (it->second);
    }
}

void Xml::traverseNode (const boost::property_tree::ptree& pt)
{
    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        if (it->first == "id")
            std::cout << it->first << ": " << 
                it->second.get_value<std::string>();
        else if (it->first == "lat")
            std::cout << " (" << it->second.get_value<std::string>();
        else if (it->first == "lon")
            std::cout << ", " << it->second.get_value<std::string>() <<
                ")" << std::endl;
        else if (it->first == "k")
            std::cout << it->second.get_value<std::string>();
        else if (it->first == "v")
            std::cout << ": " << it->second.get_value<std::string>() << std::endl;
        traverseNode (it->second);
    }
}
