/***************************************************************************
 *  Project:    osm-router
 *  File:       xml-parser.c++
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

#include "xml-parser.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    std::string city;

    Xml xml (true); // returns all nodes and ways
    std::cout << "Parsed " << xml.nodes.size () << " nodes and " << 
        xml.ways.size () << " ways" << std::endl;

    return 0;
};


