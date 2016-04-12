#include <assert.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <stdlib.h>

/***************************************************************************
 * This program extracts OSM graph weighting profiles from the configuration
 * file used in the Routino project (http://routino.org/) in order to use
 * them in the osm-routing project.
 *
 * Author: Andreas Petutschnig
 * Email: andreas@petutschnig.de
 *
 * License: GNU General Public License v3.0
 **************************************************************************/

int main (int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: extract <profiles.xml>" << std::endl;
        exit (1);
    }
    std::string profileType, profileName, osmHighway, percent_str, linetxt, xml_in = argv[1];
    std::string output_dir = "../data/weighting_profiles";
    unsigned ipos;
    float percent, percent_default = 0.01;
    std::ifstream stream_in;
    std::ofstream stream_out;
    stream_in.open (xml_in.c_str (), std::ifstream::in);
    assert (!stream_in.fail ());

    //create output directory if not present
    boost::filesystem::path dir (output_dir);
    if (!(boost::filesystem::exists (dir)))
    {
        std::cout << "Creating output directory at " << output_dir << std::endl;
        boost::filesystem::create_directory (dir);
    }

    while (getline (stream_in, linetxt, '\n'))
    {
        //next profile
        if (linetxt.find ("profile name") != std::string::npos)
        {
            //close previous output file
            stream_out.close ();

            //extract profile name and open output file
            ipos = linetxt.find ("\"");
            profileType = linetxt.substr (ipos + 1, linetxt.length ());
            ipos = profileType.find ("\"");
            profileType = profileType.substr (0, ipos);

            profileName = output_dir + "/profile_" + profileType + ".cfg";
            stream_out.open (profileName.c_str ());
        }
        if (linetxt.find ("preference highway") != std::string::npos)
        {
            //get OSM way type
            ipos = linetxt.find ("\"");
            osmHighway = linetxt.substr (ipos + 1, linetxt.length ());
            ipos = osmHighway.find ("\"");
            osmHighway = osmHighway.substr (0, ipos);

            //get percentage
            ipos = linetxt.find ("percent");
            percent_str = linetxt.substr (ipos , linetxt.length ());
            ipos = percent_str.find ("\"");
            percent_str = percent_str.substr (ipos + 1, percent_str.length ());
            ipos = percent_str.find ("\"");
            percent_str = percent_str.substr (0, ipos);
            percent = atof (percent_str.c_str ());
            percent /= 100;

            if (percent == 0)
                percent = percent_default;

            //write to file
            stream_out << osmHighway << ", " << percent << "\n";
        }
    }
}
