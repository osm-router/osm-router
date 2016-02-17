#include <assert.h>
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
    std::string profileType, profileName, osmHighway, percent_str, linetxt, xml_in = argv[1];
    unsigned ipos;
    float percent, percent_default = 0.01;
    std::ifstream stream_in;
    std::ofstream stream_out;
    stream_in.open (xml_in.c_str (), std::ifstream::in);
    assert (!stream_in.fail ());

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

            profileName = "profile_" + profileType + ".cfg";
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
