#include <iostream>
#include <cstdio>
using namespace std;

#include "src/AcerDET.h"
using namespace AcerDet;
using namespace AcerDet::conf;

int main( int argc, char **argv ) {
/*
	if( argc < 2 ) {
		cout << "Usage: " << argv[0] << " <output_file> " << endl;
		return -1;
	}
*/
	const std::string configFileName = "template.dat"; //std::string( argv[1] );
	Configuration config = Configuration::fromFile( configFileName );

/*
	Configuration config = Configuration::getDefault();
	Configuration::save( config, "template.dat" );
*/

	Configuration::save( config, "tres.dat" );
	return 0;
}
