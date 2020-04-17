#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <vector>
#include <sstream>
#include <list>
#include <deque>
#include <sys/stat.h>
#include <typeinfo>
#include <sys/resource.h>
#include <string.h>
#include <cerrno>

#include <gd.h>
#include <gdfontl.h>

#ifndef DISABLE_PNGWRITER
	#include <pngwriter.h>
#endif

#ifndef DISABLE_NETCDF4
	#include <netcdf.h>
	#include <hdf5.h>
#endif



//#include <boost/algorithm/string.hpp>
//#include <boost/lambda/lambda.hpp>
//#include <boost/lexical_cast.hpp>


#ifndef DISABLE_BOOST_SERIALIZATION
	#include <boost/serialization/list.hpp>
	#include <boost/serialization/string.hpp>
	#include <boost/serialization/vector.hpp>
	#include <boost/archive/binary_oarchive.hpp>
	#include <boost/archive/binary_iarchive.hpp>
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
	using namespace boost;
#endif

//#include <boost/mpi.hpp>

//using namespace boost::lambda;

// da prav dela error(..) - da prav displaya line number
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define FU __PRETTY_FUNCTION__
#define TT "\t"
#define tabt "\t"
#define ERRORIF(x) if (x) error(AT,FU, #x)



#define BAD_DATA_FLOAT -9999


const rlim_t kStackSize = 1000 * 1024 * 1024;   // min stack size = 16 MB

long random_number_seed=-1;

using namespace std;

//string ramdisk="/mnt/ramdisk";

#include "CU_utils.cc"
#include "CU_Archive.cc"
#include "CU_2Dunit.cc"
#include "CU_CField_class.cc"
#include "CU_C3DSupport_code.cc"
#include "CU_C3DObject_classes.cc"
#include "CU_C3DNamelist_read.cc"
#include "CU_Polygon.cc"
#include "CU_read_write.cc"
#include "CU_C3DField_class.cc"
#include "CU_3Dunit.cc"
#include "CU_Object_classes.cc"
#include "CU_Object_classes_extended_for_tropical_cyclones.cc"
#include "CU_CDomain_class.cc"
#include "CU_CHistogram_class.cc"
#include "CU_tropical_cyclones_functions.cc"
#include "CU_CFast_Convolution_object.cc"
#include "CU_Visualization.cc"

