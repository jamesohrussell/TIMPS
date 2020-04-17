
// This is all done in order to serialize the output and reading from binary files
// So far it supports the write/read of:
// - C++ native numeric type = int, long, double
// - C++ vectors of types mentioned above
// - C++ native strings
// - C++ objects that must have defined the following methods: freememory, deep_copy an a template for archive_object
// - C++ vectors of above mentioned objects
//
// To read/write a value one has to open a file using either ofstream or ifstream
// Then one simply call the Archive(filestream, value) function which does everyting
// Depending if the filestream is defined as ofstream or ifstrem the Archive will write or read
// Everthing has to be written in the same order as it will be read
// For writing of objects a special archive_object function has to be defined in which all the desired object properties are Archived using the Archive function

// This is something similar to what the BOOST serialization is doing - but more simple

// --------------------------------------------
//  An axample of archive_object method
// --------------------------------------------
//
//		 class CTropical_cyclone_NEW
//			{
//			public:
//			vector <CTropical_cyclone_frame_NEW> frame;
//			string name;
//
//			template <typename T>
//			void archive_object(T & myFile)
//				{
//				Archive(myFile, name);
//				Archive(myFile, frame);
//				}
//			};




// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//            Functions for archiving to and from binary file
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------
// Generic function which are later used by Archive function to overload for different types
// -------------------------------------------------------------------------------------------
template<typename T>
    void Archive_basic_datatype(ofstream & myFile, const T & val)
    {
	myFile.write ((char *) &val, sizeof(T));
    }
template<typename T>
    void Archive_basic_datatype(ifstream & myFile, T & val)
    {
	myFile.read ((char *) &val, sizeof(T));
    }
template<typename T>
    void Archive_string(ofstream & myFile, const T & val)
    {
	Write_append_string_to_opened_binary_file(val, myFile);
    }
template<typename T>
    void Archive_string(ifstream & myFile, T & val)
    {
	val=Read_string_from_opened_binary_file(myFile);
    }
template<typename T>
    void Archive_vector_of_nonobjects(ofstream & myFile, std::vector<T>& vec)
    {
	long count=vec.size();
	myFile.write ((char *) &count, sizeof(long));
	if (count > 0)
		myFile.write ((char *) &vec.front(), count*sizeof(T));
    }
template<typename T>
    void Archive_vector_of_nonobjects(ifstream & myFile, std::vector<T>& vec)
    {
	long count;
	myFile.read ((char *) &count, sizeof(long));
	if (count > 0)
		{
		vec.assign(count,0);
		myFile.read ((char *) &vec.front(), count*sizeof(T));
		}
    }

template<typename T>
    void Archive_object_using_archive_object_method(ofstream & myFile, T & val)
    {
	val.archive_object(myFile);
    }
template<typename T>
    void Archive_object_using_archive_object_method(ifstream & myFile, T & val)
    {
	val.freememory();
	val.archive_object(myFile);
    }

template<typename T>
    void Archive_vector_of_objects_using_archive_object_method(ofstream & myFile, std::vector<T>& vec)
    {
	long count=vec.size();
	myFile.write ((char *) &count, sizeof(long));
	for (long it=0; it < count; it++)
		Archive_object_using_archive_object_method(myFile,vec[it]);
    }
template<typename T>
    void Archive_vector_of_objects_using_archive_object_method(ifstream & myFile, std::vector<T>& vec)
    {
	vec.clear();
	long count;
	myFile.read ((char *) &count, sizeof(long));
	//cout << "--------------" << count << endl;
	T temp;
	for (long it=0; it < count; it++)
		{
		Archive_object_using_archive_object_method(myFile,temp);
		vec.push_back(temp.deep_copy());
		}
    }

// -------------------------------------------------------------------------------------------
// Overloaded variants of the general Archive function for different types
// -------------------------------------------------------------------------------------------

// single integer
template<typename fstream_T>
void Archive(fstream_T & myFile, int & val)
    {
	Archive_basic_datatype(myFile,val);
    }

// single long
template<typename fstream_T>
void Archive(fstream_T & myFile, long & val)
    {
	Archive_basic_datatype(myFile,val);
    }

// single double
template<typename fstream_T>
void Archive(fstream_T & myFile, double & val)
    {
	Archive_basic_datatype(myFile,val);
    }

// single string
template<typename fstream_T>
void Archive(fstream_T & myFile, string & val)
    {
	Archive_string(myFile,val);
    }

// vector integer
template<typename fstream_T>
void Archive(fstream_T & myFile, std::vector<int>& vec)
    {
	Archive_vector_of_nonobjects(myFile,vec);
    }

// vector long
template<typename fstream_T>
void Archive(fstream_T & myFile, std::vector<long>& vec)
    {
	Archive_vector_of_nonobjects(myFile,vec);
    }

// vector double
template<typename fstream_T>
void Archive(fstream_T & myFile, std::vector<double>& vec)
    {
	Archive_vector_of_nonobjects(myFile,vec);
    }

// vector of objects using the archive_object method
template<typename fstream_T, typename T>
void Archive(fstream_T & myFile, std::vector<T> & vec)
    {
	Archive_vector_of_objects_using_archive_object_method(myFile,vec);
    }

// single object using the archive_object method
template<typename fstream_T, typename T>
void Archive(fstream_T & myFile, T & val)
    {
    Archive_object_using_archive_object_method(myFile, val);
    }
