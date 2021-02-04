
// error function
void error(const char *location,const char *function,string txt)
	{
	cout << "ERROR: " << location << " " << function << ": " <<  txt << endl;
	exit(1);
	}

double min(double x, double y)
	{
	if (x<y) return(x);
	else return(y);
	}

double max(double x, double y)
	{
	if (x>y) return(x);
	else return(y);
	}

double sqr(double x)
	{
	return (x*x);
	}

void switch_values(double &x, double &y)
	{
	double temp=y;
	y=x;
	x=temp;
	}

string output_logic_as_string(bool a)
	{
	if (a==true) return("true");
	else return("false");
	}

bool compare_double_with_epsilon_tolerance(double x, double y, double epsilon)
	{
	if (fabs(x-y) < epsilon)
		return(true);
	return(false);
	}

bool compare_double_vector_with_epsilon_tolerance(vector <double> x, vector <double> y, double epsilon)
	{
	if (x.size() != y.size())
		error(AT,FU, "x.size() != y.size()");

	for (unsigned long il=0; il < x.size(); il++ )
		if (!compare_double_with_epsilon_tolerance(x[il],y[il], epsilon))
			return(false);

	return(true);
	}


// number to string conversion
string n2s(double x)
	{
	ostringstream s1;
	s1.str("");
	s1 << x;
	return(s1.str());
	}

template <class ForwardIterator, class T>
  ForwardIterator lower_bound_which_never_returns_an_out_of_bounds_result (ForwardIterator first, ForwardIterator last, const T& val)
{
  ForwardIterator it;
  typename iterator_traits<ForwardIterator>::difference_type count, step;
  count = distance(first,last);
  while (count>0)
  {
    it = first; step=count/2; advance (it,step);
    if (*it<val) {                 // or: if (comp(*it,val)), for version (2)
      first=++it;
      count-=step+1;
    }
    else count=step;
  }

  if (first == last)
  	return last-1;
  else
	return first;
}

template<typename T>
void multiply_with_scalar_value_vector(std::vector<T>& vec, double value)
	{
	for (long il=0; il < (long)vec.size(); il++)
		vec[il]*=value;
	}

template<typename T>
void add_scalar_value_vector(std::vector<T>& vec, double value)
	{
	for (long il=0; il < (long)vec.size(); il++)
		vec[il]+=value;
	}

template<typename T>
T sum_vector(std::vector<T>& vec)
	{
	T sum=0;
	for (long il=0; il < (long)vec.size(); il++)
		sum+=vec[il];
	return(sum);
	}

template<typename T>
T percentile_value_which_changes_vector_order(std::vector<T>& vec, double percentile)
	{
	long nth=floor((double)vec.size()*percentile);
	nth_element (vec.begin(), vec.begin()+nth, vec.end());
	return(vec[nth]);
	}


template<typename T>
T avg_vector(std::vector<T>& vec)
	{
	return(sum_vector(vec)/(double)vec.size());
	}

template<typename T>
T mediana_vector(const std::vector<T>& vec )
	{
	std::vector<T> vectemp=vec;
    std::sort(vectemp.begin(), vectemp.end());
	return(vectemp.at(floor((double)vectemp.size()/2)));
	}


template<typename T>
T max_vector(std::vector<T>& vec)
	{
	if(vec.size()==0)
		error(AT,FU, "vec.size()==0");
	T max=vec[0];
	for (long il=0; il < (long)vec.size(); il++)
		if (vec[il] > max)
			max = vec[il];
	return(max);
	}

template<typename T>
T max_vector_with_index(std::vector<T>& vec, long &index)
	{
	if(vec.size()==0)
		error(AT,FU, "vec.size()==0");
	T max=vec[0];
	for (long il=0; il < (long)vec.size(); il++)
		if (vec[il] > max)
			{
			max = vec[il];
			index=il;
			}
	return(max);
	}


template<typename T>
T min_vector(std::vector<T>& vec)
	{
	if(vec.size()==0)
		error(AT,FU, "vec.size()==0");
	T min=vec[0];
	for (long il=0; il < (long)vec.size(); il++)
		if (vec[il] < min)
			min = vec[il];
	return(min);
	}

template<typename T>
bool is_monotonic_increasing(std::vector<T>& vec)
	{
	if(vec.size()>1)
		{
		for (long il=1; il < (long)vec.size(); il++)
			if (vec[il-1] > vec[il])
				return (false);
		}
	else
		return (true);
	}

template<typename T>
bool is_monotonic_decreasing(std::vector<T>& vec)
	{
	if(vec.size()>1)
		{
		for (long il=1; il < (long)vec.size(); il++)
			if (vec[il-1] < vec[il])
				return (false);
		}
	else
		return (true);
	}


template<typename T>
string output_vector_as_string(std::vector<T>& vec, string separator)
	{
	ostringstream s1;
	for (long il=0; il < (long)vec.size(); il++)
		{
		s1 << vec[il];
		if (il< (long)vec.size() - 1)
			s1 << separator;
		}
	return(s1.str());
	}

template<typename T>
string output_vector_as_string_with_index_and_newline(std::vector<T>& vec)
	{
	ostringstream s1;
	for (long il=0; il < (long)vec.size(); il++)
		{
		s1 << il << "\t" << vec[il] << endl;
		}
	return(s1.str());
	}



// from - http://learningcppisfun.blogspot.com/2008/04/remove-duplicates-from-vector.html
template<typename T>
    void removeDuplicates(std::vector<T>& vec)
    {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }


bool convert_string_to_double(string s, double &val)
	{
	// assume: "char * mystr" is a null-terminated string
	char * e;
	errno = 0;
	val = std::strtod(s.c_str(), &e);

	if (*e != '\0' ||  // error, we didn't consume the entire string
		errno != 0 )   // error, overflow or underflow
			{
			cout << s << endl;
			return (false);
			}
	return(true);
	}

bool convert_string_to_long(string s, long &val)
	{
	// assume: "char * mystr" is a null-terminated string
	char * e;
	errno = 0;
	val = std::strtol(s.c_str(), &e, 10);

	if (*e != '\0' ||  // error, we didn't consume the entire string
		errno != 0 )   // error, overflow or underflow
			{
			cout << s << endl;
			return (false);
			}
	return(true);
	}

vector <double> convert_string_vector_to_double_vector(vector <string> s)
	{
	vector <double> vec;
	double val;
	for (unsigned long il=0; il < s.size(); il++)
		{
		if (convert_string_to_double(s[il], val))
			vec.push_back(val);
		else
			error(AT,FU, "!convert_string_to_double(s[il], val)");
		}
	return(vec);
	}

vector <long> convert_string_vector_to_long_vector(vector <string> s)
	{
	vector <long> vec;
	long val;
	for (unsigned long il=0; il < s.size(); il++)
		{
		if (convert_string_to_long(s[il], val))
			vec.push_back(val);
		else
			error(AT,FU, "!convert_string_to_long(s[il], val)");
		}
	return(vec);
	}

void readCSV_into_string_vector (string fname, char delimiter, vector< vector<string> > &output)
	{
	fstream file(fname.c_str(), ios::in);
		if(!file.is_open())

	error(AT,FU, "!file.is_open()");

	output.clear();

	string csvLine;
	// read every line from the stream
	while( getline(file, csvLine) )
	// ignore lines which dont have the delimiter
	if (csvLine.find(delimiter) !=  string::npos)
		{
		istringstream csvStream(csvLine);
		vector<string> csvColumn;
		string csvElement;
		// read every element from the line that is seperated by commas
		// and put it into the vector or strings
		//            while( getline(csvStream, csvElement, ',') )
		while( getline(csvStream, csvElement, delimiter ))
			{
			csvColumn.push_back(csvElement);
			}
		output.push_back(csvColumn);
		}
	file.close();

	}


void readCSV_into_string_vector_separate_columns (string fname, char delimiter, long ignore_first_n_lines, vector< vector<string> > &csv_columns)
	{

	vector< vector<string> > csv_lines;
	readCSV_into_string_vector ( fname,  delimiter, csv_lines);

	// check if at least one line
	if (csv_lines.size() < 1)
		error(AT,FU, "csv_lines.size() < 1");

	// check if all lines have the same number of elements
	for (unsigned long il=0; il < csv_lines.size(); il++)
		if (csv_lines[il].size() != csv_lines[0].size())
			error(AT,FU, "csv_lines[il].size() != csv_lines[0].size()");

	csv_columns.assign(csv_lines[0].size(), vector <string> (0,"") );

	for (unsigned long il=ignore_first_n_lines; il < csv_lines.size(); il++)
	for (unsigned long ic=0; ic < csv_lines[il].size(); ic++)
		csv_columns[ic].push_back(csv_lines[il][ic]);

	//cout << csv_lines.size() << endl;
	//cout << csv_lines[0].size() << endl;
	//cout << csv_columns.size() << endl;
	//cout << csv_columns[0].size() << endl;


	}



// is value inside interval test
bool is_inside_interval(double lower_limit, double upper_limit, double value)
	{
	if (value >= lower_limit && value <= upper_limit) return(true);
	return(false);
	}


string extract_value_from_a_string(string s)
	{
	// remove what is in front of = sign
	s=s.substr(s.find("=")+1,s.size());
	// remove all spaces
	while (s.find(" ") != s.npos)
		s=s.erase( s.find(" "), 1 );
	//cout << "'" << s << "'" << endl;
	return(s);
	}


template<char Remove> bool BothAre(char lhs, char rhs) { return lhs == rhs && lhs == Remove; }

void Remove_repeating_whitespace_characters_from_string(string &str)
	{
	str.erase(std::unique(str.begin(), str.end(), BothAre<' '>), str.end());
	}


/*
void pause()
{
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max());   //Clear whatever's still in the buffer
    std::cout << "Press Enter to continue . . .\n";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}
*/

double deg2rad(double kot)
	{
	return (kot*M_PI/180.0);
	}

double rad2deg(double kot)
	{
	return (180.0*kot/M_PI);
	}


// Great-circle distance - from http://en.wikipedia.org/wiki/Great-circle_distance
// the The haversine formula1: http://www.movable-type.co.uk/scripts/latlong.html
// Has to be MULTIPLIED by earths radius to get distance on earth!!!!!
// angles have to be in RADIANS
double distance_on_a_sphere_specified_by_radians(double lat1, double lon1, double lat2, double lon2)
	{
	double dLat,dLon,a,c;
	dLat = (lat2-lat1);
	dLon = (lon2-lon1);
	a = sin(dLat/2) * sin(dLat/2) + sin(dLon/2) * sin(dLon/2) * cos(lat1) * cos(lat2);
	c = 2 * atan2(sqrt(a), sqrt(1-a));
	return(c);
	}

// Great-circle distance - from http://en.wikipedia.org/wiki/Great-circle_distance
// angles have to be in DEGREES
double distance_on_earth_in_km_specified_by_deg(double lat1, double lon1, double lat2, double lon2)
	{
	return(6371*distance_on_a_sphere_specified_by_radians(deg2rad(lat1),deg2rad(lon1), deg2rad(lat2), deg2rad(lon2)));
	}

// Destination point given distance and bearing from start point
// Given a start point, initial bearing, and distance, this will calculate the destination point and final bearing travelling along a (shortest distance) great circle arc.
//? is latitude, ? is longitude, ? is the bearing (clockwise from north), ? is the angular distance d/R; d being the distance travelled, R the earth’s radius
// all angles in radians
// http://www.movable-type.co.uk/scripts/latlong.html
void destination_point_given_distance_and_bearing_from_start_point(double lon1, double lat1, double bearing, double d, double r, double &lon2, double &lat2)
	{
	lat2 = asin(sin(lat1)*cos(d/r) + cos(lat1)*sin(d/r)*cos(bearing));
	lon2 = lon1 + atan2(sin(bearing)*sin(d/r)*cos(lat1), cos(d/r)-sin(lat1)*sin(lat2));

	}



double squared_euclidian_distance(double x1, double y1, double x2, double y2)
	{
	return((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	}

double euclidian_distance(double x1, double y1, double x2, double y2)
	{
	return(sqrt(squared_euclidian_distance(x1,y1,x2,y2)));
	}


double euclidian_distance_domain_periodic_in_x_dimension(double x1, double y1, double x2, double y2, int dimx)
	{
	double dx=fabs(x1-x2);
	if (dx > (double)dimx/2)
		dx=(double)dimx - dx;
	return(sqrt(dx*dx + (y2-y1)*(y2-y1)));
	}


double squared_euclidian_distance_3D(double x1, double y1, double z1, double x2, double y2, double z2)
	{
	return((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	}


bool is_location_inside_rotated_elipse(int x,int y, int x1, int y1, int r1, double e1, double fi1 )
	{
	double xx,yy;

	// coordinates in rotated system
	xx=(double)(x-x1) * cos(fi1) + (double)(y-y1)* sin(fi1);
	yy=(double)(y-y1) * cos(fi1) - (double)(x-x1)* sin(fi1);
	//if (sqr(xx-(double)x1)/sqr((double)r1) + sqr(yy-(double)y1)/sqr((double)r1*e1) <= 1.0)
	if (sqr(xx-(double)0)/sqr((double)r1) + sqr(yy-(double)0)/sqr((double)r1*e1) <= 1.0)
		return(true);
	return(false);
	}


// -----------------------------------------------------------------------------------
// sort function for sorting class with two double values by first value as index
// -----------------------------------------------------------------------------------
class Csingle_value
	{
	public:
	double x1,x2;

	Csingle_value()
		{x1=BAD_DATA_FLOAT;x2=BAD_DATA_FLOAT;}

	bool operator< (const Csingle_value &rhs) const
		{
		if (x1 < rhs.x1)
			return true;
		else
			return false;
		}
	};
class Clist_of_values_to_sort
	{
	public:
	vector <Csingle_value> single_value;

	void sort_two_vectors_of_doubles_by_first_index(vector <double> &first, vector <double> &second)
		{
		if (first.size() != second.size())
			error(AT,FU, "first.size() != second.size()");

		long il;
		single_value.clear();
		Csingle_value temp;
		if  (first.size() != 0)
			{
			// populate
			for (il=0; il < (long)first.size(); il++)
				{
				temp.x1=first[il];
				temp.x2=second[il];
				single_value.push_back(temp);
				}

			// sort
			sort(single_value.begin(), single_value.end());

			// get sorted values
			for (il=0; il < (long)first.size(); il++)
				{
				first[il]=single_value[il].x1;
				second[il]=single_value[il].x2;
				}
			}
		}
	};


void Write_append_string_to_opened_binary_file(string str, ofstream &myFile)
	{
	long c3;
	char *cstr;
	if (!myFile.is_open())
		error(AT,FU, "!myFile.is_open()");

	cstr = new char [str.size()+1];
	strcpy (cstr, str.c_str());
	c3=strlen(cstr);
	myFile.write ((char *) &c3, sizeof(long));
	myFile.write ((char *) cstr, sizeof(char)*(c3+1));
	delete[] cstr;
	}

string Read_string_from_opened_binary_file(ifstream &myFile)
	{
	string str;
	long c3;
	char *cstr;
	if (!myFile.is_open())
		error(AT,FU, "!myFile.is_open()");

	myFile.read ((char *) &c3, sizeof(long));
	cstr = new char [c3+1];
	myFile.read ((char *) cstr, sizeof(char)*(c3+1));
	str.assign(cstr);
	delete[] cstr;
	return(str);
	}


// checkig wether file exists
bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }

  return(blnReturn);
}

// creates directory (including all needed parent directories)
// checks first if directory allready exists
void mkpath(string path)
{
string command;
if (!FileExists(path))
	{
	command = "mkdir -p " + path;
	//cout << path << endl;
	system(command.c_str());
	}
}

// similar to mkpath except the "path" mush inlude the filename
void mkpath_from_path_with_filenane(string path)
{
	// find the last occurence of slash
	size_t found;
	found=path.find_last_of("/\\");

	if (found==string::npos)
		error(AT,FU, path + " does not include / !");

	mkpath(path.substr(0,found));
}


void write_string_to_new_file(string filename, string s)
	{
	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << s;
	myfile.close();
	}

string read_file_to_a_string(string filename)
	{
	std::ifstream t(filename.c_str());
	std::stringstream buffer;
	buffer << t.rdbuf();
	string s;
	s=buffer.str();
	return(s);
	}


void write_string_to_new_file_and_create_dir_if_necessary(string filename, string s)
	{
	// find the last occurence of slash
	size_t found;
	found=filename.find_last_of("/\\");
	if (found!=string::npos)
		mkpath_from_path_with_filenane(filename);

	write_string_to_new_file(filename, s);
	}


string basename_from_filename(string filename)
{
	// find the last occurence of slash
	size_t found;
	found=filename.find_last_of("/\\");

	if (found==string::npos)
		error(AT,FU, filename + " does not include / !");

	return(filename.substr(found+1));
}

string dirname_from_filename(string filename)
{
	// find the last occurence of slash
	size_t found;
	found=filename.find_last_of("/\\");

	if (found==string::npos)
		error(AT,FU, filename + " does not include / !");

	return(filename.substr(0,found+1));
}




// replace string with string  - version using boost
string replace_string_with_string(std::string s , const std::string& from, const std::string& to)
	{
	if ( from != to || from != "")
    	for(size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
        	s.replace(pos, from.size(), to);
	return (s);
	}

// replace character with string - using bool
string replace_string(string haytack, const string &needele, const string &replacement )
		{
		return(replace_string_with_string(haytack, needele, replacement));
		}


// replace string with string  - version using boost
/*string replace_string_with_string(string haytack, string needele, string replacement )
	{
    assert( needele != replacement );

    string::size_type pos = 0;
    while ( (pos = haytack.find(needele, pos)) != string::npos ) {
        haytack.replace( pos, needele.size(), replacement );
        //pos++;
        pos+=replacement.size();
    }
	    return(haytack);
	}
*/


bool is_string_an_Integer_test(string s)
	{
	char* p;
	strtol(s.c_str(), &p, 10);
	if (*p)
		return(false);
	else
		return(true);
	}

vector<string> StringExplode(string str, string separator)
	{
    unsigned long found;

    vector<string> results;

    found = str.find_first_of(separator);
    while(found != string::npos){
        if(found > 0){
            results.push_back(str.substr(0,found));
        }
        else
            results.push_back("");
        str = str.substr(found+1);
        found = str.find_first_of(separator);
    }
    if(str.length() > 0){
        results.push_back(str);
    }

	return(results);

}

// read lines of a txt file and make them string vector
vector <string> get_each_line_as_string_from_txt_file(string filename)
	{
	vector <string> out;
	string line;

	if (!FileExists(filename))
		error(AT,FU, filename + " does not exist" );

	ifstream myfile(filename.c_str());
	if (!myfile.is_open())
		error(AT,FU, "could not open file: " + filename);

	while (! myfile.eof() )
		{
		getline (myfile,line);

		// remove carrigae return if present - usualy if file is edited with windows
		line.erase( std::remove(line.begin(), line.end(), '\r'), line.end() );
		//if (line.find_first_of("\r"))
		//	line=line.substr(0,line.find_first_of("\r"));

		out.push_back(line);
		}

	myfile.close();

	return(out);
	}



template<typename T>
string merge_numeric_vector_into_string_and_put_a_separator_between(const vector <T> &vec, string separator)
	{
	ostringstream s1;
	s1.str("");
	for (long il=0; il < (long)vec.size(); il++)
		{
		s1 << vec[il];

		if (il < (long)vec.size() - 1 )
			s1 << separator;
		//cout << s1.str() << endl;
		}

	//cout << "aaa" << endl;
	//cout << s1.str() << endl;
	return(s1.str());
	}


string merge_string_vector_and_put_a_separator_between(const vector <string> &lines, string separator)
	{
	ostringstream s1;
	s1.str("");
	for (long il=0; il < (long)lines.size(); il++)
		{
		s1 << lines[il];

		if (il < (long)lines.size() - 1 )
			s1 << separator;
		//cout << s1.str() << endl;
		}

	//cout << "aaa" << endl;
	//cout << s1.str() << endl;
	return(s1.str());
	}

// returns filename withouth extension
string remove_extension(string full_filename)
		{
		string result;
	    size_t found=0;

	    found = full_filename.find_last_of(".");
		if (found == string::npos)
		{cout << "ERROR: display_withouth_extension: Could not find the extension!" << endl; exit(1);}

		result=full_filename.substr(0,found);
		return(result);
		}





// Generate a random character string  -  http://php.net/manual/en/function.rand.php
string rand_str(int length)
{
	float ran2(long *idum);

	string out;
    string chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

	int ix;

    // Generate random string
    for (ix = 0; ix < length; ix++)
	    {
        // Grab a random character from our list
        out.push_back(chars[(int)(ran2(&random_number_seed)*(double)chars.size())]);
	    }

 	// Return the string
    return (out);
}

// Generate a random character string  -  http://php.net/manual/en/function.rand.php
string rand_str_randomize_seed_by_timestamp_the_firsttime(int length)
{
	float ran2(long *idum);

	if (random_number_seed == -1)
		random_number_seed=-time(0);

	//cout << random_number_seed << endl;
	//cout << ran2(&random_number_seed) << endl;
	//cout << ran2(&random_number_seed) << endl;

	string out;
    string chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

	int ix;

    // Generate random string
    for (ix = 0; ix < length; ix++)
	    {
        // Grab a random character from our list
        out.push_back(chars[(int)(ran2(&random_number_seed)*(double)chars.size())]);
	    }

 	// Return the string
    return (out);
}



// rounds a double to digits
double round_to_digits(double x, int digits)
	{
	return(round(pow(10,(double)digits)*x)/pow(10,(double)digits));
	}

// rounds a double to digits
string convert_double_to_string(double x)
	{
	ostringstream s1;
	s1 << x;
	return(s1.str());
	}




// for long number il it adds 0 in front so the number in mest chaacters long
string output_leading_zero_string(long il,int mest)
	{
	int ix;

	ostringstream s1;
	s1.str("");
	int stevke;
	string s;

	if (il==0) stevke=1;
	else stevke=(int)floor(log10((double)il))+1;

	s1 << il;
	s=s1.str();
	s1.str("");
	for (ix=0;ix< mest-stevke;ix++)
		{
		s1 << "0" << s;
		s=s1.str();
		s1.str("");
		}
	return(s);
	}




/*// distance between two points on earth surface - has to be multiplied by earth radius to get the result
double distance_between_two_points_on_earth_surface(double dlon, double dlat)
	{
	return(2 * asin(sqrt(sin(dlat/2)*sin(dlat/2)             )));
	}
*/



double average_from_vector( const vector <double> &x)
	{
	long il;
	double ave;

	ave=0;
	for (il=0; il < (long)x.size(); il++)
		ave+=x.at(il);
	return(ave/(double)x.size());
	}

double std_deviation_from_vector( const vector <double> &x)
	{
	long il;
	double ave;
	double std_dev;

	ave=average_from_vector(x);
	std_dev=0;

	for (il=0; il < (long)x.size(); il++)
		std_dev+=(x.at(il)-ave)*(x.at(il)-ave);

	return(sqrt(std_dev/(double)x.size()));
	}

double pearsons_correlation_from_vector( const vector <double> &x, const vector <double> &y)
	{
	long il;
	double r,ax,ay,sx,sy;

	ax=average_from_vector(x);
	ay=average_from_vector(y);
	sx=std_deviation_from_vector(x);
	sy=std_deviation_from_vector(y);
	r=0;
	for (il=0; il < (long)x.size(); il++)
		r+=(x.at(il)-ax)*(y.at(il)-ay);


	if (x.size()-1 == 0)
		error(AT,FU, "x.size()-1  == 0");

	if (sx==0 || sy==0)
		return(0);

	return(r/((double)x.size()-1)/sx/sy);
	}

/*string get_product_string_from_number_CCSM(int product_number)
	{
	string out;

	if (product_number==0) out="TRMM";
	else if (product_number==1) out="CCSM";
	else
	{cout << "ERROR: get_product_string_from_number: product_number > 1 ""\n";exit(1);}

	return out;
	}
*/

string get_product_string_from_number(int product_number)
	{
	string out;

	if (product_number==0) out="TRMM";
	else if (product_number==1) out="WRF";
	else
	{cout << "ERROR: get_product_string_from_number: product_number > 1 ""\n";exit(1);}

	return out;
	}


// Class for easy conversion of dates. For example: adding a day...
// the original C++ date and time functions are too clumsy to use. The main reason being the timezones and stupid use of year (year is defined as year-1900) and month (January is stupidly = 0)
// to add a day, add 24*3600 to timestamp and convert back to date
class CMy_time_format
	{
	public:
	int year;
	int month;
	int day;
	int hour;
	int minute;
	int second;
	time_t timestamp;


	int season; // DJF=1, MAM=2, JJA=3, SON=4
	string season_string;

	string month_leading_zero;
	string day_leading_zero;
	string hour_leading_zero;

	// constructor timestamp
	CMy_time_format(time_t Ttimestamp)
		{
		timestamp=Ttimestamp;
		timestamp_to_date();
		};

	// constructor date
	CMy_time_format(int Tyear,int Tmonth,int Tday,int Thour,int Tminute,int Tsecond)
		{
		year=Tyear;
		month=Tmonth;
		day=Tday;
		hour=Thour;
		minute=Tminute;
		second=Tsecond;
		if (month < 10) month_leading_zero="0"; else  month_leading_zero="";
		if (day < 10) day_leading_zero="0"; else  day_leading_zero="";
		if (hour < 10) hour_leading_zero="0"; else  hour_leading_zero="";
		date_to_timestamp();
		set_season();
		};

	// constructor empty
	CMy_time_format()
		{
		timestamp=0;
		}

	bool compare_are_they_the_same_by_comparing_timestamp(const CMy_time_format &tt) const
		{
		if (timestamp != tt.timestamp) return(false);
		return(true);
		}

	// overload operators == nad !=
	bool operator== (const CMy_time_format& tt) const
		{return(compare_are_they_the_same_by_comparing_timestamp(tt));}
	bool operator!= (const CMy_time_format& tt) const
		{return(!compare_are_they_the_same_by_comparing_timestamp(tt));}
	bool operator< (const CMy_time_format& tt) const
		{
		if (timestamp < tt.timestamp) return(true);
		return(false);
		}
	bool operator> (const CMy_time_format& tt) const
		{
		if (*this < tt || *this==tt) return(false);
		return(true);
		}
	bool operator>= (const CMy_time_format& tt) const
		{
		if (*this < tt) return(false);
		return(true);
		}
	bool operator<= (const CMy_time_format& tt) const
		{
		if (*this > tt) return(false);
		return(true);
		}

	CMy_time_format& operator= (const CMy_time_format &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			timestamp=rhs.timestamp;
			timestamp_to_date();
			}
		return *this;    // Return ref for multiple assignment
		}

	CMy_time_format& operator+= (const long &rhs)
		{
 		timestamp+=rhs;
		timestamp_to_date();
		return *this;    // Return ref for multiple assignment
		}

	CMy_time_format& operator-= (const long &rhs)
		{
 		timestamp-=rhs;
		timestamp_to_date();
		return *this;    // Return ref for multiple assignment
		}

	const CMy_time_format operator+ (const long &rhs) const
		{
	    CMy_time_format result = *this;     // Make a copy of myself.  Same as MyClass result(*this);
	    result += rhs;            // Use += to add other to the copy.
    	return result;              // All done!
		}


	// set date
	void set_by_date(int Tyear,int Tmonth,int Tday,int Thour,int Tminute,int Tsecond)
		{
		year=Tyear;
		month=Tmonth;
		day=Tday;
		hour=Thour;
		minute=Tminute;
		second=Tsecond;
		if (month < 10) month_leading_zero="0"; else  month_leading_zero="";
		if (day < 10) day_leading_zero="0"; else  day_leading_zero="";
		if (hour < 10) hour_leading_zero="0"; else  hour_leading_zero="";
		date_to_timestamp();
		set_season();
		}


	// set date
	void set_by_timestamp(long temp_timestamp)
		{
		timestamp=temp_timestamp;
		timestamp_to_date();
		}

	// set by_hours_since_1900_01_01 - this is the ECMWF format for data
	void set_by_hours_since_1900_01_01(double hours_since_1900_01_01)
		{
		// time_in_data is in hours since 1900-01-01 00:00:0.0
		CMy_time_format ttemp(1970,1,1,0,0,0);
		// substract 613608 hours - hours betwwen 1900 and 1970
		double time_in_data=hours_since_1900_01_01-613608;
		set_by_timestamp(ttemp.timestamp +  time_in_data*3600);
		}


	// converts the intigers representing the date to timestamp
	void date_to_timestamp()
		{
		tm now;
		now.tm_year=year-1900;
		now.tm_mon=month-1;
		now.tm_mday=day;
		now.tm_hour=hour;
		now.tm_min=minute;
    	now.tm_sec=second;
	    timestamp=timegm(&now);
		}

	// converts timestamp to date
	void timestamp_to_date()
		{
		tm *now = gmtime(&timestamp);
		year=now->tm_year+1900;
		month=now->tm_mon+1;
		day=now->tm_mday;
		hour=now->tm_hour;
		minute=now->tm_min;
		second=now->tm_sec;
		if (month < 10) month_leading_zero="0"; else  month_leading_zero="";
		if (day < 10) day_leading_zero="0"; else  day_leading_zero="";
		if (hour < 10) hour_leading_zero="0"; else  hour_leading_zero="";
		set_season();
		}



	string output_date_string(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "_" << output_leading_zero_string(month,2) << "_" << output_leading_zero_string(day,2) << "_" << output_leading_zero_string(hour,2) << ":" << output_leading_zero_string(minute,2) << ":" << output_leading_zero_string(second,2);
		return (s1.str());
		}

	string output_date_string_with_custom_seperator(string sep) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "_" << output_leading_zero_string(month,2) << "_" << output_leading_zero_string(day,2) << "_" << output_leading_zero_string(hour,2) << "_" << output_leading_zero_string(minute,2) << "_" << output_leading_zero_string(second,2);
		return (s1.str());
		}

	string output_date_string_to_minute(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "_" << output_leading_zero_string(month,2) << "_" << output_leading_zero_string(day,2) << "_" << output_leading_zero_string(hour,2) << "_" << output_leading_zero_string(minute,2);
		return (s1.str());
		}

	string output_date_string_to_hour(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "_" << output_leading_zero_string(month,2) << "_" << output_leading_zero_string(day,2) << "_" << output_leading_zero_string(hour,2) ;
		return (s1.str());
		}


	string output_date_string_to_day(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "_" << output_leading_zero_string(month,2) << "_" << output_leading_zero_string(day,2);
		return (s1.str());
		}

	string yyyy()
		{return(output_leading_zero_string(year,4));}
	string mm()
		{return(output_leading_zero_string(month,2));}
	string dd()
		{return(output_leading_zero_string(day,2));}
	string hh()
		{return(output_leading_zero_string(hour,2));}
	string minmin()
		{return(output_leading_zero_string(minute,2));}


	string output_date_string_to_hour_no_leading_zeroes(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << year << "_" << month << "_" << day << "_" <<hour ;
		return (s1.str());
		}


	string output_date_string_to_day_with_minuses(void) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << "-" << output_leading_zero_string(month,2) << "-" << output_leading_zero_string(day,2);
		return (s1.str());
		}

	string output_date_string_to_day_with_custom_separator(string sep) const
		{
		ostringstream s1;
		s1.str("");
		s1 << output_leading_zero_string(year,4) << sep << output_leading_zero_string(month,2) << sep << output_leading_zero_string(day,2);
		return (s1.str());
		}

	void set_by_julian_day(double julianday)
		{
		// substract linux epoch expressed in julianday notantion = 2440587.5;
		timestamp=(long)((julianday - 2440587.5)*3600*24);
		timestamp_to_date();
		}

	void set_by_modified_julian_day(double modified_julianday)
		{
		// http://en.wikipedia.org/wiki/Julian_day
		set_by_julian_day(modified_julianday + 2400000.5);
		}




	void set_season()
		{
		if (month==12 || month == 1 || month == 2) {season=1; season_string="DJF";}
		if (month== 3 || month == 4 || month == 5) {season=2; season_string="MAM";}
		if (month== 6 || month == 7 || month == 8) {season=3; season_string="JJA";}
		if (month== 9 || month == 10 || month == 11) {season=4; season_string="SON";}
		if (month== 0) {season=0; season_string="annual";}
		}

	string get_month_string_from_month_number(int mnumber) const
	    {
		string out;
		     if (mnumber==1) out="JAN";
		else if (mnumber==2) out="FEB";
		else if (mnumber==3) out="MAR";
		else if (mnumber==4) out="APR";
		else if (mnumber==5) out="MAY";
		else if (mnumber==6) out="JUN";
		else if (mnumber==7) out="JUL";
		else if (mnumber==8) out="AUG";
		else if (mnumber==9) out="SEP";
		else if (mnumber==10) out="OCT";
		else if (mnumber==11) out="NOV";
		else if (mnumber==12) out="DEC";
		else
			error(AT,FU, "Wrong mnumber!");

	    return (out);
	    }

	  string get_month_string_from_month_number_include_annual_for_zero(int mnumber) const
	    {
		string out;
	    if (mnumber==0) out="Annual";
		else out=get_month_string_from_month_number(mnumber);
	    return (out);
	    }


	  string get_season_string_from_month_number(int mnumber) const
	    {
	    CMy_time_format begin(2000,mnumber,1,3,0,0);
	    return (begin.season_string);
	    }

	  int get_season_number_from_month_number(int mnumber) const
	    {
	      CMy_time_format begin(2000,mnumber,1,3,0,0);
	      return (begin.season);
	    }

	  string get_season_sting_from_season_number(int snumber) const
	     {
	       CMy_time_format begin(2000,1,1,3,0,0);

	       int ix=1;
	       while (begin.season != snumber && ix <= 12)
		 {
		   begin.set_by_date(2000,ix,1,3,0,0);
		   ix++;
		 }

	       if (snumber==0) begin.season_string="Annual";

	       return (begin.season_string);
	     }



	};


string output_date_string_from_timestamp(long timestamp)
	{
	CMy_time_format tt;
	tt.set_by_timestamp(timestamp);
	return (tt.output_date_string());
	}


string output_date_string_to_hour_from_timestamp(long timestamp)
	{
	CMy_time_format tt;
	tt.set_by_timestamp(timestamp);
	return (tt.output_date_string_to_hour());
	}


int number_of_days_in_a_month(int year, int month)
	{
	CMy_time_format t1(year,month,1,0,0,0);
	CMy_time_format t2,t3;
	t2.set_by_timestamp(t1.timestamp + 3600*24*35);
	t3.set_by_date(t2.year,t2.month,1,0,0,0);
	return(int((t3.timestamp-t1.timestamp)/3600/24));
	}

int day_of_the_year(long timestamp)
	{
	CMy_time_format t1,t2;
	t1.set_by_timestamp(timestamp);
	t2.set_by_date(t1.year,1,1,0,0,0);
	return(int((t1.timestamp-t2.timestamp)/3600/24));
	}





// Calculates how many seconds elapsed in a specific month in a time period. For example, how many January seconds passed in period 2000.3 - 2005.12 - there were 5 Januarys in this period so the answer is 31*5days. It is more complicated for February since days wary from year to year.
long Seconds_in_a_month(int start_year, int start_month, int end_year, int end_month, int month)
	{
	CMy_time_format temp, next_month,first_day,last_day,after_last_month;
	long sum=0;


	temp.set_by_date(end_year,end_month,15,0,0,0);
	after_last_month.timestamp=temp.timestamp + 30*24*3600;
	after_last_month.timestamp_to_date();

	temp.set_by_date(start_year,start_month,15,0,0,0);
	while (temp.year != after_last_month.year || temp.month != after_last_month.month )
		{
		next_month.timestamp=temp.timestamp + 30*24*3600;
		next_month.timestamp_to_date();

		first_day.set_by_date(temp.year,temp.month,1,0,0,0);
		last_day.set_by_date(next_month.year,next_month.month,1,0,0,0);

		if (temp.month==month)
			sum+=(last_day.timestamp-first_day.timestamp);

//		cout << temp.year << " " << temp.month << " " << sum  <<  endl;

		temp.set_by_date(next_month.year,next_month.month,15,0,0,0);
		}

	return (sum);
}


// display start and end times of program run
void timer(CMy_time_format &timer_start)
	{

	// start timer if timestamp not set and display time
	if (timer_start.timestamp == 0)
		{
		timer_start.timestamp = time(NULL);
		timer_start.timestamp_to_date();

		cout << "Timer START: ";
		if (timer_start.hour < 10) cout << "0";
		cout << timer_start.hour << ":";
		if (timer_start.minute < 10) cout << "0";
		cout << timer_start.minute << ":";
		if (timer_start.second < 10) cout << "0";
		cout << timer_start.second << endl;
		}

	else
		{
		CMy_time_format tt;
		tt.timestamp = time(NULL);
		tt.timestamp_to_date();

		cout << "Timer TIME:  ";
		if (tt.hour < 10) cout << "0";
		cout << tt.hour << ":";
		if (tt.minute < 10) cout << "0";
		cout << tt.minute << ":";
		if (tt.second < 10) cout << "0";
		cout << tt.second;

		long all_seconds = tt.timestamp - timer_start.timestamp;
		long hours = all_seconds / 3600;
		long minutes = (all_seconds  - hours * 3600)/60;
		long seconds = all_seconds  - hours * 3600 - minutes * 60;

		double all_minutes = (double)all_seconds/60;

		cout << "     " << hours << "h " << minutes << "min " << seconds << "sec                 " << all_minutes << "m"<< endl;
		}
	}

// display start and end times of program run
void timer2(CMy_time_format &timer_start, CMy_time_format &timer_last_checkpoint, string text)
	{
	// start timer if timestamp not set and display time
	if (timer_start.timestamp == 0)
		{
		timer_start.timestamp = time(NULL);
		timer_start.timestamp_to_date();
		timer_last_checkpoint= timer_start;

		cout << "Timer START: ";
		if (timer_start.hour < 10) cout << "0";
		cout << timer_start.hour << ":";
		if (timer_start.minute < 10) cout << "0";
		cout << timer_start.minute << ":";
		if (timer_start.second < 10) cout << "0";
		cout << timer_start.second << endl;
		}

	else
		{
		CMy_time_format tt;
		tt.timestamp = time(NULL);
		tt.timestamp_to_date();

		cout << "Timer TIME:  ";
		if (tt.hour < 10) cout << "0";
		cout << tt.hour << ":";
		if (tt.minute < 10) cout << "0";
		cout << tt.minute << ":";
		if (tt.second < 10) cout << "0";
		cout << tt.second;

		long all_seconds = tt.timestamp - timer_start.timestamp;
		long hours = all_seconds / 3600;
		long minutes = (all_seconds  - hours * 3600)/60;
		long seconds = all_seconds  - hours * 3600 - minutes * 60;

		double all_minutes = (double)all_seconds/60;
		double all_minutes_from_last_checkup = (double)(tt.timestamp - timer_last_checkpoint.timestamp)/60;
		//double all_seconds_from_last_checkup = tt.timestamp - timer_last_checkpoint.timestamp;

		cout << "     " << hours << "h " << minutes << "min " << seconds << "sec                 " << all_minutes << "m    LAST:"<< all_minutes_from_last_checkup <<"m OR " << tt.timestamp - timer_last_checkpoint.timestamp << "s" << endl;

		timer_last_checkpoint=tt;
		}
	}

// Timer to measure the execution time in miliseconds
// by Andreas Bonini
// http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
typedef long long int64;
typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64()
{
#ifdef _WIN32
 /* Windows */
 FILETIME ft;
 LARGE_INTEGER li;

 /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
  * to a LARGE_INTEGER structure. */
 GetSystemTimeAsFileTime(&ft);
 li.LowPart = ft.dwLowDateTime;
 li.HighPart = ft.dwHighDateTime;

 uint64 ret = li.QuadPart;
 ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
 ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

 return ret;
#else
 /* Linux */
 struct timeval tv;

 gettimeofday(&tv, NULL);

 uint64 ret = tv.tv_usec;
 /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
 ret /= 1000;

 /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
 ret += (tv.tv_sec * 1000);

 return ret;
#endif
}


///-------------------------------------------
///  NRC Functions
///
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

//static double dsqrarg;
//#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
//static double dmaxarg1,dmaxarg2;
//#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
//
//static double dminarg1,dminarg2;
//#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?  (dminarg1) : (dminarg2))
//
//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
//
//static float minarg1,minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?  (minarg1) : (minarg2))
//
//static long lmaxarg1,lmaxarg2;
//#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?   (lmaxarg1) : (lmaxarg2))
//
//static long lminarg1,lminarg2;
//#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?   (lminarg1) : (lminarg2))
//
//static int imaxarg1,imaxarg2;
//#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?   (imaxarg1) : (imaxarg2))
//
//static int iminarg1,iminarg2;
//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?   (iminarg1) : (iminarg2))
//
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran2(long *idum)
/*Long period (> 2  1018) random number generator of L'Ecuyer with Bays-Durham shue
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.*/
{
int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0)
{
    // Initialize.
    if (-(*idum) < 1) *idum=1;
    // Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;
    j>=0;
    j--)
    {
        // Load the shue table (after 8 warm-ups).
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
}
k=(*idum)/IQ1;
/// Start here when not initializing.
*idum=IA1*(*idum-k*IQ1)-k*IR1;
//Compute idum=(IA1*idum) % IM1 without
if (*idum < 0) *idum += IM1;
//overflows by Schrage's method.
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
//Compute idum2=(IA2*idum) % IM2 likewise.
if (idum2 < 0) idum2 += IM2;
j=iy/NDIV;
// Will be in the range 0..NTAB-1.
iy=iv[j]-idum2;
// Here idum is shued, idum and idum2 are
iv[j] = *idum;
// combined to generate output.
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;
// Because users don't expect endpoint values.
else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
// no_static_variables_for_use_with_openMP
float ran2_no_static_variables_for_use_with_openMP(long *idum, long &idum2, long &iy, vector <long> &iv)
{
int j;
long k;
//static long idum2=123456789;
//static long iy=0;
//static long iv[NTAB];
float temp;
if (*idum <= 0)
{
    // Initialize.
    if (-(*idum) < 1) *idum=1;
    // Be sure to prevent idum = 0.
    else *idum = -(*idum);
	iv.assign(NTAB,0);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--)
    {
        // Load the shue table (after 8 warm-ups).
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
}
k=(*idum)/IQ1;
/// Start here when not initializing.
*idum=IA1*(*idum-k*IQ1)-k*IR1;
//Compute idum=(IA1*idum) % IM1 without
if (*idum < 0) *idum += IM1;
//overflows by Schrage's method.
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
//Compute idum2=(IA2*idum) % IM2 likewise.
if (idum2 < 0) idum2 += IM2;
j=iy/NDIV;
// Will be in the range 0..NTAB-1.
iy=iv[j]-idum2;
// Here idum is shued, idum and idum2 are
iv[j] = *idum;
// combined to generate output.
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;
// Because users don't expect endpoint values.
else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#include <math.h>
#define NRANSI

void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
	float *b, float *siga, float *sigb, float *chi2, float *q)
{
	float gammq(float a, float x);
	int i;
	float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	if (mwt == 0) {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=1;i<=ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		*q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}
#undef NRANSI

float gammq(float a, float x)
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;

	char msg[]="Invalid arguments in routine gammq";
	if (x < 0.0 || a <= 0.0) nrerror(msg);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(float *gammcf, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}

	char msg[]="a too large, ITMAX too small in gcf";
	if (i > ITMAX) nrerror(msg);
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

#include <math.h>

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		char msg1[]="x less than 0 in routine gser";
		if (x < 0.0) nrerror(msg1);
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		char msg2[]="a too large, ITMAX too small in routine gser";
		if (x < 0.0) nrerror(msg2);
		return;
	}
}
#undef ITMAX
#undef EPS

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


void generate_random_color_based_on_seed(long idum, double &r, double &g,  double &b)
	{
	r=b=g=0;
	while (r+g+b < 0.5)
		{
		r=ran2(&idum);
		g=ran2(&idum);
		b=ran2(&idum);
		}
	}



void generate_random_binary_vector_with_fixed_number_of_1_values(long p,long N,vector <double> &v, long &seed)
	{
	v.clear();
	if (N < p)
		error(AT,FU, "N < p");

	long count=0;
	long pos;

	double to_set;
	double to_fill_on_start;
	long how_many_to_add;

	if ((double)p/(double)N < 0.5)
		{
		to_set=1;
		to_fill_on_start=0;
		how_many_to_add=p;
		}
	else
		{
		to_set=0;
		to_fill_on_start=1;
		how_many_to_add=N-p;
		}

	v.assign(N,to_fill_on_start);

	while (count < how_many_to_add)
			{
			pos = floor(ran2(&seed)*(double)N);
			if (pos > N - 1)
				error(AT,FU, "pos > N - 1");

			if (v[pos]==to_fill_on_start)
				{
				v[pos]=to_set;
				count++;
				}
			}

	// make a final check
	if (sum_vector(v) != p)
		error(AT,FU, "sum_vector(v) != p");
	}


void generate_random_binary_vector_with_fixed_number_of_1_values_openMP_friendly(long p,long N,vector <double> &v, long &idum, long &idum2, long &iy, vector <long> &iv)
	{
	v.clear();
	if (N < p)
		error(AT,FU, "N < p");

	long count=0;
	long pos;

	double to_set;
	double to_fill_on_start;
	long how_many_to_add;

	if ((double)p/(double)N < 0.5)
		{
		to_set=1;
		to_fill_on_start=0;
		how_many_to_add=p;
		}
	else
		{
		to_set=0;
		to_fill_on_start=1;
		how_many_to_add=N-p;
		}

	v.assign(N,to_fill_on_start);

	while (count < how_many_to_add)
			{
			pos = floor(ran2_no_static_variables_for_use_with_openMP(&idum, idum2, iy, iv)*(double)N);
			if (pos > N - 1)
				error(AT,FU, "pos > N - 1");

			if (v[pos]==to_fill_on_start)
				{
				v[pos]=to_set;
				count++;
				}
			}

	// make a final check
	if (sum_vector(v) != p)
		error(AT,FU, "sum_vector(v) != p");
	}


void generate_random_binary_vector(double f,long n,vector <double> &v, long &seed)
	{
	v.clear();
	for (long il=0; il<n; il++)
		{
		if (ran2(&seed) < f)
			v.push_back(1);
		else
			v.push_back(0);
		}
	}

void generate_random_binary_vector_openMP_friendly(double f,long n,vector <double> &v, long &idum, long &idum2, long &iy, vector <long> &iv)
	{
	v.clear();
	for (long il=0; il<n; il++)
		{
		if (ran2_no_static_variables_for_use_with_openMP(&idum, idum2, iy, iv) < f)
			v.push_back(1);
		else
			v.push_back(0);
		}
	}

// this function draws the pixel with gd and you dont have to worry about assigning new colors
// the figure is inverted in y dimension - so it is the same as in pngwriter
// r,g,b components should be in the range 0-1
void gd_draw_pixel(gdImagePtr image, int x, int y, double r, double g, double b)
	{
	int color;
	color=gdImageColorExact(image,(int)(r*255.),(int)(g*255.),(int)(b*255.));
	// if color not defined previously define it now
	if (color == -1)
	 	color=gdImageColorAllocate(image, (int)(r*255.),(int)(g*255.),(int)(b*255.));
	gdImageSetPixel( image , x, image->sy - y - 1, color );
	}

// write gd image to jpg or png file
void Write_gd_image_to_file(gdImagePtr image, string filename, string file_type)
	{
	FILE *out;
	out = fopen(filename.c_str(), "wb");

	if (file_type == "jpg")
		gdImageJpeg(image, out, 80);

	if (file_type == "png")
		gdImagePng(image, out);

	fclose(out);
	}

// write gd image to jpg or png file
void Write_gd_image_to_png(gdImagePtr image, string filename)
	{
	FILE *out;
	out = fopen(filename.c_str(), "wb");
	gdImagePng(image, out);
	fclose(out);
	}

void Write_y_inverted_gd_image_to_file(gdImagePtr image, string filename, string file_type)
	{
	int ix, iy, c;
	gdImagePtr im_out;
	im_out = gdImageCreateTrueColor(image->sx, image->sy);
	for (ix=0; ix < image->sx; ix++)
		for (iy=0; iy < image->sy; iy++)
			{
			c = gdImageGetPixel(image,ix, iy);
			gdImageSetPixel  ( im_out , ix, image->sy - iy -1, c );
			}
	Write_gd_image_to_file(im_out, filename, file_type);
	gdImageDestroy(im_out);
	}


