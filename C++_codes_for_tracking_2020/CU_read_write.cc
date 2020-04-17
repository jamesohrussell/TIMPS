void WriteField_to_CSV_file_and_create_dir_if_necessary(string fname,  const CField &f )
	{
	ostringstream s1,s2,s3,s4,s5,s6,s7,dataset;

	if (!f.check_if_initilized())
		error(AT,FU, "f.check_if_initilized()");

	s2.str("");
	for (long iy=0; iy < f.dimy; iy++)
		{
		for (long ix=0; ix < f.dimx; ix++)
			{
			s2 << f.get(ix,iy);
			if (ix < f.dimx - 1)
				s2 << ";";
			}
		s2 << endl;
		}
	write_string_to_new_file_and_create_dir_if_necessary(fname, s2.str());
	}


CField ReadField_from_CSV_file_no_lat_lon_data(string fname)
	{
	vector< vector<string> > data_str;
	readCSV_into_string_vector (fname, ';', data_str);

	CField f;
	f.dimy=data_str.size();
	f.dimx=data_str[0].size();
	f.allocatememory();
	f.allocatememory_lat_lon_array();

	for (int ix=0; ix < f.dimx; ix++)
	for (int iy=0; iy < f.dimy; iy++)
		f.set(ix,iy,atof(data_str[iy][ix].c_str()));

	return(f);
	}

CField ReadField_from_CSV_file_no_lat_lon_data_no_newline(string fname, int dimx, int dimy)
	{
	vector< vector<string> > data_str;
	readCSV_into_string_vector (fname, ';', data_str);

	if ((long)data_str[0].size() != (long)dimx*(long)dimy)
		error(AT,FU, "(long)data_str[0].size() != (long)dimx*(long)dimy");

	CField f;
	f.dimy=dimy;
	f.dimx=dimx;
	f.allocatememory();
	f.allocatememory_lat_lon_array();

	for (long il=0; il < f.size(); il++)
		f.data[il]=atof(data_str[0][il].c_str());

	return(f);
	}



#ifndef DISABLE_NETCDF4
void WriteField_compressed_with_deflate_level(const string fname, const CField &f , int deflate_level)
	{
    int ncid, x_dimid, y_dimid, varid1, varid2, varid3;
    int dimids[3];
    int retval;
	//int ix;

	//cout << f.lat_data.size() << " " << f.lon_data.size() << " " << f.data.size() << endl;
	// check the basics
	if (!f.check_if_initilized())
		error(AT,FU, "f.check_if_initilized()");

	// smaller than zero longitudes check
	//for (ix=0; ix < f.dimx; ix++)
	//	if (f.lon_data[ix] < 0)
	//		error(AT,FU, "f.lon_data[ix] < 0 ");

    // Create the file. The NC_NETCDF4 flag tells netCDF to create a netCDF-4/HDF5 file.
    if ((retval = nc_create(fname.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid)))
    	{
		cout << "error: " << fname << endl;
		error(AT,FU, "retval = nc_create(fname.c_str(), NC_NOCLOBBER|NC_NETCDF4, &ncid)");

		}

	// Define the dimensions. NetCDF will hand back an ID for each.
    if ((retval = nc_def_dim(ncid, "lat", f.dimy, &y_dimid)))
		error(AT,FU, "retval");
    if ((retval = nc_def_dim(ncid, "lon", f.dimx, &x_dimid)))
		error(AT,FU, "retval");


     //The dimids array is used to pass the IDs of the dimensions of the variable.
	dimids[0] = y_dimid;
	dimids[1] = x_dimid;

	//Define the variable. The type of the variable in this case is NC_DOUBLE
	if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimids[0], &varid1)))
		error(AT,FU, "retval");
	if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimids[1], &varid2)))
		error(AT,FU, "retval");
	if ((retval = nc_def_var(ncid, "precip", NC_DOUBLE, 2, dimids, &varid3)))
		error(AT,FU, "retval");

	double missing=BAD_DATA_FLOAT;

	// Assign units attributes to coordinate variables.
	if ((retval = nc_put_att_text(ncid, varid1, "units", strlen("degrees_north"), "degrees_north")))
		error(AT,FU, "retval");
	if ((retval = nc_put_att_double(ncid, varid1, "_FillValue", NC_DOUBLE, 1,&missing)))
		error(AT,FU, "retval");
	if ((retval = nc_put_att_text(ncid, varid2, "units", strlen("degrees_east"), "degrees_east")))
		error(AT,FU, "retval");
	if ((retval = nc_put_att_double(ncid, varid2, "_FillValue", NC_DOUBLE, 1,&missing)))
		error(AT,FU, "retval");
	if ((retval = nc_put_att_double(ncid, varid3, "_FillValue", NC_DOUBLE, 1,&missing)))
		error(AT,FU, "retval");

	// compress variable
	if (deflate_level != 0)
		nc_def_var_deflate(ncid, varid3, 0, 1, deflate_level);

	//End define mode. This tells netCDF we are done defining metadata.
	if ((retval = nc_enddef(ncid)))
		error(AT,FU, "retval");

	// Write the data to the file. Although netCDF supports reading and writing subsets of data, in this case we write all the data in one operation.
	if ((retval = nc_put_var_double(ncid, varid1, &f.lat_data.at(0) )))
		error(AT,FU, "retval");
	if ((retval = nc_put_var_double(ncid, varid2, &f.lon_data.at(0) )))
		error(AT,FU, "retval");
	if ((retval = nc_put_var_double(ncid, varid3, &f.data.at(0) )))
		error(AT,FU, "retval");

	// Close the file. This frees up any internal netCDF resources associated with the file, and flushes any buffers.
	if ((retval = nc_close(ncid)))
		error(AT,FU, "retval");
	}
#endif

#ifndef DISABLE_NETCDF4
void WriteField(const string fname, const CField &f)
	{
	WriteField_compressed_with_deflate_level(fname, f , 0);
	}
#endif

#ifndef DISABLE_NETCDF4
void WriteField_compressed(const string fname, const CField &f)
	{
	WriteField_compressed_with_deflate_level(fname, f , 3);
	}
#endif

// same as WriteField but also create the necessary path if it does not allready exists
#ifndef DISABLE_NETCDF4
void WriteField_and_create_path(const string fname, const CField &f )
	{
	mkpath_from_path_with_filenane(fname);
	WriteField(fname,f);
	}
#endif

// same as WriteField_and_create_path creates gziped version of the file
#ifndef DISABLE_NETCDF4
void WriteField_and_create_path_and_gzip(const string fname, const  CField &f )
	{
	mkpath_from_path_with_filenane(fname);
	WriteField(fname,f);
	string command;
	command = "gzip -f " + fname;
	system(command.c_str());
	}
#endif

// same as WriteField but also create the necessary path if it does not allready exists
#ifndef DISABLE_NETCDF4
void WriteField_compressed_and_create_path(const string fname, const CField &f )
	{
	mkpath_from_path_with_filenane(fname);
	WriteField_compressed(fname,f);
	}
#endif

// ******************************************************
// ******************************************************
// BRANJE iz NetCDF fajla, ki se lahko dajo v MODE tool (notri je le eno polje)
//

#ifndef DISABLE_NETCDF4
string get_filename_path_from_opened_netcdf_file(const int ncid)
	{
	unsigned long MAX_FILENAME_LENGTH=1024;

	int  status;
	// get filename
	size_t pathlen;;
	status = nc_inq_path (ncid, &pathlen, NULL);
	if (status != NC_NOERR)
		error(AT,FU, "Can't read file path for netcdf file with specified ncid");

	char path[MAX_FILENAME_LENGTH];
	if (pathlen > MAX_FILENAME_LENGTH -1)
		error(AT,FU, "netcdf filename/path too long");

	status = nc_inq_path (ncid, &pathlen, path);
	if (status != NC_NOERR)
		error(AT,FU, "Can't read file path for netcdf file with specified ncid");

	return((string)path);
	}
#endif

#ifndef DISABLE_NETCDF4
void inquire_about_netcdf_variable_from_opened_netcdf_file(const int ncid, const string var_name, vector <long> &var_dim_sizes, vector <string> &var_dim_names, long &total_variable_size)
	{
	int  status;

	string fname=get_filename_path_from_opened_netcdf_file(ncid);

	// get varible id
	int var_id;
	status = nc_inq_varid (ncid, var_name.c_str(), &var_id);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

	// get info about a variable
	nc_type var_type;
	int var_dimids[NC_MAX_VAR_DIMS];
	int var_ndims;
	int var_natts;
	status = nc_inq_var (ncid, var_id, 0, &var_type, &var_ndims, var_dimids, &var_natts);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

  	char recname[NC_MAX_NAME+1];
	size_t length;

	var_dim_sizes.clear();
	var_dim_names.clear();

	// loop over all variable dimensions
	for (long id=0; id < var_ndims; id++)
		{
		status = nc_inq_dim( ncid, var_dimids[id], recname, &length);
		if (status != NC_NOERR)
			error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

		var_dim_sizes.push_back(length);
		var_dim_names.push_back((string)recname);
		}

	// get total variable size
	total_variable_size=1;
	for (unsigned long id=0; id < var_dim_sizes.size(); id++)
		total_variable_size*=var_dim_sizes[id];
	}
#endif


#ifndef DISABLE_NETCDF4
void Read_whole_netcdf_double_variable_from_opened_netcdf_file(const int ncid, const string var_name,  vector <double> &var_data, vector <long> &var_dim_sizes, vector <string> &var_dim_names)
	{
	int  status;

	string fname=get_filename_path_from_opened_netcdf_file(ncid);

	// get dimension info about a variable
	long total_variable_size;
	inquire_about_netcdf_variable_from_opened_netcdf_file(ncid, var_name , var_dim_sizes, var_dim_names, total_variable_size);

	// get varible id
	int var_id;
	status = nc_inq_varid (ncid, var_name.c_str(), &var_id);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

	var_data.assign(total_variable_size,0);
	status = nc_get_var_double(ncid, var_id, &var_data[0]);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);
	}
#endif

#ifndef DISABLE_NETCDF4
void Read_subset_of_netcdf_double_variable_from_opened_netcdf_file(const int ncid, const string var_name,  vector <double> &var_data, vector <long> &var_dim_sizes, vector <string> &var_dim_names, const vector <long> start, const vector <long> length)
	{
	int  status;

	string fname=get_filename_path_from_opened_netcdf_file(ncid);

	// get dimension info about a variable
	long total_variable_size;
	inquire_about_netcdf_variable_from_opened_netcdf_file(ncid, var_name , var_dim_sizes, var_dim_names, total_variable_size);

	// get varible id
	int var_id;
	status = nc_inq_varid (ncid, var_name.c_str(), &var_id);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

	if (start.size() == 0 || length.size() == 0)
		error(AT,FU, "start.size() == 0 || length.size() == 0");

	if (start.size() != length.size())
		error(AT,FU, "start.size() != length.size()");

	if (start.size() != var_dim_sizes.size())
		error(AT,FU, "start.size() != var_dim_sizes.size()");


	vector <size_t> start_size_t;
	for (unsigned long il=0; il < start.size(); il++)
		start_size_t.push_back(start[il]);

	vector <size_t> length_size_t;
	for (unsigned long il=0; il < length.size(); il++)
		{
		if (length[il] < 0)
			length_size_t.push_back(var_dim_sizes[il]);
		else
			length_size_t.push_back(length[il]);
		}

	// get total variable size
	total_variable_size=1;
	for (unsigned long id=0; id < length_size_t.size(); id++)
		total_variable_size*=length_size_t[id];
	if (total_variable_size < 1)
		error(AT,FU, "total_variable_size < 1");
	var_data.assign(total_variable_size,0);

	//cout << output_vector_as_string_with_index_and_newline(start_size_t);
	//cout << output_vector_as_string_with_index_and_newline(length_size_t);
	//cout << total_variable_size << endl;

	status = nc_get_vara_double(ncid, var_id, &start_size_t[0], &length_size_t[0], &var_data[0]);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

	}
#endif


#ifndef DISABLE_NETCDF4
CField ReadField(const string fname)
	{
	int  status;
	int ncid;

	status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't open input file: " + fname);

	vector <long> var_dim_sizes;
	vector <string> var_dim_names;
	vector <double> var_data;
	vector <double> lat_data;
	vector <double> lon_data;

	// read variable precip
	Read_whole_netcdf_double_variable_from_opened_netcdf_file( ncid, "precip",  var_data, var_dim_sizes, var_dim_names);

	if (var_dim_sizes.size() != 2)
		error(AT,FU, "Error openning input file: " + fname + ". var_dim_sizes.size() != 2");

	CField f;
	f.dimy=var_dim_sizes[0];
	f.dimx=var_dim_sizes[1];
	f.data=var_data;

	// read lat lon
	Read_whole_netcdf_double_variable_from_opened_netcdf_file( ncid, "lat",  lat_data, var_dim_sizes, var_dim_names);
	if (var_dim_sizes.size() != 1)
		error(AT,FU, "Error openning input file: " + fname + ". var_dim_sizes.size() != 1");
	Read_whole_netcdf_double_variable_from_opened_netcdf_file( ncid, "lon",  lon_data, var_dim_sizes, var_dim_names);
	if (var_dim_sizes.size() != 1)
		error(AT,FU, "Error openning input file: " + fname + ". var_dim_sizes.size() != 1");

	f.lat_data=lat_data;
	f.lon_data=lon_data;

	status = nc_close(ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't close input file: " + fname);

	return (f);
	}
#endif


#ifndef DISABLE_NETCDF4
double Read_attribute_as_double(string fname, string var_name, string attribute_name)
	{
	int  status;
	int ncid;

	status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't open input file: " + fname);

	// get varible id
	int var_id;
	status = nc_inq_varid (ncid, var_name.c_str(), &var_id);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ var_name +" variable from netcdf file: " + fname);

	double attribute_value;
	status = nc_get_att_double(ncid, var_id, attribute_name.c_str(), &attribute_value);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't read "+ attribute_name +" attribute from netcdf file: " + fname);

	status = nc_close(ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't close input file: " + fname);

	return(attribute_value);
	}
#endif

#ifndef DISABLE_NETCDF4
void Read_double_vector_from_netcdf_file(string fname, string variable_name, vector <double> &out)
	{
	int  status;
	int ncid;

	status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't open input file: " + fname);

	vector <long> var_dim_sizes;
	vector <string> var_dim_names;

	// read variable precip
	Read_whole_netcdf_double_variable_from_opened_netcdf_file( ncid, variable_name,  out, var_dim_sizes, var_dim_names);

	if (var_dim_sizes.size() != 1)
		error(AT,FU, "Error openning input file: " + fname + ". var_dim_sizes.size() != 1");

	status = nc_close(ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't close input file: " + fname);
	}
#endif

#ifndef DISABLE_NETCDF4
CField Read_general_2D_variable_netcdf(string fname, int timenum, string var_name)
	{
	int  status;
	int ncid;

	status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't open input file: " + fname);

	vector <long> var_dim_sizes;
	vector <string> var_dim_names;
	vector <double> var_data;

	vector <long> start;
	start.push_back(timenum);
	start.push_back(0);
	start.push_back(0);

	vector <long> length;
	length.push_back(1);
	length.push_back(-1);
	length.push_back(-1);

	Read_subset_of_netcdf_double_variable_from_opened_netcdf_file( ncid, var_name,  var_data, var_dim_sizes, var_dim_names, start, length);

	if (var_dim_sizes.size() != 3)
		error(AT,FU, "var_dim_sizes.size() != 3");

	if (var_dim_sizes[1]*var_dim_sizes[2] != (long)var_data.size())
		error(AT,FU, "var_dim_sizes[1] *var_dim_sizes[2] != var_data.size()");

	CField f;
	f.dimy=var_dim_sizes[1];
	f.dimx=var_dim_sizes[2];
	f.data=var_data;

	status = nc_close(ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't close input file: " + fname);

	return (f);
	}
#endif


#ifndef DISABLE_NETCDF4
CField Read_general_netcdf_file(string fname, string lat_var_name, string lon_var_name, string var_name, int timenum)
	{
	CField f= Read_general_2D_variable_netcdf(fname, timenum, var_name);
	Read_double_vector_from_netcdf_file(fname, lat_var_name, f.lat_data);
	Read_double_vector_from_netcdf_file(fname, lon_var_name, f.lon_data);
	return (f);
	}
#endif

#ifndef DISABLE_NETCDF4
CField Read_general_2D_variable_netcdf_no_time_dimension(string fname, string var_name)
	{
	//error(AT,FU, "Treba se preveriti ali ta funcija res OK dela ko sem jo prepisal da deluje brez netcdf-C++ interface");

	int  status;
	int ncid;

	status = nc_open(fname.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't open input file: " + fname);

	vector <long> var_dim_sizes;
	vector <string> var_dim_names;
	vector <double> var_data;
	Read_whole_netcdf_double_variable_from_opened_netcdf_file( ncid, var_name,  var_data, var_dim_sizes, var_dim_names);

	if (var_dim_sizes.size() != 2)
		error(AT,FU, "var_dim_sizes.size() != 2");

	if (var_dim_sizes[0] *var_dim_sizes[1] != (long)var_data.size())
		error(AT,FU, "var_dim_sizes[0] *var_dim_sizes[1] != var_data.size()");

	CField f;
	f.dimy=var_dim_sizes[0];
	f.dimx=var_dim_sizes[1];
	f.data=var_data;

	status = nc_close(ncid);
	if (status != NC_NOERR)
		error(AT,FU, (string)"Can't close input file: " + fname);

	return (f);
	}
#endif


#ifndef DISABLE_NETCDF4
void convert_scaled_and_offset_netcdf_variable(string fname, string var_name, vector <double> &data)
	{
	// read attributes for scale and offset parameters
	double scale_factor=Read_attribute_as_double(fname, var_name, "scale_factor");
	double offset=Read_attribute_as_double(fname, var_name, "add_offset");
	double missing_value=Read_attribute_as_double(fname, var_name, "missing_value");

	// apply scale and offset factors
	for (unsigned long il=0; il < data.size(); il++)
		{
		if (data[il] == missing_value)
			data[il] = BAD_DATA_FLOAT;
		else
			data[il]=data[il]*scale_factor + offset;
		}
	}
#endif


// *******************************************************************************************
// *******************************************************************************************
// ECWMF data
// ***********************************

// ******************************************************
// Reads ECMWF netcdf file -- NEW
//
#ifndef DISABLE_NETCDF4
CField Read_ECMWF_2Dfield(string fname, int timenum, string var_name)
	{
	// read unscaled variable
	CField f=Read_general_2D_variable_netcdf(fname,  timenum, var_name);
	// read lon and lat
	Read_double_vector_from_netcdf_file(fname, "latitude", f.lat_data);
	Read_double_vector_from_netcdf_file(fname, "longitude", f.lon_data);

	convert_scaled_and_offset_netcdf_variable(fname, var_name, f.data);

	// read time
	vector <double> times;
	Read_double_vector_from_netcdf_file(fname, "time", times);
	CMy_time_format tt;
	tt.set_by_hours_since_1900_01_01(times[timenum]);
	f.unique_timestamp=tt.timestamp;

	return(f);
	}
#endif




