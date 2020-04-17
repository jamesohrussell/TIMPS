
class CField;
void WriteField(string ,  const CField &);
void WriteField_PNG(string ,  const CField &);
void draw_pngimage_of_binary_field(string , const CField & );


// ******************************************************
// ******************************************************
// Class CField. Stores one 2-D field with dimension atributes and raw data.
//
class CField {
	public:
	int dimx, dimy;
	vector <double> data;
	vector <double> lat_data;
	vector <double> lon_data;
	// used optionally
	long unique_timestamp;

	// constructor
	CField() {freememory();}
	~CField() {freememory();}

	// frees memory
	void freememory()
		{
		data.clear();
		lat_data.clear();
		lon_data.clear();
		unique_timestamp=0;
		dimx=0;
		dimy=0;
		}

	string output()
		{
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		ostringstream s1;
		s1.str("");

		long ix,iy;
		for (iy=0; iy < dimy; iy++)
			{
			s1 << iy << ": ";
			for (ix=0; ix < dimx; ix++)
				s1 << get(ix,iy) << " ";
			s1 << endl;
			}
		return(s1.str());
		}


	string output_summary()
		{
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		ostringstream s1;
		s1.str("");
		s1 << dimx << "x" << dimy << " data.size():"<< data.size() << " lat_data.size():"<< lat_data.size() << " lon_data.size():"<< lon_data.size() << " " << output_date_string_from_timestamp(unique_timestamp);
		return(s1.str());
		}


	bool check_if_initilized() const
		{
		if (lat_data.size()==0 || lon_data.size()==0 || data.size()==0) return (false);
		return(true);
		}

	bool check_if_data_is_initilized() const
		{
		if (data.size()==0 ) return (false);
		return(true);
		}

	bool compare_have_they_the_same_dimensions(const CField &f) const
		{
		if (dimx != f.dimx) return(false);
		if (dimy != f.dimy) return(false);
		return(true);
		}

	bool compare_have_they_the_same_dimensions_and_coordinates(const CField &f) const
		{
		if (lat_data != f.lat_data) return(false);
		if (lon_data != f.lon_data) return(false);
		if (dimx != f.dimx) return(false);
		if (dimy != f.dimy) return(false);
		return(true);
		}


	bool compare_are_they_the_same(const CField &f) const
		{
		if (unique_timestamp != f.unique_timestamp) return(false);
		if (lat_data != f.lat_data) return(false);
		if (lon_data != f.lon_data) return(false);
		if (data != f.data) return(false);
		return(true);
		}

	// overload operators == nad !=
	bool operator== (const CField& f) const
		{return(compare_are_they_the_same(f));}
	bool operator!= (const CField& f) const
		{return(!compare_are_they_the_same(f));}

	// area size of field
	long size() const
		{
		return((long)dimx*(long)dimy);
		}

	// gets value at x,y
	double get(int x, int y) const
		{
		return(data[(long)y*(long)dimx + (long)x]);
		}

	// gets value at x,y in periodic domain
	double get_periodic_x_y(int x, int y) const
		{
		int ix,iy;
		if (x<0)
			ix=dimx-(-x)%dimx;
		else
			ix=x%dimx;

		if (y<0)
			iy=dimy-(-y)%dimy;
		else
			iy=y%dimy;

		return(get(ix,iy));
		}

	// gets value at x,y in periodic domain
	double get_periodic_x(int x, int y) const
		{
		int ix;
		if (x<0)
			ix=dimx-(-x)%dimx;
		else
			ix=x%dimx;

		return(get(ix,y));
		}

	// gets value at x,y in periodic sphere - for y direction only crossing of border for one grid point is supported
	double get_periodic_sphere(int x, int y) const
		{
		int ix,iy;
		if (x<0)
			ix=dimx-(-x)%dimx;
		else
			ix=x%dimx;

		if (y < -1 || y > dimy)
			error(AT,FU, "y < -1 || y > dimy");

		if (y<0)
			{
			iy=-1-y;
			ix=(ix+dimx/2)%dimx;
			}
		else if (y > dimy-1)
			{
			iy=dimy-1-(y-dimy);
			ix=(ix+dimx/2)%dimx;
			}
		else
			iy=y;

		return(get(ix,iy));
		}

	// value at x,y in periodic sphere - for y direction only crossing of border for one grid point is supported
	void set_periodic_sphere(int x, int y, double value)
		{
		int ix;
		if (x<0)
			ix=dimx-(-x)%dimx;
		else
			ix=x%dimx;

		if (y < -1 || y > dimy)
			error(AT,FU, "y < -1 || y > dimy");

		// sets all points at lat=-90
		if (y==-1)
			for (long ix2=0; ix2 < dimx; ix2++)
				set(ix2,0,value);
		else if (y==dimy)
			for (long ix2=0; ix2 < dimx; ix2++)
				set(ix2,dimy-1,value);
		else
			set(ix,y,value);
		}

	bool check_coordianes(int x, int y)  const
		{
		if (x >= 0 && x < dimx && y >= 0 && y < dimy) return(true);
		else return(false);
		}


	// gets values at x,y and checks if x and y are valid
	double get_but_check_coordianes(int x, int y) const
		{
		if (check_coordianes(x,y))
			return(get(x,y));
		else return(BAD_DATA_FLOAT);
		}

	// gets values at x,y and checks if x and y are valid
	double get_but_check_coordinates_error(int x, int y) const
		{
		double out;

		out=get_but_check_coordianes( x,  y);

		if (out == BAD_DATA_FLOAT)
			{
			ostringstream s1;
			s1.str("");
			s1 << x << " " << y << " is invalid!";
			error(AT,FU, s1.str());
			}
		return(out);
		}

	// gets values at x,y and checks if x and y are valid
	double get_periodic_x_and_y(int x, int y) const
		{
		int x2;
		int y2;

		if (x >= 0) x2=x%dimx;
		else x2=dimx-((-x)%dimx);

		if (y >= 0) y2=y%dimy;
		else y2=dimy-((-y)%dimy);

		return(get_but_check_coordinates_error(x2,y2));
		}

	// gets values at x,y and checks if x and y are valid
	double get_but_outside_of_domain_zero_value(int x, int y) const
		{
		if (check_coordianes(x,y))
			return(get(x,y));
		else return(0);
		}


	CField& operator= (const CField &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			if (!rhs.check_if_initilized())
				error(AT,FU, "!rhs.check_if_initilized()");
			dimx=rhs.dimx;
			dimy=rhs.dimy;
			unique_timestamp=rhs.unique_timestamp;
			lon_data=rhs.lon_data;
			lat_data=rhs.lat_data;
			data=rhs.data;
			}
		return *this;    // Return ref for multiple assignment
		}

	CField& operator+= (const CField &rhs)
		{
		if (!rhs.check_if_initilized())
			error(AT,FU, "!rhs.check_if_initilized()");
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");
  		if (!compare_have_they_the_same_dimensions_and_coordinates(rhs))
			error(AT,FU, "!compare_have_they_the_same_dimensions_and_coordinates(f)");

		for (long il=0; il < rhs.size(); il++)
			{
  			if (rhs.data[il] == BAD_DATA_FLOAT || data[il] == BAD_DATA_FLOAT)
  				data[il] = BAD_DATA_FLOAT;
			else
				data[il]+=rhs.data[il];
			}
		return *this;    // Return ref for multiple assignment
		}


	CField& operator-= (const CField &rhs)
		{
		if (!rhs.check_if_initilized())
			error(AT,FU, "!rhs.check_if_initilized()");
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");
  		if (!compare_have_they_the_same_dimensions_and_coordinates(rhs))
			error(AT,FU, "!compare_have_they_the_same_dimensions_and_coordinates(f)");

		for (long il=0; il < rhs.size(); il++)
			{
   			if (rhs.data[il] == BAD_DATA_FLOAT || data[il] == BAD_DATA_FLOAT)
 				data[il] = BAD_DATA_FLOAT;
			else
				data[il]-=rhs.data[il];
			}
		return *this;    // Return ref for multiple assignment
		}


	const CField operator+ (const CField &rhs) const
		{
	    CField result = *this;     // Make a copy of myself.  Same as MyClass result(*this);
	    result += rhs;            // Use += to add other to the copy.
    	return result;              // All done!
		}

	const CField operator- (const CField &rhs) const
		{
	    CField result = *this;     // Make a copy of myself.  Same as MyClass result(*this);
	    result -= rhs;            // Use += to add other to the copy.
    	return result;              // All done!
		}

	void multiply_by_double(double x)
		{
		for (long il=0; il < size(); il++)
  			if (data[il] != BAD_DATA_FLOAT)
				data[il]*=x;
		}

	void add_double(double x)
		{
		for (long il=0; il < size(); il++)
  			if (data[il] != BAD_DATA_FLOAT)
				data[il]+=x;
		}


	double get_sum() const
		{
		double sum=0;
		for (long il=0; il < size(); il++)
			if (data[il] != BAD_DATA_FLOAT)
				sum+=data[il];
		return(sum);
		}

	double get_sum_row(long iy) const
		{
		double sum=0;
		for (long ix=0; ix < dimx; ix++)
			if (get(ix,iy) != BAD_DATA_FLOAT)
				sum+=get(ix,iy);
		return(sum);
		}

	double get_sum_column(long ix) const
		{
		double sum=0;
		for (long iy=0; iy < dimy; iy++)
			if (get(ix,iy) != BAD_DATA_FLOAT)
				sum+=get(ix,iy);
		return(sum);
		}

	long get_count_value_in_subdomain(double value, int xmin, int xmax, int ymin, int ymax) const
		{
		long count=0;
		for (long ix=xmin; ix <= xmax; ix++)
		for (long iy=ymin; iy <= ymax; iy++)
			if (get(ix,iy)==value)
				count++;
		return(count);
		}

	double get_maximum() const
		{
		double max=BAD_DATA_FLOAT;
		for (long il=0; il < size(); il++)
			if (data[il] > max)
				max=data[il];
		return(max);
		}

	double get_minimum() const
		{
		double min=1E100;
		for (long il=0; il < size(); il++)
			if (data[il] < min)
				min=data[il];
		return(min);
		}

	bool check_if_binary_field_no_missing_values() const
		{
		for (long il=0; il < size(); il++)
			if (data[il] == BAD_DATA_FLOAT)
				return(false);

		return(true);
		}

	bool check_if_field_is_binary_allow_missing_values() const
		{
		//cout << "bbb" << endl;
		for (long il=0; il < size(); il++)
			if (data[il] != 0 && data[il] != 1 && data[il] != BAD_DATA_FLOAT)
				{
				//cout << "aaa" << endl;
				cout << data[il] << endl;
				return(false);
				}

		return(true);
		}
 	bool check_if_field_has_no_missing_values() const
		{
		for (long il=0; il < size(); il++)
			if (data[il] == BAD_DATA_FLOAT)
				return(false);

		return(true);
		}

	void get_max_and_min_in_a_cricle(int x0, int y0, double r, double &max, double &min) const
		{
		min=1E20;
		max=1E-20;

		int dx,dy;
		for (dx=-(int)r; dx<=(int)r; dx++)
			for (dy=-(int)sqrt(r*r-(double)dx*(double)dx); dy<=(int)sqrt(r*r-(double)dx*(double)dx); dy++)
				if (check_coordianes(x0+dx, y0+dy))
					{
					if (get(x0+dx, y0+dy) > max) max=get(x0+dx, y0+dy);
					if (get(x0+dx, y0+dy) < min) min=get(x0+dx, y0+dy);
					}
		}


	string output_lat_lon() const
		{
		ostringstream s1;
		s1.str("");
		long il;
		for (il=0; il < dimx; il++)
			s1 << lon_data[il] << " ";
		s1 << endl;
		for (il=0; il < dimy; il++)
			s1 << lat_data[il] << " ";
		s1 << endl;
		return(s1.str());
		}

	string output_data() const
		{
		ostringstream s1;
		s1.str("");
		long il;
		for (il=0; il < size(); il++)
			s1 << data[il] << " ";
		s1 << endl;
		return(s1.str());
		}

	// reduces the resultion of the field by an integer factor using aggregation
	// is not verified that it works properly
	CField aggregate(int factor) const
		{
		double oldv;

		// dimenstions should be dividable ba factor
		if (dimx%factor!=0 || dimy%factor!=0)
			error(AT,FU, "dimx%factor!=0 || dimy%factor!=0");

		// there should be no missing data
		if (Check_missing_data(BAD_DATA_FLOAT) > 0)
			error(AT,FU, "Check_missing_data(BAD_DATA_FLOAT) > 0");

		CField f;
		f.allocatememory_and_set_dimx_dimy(dimx/factor, dimy/factor);
		f.allocatememory_lat_lon_array();
		f.Set_all_values(0);

		for (int ix=0; ix < dimx; ix++)
		for (int iy=0; iy < dimy; iy++)
			{
			oldv=f.get(ix/factor,iy/factor);
			f.set(ix/factor,iy/factor,oldv+get(ix,iy));
			}

		// normalize by numer of old points vs new points
		f.multiply_by_double(1.0/((double)factor*(double)factor));

		return(f);

		}


	// reduces the resultion of the field by an integer factor - it just skips points
	CField resample(int factor) const
		{

		int ix,iy;

		// dimenstions should be dividable by factor
		if (dimx%factor!=0 || dimy%factor!=0)
			error(AT,FU, "dimx%factor!=0 || dimy%factor!=0");

		CField f;
		f.allocatememory_and_set_dimx_dimy(dimx/factor, dimy/factor);
		f.allocatememory_lat_lon_array();

		for (ix=0; ix < dimx; ix+=factor)
		for (iy=0; iy < dimy; iy+=factor)
			{
			f.set(ix/factor,iy/factor,get(ix,iy));
			f.lon_data[ix/factor]=lon_data[ix];
			f.lat_data[iy/factor]=lat_data[iy];
			}
		return(f);
		}


	// reduces the resultion of the field by an integer factor - it just skips points
	CField resample2(int factor) const
		{

		int ix,iy;

		long newdimx=dimx/factor;
		long newdimy=dimy/factor;

		CField f;
		f.allocatememory_and_set_dimx_dimy(newdimx, newdimy);
		f.allocatememory_lat_lon_array();

		for (ix=0; ix < newdimx; ix++)
		for (iy=0; iy < newdimy; iy++)
			{
			f.set(ix,iy,get(ix*factor,iy*factor));
			f.lon_data[ix]=lon_data[ix*factor];
			f.lat_data[iy]=lat_data[iy*factor];
			}
		return(f);
		}

	// reduces the resultion of the field by an integer factor - it just skips points
	CField resize_no_anti_aliasing(double scale) const
		{
		int ix,iy;

		long newdimx=scale*(double)dimx;
		long newdimy=scale*(double)dimy;

		CField f;
		f.allocatememory_and_set_dimx_dimy(newdimx, newdimy);
		f.allocatememory_lat_lon_array();

		for (ix=0; ix < newdimx; ix++)
		for (iy=0; iy < newdimy; iy++)
			f.set(ix,iy,get(round((double)ix/scale),round((double)iy/scale)));
		return(f);
		}

	// reduces the resultion of the field by an integer factor - it just skips points
	void invert_grayscale_field_between_0_and_1()
		{
		ERRORIF(!check_if_binary_field_no_missing_values());
		for (long il=0; il < size(); il++)
			data[il]=1-data[il];

		}


	// reduces the resultion of the field
	CField resize_with_anti_aliasing(double scale) const
		{
		int ix,iy;

		long newdimx=scale*(double)dimx;
		long newdimy=scale*(double)dimy;

		CField f;
		f.allocatememory_and_set_dimx_dimy(newdimx, newdimy);
		f.allocatememory_lat_lon_array();

		// calculate fractions
		ERRORIF(!check_if_binary_field_no_missing_values());

		long n=round(1.0/scale);
		CField ffrac=calculate_fractions_from_binary_field(n);

		#pragma omp parallel for
		for (ix=0; ix < newdimx; ix++)
		for (iy=0; iy < newdimy; iy++)
			f.set(ix,iy,ffrac.get(round((double)ix/scale),round((double)iy/scale)));
		return(f);
		}

	// sets value at x,y
	void set(int x, int y, double value)
		{
		data[(long)y*(long)dimx + (long)x]=value;
		}

	// sets value at x,y but makes a check first that x,y are inside domain
	void set_but_check_first(int x, int y, double value)
		{
		if (x >= 0 && x < dimx && y >= 0 && y < dimy)
			set(x,y,value);
		}

	// grows domain on all sides - the new area have the value of border_value
	CField grow_domain_by_n_pixels(long n, double border_value) const
		{
		CField f;
		f.allocatememory_and_set_dimx_dimy(dimx + 2*n, dimy + 2*n);
		f.Set_all_values(border_value);
		f.allocatememory_lat_lon_array();


		for (int ix=0; ix < dimx; ix++)
		for (int iy=0; iy < dimy; iy++)
			f.set(ix+n,iy+n,get(ix,iy));
		return(f);
		}

	// grows domain on all sides - the new area have the value of border_value
	CField grow_domain_by_n_pixels_only_vertical(long ntop, long nbottom, double border_value) const
		{
		CField f;
		f.allocatememory_and_set_dimx_dimy(dimx, dimy + ntop + nbottom);
		f.Set_all_values(border_value);
		f.allocatememory_lat_lon_array();


		for (int ix=0; ix < dimx; ix++)
		for (int iy=0; iy < dimy; iy++)
			f.set(ix,iy+nbottom,get(ix,iy));
		return(f);
		}

	// gets lat/lon location from t x,y
	void XY2LatLon(int x, int y, double& lon, double& lat) const
		{
		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->XY2LatLon ERROR: lat_data.size()==0 || lon_data.size()==0 \n";exit(1);}

		if (x < 0 || y < 0 || x > dimx-1 || y > dimy - 1)
			{cout << "ERROR: CField->XY2LatLon ERROR: x < 0 || y < 0 || x > dimx-1 || y > dimy - 1 \n";exit(1);}

		lon=lon_data[x];
		lat=lat_data[y];
		}

	// gets x,y location from lat/lon
	void LatLon2XY(double lon, double lat, int &x, int &y ) const
		{
		int ix;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->LatLon2XY ERROR: lat_data.size()==0 || lon_data.size()==0 \n";exit(1);}

		x=BAD_DATA_FLOAT;
		y=BAD_DATA_FLOAT;
		for (ix = 0; ix < dimx-1; ix++)
			if (lon_data[ix] <= lon && lon_data[ix+1] >= lon)
				{
				if (fabs(lon_data[ix] - lon) > fabs(lon_data[ix+1] - lon))
					x=ix+1;
				else
					x=ix;
				}

		for (ix = 0; ix < dimy-1; ix++)
			if (lat_data[ix] <= lat && lat_data[ix+1] >= lat)
				{
				if (fabs(lat_data[ix] - lat) > fabs(lat_data[ix+1] - lat))
					y=ix+1;
				else
					y=ix;
				}

		if (x== BAD_DATA_FLOAT || y== BAD_DATA_FLOAT)
			{cout << "ERROR: CField->LatLon2XY ERROR: x== BAD_DATA_FLOAT || y== BAD_DATA_FLOAT \n"; cout << lon << " " << lat << " " << lon_data[0] << endl; exit(1);}
		}


	// gets x,y location from lat/lon
	void LatLon2XY_new(double lon, double lat, int &x, int &y ) const
		{
		int ix,iy;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->LatLon2XY ERROR: lat_data.size()==0 || lon_data.size()==0 \n";exit(1);}

		x=BAD_DATA_FLOAT;
		y=BAD_DATA_FLOAT;

		double min;
		min=1000000;
		// find colosest distance to the lat
		for (ix = 0; ix < dimx; ix++)
			if (fabs(lon_data[ix]-lon) < min)
				{
				x=ix;
				min=fabs(lon_data[ix]-lon);
				}
		// if x= border value - check size of min
		if (x == 0 && min > fabs(lon_data[0]-lon_data[1])) x=BAD_DATA_FLOAT;
		if (x == dimx-1 && min > fabs(lon_data[dimx-1]-lon_data[dimx-2])) x=BAD_DATA_FLOAT;


		min=1000000;
		for (iy = 0; iy < dimy; iy++)
			if (fabs(lat_data[iy]-lat) < min)
				{
				y=iy;
				min=fabs(lat_data[iy]-lat);
				}

		// if y= border value - check size of min
		if (y == 0 && min > fabs(lat_data[0]-lat_data[1])) y=BAD_DATA_FLOAT;
		if (y == dimy-1 && min > fabs(lat_data[dimy-1]-lat_data[dimy-2])) y=BAD_DATA_FLOAT;

		if (x== BAD_DATA_FLOAT || y== BAD_DATA_FLOAT)
			{
			cout << lon << " " << lat << " " << lon_data[0] << " " << lon_data[dimx-1] << " " << lat_data[0] << " " << lat_data[dimy-1] << endl;
			error(AT,FU, "x== BAD_DATA_FLOAT || y== BAD_DATA_FLOAT");
			}
		}


	// gets x,y location from lat/lon
	double get_by_lon_and_lat(double lon, double lat) const
		{
		int x,y;
		//cout << lon << "\t" << lat << endl;
		LatLon2XY_new(lon, lat, x, y );
		return(get_but_check_coordinates_error(x,y));
		}

	// gets x,y location from lat/lon
	double get_by_lon_and_lat_v2(double lon, double lat) const
		{
		int x,y;
		LatLon2XY_new(lon, lat, x, y );
		//cout << x << "\t" << y << "\t" << get(x,y) << endl;
		return(get_but_check_coordianes(x,y));
		}

	// set x,y location from lat/lon
	void set_by_lon_and_lat(double lon, double lat, double value)
		{
		int x,y;
		LatLon2XY_new(lon, lat, x, y );
		if (!check_coordianes(x,y))
			error(AT,FU, "!check_coordianes(x,y)");
		set(x,y,value);
		}

	// calculates x,y location from the array offset
	void DataLoc2XY(long ix, int &x, int &y) const
		{
		x=ix%dimx;
		y=ix/dimx;
		}

	// calculates x location from the array offset
	int DataLoc2X(long ix) const
		{
		int x,y;
		DataLoc2XY(ix,x,y);
		return(x);
		}

	// calculates y location from the array offset
	int DataLoc2Y(long ix) const
		{
		int x,y;
		DataLoc2XY(ix,x,y);
		return(y);
		}


	// calculates array offset from x,y
	long XY2DataLoc(int x, int y) const
		{
		return((long)y*(long)dimx + (long)x);
		}

	void Point_XY2LatLon(Point &p) const
		{
		XY2LatLon((int)floor(p.x), (int)floor(p.y), p.lon, p.lat);
		}

	void Point_LatLon2XY(Point &p) const
		{
		int x,y;
		LatLon2XY(p.lon, p.lat, x, y);
		p.x=(double)x;
		p.y=(double)y;
		}

	// calculates x,y location from the array offset
	void Point_DataLoc2XY(long ix, Point &p) const
		{
		p.x=ix%dimx;
		p.y=ix/dimx;
		}

	double area_of_grid_square_in_km2(int ix, int iy) const
		{
		int dx,dy;
		double area_in_km2;
		double Rz=6371;
		double lon,lat,lon2,lat2;

		// get longitude and latitude of point
		XY2LatLon(ix,iy,lon,lat);

		// get longitude and latitude of nearby point
		if (ix != dimx - 1) dx=1;
		else dx=-1;
		if (iy != dimy - 1) dy=1;
		else dy=-1;
		XY2LatLon(ix+dx,iy+dy,lon2,lat2);

		if ( fabs(lat-lat2) > 5 || fabs(lon-lon2) > 5)
			error(AT,FU, "fabs(lat-lat2) > 5 || fabs(lon-lon2) > 5");

		area_in_km2=fabs(deg2rad(lat-lat2)*Rz*deg2rad(lon-lon2)*Rz*cos(deg2rad(lat)));

		return(area_in_km2);
		//area_cos+=fabs(cos(lat*M_PI/180.0)*(lat-lat2)*(lon-lon2));
		}


	void Set_all_values_to_grid_square_area_in_km2()
		{
		int ix,iy;
		for (ix = 0; ix < dimx; ix++)
			for (iy=0; iy < dimy ; iy++)
				{
				set(ix,iy,area_of_grid_square_in_km2(ix,iy));
				//cout << area_of_grid_square_in_km2(ix,iy) << endl;
				}
		}

	double calculate_sqr_distance_between_two_grid_points(int x1, int y1, int x2, int y2) const
		{
		return((double)(x1-x2)*(double)(x1-x2) + (double)(y1-y2)*(double)(y1-y2));
		}

	double calculate_sqr_distance_between_two_grid_points_specified_by_DataLoc(long il1, long il2) const
		{
		return(calculate_sqr_distance_between_two_grid_points(DataLoc2X(il1), DataLoc2Y(il1), DataLoc2X(il2), DataLoc2Y(il2)));
		}


	double calculate_distance_on_earth_in_km_between_two_grid_points(int x1, int y1, int x2, int y2) const
		{
		double lat1,lon1,lat2,lon2;
		XY2LatLon(x1, y1, lon1, lat1);
		XY2LatLon(x2, y2, lon2, lat2);
		return(distance_on_earth_in_km_specified_by_deg(lat1, lon1, lat2, lon2));
		}


	// allocates memory to data pointer using dimx, dimy
	void allocatememory()
		{
		data.assign((long)dimx*(long)dimy,0);
		}

	// allocates memory to data pointer using dimx, dimy
	void allocatememory_and_set_dimx_dimy(int dimx_, int dimy_)
		{
		dimx=dimx_;
		dimy=dimy_;
		allocatememory();
		}


	// allocates memory for lat and lon array using dimx, dimy
	void allocatememory_lat_lon_array()
		{
		lat_data.assign((long)dimy,0);
		lon_data.assign((long)dimx,0);
		}

	void set_lat_lon_to_zero()
		{
		long il;
		for (il=0; il < dimx; il++) lon_data[il]=0;
		for (il=0; il < dimy; il++) lat_data[il]=0;
		}

	void set_epsilon_zero(double epsilon)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			if (fabs(data[ixl]) < epsilon)
				data[ixl]=0;
		}

	void set_epsilon_zero_for_negative_values_only(double epsilon)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			if (fabs(data[ixl]) < epsilon && data[ixl] < 0)
				data[ixl]=0;
		}

	// sets all values of a field to value
	void Set_all_values(double value)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			data[ixl]=value;
		}

	// sets all values of a field to value
	void Set_all_values_in_value_range(double value, double lower_limit, double upper_limit)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			if (data[ixl] >= lower_limit && data[ixl] <= upper_limit)
				data[ixl]=value;
		}

	void set_to_one_above_threshold_else_zero_replace_also_missing_data(double thrs)
		{
		for (long il=0; il<size(); il++)
			{
			if (data[il] > thrs)
				data[il]=1;
			else
				data[il]=0;
			}
		}


	bool are_there_different_values_around(int x, int y, int radius, double value)
		{
		for (int ix=x-radius; ix <= x+radius; ix++)
		for (int iy=y-radius; iy <= y+radius; iy++)
			if (check_coordianes(ix,iy))
				if (get(ix,iy) != value)
					return(true);
		return(false);
		}

	bool are_there_different_values_around_or_domain_border(int x, int y, int radius, double value)
		{
		for (int ix=x-radius; ix <= x+radius; ix++)
		for (int iy=y-radius; iy <= y+radius; iy++)
			if (get_but_check_coordianes(ix,iy) != value)
				return(true);
		return(false);
		}



	// replace all values oldV with newV
	void Replace_values(double oldV, double newV)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			if (data[ixl]==oldV) data[ixl]=newV;
		}

	// replace all values oldV with newV
	void Replace_values_in_certain_range(double oldVmin, double oldVmax, double newV)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			if (data[ixl]>=oldVmin && data[ixl] <=oldVmax ) data[ixl]=newV;
		}


	// replace all values oldV with newV
	void Replace_values_only_to_certain_line(double oldV, double newV, int maxy)
		{
		int ix,iy;

		for (iy=0;iy<=maxy;iy++)
			for (ix=0;ix<dimx;ix++)
				if (get(ix,iy) == oldV)
					set(ix,iy,newV);
		}

	// replace all values oldV with newV
	void Replace_values_only_in_certain_area(double oldV, double newV, int xmin, int xmax, int ymin, int ymax)
		{
		if ( !check_coordianes(xmin, ymin) || !check_coordianes(xmax, ymax))
			{
			ostringstream s1;
			s1.str("");
			s1 << xmin << " " << ymin << " " << xmax << " " << ymax << " is invalid!";
			error(AT,FU, s1.str());
			}

		int ix,iy;

		for (iy=ymin;iy<=ymax;iy++)
			for (ix=xmin;ix<=xmax;ix++)
				if (get(ix,iy) == oldV)
					set(ix,iy,newV);
		}

	void Replace_list_of_values_with_a_single_value ( vector <double> old_values_list, double new_value)
		{
		sort(old_values_list.begin(), old_values_list.end());
		for (long il=0; il < size(); il++)
			if (binary_search (old_values_list.begin(), old_values_list.end(), data[il]))
				data[il]=new_value;
		}

	void Replace_list_of_values_with_a_list_of_new_values (vector <double> old_values_list, vector <double> new_values_list)
		{
		if (old_values_list.size() != new_values_list.size())
			error(AT,FU, "old_values_list.size() != new_values_list.size()");


		vector <double> sorted_old_values_list=old_values_list;
		sort(sorted_old_values_list.begin(), sorted_old_values_list.end());
		for (long il=0; il < size(); il++)
			// first check that the valuje is in the old_values_list - use binary_search on a sorted array which is very fast
			if (binary_search (sorted_old_values_list.begin(), sorted_old_values_list.end(), data[il]))
				// if found then go over list and select the correspoding new value
				for (long ix=0; ix < (long)old_values_list.size(); ix++)
					if (old_values_list[ix]==data[il])
						data[il]=new_values_list[ix];
		}

	// Identify_representative_axis_of_objects - calculated orientation and symetry of region vith values 1 and 0 elsewhere
	void calculate_object_orientation_and_symetry_using_image_algebra_book_algorithm(int xmin, int xmax, int ymin, int ymax, double &orientation, double &symetrydouble) //const
		{
		//*this->allocatememory_lat_lon_array();
		int ix,iy;

		/*Set_all_values(0);
		for (ix=100; ix < 105; ix++)
		for (iy=100; iy < 150; iy++)
			set(ix,iy,1);
		Set_all_values(0);
		for (ix=0; ix < dimx; ix++)
		for (iy=0; iy < dimy; iy++)
			if (
				pow((double)(ix-150)*cos(M_PI/4)+(double)(iy-150)*sin(M_PI/4),2)/1000 +
				pow((double)(ix-150)*sin(M_PI/4)-(double)(iy-150)*cos(M_PI/4),2)/300
				< 1)
				set(ix,iy,1);
	*/

		allocatememory_lat_lon_array();
		WriteField("aaa1.nc",*this);
		WriteField_PNG("aaa1.png",*this);

		double xc,yc,m20,m11,m02,objsize=0;

		// center
		for (ix=0; ix < dimx; ix++)
		for (iy=0; iy < dimy; iy++)
			{
			xc+=get(ix,iy)*(double)ix;
			yc+=get(ix,iy)*(double)iy;
			objsize+=get(ix,iy);
			}
		xc/=objsize;
		yc/=objsize;

		// moments
		for (ix=0; ix < dimx; ix++)
		for (iy=0; iy < dimy; iy++)
			{
			m20+=get(ix,iy)*((double)ix-xc)*((double)ix-xc);
			m11+=get(ix,iy)*((double)ix-xc)*((double)iy-yc);
			m02+=get(ix,iy)*((double)iy-yc)*((double)iy-yc);
			}

		// orientation, symetry
		double a,fi1,fi2,mmm1,mmm2,symetryparameter,fimax,fimin;
		if (m11 != 0)
			{
			a=(m20-m02)/m11;
			fi1=atan(0.5* (-a - sqrt(4 + a*a)));
			fi2=atan(0.5* (-a + sqrt(4 + a*a)));
			}
		else
			{
			fi1=0;
			fi2=M_PI/2;
			}

		mmm1=m20*sin(fi1)*sin(fi1)-2*m11*sin(fi1)*cos(fi1)+m02*cos(fi1)*cos(fi1);
		mmm2=m20*sin(fi2)*sin(fi2)-2*m11*sin(fi2)*cos(fi2)+m02*cos(fi2)*cos(fi2);


		if (mmm1 > mmm2)
			{
			fimax=fi1;
			fimin=fi2;
			symetryparameter=mmm2/mmm1;
			}
		else
			{
			fimax=fi2;
			fimin=fi1;
			symetryparameter=mmm1/mmm2;
			}

		cout << objsize << endl;
		cout << xc << endl;
		cout << yc << endl;
		cout << m20 << endl;
		cout << m11 << endl;
		cout << m02 << endl;
		cout << fi1 *180/M_PI << endl;
		cout << fi2 *180/M_PI << endl;
		cout << mmm1 << endl;
		cout << mmm2 << endl;
		cout << fimax *180/M_PI << endl;
		cout << fimin *180/M_PI << endl;
		cout << symetryparameter << endl;
		exit(1);

		}


	// draws a filled circle with value. Circle center is x0,y0 and radius
	void draw_filled_circle(int x0, int y0, double radius, double value)
		{
		int dx,dy;
		for (dx=-(int)radius; dx<=(int)radius; dx++)
			for (dy=-(int)sqrt(radius*radius-(double)dx*(double)dx); dy<=(int)sqrt(radius*radius-(double)dx*(double)dx); dy++)
				set_but_check_first(x0+dx, y0+dy, value);
		}

	// draws a filled circle with value. Circle center is x0,y0 and radius
	void draw_gauss(int x0, int y0, double sigma, double amplitude, double bias)
		{
		for (long ix=0; ix < dimx; ix++)
		for (long iy=0; iy < dimy; iy++)
			set(ix,iy,bias + amplitude*exp(-(
					(double)(ix-x0)*(double)(ix-x0)/(2*sigma*sigma) +
					(double)(iy-y0)*(double)(iy-y0)/(2*sigma*sigma))
					));
		}


	// draws a gaussian function inside circle
	void draw_gauss_inside_circle_add_value(int x0, int y0, double sigma, double amplitude, double bias, double radius)
		{
		#pragma omp parallel for
		for (long dx=-(int)radius; dx<=(int)radius; dx++)
			for (long dy=-(int)sqrt(radius*radius-(double)dx*(double)dx); dy<=(int)sqrt(radius*radius-(double)dx*(double)dx); dy++)
				{
				long ix=x0+dx;
				long iy=y0+dy;
				if (check_coordianes(ix,iy))
					{
					double old_value=get(ix,iy);
					if (old_value != BAD_DATA_FLOAT)
						{
						set(ix,iy,old_value+bias + amplitude*exp(-(
						(double)(ix-x0)*(double)(ix-x0)/(2*sigma*sigma) +
						(double)(iy-y0)*(double)(iy-y0)/(2*sigma*sigma))
						));
						}
					}
				}
		}


// draws a filled roteted ellipse. e1 is the ration of rmin/rmax. e1 should be less than 1 so we can optimize the drawing - by drawing just inside the circle which encompasses the ellipse
void draw_filled_rotated_ellipse(int x1, int y1, int rmax, double e1, double fi1, double value)
		{
		if (e1 > 1 || e1 <= 0)
			error(AT,FU, "e1 > 1 || e1 <= 0)");

		int dx,dy;
		for (dx=-(int)rmax; dx<=(int)rmax; dx++)
		for (dy=-(int)sqrt(rmax*rmax-(double)dx*(double)dx); dy<=(int)sqrt(rmax*rmax-(double)dx*(double)dx); dy++)
			if (is_location_inside_rotated_elipse(x1+dx,y1+dy, x1, y1, rmax, e1, fi1 ))
				set_but_check_first(x1+dx, y1+dy, value);
		}


	// draws a filled rectangle
	void draw_filled_rectangle(int x0, int y0, int x1, int y1, double value)
		{
		if (x0 > x1 || y0 > y1)
			error(AT,FU, "x0 > x1 || y0 > y1");

		for (long ix=x0; ix <= x1; ix++)
		for (long iy=y0; iy <= y1; iy++)
			if (check_coordianes(ix,iy))
				set(ix,iy,value);
		}

	// checks is all points in a CField are missing values. Returns 0 if vhole field is empty or "count" which represents how many points are non-missing
	long Test_if_whole_field_is_missing_values() const
		{
		long ixl;
		long count=0;
		for (ixl=0;ixl<size();ixl++)
			if (data[ixl]!=BAD_DATA_FLOAT) count++;
		return(count);
		}

	// checks is all non-missing points have a specific value - good for testing if all points have value 0. Returns the number of points not havening a value of "value".
	long Test_if_whole_field_is_a_value(double value) const
		{
		long ixl;
		long count=0;
		for (ixl=0;ixl<size();ixl++)
			if (data[ixl] != BAD_DATA_FLOAT && data[ixl] != value) count++;
		return(count);
		}

	// draw a line between two points with value
	// the algorithm from http://cs.unc.edu/~mcmillan/comp136/Lecture6/Lines.html
    void lineImproved(int x0, int y0, int x1, int y1, double value, int line_width)
    {
        int dx = x1 - x0;
        int dy = y1 - y0;
        int width;

		// check if at least one point is inside the domain
		int inside_domain=0;
		if (x0 >= 0 && x0 < dimx && y0 >= 0 && y0 < dimy) inside_domain=1;
		if (x1 >= 0 && x1 < dimx && y1 >= 0 && y1 < dimy) inside_domain=1;


		// draw only if at least one point is inside the domain
		if 	(inside_domain==1)
			{
			set_but_check_first(x0, y0,value);
			if (abs(dx) > abs(dy))          // slope < 1
				{
				double m = (double) dy / (double) dx;      // compute slope
				double b = (double)y0 - m*(double)x0;
				dx = (dx < 0) ? -1 : 1;
				while (x0 != x1)
					{
					x0 += dx;
					set_but_check_first(x0, (int)round(m*(double)x0 + b), value);
					for (width=2; width <= line_width; width++)
						set_but_check_first(x0, (int)round(m*(double)x0 + b) + width - 1, value);
					}
				}
			else
			if (dy != 0)
				{                              // slope >= 1
				double m = (double) dx / (double) dy;      // compute slope
				double b = (double)x0 - m*(double)y0;
				dy = (dy < 0) ? -1 : 1;
				while (y0 != y1)
					{
					y0 += dy;
					set_but_check_first((int)round(m*(double)y0 + b), y0,value);
					for (width=2; width <= line_width; width++)
						set_but_check_first((int)round(m*(double)y0 + b) + width - 1, y0,value);
					}
				}
			}
	}

	double get_max_value() const
		{
		double max=1E-100;
		for (long il=0; il < size(); il++)
			if (data[il] > max)
				max=data[il];

		if (max==1E-100)
			error(AT,FU, "max==1E-100");

		return(max);
		}

	double calcualte_RMSE(const CField &f1) const
		{
		double rmse=0;

		if (!compare_have_they_the_same_dimensions(f1))
			error(AT,FU, "compare_have_they_the_same_dimensions(f1)");

		for (long il=0; il < size(); il++)
			rmse+=(f1.data[il]-data[il])*(f1.data[il]-data[il]);

		return(sqrt(rmse/(double)size()));
		}

	double calcualte_RMSE_in_domain(const CField &f1,const CField &fdomain) const
		{
		double rmse=0;

		if (!compare_have_they_the_same_dimensions(f1))
			error(AT,FU, "compare_have_they_the_same_dimensions(f1)");

		double counter=0;

		for (long il=0; il < size(); il++)
			if (fdomain.data[il]==1)
				{
				rmse+=(f1.data[il]-data[il])*(f1.data[il]-data[il]);
				counter++;
				}

		return(sqrt(rmse/counter));
		}

	double calcualte_correlation_in_domain(const CField &f1,const CField &fdomain) const
		{
		if (!compare_have_they_the_same_dimensions(f1))
			error(AT,FU, "compare_have_they_the_same_dimensions(f1)");

		if (!compare_have_they_the_same_dimensions(fdomain))
			error(AT,FU, "compare_have_they_the_same_dimensions(fdomain)");

		if (!check_if_initilized() || !f1.check_if_initilized() || !fdomain.check_if_initilized())
			error(AT,FU, "!check_if_initilized() || !f1.check_if_initilized() || !fdomain.check_if_initilized()");

		vector <double> x;
		vector <double> y;
		for (long il=0; il < size(); il++)
			if (fdomain.data[il]==1 && data[il] != BAD_DATA_FLOAT && f1.data[il] != BAD_DATA_FLOAT)
				{
				x.push_back(data[il]);
				y.push_back(f1.data[il]);
				}

		//if (x.size()==0 || y.size()==0 )
		//	error(AT,FU, "x.size()==0 || y.size()==0");

		if (x.size() < 2 )
			return(0);

		return(pearsons_correlation_from_vector(x,y));
		}



	// Grow areas between lower and upper limit by pixels
	void grow_area(double lower_limit, double upper_limit, int pixels)
		{
		int ix, iy, ip;
		CField f;
		long il;

		for (ip=0;ip<pixels;ip++)
			{
			f=create_a_copy_with_data();
			for (iy=0;iy<dimy;iy++)
			for (ix=0;ix<dimx;ix++)
			if (get(ix,iy) != BAD_DATA_FLOAT)
			if (get(ix,iy) < lower_limit || get(ix,iy) > upper_limit )
				{
				if (check_coordianes(ix-1,iy)) if (get(ix-1,iy) >= lower_limit && get(ix-1,iy) <= upper_limit) f.set(ix,iy,get(ix-1,iy));
				if (check_coordianes(ix+1,iy)) if (get(ix+1,iy) >= lower_limit && get(ix+1,iy) <= upper_limit) f.set(ix,iy,get(ix+1,iy));
				if (check_coordianes(ix,iy-1)) if (get(ix,iy-1) >= lower_limit && get(ix,iy-1) <= upper_limit) f.set(ix,iy,get(ix,iy-1));
				if (check_coordianes(ix,iy+1)) if (get(ix,iy+1) >= lower_limit && get(ix,iy+1) <= upper_limit) f.set(ix,iy,get(ix,iy+1));
				}

			for (il=0; il < f.size(); il++)
				data[il]=f.data[il];
			f.freememory();
			}


		}

	CField get_border_of_area(double min, double max, int border_width)
		{
		int ip;
		long il;
		CField f,f2;
		f2=create_a_copy_with_data();
		f2.Set_all_values(0);
		f=create_a_copy_with_data();
		for (ip=0;ip<border_width;ip++)
			{
			// enlarge areas with values above max or below min
			f.grow_area(max+0.00001,1E99,1);
			f.grow_area(-1E99, min-0.00001, 1);
			}
		// define border as area of difference between original and f
		for (il=0;il<f.size();il++)
			if (f.data[il] != data[il]) f2.data[il]=1;

		f.freememory();
		return(f2);
		}


	// defines border for all same value areas between min and max
	CField get_border_of_area_variant_2(double min, double max, int border_width)
		{
		CField f=*this;
		f.Set_all_values(0);
		for (int ix=0; ix < dimx;ix++)
		for (int iy=0; iy < dimy;iy++)
			if (get(ix,iy)>= min && get(ix,iy)<= max)
				if( are_there_different_values_around_or_domain_border(ix,iy,border_width, get(ix,iy)))
					f.set(ix,iy,get(ix,iy));

		return(f);
		}

	CField get_contors_of_scalar_field(vector <double> levels) const
		{
		double l,t,c;
		CField f;
		f=*this;
		f.Set_all_values(0);
		for (int ix=1; ix < dimx-1;ix++)
		for (int iy=1; iy < dimy-1;iy++)
			{
			l=get(ix-1,iy);
			t=get(ix,iy-1);
			c=get(ix,iy);

			for (long il=0; il < (long)levels.size(); il++)
				if ( (l <= levels[il] && c > levels[il]) ||
				     (l >= levels[il] && c < levels[il]) ||
				     (t <= levels[il] && c > levels[il]) ||
				     (t >= levels[il] && c < levels[il])  )
					f.set(ix,iy,1);
			}

		return(f);
		}

	CField get_contours_of_scalar_field_equaly_spaced_levels(double dlev) const
		{
		double l,t,c;
		CField f,f1;
		f=*this;
		f1=*this;
		f.Set_all_values(0);

		// floor data/dlev for the whole field
		for (long il=0; il < f1.size(); il++)
			f1.data[il]=floor(data[il]/dlev);

		for (int ix=1; ix < dimx-1;ix++)
		for (int iy=1; iy < dimy-1;iy++)
			{
			l=f1.get(ix-1,iy);
			t=f1.get(ix,iy-1);
			c=f1.get(ix,iy);

			if (l != c || t != c)
				f.set(ix,iy,1);
			}
		return(f);
		}

CField calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation() const
	{
	CField fsum=*this;

	// replace missing values with zeros
	fsum.Replace_values(BAD_DATA_FLOAT,0);

	// calculate the summed fields
	// sum by rows first - we can skip the left row since it allready contains the valeus from f1
	#pragma omp parallel for
	for (int iy=0; iy < dimy; iy++)
	for (int ix=1; ix < dimx; ix++)
		fsum.set(ix,iy,fsum.get(ix-1,iy) + fsum.get(ix,iy));

	// sum by columns next - we can skip the bottom row since it allready contains the sumemd valeus from f1sum which dont change
	#pragma omp parallel for
	for (int ix=0; ix < dimx; ix++)
	for (int iy=1; iy < dimy; iy++)
		fsum.set(ix,iy,fsum.get(ix,iy) + fsum.get(ix,iy-1));

	return(fsum);
	}

double calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(const CField &f2sum, int n) const
	{
	if (!compare_have_they_the_same_dimensions(f2sum))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f2sum)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;

	// calculate the fractions and FBS and FBSref
	double FBS=0;
	double FBSref=0;
	// the parallelized loop actualy runs slower on two pocessor than non-parralelized - no idea why
	//#pragma omp parallel for default(shared) reduction(+:FBS,FBSref)
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		double frac1;
		double frac2;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r - 1;
		ix + r >= dimx   ? x2 = dimx - 1 : x2 = ix + r;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r - 1;
		iy + r >= dimy   ? y2 = dimy - 1 : y2 = iy + r;

		frac1=get(x2,y2);
		frac2=f2sum.get(x2,y2);

		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			frac2+=-f2sum.get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			frac2+=-f2sum.get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			frac2+=f2sum.get(x1,y1);
			}
		//double FBS,FBSref;
		FBS=FBS+((frac1-frac2)*(frac1-frac2));
		FBSref=FBSref+frac1*frac1+frac2*frac2;
		}

	if (FBSref != 0)
		return(1.0-FBS/FBSref);
	else
		return(0);
	}

// Similar to Standard FSS but the neigberhood is cropped at the borders - the fractions are not divided by n^2 but with the number of neigberhood points included in the domain
double calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation_alternative_border_treatment(const CField &f2sum, int n) const
	{
	if (!compare_have_they_the_same_dimensions(f2sum))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f2sum)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;

	// calculate the fractions and FBS and FBSref
	double FBS=0;
	double FBSref=0;
	// the parallelized loop actualy runs slower on two pocessor than non-parralelized - no idea why
	//#pragma omp parallel for default(shared) reduction(+:FBS,FBSref)
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		double frac1;
		double frac2;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r - 1;
		ix + r >= dimx   ? x2 = dimx - 1 : x2 = ix + r;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r - 1;
		iy + r >= dimy   ? y2 = dimy - 1 : y2 = iy + r;

		frac1=get(x2,y2);
		frac2=f2sum.get(x2,y2);

		double dx,dy;
		x1 + 1 < 0   ? dx = x2 + 1: dx = x2 - x1;
		y1 + 1 < 0   ? dy = y2 + 1: dy = y2 - y1;

		//cout << ix << "\t" << iy << "\t" << dx << "\t" << dy << "\tx1:" << x1 << " x2:" << x2 << " y1:" << y1 << " y2:" << y2 << endl;


		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			frac2+=-f2sum.get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			frac2+=-f2sum.get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			frac2+=f2sum.get(x1,y1);
			}

		frac1*=1.0/(dx*dy);
		frac2*=1.0/(dx*dy);

		//if (!compare_double_with_epsilon_tolerance(frac1,1,0.00001) || !compare_double_with_epsilon_tolerance(frac2,1,0.00001))
		//	{
		//	cout << ix << "\t" << iy << "\t" << frac1 << "\t" << frac2 << endl;
		//	error(AT,FU, "frac1 != 1 && frac2 != 1");
		//	}

		//double FBS,FBSref;
		FBS=FBS+((frac1-frac2)*(frac1-frac2));
		FBSref=FBSref+frac1*frac1+frac2*frac2;
		}

	if (FBSref != 0)
		return(1.0-FBS/FBSref);
	else
		return(0);
	}


double calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation_in_subdomain(const CField &f2sum, int n, const CField &fdomain) const
	{
	if (!compare_have_they_the_same_dimensions(f2sum))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f2sum)");

	if (!compare_have_they_the_same_dimensions(fdomain))
		error(AT,FU, "!compare_have_they_the_same_dimensions(fdomain)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;

	// calculate the fractions and FBS and FBSref
	double FBS=0;
	double FBSref=0;
	// the parallelized loop actualy runs slower on two pocessor than non-parralelized - no idea why
	//#pragma omp parallel for default(shared) reduction(+:FBS,FBSref)
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
	if (fdomain.get(ix,iy) == 1)
		{
		double frac1;
		double frac2;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r - 1;
		ix + r >= dimx   ? x2 = dimx - 1 : x2 = ix + r;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r - 1;
		iy + r >= dimy   ? y2 = dimy - 1 : y2 = iy + r;

		frac1=get(x2,y2);
		frac2=f2sum.get(x2,y2);

		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			frac2+=-f2sum.get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			frac2+=-f2sum.get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			frac2+=f2sum.get(x1,y1);
			}
		//double FBS,FBSref;
		FBS=FBS+((frac1-frac2)*(frac1-frac2));
		FBSref=FBSref+frac1*frac1+frac2*frac2;
		}

	if (FBSref != 0)
		return(1.0-FBS/FBSref);
	else
		return(0);
	}

// fast calculation of FSS for a square neigberhood. Missing values are treated as zero preciaption
// the code is parralelized with openMP
double calculate_FSS_fast_by_Roberts(const CField &f, int n) const
	{
	// calculate the summed fields
	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	return(f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(f2sum,n));
	}

double calculate_FSS_fast_by_Roberts_alternative_border_treatment(const CField &f, int n) const
	{
	// calculate the summed fields
	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	return(f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation_alternative_border_treatment(f2sum,n));
	}

// fast calculation of FSS for a square neigberhood. Missing values are treated as zero preciaption
// the code is parralelized with openMP. Multiple neigberhoods are calculated simultaniusly
// po testih ta stvar dela za 1 n precej pocasneje kot verzija brez multiple_neighbourhood - priblizno enako pa je pri 2-h ali 3-h neigberhoodsih in takrat postane ta boljša
vector <double> calculate_FSS_fast_by_Roberts_multiple_neighbourhoods(const CField &f, const vector <int> &nvec) const
	{
	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();

	if (nvec.size() < 1)
		error(AT,FU, "nvec.size()");

	vector <double> FSSvec (nvec.size(),0);
	#pragma omp parallel for
	for (long in=0; in < (long)nvec.size(); in++)
		FSSvec[in]=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(f2sum,nvec[in]);
	return(FSSvec);
	}


vector <double> calculate_FSS_fast_by_Roberts_multiple_neighbourhoods_alternative_border_treatment(const CField &f, const vector <int> &nvec) const
	{
	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();

	if (nvec.size() < 1)
		error(AT,FU, "nvec.size()");

	vector <double> FSSvec (nvec.size(),0);
	#pragma omp parallel for
	for (long in=0; in < (long)nvec.size(); in++)
		FSSvec[in]=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation_alternative_border_treatment(f2sum,nvec[in]);
	return(FSSvec);
	}


// fast calculation of FSS for a square neigberhood. Missing values are treated as zero preciaption
// the code is parralelized with openMP. Multiple neigberhoods are calculated simultaniusly
// po testih ta stvar dela za 1 n precej pocasneje kot verzija brez multiple_neighbourhood - priblizno enako pa je pri 2-h ali 3-h neigberhoodsih in takrat postane ta boljša
vector <double> calculate_FSS_fast_by_Roberts_multiple_neighbourhoods_in_subdomain(const CField &f, const vector <int> &nvec, const CField &fdomain) const
	{
	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();

	if (nvec.size() < 1)
		error(AT,FU, "nvec.size()");

	vector <double> FSSvec (nvec.size(),0);
	#pragma omp parallel for
	for (long in=0; in < (long)nvec.size(); in++)
		FSSvec[in]=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation_in_subdomain(f2sum,nvec[in],fdomain);
	return(FSSvec);
	}

// fast calculation of neigberhood size where FSS=0.5 - using bisection method - assuming a montonic increase of FSS with n
// the returned value is the first n (assuming monotonicity) that has a FSS value LARGER than 0.5
long calculate_usefull_neighborhood_size_using_bisection(const CField &f) const
	{
	ERRORIF(!check_if_field_is_binary_allow_missing_values());
	ERRORIF(!f.check_if_field_is_binary_allow_missing_values());

	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();


    long x1=1;
    long x2=2*max(f1sum.dimx,f1sum.dimy) + 1;

	// domain size needs to be at least 2
	if (x2 < 5)
		error(AT,FU, "domain size needs to be at least 2");

	double FSS1=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(f2sum,x1);
	double FSS2=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(f2sum,x2);

	// special case when FSS>0.5 at n=1
	if (FSS1>0.5)
		return(1);

	// special case when FSS does never reach value 0.5
	if (FSS2 <= 0.5)
		error(AT,FU, "FSS does never reach value 0.5!");

    do
	    {
        // select middle point
        long xnew=(x1 + x2)/2;
        // if xnew is even add 1
        if (xnew%2==0) xnew++;
        // calulcate FSS vale at xnew
        double FSSnew=f1sum.calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(f2sum,xnew);

        // move x1 or x2 to the middle point according to FSS value at xnew
        if (FSSnew > 0.5)
        	{
			x2=xnew;
			FSS2=FSSnew;
			}
		else
			{
			x1=xnew;
			FSS1=FSSnew;
			}
	    }
    while ( x2 - x1 > 2 );

	return(x2);
	}


double calculate_dFSS_using_bisection_return_minus_one_on_large_bias(const CField &f) const
	{
	if (get_sum() == 0 || f.get_sum() == 0)
		return(-1);

	double bias1=max(get_sum()/f.get_sum(),f.get_sum()/get_sum());

	if ( bias1 > 2)
		return(-1);

	double FSSn1=calculate_FSS_fast_by_Roberts(f, 1) ;

	double dFSS=0;
	if (FSSn1 == 1)
		{
		dFSS=0;
		}
	else
		{
		CField f1bin=*this;
		CField f2bin=f;

		// remove overlapping areas
		f1bin.remove_overlapping_areas_in_binary_fields(f2bin);

		if (f1bin.get_sum() == 0 || f2bin.get_sum() == 0)
			return(-1);

		double bias2=max(f1bin.get_sum()/f2bin.get_sum(),f2bin.get_sum()/f1bin.get_sum());
		if ( bias2 > 1.5)
			return(-1);

		double dFSS05=floor((double)f1bin.calculate_usefull_neighborhood_size_using_bisection(f2bin)/2.);

		dFSS=(1-FSSn1)*dFSS05;
		}
	return(dFSS);
	}


// calculation of neigberhood size at FSS=0.5  - above is a faster BISECTION alghorithm that usually works faster
int calculate_usefull_neigberhood_size_brute_force_fast(const CField &f) const
	{
	vector <int> nvec;
	vector <double> FSSvec;
	int nmax=2*max(dimx,dimy)+1;
	int dn=nmax/40/2;
	int counter=0;

	//cout << nmax << "\t" << dn << endl;

	// calculate is steps of 1/20 od the possible n values
	for (int n=1; n < nmax; n+=2)
		{
		nvec.push_back(n);

		// try to find if the n(FSS=0.5) was allready passed
		if (counter > dn)
			{
			counter=0;
			FSSvec=calculate_FSS_fast_by_Roberts_multiple_neighbourhoods(f, nvec);

			for (long in=0; in < (long)nvec.size(); in++)
				{
				//cout << nvec[in] << "\t" << FSSvec[in] << endl;
				if (FSSvec[in]>=0.5) return(nvec[in]);
				}
			nvec.clear();
			}

		counter++;
		}

	// if no FSS>0.5 value has been found return error
	error(AT,FU, "no FSS>0.5 value found for any nighberhood size");
	return(0);
	}



double calculate_FSS_brute_force_outside_of_domain_is_value_zero_sum_is_done_only_inside_domain_with_square_neigberhood_mask(CField &f, int n) const
	{
	int ix,iy,x,y;
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)n*(double)n;

	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;


	//double count;
	for (ix=0; ix < dimx; ix++)
	for (iy=0; iy < dimy; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (x=ix-r; x <= ix+r; x++)
		for (y=iy-r; y <= iy+r; y++)
			{
			if (check_coordianes(x,y))
				{
				//count++;
				fr1+=get(x,y);
				fr2+=f.get(x,y);
				}
			}
		//fr1*=1.0/count;
		//fr2*=1.0/count;
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	if (MSEref != 0)
		return(1.0-MSE/MSEref);
	else
		return(0);
	}

double calculate_FSS_brute_force_outside_of_domain_is_value_zero_sum_is_done_only_inside_domain_with_square_neigberhood_mask_return_minusone_if_MSRef_zero(CField &f, int n)
	{
	int ix,iy,x,y;
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)n*(double)n;

	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;


	//double count;
	for (ix=0; ix < dimx; ix++)
	for (iy=0; iy < dimy; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (x=ix-r; x <= ix+r; x++)
		for (y=iy-r; y <= iy+r; y++)
			{
			if (check_coordianes(x,y))
				{
				//count++;
				fr1+=get(x,y);
				fr2+=f.get(x,y);
				}
			}
		//fr1*=1.0/count;
		//fr2*=1.0/count;
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	if (MSEref != 0)
		return(1.0-MSE/MSEref);
	else
		return(-1);
	}



double calculate_FSS_brute_periodic_domain_x_any_y_square_neigberhood_mask(CField &f, int n)
	{
	int ix,iy,x,y;
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)n*(double)n;

	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;


	//double count;
	for (ix=0; ix < dimx; ix++)
	for (iy=0; iy < dimy; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (x=ix-r; x <= ix+r; x++)
		for (y=iy-r; y <= iy+r; y++)
			{
				fr1+=get_periodic_x_and_y(x,y);
				fr2+=f.get_periodic_x_and_y(x,y);
			}
		//fr1*=1.0/count;
		//fr2*=1.0/count;
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	return(1.0-MSE/MSEref);
	}

double calculate_FSS_brute_force_on_the_subdomain_square_neighborhood(CField &f, int n, int xmin, int ymin, int xmax, int ymax)
	{
	int ix,iy,x,y;
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)n*(double)n;

	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;


	//double count;
	for (ix=xmin; ix <= xmax; ix++)
	for (iy=ymin; iy <= ymax; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (x=ix-r; x <= ix+r; x++)
		for (y=iy-r; y <= iy+r; y++)
			{
			if (check_coordianes(x,y))
				{
				//count++;
				fr1+=get(x,y);
				fr2+=f.get(x,y);
				}
			}
		//fr1*=1.0/count;
		//fr2*=1.0/count;
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	double inside_ares_size=(double)(xmax-xmin+1)*(double)(ymax-ymin+1);
	MSE*=1.0/inside_ares_size;
	MSEref*=1.0/inside_ares_size;
	if (MSEref != 0)
		return(1.0-MSE/MSEref);
	else
		return(0);
	//	error(AT,FU, "MSEref == 0");
	}

double calculate_FSS_circular_neigberhood_equdistant_grid_using_summed_fields(CField &f, int n)
	{
	double MSE=0;
	double MSEref=0;

	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;

	CField fconv1=convolution_fast_using_summed_fields(r);
	CField fconv2=f.convolution_fast_using_summed_fields(r);

	for (long il=0; il < fconv1.size(); il++)
		{
		double fr1=fconv1.data[il];
		double fr2=fconv2.data[il];
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	return(1.0-MSE/MSEref);
	}

// this code might not be entirely correct - needs to be rechecked
double calculate_FSS_faster_circular_neigberhood_equdistant_grid_not_completely_finished(CField &f, int n) const
	{
	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	CField fsum1=*this;
	CField fsum2=f;

	// replace missing values with zeros
	fsum1.Replace_values(BAD_DATA_FLOAT,0);
	fsum2.Replace_values(BAD_DATA_FLOAT,0);

	// calculate the summed fields - rows only
	#pragma omp parallel for
	for (int iy=0; iy < dimy; iy++)
	for (int ix=1; ix < dimx; ix++)
		{
		fsum1.set(ix,iy,fsum1.get(ix-1,iy) + fsum1.get(ix,iy));
		fsum2.set(ix,iy,fsum2.get(ix-1,iy) + fsum2.get(ix,iy));
		}

	int r;
	r=(n-1)/2+1;

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	// caluclate fractions
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)r*(double)r*M_PI;
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (int y=iy-r; y <= iy+r; y++)
		if ( y < dimy && y >= 0)
			{
			double dy=y-iy;
			double rr=r;
			int x1=floor(sqrt(rr*rr-dy*dy))-1;
			int x2=floor(sqrt(rr*rr+dy*dy));

			if (x2 >= dimx) x2=dimx-1;

			if (x1<0)
				{
				fr1+=fsum1.get(x2,y);
				fr2+=fsum2.get(x2,y);
				}
			else
				{
				fr1+=fsum1.get(x2,y)-fsum1.get(x1,y);
				fr2+=fsum2.get(x2,y)-fsum2.get(x1,y);;
				}
			}
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	if (MSEref != 0)
		return(1.0-MSE/MSEref);
	else
		return(0);
	}

// this code might not be entirely correct - needs to be rechecked
double calculate_FSS_faster_square_neigberhood_equdistant_grid_not_completely_finished(CField &f, int n) const
	{
	if (!compare_have_they_the_same_dimensions(f))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f)");

	CField fsum1=*this;
	CField fsum2=f;

	// replace missing values with zeros
	fsum1.Replace_values(BAD_DATA_FLOAT,0);
	fsum2.Replace_values(BAD_DATA_FLOAT,0);

	// calculate the summed fields - rows only
	#pragma omp parallel for
	for (int iy=0; iy < dimy; iy++)
	for (int ix=1; ix < dimx; ix++)
		{
		fsum1.set(ix,iy,fsum1.get(ix-1,iy) + fsum1.get(ix,iy));
		fsum2.set(ix,iy,fsum2.get(ix-1,iy) + fsum2.get(ix,iy));
		}

	int r;
	r=(n+1)/2-1;

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	// caluclate fractions
	double MSE=0;
	double MSEref=0;
	double fr1,fr2;
	double neighborhood_area;
	neighborhood_area=(double)n*(double)n;
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		//count=0;
		fr1=fr2=0;
		for (int y=iy-r; y <= iy+r; y++)
		if ( y < dimy && y >= 0)
			{
			int x1=ix-r-1;
			int x2=ix+r;

			if (x2 >= dimx) x2=dimx-1;

			if (x1<0)
				{
				fr1+=fsum1.get(x2,y);
				fr2+=fsum2.get(x2,y);
				}
			else
				{
				fr1+=fsum1.get(x2,y)-fsum1.get(x1,y);
				fr2+=fsum2.get(x2,y)-fsum2.get(x1,y);;
				}
			}
		fr1*=1.0/neighborhood_area;
		fr2*=1.0/neighborhood_area;
		MSE+=(fr1-fr2)*(fr1-fr2);
		MSEref+=fr1*fr1+fr2*fr2;
		}
	MSE*=1.0/(double)size();
	MSEref*=1.0/(double)size();
	if (MSEref != 0)
		return(1.0-MSE/MSEref);
	else
		return(0);
	}


// calculates the FSSwind for a single neigberhood size. The input fields should be class fields containing only integer values with the first class haveing a value of 1
double calculate_the_FSSwind_fast(const CField &f2, int n, const long  number_of_classes) const
	{
	if (!compare_have_they_the_same_dimensions(f2))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f2sum)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	// Get max index
	double max_index=max(get_maximum(),f2.get_maximum());
	ERRORIF(number_of_classes + 1 < max_index);

	// Decompose into binary fields by classes
	vector <CField> f1bin;
	vector <CField> f2bin;
	CField ftemp=f2;
	ftemp.Set_all_values(0);

	for (long il=0; il <= number_of_classes + 1; il++)
		{
		f1bin.push_back(ftemp.create_a_copy_with_data());
		f2bin.push_back(ftemp.create_a_copy_with_data());
		}


	for (long il=0; il < size(); il++)
		{
		// the 0 and smaller values are ignored - for example BAD_DATA_FLOAT
		if (data[il] > 0) f1bin[data[il]].data[il]=1;
		if (f2.data[il] > 0) f2bin[f2.data[il]].data[il]=1;
		}

	/*for (long il=0; il <= number_of_classes + 1; il++)
		cout <<  "\t" << f1bin[il].get_sum();
	cout << endl;
	for (long il=0; il <= number_of_classes + 1; il++)
		cout <<  "\t" << f2bin[il].get_sum();
	cout << endl;*/


    // draw image and netcdf of binary fields
	/*
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		ostringstream s1;
		s1.str("");
		s1 << "f1bin_" << il;
		WriteField(s1.str()+".nc",f1bin[il]);
		draw_pngimage_of_binary_field(s1.str()+".png", f1bin[il]);

		s1.str("");
		s1 << "f2bin_" << il;
		WriteField(s1.str()+".nc",f2bin[il]);
		draw_pngimage_of_binary_field(s1.str()+".png", f2bin[il]);
		}
		*/

	// Calculate summed fields for all classes
	vector <CField> f1sum;
	vector <CField> f2sum;
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		f1sum.push_back((f1bin[il].calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation()).create_a_copy_with_data());
		f2sum.push_back((f2bin[il].calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation()).create_a_copy_with_data());
		}

	// Calculate fraction fields for all classes
	vector <CField> f1frac;
	vector <CField> f2frac;
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		f1frac.push_back((f1sum[il].calculate_fractions_from_the_summed_field(n)).create_a_copy_with_data());
		f2frac.push_back((f2sum[il].calculate_fractions_from_the_summed_field(n)).create_a_copy_with_data());
		}

    // write netcdf of fraction fields
	/*
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		ostringstream s1;
		s1.str("");
		s1 << "f1frac_" << il;
		WriteField(s1.str()+".nc",f1frac[il]);

		s1.str("");
		s1 << "f2frac_" << il;
		WriteField(s1.str()+".nc",f2frac[il]);
		}
	*/
	// Calculate FSSwind
	double sumOM2=0;
	double sumO2=0;
	double sumM2=0;
	for (long ic=0; ic <= number_of_classes + 1; ic++)
	for (long il=0; il < f1frac[ic].size(); il++)
		{
		sumOM2+=sqr(f1frac[ic].data[il]-f2frac[ic].data[il]);
		sumO2+=sqr(f1frac[ic].data[il]);
		sumM2+=sqr(f2frac[ic].data[il]);
		}

	// sumO2+sumM2 should be larger than zero since FSSwind will have classes defined at all locations in the domain
	if (sumO2+sumM2 <= 0)
		error(AT,FU, "sumO2+sumM2 <= 0");

	return(1-sumOM2/(sumO2+sumM2));
	}



// calculates the FSSwind for multiple neigberhood sizes - the calulacion is not optimized since each core repeats the calculation of the binary and summed fields. The input fields should be class fields containing only integer values with the first class haveing a value of 1
vector <double> calculate_the_FSSwind_fast_multiple_neigberhood_not_optimized(const CField &f2, const vector <long> &nvec, const long  number_of_classes) const
	{
	vector <double> FSSwindvec (nvec.size(),0);

	#pragma omp parallel for
	for (long in=0; in < (long)nvec.size(); in++)
		FSSwindvec[in]=calculate_the_FSSwind_fast(f2, nvec[in], number_of_classes);

	return(FSSwindvec);
	}

// calculates the FSSwind for multiple neigberhood sizes. The input fields should be class fields containing only integer values with the first class haveing a value of 1
vector <double> calculate_the_FSSwind_fast_multiple_neigberhoods_optimized(const CField &f2, const vector <long> &nvec, const long  number_of_classes) const
	{
	ERRORIF(!compare_have_they_the_same_dimensions(f2));

	// Get max index
	double max_index=max(get_maximum(),f2.get_maximum());
	ERRORIF(number_of_classes + 1 < max_index);


	// Decompose into binary fields by classes
	vector <CField> f1bin;
	vector <CField> f2bin;
	CField ftemp=f2;
	ftemp.Set_all_values(0);

	for (long il=0; il <= number_of_classes + 1; il++)
		{
		f1bin.push_back(ftemp.create_a_copy_with_data());
		f2bin.push_back(ftemp.create_a_copy_with_data());
		}

	for (long il=0; il < size(); il++)
		{
		// the 0 and smaller values are ignored - for example BAD_DATA_FLOAT
		if (data[il] > 0) f1bin[data[il]].data[il]=1;
		if (f2.data[il] > 0) f2bin[f2.data[il]].data[il]=1;
		}

	// Calculate summed fields for all classes
	vector <CField> f1sum;
	vector <CField> f2sum;
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		f1sum.push_back((f1bin[il].calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation()).create_a_copy_with_data());
		f2sum.push_back((f2bin[il].calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation()).create_a_copy_with_data());
		}

	vector <double> FSSwindvec (nvec.size(),0);

	#pragma omp parallel for
	for (long in=0; in < (long)nvec.size(); in++)
		FSSwindvec[in]=calculate_the_FSSwind_fast_from_summed_field(f1sum, f2sum, nvec[in]);

	return(FSSwindvec);
	}


double calculate_the_FSSwind_fast_from_summed_field(const vector <CField> &f1sum, const vector <CField> &f2sum, int n) const
	{
	ERRORIF(n%2 != 1);

	// Calculate fraction fields for all classes
	vector <CField> f1frac;
	vector <CField> f2frac;
	for (long il=0; il < (long)f1sum.size(); il++)
		{
		f1frac.push_back((f1sum[il].calculate_fractions_from_the_summed_field(n)).create_a_copy_with_data());
		f2frac.push_back((f2sum[il].calculate_fractions_from_the_summed_field(n)).create_a_copy_with_data());
		}

    // write netcdf of fraction fields
	/*
	for (long il=0; il <= number_of_classes + 1; il++)
		{
		ostringstream s1;
		s1.str("");
		s1 << "f1frac_" << il;
		WriteField(s1.str()+".nc",f1frac[il]);

		s1.str("");
		s1 << "f2frac_" << il;
		WriteField(s1.str()+".nc",f2frac[il]);
		}
	*/
	// Calculate FSSwind
	double sumOM2=0;
	double sumO2=0;
	double sumM2=0;
	for (long ic=0; ic < (long)f1sum.size(); ic++)
	for (long il=0; il < f1frac[ic].size(); il++)
		{
		sumOM2+=sqr(f1frac[ic].data[il]-f2frac[ic].data[il]);
		sumO2+=sqr(f1frac[ic].data[il]);
		sumM2+=sqr(f2frac[ic].data[il]);
		}

	// sumO2+sumM2 should be larger than zero since FSSwind will have classes defined at all locations in the domain
	if (sumO2+sumM2 <= 0)
		error(AT,FU, "sumO2+sumM2 <= 0");

	return(1-sumOM2/(sumO2+sumM2));
	}



void remove_overlapping_areas_in_binary_fields(CField &f2)
	{
	ERRORIF(!compare_have_they_the_same_dimensions(f2));

	for (long il=0; il < size(); il++)
		if (data[il]==1 && f2.data[il]==1)
			{
			data[il]=0;
			f2.data[il]=0;
			}
	}


CField calculate_fractions_from_the_summed_field(int n) const
	{
	CField f_frac;
	f_frac=*this;

	int r;
	r=(n+1)/2-1;

	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		double frac1;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r - 1;
		ix + r >= dimx   ? x2 = dimx - 1 : x2 = ix + r;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r - 1;
		iy + r >= dimy   ? y2 = dimy - 1 : y2 = iy + r;

		frac1=get(x2,y2);

		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			}
		f_frac.set(ix,iy,frac1/((double)n*(double)n));
		}

	return(f_frac);
	}

CField calculate_fractions_from_the_summed_field_n_can_be_even(int n) const
	{
	CField f_frac;
	f_frac=*this;

	int r1=n/2;
	int r2=n-n/2-1;

	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		double frac1;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r1 - 1;
		ix + r2 >= dimx   ? x2 = dimx - 1 : x2 = ix + r2;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r1 - 1;
		iy + r2 >= dimy   ? y2 = dimy - 1 : y2 = iy + r2;

		frac1=get(x2,y2);

		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			}
		f_frac.set(ix,iy,frac1/((double)n*(double)n));
		}

	return(f_frac);
	}

CField calculate_fractions_from_binary_field(int n) const
	{
	//cout << n << endl;
	CField fsum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField ffrac=fsum.calculate_fractions_from_the_summed_field_n_can_be_even(n);

	//ffrac.multiply_by_double(1.0/((double)n*(double)n));
	return (ffrac);
	}

CField smooth_using_square_kernel(int n) const
	{
	//cout << n << endl;

	ERRORIF(Check_missing_data(BAD_DATA_FLOAT) > 0);

	CField fsum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f_smoothed=fsum.calculate_fractions_from_the_summed_field_n_can_be_even(n);

	//ffrac.multiply_by_double(1.0/((double)n*(double)n));
	return (f_smoothed);
	}

// calculate precipitation neighborhood score value for specific n
double calculate_PNS(const CField &f2, long n, double p) const
	{
	ERRORIF(!compare_have_they_the_same_dimensions(f2));

	ERRORIF(n%2 != 1);

	CField f1_smooth=smooth_using_square_kernel(n);
	CField f2_smooth=f2.smooth_using_square_kernel(n);


	double top=0;
	double bottom=0;
	for (long il=0; il < f1_smooth.size(); il++)
		{
		double v1=f1_smooth.data[il];
		double v2=f2_smooth.data[il];
		top+=pow(fabs(v1-v2),p);
		bottom+=pow(fabs(v1),p)+pow(fabs(v2),p);
		}

	if ( bottom != 0)
		return(1.0-top/bottom);
	else
		{
		ERRORIF(bottom == 0 );
		return(0);
		}
	}

// calculate precipitation neighborhood score value for specific n
double calculate_PNS_from_summed_fields(const CField &f2sum, long n) const
	{
	ERRORIF(!compare_have_they_the_same_dimensions(f2sum));

	ERRORIF(n%2 != 1);

	CField f1_smooth=calculate_fractions_from_the_summed_field_n_can_be_even(n);
	CField f2_smooth=f2sum.calculate_fractions_from_the_summed_field_n_can_be_even(n);

	double top=0;
	double bottom=0;
	for (long il=0; il < f1_smooth.size(); il++)
		{
		double v1=f1_smooth.data[il];
		double v2=f2_smooth.data[il];
		top+=fabs(v1-v2);
		bottom+=fabs(v1)+fabs(v2);
		}

	if ( bottom != 0)
		return(1.0-top/bottom);
	else
		{
		ERRORIF(bottom == 0 );
		return(0);
		}
	}


// fast calculation of neigberhood size where PNS=0.5 - using bisection method - assuming a montonic increase of PNS with n
// the returned value is the first n (assuming monotonicity) that has a FSS value LARGER than 0.5
double calculate_dPNS_using_bisection_from_non_overlapping_fields(const CField &f, double Q) const
	{
	// check missing data
	ERRORIF(!check_if_field_has_no_missing_values());
	ERRORIF(!f.check_if_field_has_no_missing_values());

	// check if the fields really contain only non-overlaping precipitation
	for (long il=0; il < size(); il++)
		ERRORIF(data[il] != 0 && f.data[il] != 0);

	// check if the fields are not empty - if so the displacment is zero since the original fields were identical
	if (get_sum() == 0)
		return(0);

	CField f1sum=calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();

    long x1=1;
    long x2=2*max(f1sum.dimx,f1sum.dimy) + 1;

	// domain size needs to be at least 2
	if (x2 < 5)
		error(AT,FU, "domain size needs to be at least 2");
	double PNS1=f1sum.calculate_PNS_from_summed_fields(f2sum,x1);
	double PNS2=f1sum.calculate_PNS_from_summed_fields(f2sum,x2);

	// special case when PNS>0.5 at n=1
	if (PNS1>0.5)
		return(1);

	// special case when PNS does never reach value 0.5
	if (PNS2 <= 0.5)
		error(AT,FU, "PNS does never reach value 0.5!");

    do
	    {
        // select middle point
        long xnew=(x1 + x2)/2;
        // if xnew is even add 1
        if (xnew%2==0) xnew++;
        // calulcate FSS vale at xnew
        double PNSnew=f1sum.calculate_PNS_from_summed_fields(f2sum,xnew);

        // move x1 or x2 to the middle point according to FSS value at xnew
        if (PNSnew > 0.5)
        	{
			x2=xnew;
			PNS2=PNSnew;
			}
		else
			{
			x1=xnew;
			PNS1=PNSnew;
			}
	    }
    while ( x2 - x1 > 2 );

	//cout << x2 << "\t" << floor((double)x2/2.0*(1.0-Q)) << endl;
	//cout << x2 << endl;

	return(((double)(x2-1)/2.0*(1.0-Q)));
	}

// fast calculation of dPNS
double calculate_dPNS(const CField &f, double &Q, double &R) const
	{
	// check missing data
	ERRORIF(!check_if_field_has_no_missing_values());
	ERRORIF(!f.check_if_field_has_no_missing_values());

	// Calulate and remove bias
	double sumf1=get_sum();
	double sumf2=f.get_sum();

	ERRORIF(sumf1 == 0 || sumf2 == 0);

	R=sumf1/sumf2;

	CField f1unbiased=*this;
	CField f2unbiased=f;
	f2unbiased.multiply_by_double(R);
	//cout << bias << endl;

	double subtracted_precipitation_sum;
	CField f1subtracted=f1unbiased;
	CField f2subtracted=f2unbiased;
	f1subtracted.generate_subtracted_precipation_fields_for_PNS_calculation(f2subtracted,subtracted_precipitation_sum);
	Q=subtracted_precipitation_sum/sumf1;
	return(f1subtracted.calculate_dPNS_using_bisection_from_non_overlapping_fields(f2subtracted, Q));
	}

// fast calculation of dPNS
double calculate_dPNS_with_enlarged_domain(const CField &f, double &Q, double &R, long grow_domain) const
	{
	CField f1enlarged=grow_domain_by_n_pixels(grow_domain, 0);
	CField f2enlarged=f.grow_domain_by_n_pixels(grow_domain, 0);

	return(f1enlarged.calculate_dPNS(f2enlarged, Q, R));
	}

// fast calculation of dPNS
vector <double> calculate_PNS_for_multiple_neighborhoods(const CField &f, double &R, const vector <long> nvec) const
	{
	// check missing data
	ERRORIF(!check_if_field_has_no_missing_values());
	ERRORIF(!f.check_if_field_has_no_missing_values());

	// Calulate and remove bias
	double sumf1=get_sum();
	double sumf2=f.get_sum();

	ERRORIF(sumf1 == 0 || sumf2 == 0);

	R=sumf1/sumf2;

	CField f1unbiased=*this;
	CField f2unbiased=f;
	f2unbiased.multiply_by_double(R);

	CField f1sum=f1unbiased.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();
	CField f2sum=f2unbiased.calculate_the_summed_field_in_both_directions_used_by_fast_FSS_calculation();

	vector <double> PNSvec;
	for (unsigned long in=0; in < nvec.size(); in++)
		PNSvec.push_back(f1sum.calculate_PNS_from_summed_fields(f2sum, nvec[in]));

	return(PNSvec);
	}


void generate_subtracted_precipation_fields_for_PNS_calculation(CField &f2, double &subtracted_precipitation_sum)
	{
	subtracted_precipitation_sum=0;
	for (long il=0; il < size(); il++)
		{
		// subtract the lower value from both fields
		double val=min(data[il],f2.data[il]);
		data[il]-=val;
		f2.data[il]-=val;
		subtracted_precipitation_sum+=val;
		}
	}

/*

double calculate_the_FSS_value_from_the_summed_field_used_by_fast_FSS_calculation(const CField &f2sum, int n) const
	{
	if (!compare_have_they_the_same_dimensions(f2sum))
		error(AT,FU, "!compare_have_they_the_same_dimensions(f2sum)");

	if (n%2 != 1)
		error(AT,FU, "n%2 != 1");

	int r;
	r=(n+1)/2-1;

	// calculate the fractions and FBS and FBSref
	double FBS=0;
	double FBSref=0;
	// the parallelized loop actualy runs slower on two pocessor than non-parralelized - no idea why
	//#pragma omp parallel for default(shared) reduction(+:FBS,FBSref)
	for (int ix=0; ix < dimx; ix++)
	for (int iy=0; iy < dimy; iy++)
		{
		double frac1;
		double frac2;
		double x1,x2,y1,y2;
		//ix - r < 0       ? x1 = 0        : x1 = ix - r;
		x1 = ix - r - 1;
		ix + r >= dimx   ? x2 = dimx - 1 : x2 = ix + r;
		//iy - r < 0       ? y1 = 0        : y1 = iy - r;
		y1 = iy - r - 1;
		iy + r >= dimy   ? y2 = dimy - 1 : y2 = iy + r;

		frac1=get(x2,y2);
		frac2=f2sum.get(x2,y2);

		// subctract the left part if needed
		if (x1 >= 0)
			{
			frac1+=-get(x1,y2);
			frac2+=-f2sum.get(x1,y2);
			}
		// subctract the bottom part if needed
		if (y1 >= 0)
			{
			frac1+=-get(x2,y1);
			frac2+=-f2sum.get(x2,y1);
			}

		// put back in the left bottom corner part if needed
		if (x1 >= 0 && y1 >= 0)
			{
			frac1+=get(x1,y1);
			frac2+=f2sum.get(x1,y1);
			}
		//double FBS,FBSref;
		FBS=FBS+((frac1-frac2)*(frac1-frac2));
		FBSref=FBSref+frac1*frac1+frac2*frac2;
		}

	if (FBSref != 0)
		return(1.0-FBS/FBSref);
	else
		return(0);
	}
*/

	// calculate threshold by percentiles - the top 5% should be specified as percentile=0.95
	double calculate_threshold_by_percentiles(double percentile)
		{
		vector <double> temp=data;
		long nth=floor((double)temp.size()*percentile);
		nth_element (temp.begin(), temp.begin()+nth, temp.end());
		return(temp[nth]);
		}

	// calculate Mean-error distance  - the code is not thourogly checked and might be wrong
	double calculate_MED(const CField f) const
		{
		CField f1=*this;
		CField f2=f;

		ERRORIF(!f1.check_if_initilized());
		ERRORIF(!f2.check_if_initilized());

		ERRORIF(!f1.check_if_field_is_binary_allow_missing_values());
		ERRORIF(!f2.check_if_field_is_binary_allow_missing_values());

		ERRORIF(!f1.check_if_field_has_no_missing_values());
		ERRORIF(!f2.check_if_field_has_no_missing_values());

		ERRORIF(!f1.compare_have_they_the_same_dimensions_and_coordinates(f2));

		if (f1.get_sum() == 0 || f2.get_sum() == 0)
			return(-1);

		// initialize distance field
		CField f_dist=f1;
		f_dist.Set_all_values(BAD_DATA_FLOAT);

		// generate border areas from the other field
		//CField f_border1=f1.get_border_of_area(0.5, 1.5, 1);
		CField f_border2=f2.get_border_of_area(0.5, 1.5, 1);

		//vector <long> border_list1;
		vector <long> border_list2;
		for (long il=0; il < f1.size(); il++)
			{
			//if (f_border1.data[il] == 1)
			//	border_list1.push_back(il);
			if (f_border2.data[il] == 1)
				border_list2.push_back(il);
			}

		// go over all the non-ovelapping pixels in f1 and calcualte the closest border pixel in f2
		for (long il=0; il < f1.size(); il++)
			{
			if (f1.data[il] == 1)
				{
				//all points that overlap have distance zero
				if (f2.data[il] == 1)
					f_dist.data[il]=0;
				else
					{
					int x1,y1,x2,y2;
					DataLoc2XY(il, x1, y1);
					// loop over all border points in f2
					double min_dst=1E100;
					for (unsigned long ib=0; ib < border_list2.size(); ib++)
						{
						DataLoc2XY(border_list2[ib], x2, y2);
						double dst=squared_euclidian_distance(x1, y1, x2, y2);
						if (min_dst > dst) min_dst=dst;
						}
					f_dist.data[il]=sqrt(min_dst);
					}
				}
			}

		/*WriteField("aaa1.nc",f2);
		WriteField("aaa2.nc",f_border2);

		CField ftemp=f_border2;
		for (long il=0; il < f1.size(); il++)
			ftemp.data[il]+=f2.data[il];

		WriteField("aaa3.nc",ftemp);
		WriteField("aaa4.nc",f_dist);
		*/

		double MED=0;
		for (long il=0; il < f1.size(); il++)
			if (f_dist.data[il] > 0)
				MED+=f_dist.data[il];
		MED=MED/(double)f1.get_sum();

		return(MED);
		}

	// gest the locations of local minimums
	void calculate_locations_of_local_extreems_in_a_scalar_field(vector <long> &index_of_minimums, vector <long> &index_of_maximums) const
		{
		double l,r,t,b,c;

		index_of_minimums.clear();
		index_of_maximums.clear();

		for (int ix=1; ix < dimx-1;ix++)
		for (int iy=1; iy < dimy-1;iy++)
			{
			l=get(ix-1,iy);
			r=get(ix+1,iy);
			t=get(ix,iy-1);
			b=get(ix,iy+1);
			c=get(ix,iy);

			if ( c < l && c < r && c < t && c < b)
				index_of_minimums.push_back(XY2DataLoc(ix, iy));

			if ( c > l && c > r && c > t && c > b)
				index_of_maximums.push_back(XY2DataLoc(ix, iy));
			}
		}




	// gest the locations of local minimums which represent smallest values in a specified radius
	void calculate_locations_of_local_extreems_in_a_scalar_field_using_a_radius_check(double radius, vector <long> &index_of_minimums, vector <long> &index_of_maximums) const
		{
		calculate_locations_of_local_extreems_in_a_scalar_field(index_of_minimums, index_of_maximums);

		vector <long> index_of_minimums_temp=index_of_minimums;
		index_of_minimums.clear();

		double radius_sqr;
		radius_sqr=radius*radius;


		while (index_of_minimums_temp.size() > 0)
			{
			bool all_ok=true;
			// check all the other mimumuns agains the first one
			for (long il=1; il < (long)index_of_minimums_temp.size(); il++)
				// check if they are closer than radius
				if (calculate_sqr_distance_between_two_grid_points_specified_by_DataLoc(index_of_minimums_temp.front(), index_of_minimums_temp[il]) < radius_sqr)
					{
					all_ok=false;
					if (data[index_of_minimums_temp.front()] < data[index_of_minimums_temp[il]])
						index_of_minimums_temp.erase(index_of_minimums_temp.begin()+il);
					else
						index_of_minimums_temp.erase(index_of_minimums_temp.begin());
					break;
					}
			if (all_ok)
				{
				index_of_minimums.push_back(index_of_minimums_temp.front());
				index_of_minimums_temp.erase(index_of_minimums_temp.begin());
				}

			}
		}

	// gest the locations of local minimums which represent smallest values in a specified radius
	void calculate_locations_of_local_extreems_in_a_scalar_field_using_a_radius_check_where_you_demand_a_span_of_values(double radius, double reqired_deltavalue, vector <long> &index_of_minimums, vector <long> &index_of_maximums) const
		{
		calculate_locations_of_local_extreems_in_a_scalar_field(index_of_minimums, index_of_maximums);

		double max,min;

		vector <long> index_of_minimums_old=index_of_minimums;
		index_of_minimums.clear();
		for (long il=0; il < (long)index_of_minimums_old.size(); il++)
			{
			get_max_and_min_in_a_cricle(DataLoc2X(index_of_minimums_old[il]),DataLoc2Y(index_of_minimums_old[il]),radius, max, min);
			if (min == data[index_of_minimums_old[il]] && max-min >= reqired_deltavalue)
				index_of_minimums.push_back(index_of_minimums_old[il]);
			}

		vector <long> index_of_maximums_old=index_of_maximums;
		index_of_maximums.clear();
		for (long il=0; il < (long)index_of_maximums_old.size(); il++)
			{
			get_max_and_min_in_a_cricle(DataLoc2X(index_of_maximums_old[il]),DataLoc2Y(index_of_maximums_old[il]),radius, max, max);
			if (max == data[index_of_maximums_old[il]] && max-max >= reqired_deltavalue)
				index_of_maximums.push_back(index_of_maximums_old[il]);
			}



		}



	// This crop only works on grids with parralel axis
	CField crop(int x1,int y1,int x2,int y2) const
		{
		CField Field2;
		int ix,iy;
		double lat1,lon1;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->crop ERROR: lat_data.size()==0 || lon_data.size()==0 ""\n";exit(1);}

		// dimenzije. Uposteva ce gre cez datumsko mejo
		if (x2 >= x1) Field2.dimx=x2-x1+1;
		else Field2.dimx=dimx-x1+x2+1;
		Field2.dimy=y2-y1+1;

		Field2.allocatememory();
		Field2.allocatememory_lat_lon_array();

		for (ix=0 ; ix < Field2.dimx; ix++)
			for (iy=0 ; iy < Field2.dimy; iy++)
				{
				if (ix+x1<dimx)
					{
					Field2.set(ix,iy,get(ix+x1,iy+y1));
					Field2.lon_data[ix]=lon_data[ix+x1];
					Field2.lat_data[iy]=lat_data[iy+y1];
					}
				else
					{
					Field2.set(ix,iy,get(ix+x1-dimx,iy+y1));
					Field2.lon_data[ix]=lon_data[ix+x1-dimx];
					Field2.lat_data[iy]=lat_data[iy+y1];
					}
				}


		// doloci se nove lat/lon konstante
		XY2LatLon(x1,y1, lon1, lat1);

		return(Field2);
		}

	// horizontal merge - the *this file is on the left side
	CField merge_horizontal(const CField &f) const
		{
		if (dimy!=f.dimy)
			error(AT,FU, "dimy!=f.dimy");

		CField fout;
		fout.allocatememory_and_set_dimx_dimy(dimx+f.dimx,dimy);
		// copy data
		for (int ix=0; ix < dimx+f.dimx; ix++)
			for (int iy=0; iy < dimy; iy++)
				{
				if (ix < dimx)
					fout.set(ix,iy,get(ix,iy));
				else
					fout.set(ix,iy,f.get(ix-dimx,iy));
				}

		// copy lon nad lat
		fout.lat_data=lat_data;
		for (int ix=0; ix < dimx+f.dimx; ix++)
			{
			if (ix < dimx)
				fout.lon_data.push_back(lon_data[ix]);
			else
				fout.lon_data.push_back(f.lon_data[ix-dimx]);
			}
		return(fout);
		}

	// horizontal merge - the *this file is on the top side
	// DOES NOT WORK
	CField merge_vertical(const CField &f) const
		{
		ERRORIF(dimx!=f.dimx);

		//sdasds
		CField fout;
		fout.allocatememory_and_set_dimx_dimy(dimx,dimy+f.dimy);
		// copy data
		for (int ix=0; ix < dimx; ix++)
			for (int iy=0; iy < dimy +f.dimy; iy++)
				{
				if (iy < dimy)
					fout.set(ix,iy,get(ix,iy));
				else
					fout.set(ix,iy,f.get(ix,iy-dimy));
				}

		// copy lon nad lat
		fout.lon_data=lon_data;
		for (int iy=0; iy < dimy+f.dimy; iy++)
			{
			if (iy < dimy)
				fout.lat_data.push_back(lat_data[iy]);
			else
				fout.lat_data.push_back(f.lat_data[iy-dimy]);
			}
		return(fout);
		}

	// merge a vector of CFields
	void merge_horizontal_vector_of_fields(const vector <CField> &f)
		{
		ERRORIF(f.size() < 1);
		*this=f[0];
		CField ff2;
		for (long il=1; il < (long)f.size(); il++)
			{
			ff2=this->merge_horizontal(f[il]);
			*this=ff2;
			}
		}

	// merge a vector of CFields -> from top to bottom
	void merge_vertical_vector_of_fields(const vector <CField> &f)
		{
		ERRORIF(f.size() < 1);

		*this=f.back();
		CField ff2;
		for (long il=f.size()-2; il >= 0; il--)
			{
			ff2=this->merge_vertical(f[il]);
			*this=ff2;
			}
		}

	// resize_canvas old Field is in the botom left corner
	void resize_canvas(long new_dimx,long new_dimy, double value)
		{
		ERRORIF(new_dimx < dimx || new_dimy < dimy );

		CField fout;
		fout.allocatememory_and_set_dimx_dimy(new_dimx,new_dimy);
	    fout.allocatememory_lat_lon_array();
		fout.Set_all_values(value);

		for (int ix=0; ix < dimx; ix++)
			for (int iy=0; iy < dimy; iy++)
			 fout.set(ix,iy,get(ix,iy));
		*this=fout;
		}


	// trim on the left and right sides horizontal
	void trim_horizontal(double trim_value, long leave_space)
		{
		long right_trim=0;
		for (int ix=0; ix < dimx; ix++)
			for (int iy=0; iy < dimy; iy++)
				if (get(ix,iy) != trim_value)
					right_trim=ix;

		long left_trim=0;
		for (int ix=dimx-1; ix >= 0; ix--)
			for (int iy=0; iy < dimy; iy++)
				if (get(ix,iy) != trim_value)
					left_trim=ix;

		// leave some space on the sides
		left_trim=max(0,left_trim-leave_space);
		right_trim=min(dimy-1,right_trim+leave_space);
		//cout << left_trim << "\t" << right_trim << endl;

		CField f_trimmed=crop(left_trim,0,right_trim,dimy-1);
		*this=f_trimmed;
		}

	void mirror_vertical()
		{
		CField f;
		f=*this;
		for (int ix=0; ix < dimx; ix++)
			for (int iy=0; iy < dimy; iy++)
				set(ix,iy,f.get(ix,dimy-1-iy));

		// mirror lat array
		for (int iy=0; iy < dimy; iy++)
			lat_data[iy]=f.lat_data[dimy-1-iy];
		}

	CField crop_by_lat_and_lon(double lon_min, double lon_max, double lat_min, double lat_max)
		{
		int xmin,xmax,ymin,ymax;

		if (!check_if_initilized())
			error(AT,FU, " check_if_initilized()");

		if (lon_min >= lon_max || lat_min >= lat_max)
			error(AT,FU, " lon_min >= lon_max || lat_min >= lat_max");

		LatLon2XY_new(lon_min, lat_min, xmin, ymin);
		LatLon2XY_new(lon_max, lat_max, xmax, ymax);

		if (ymin == ymax || xmin == xmax)
			error(AT,FU, "ymin == ymax || xmin == xmax");


		return(crop(xmin,ymin,xmax, ymax));
		}

	// creates a o copy of CField. Copies only atributes, does not alocatememory or copy data
	/*CField create_a_copy()
		{
		CField f;

		f.dimx=dimx;
		f.dimy=dimy;
		f.unique_timestamp=unique_timestamp;

		int ix;

		if (lat_data.size()==0 || lon_data.size()==0)
			error(AT,FU, " lat_data.size()==0 || lon_data.size()==0");

		// copy lat lon data if available
		if (lat_data.size() != 0 && lon_data.size() != 0)
		  {
		    f.allocatememory_lat_lon_array();

		    for (ix=0; ix < f.dimx; ix++)
		      f.lon_data[ix]=lon_data[ix];
		    for (ix=0; ix < f.dimy; ix++)
		      f.lat_data[ix]=lat_data[ix];
		  }

		return(f);
		}
	*/

	// creates a o copy of CField. Also alocates memory and copies data
	CField create_a_copy_with_data() const
		{
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		CField f;
		f.dimx=dimx;
		f.dimy=dimy;
		f.unique_timestamp=unique_timestamp;
		f.lon_data=lon_data;
		f.lat_data=lat_data;
		f.data=data;
		return(f);
		}


	bool check_compare_dimensions(CField &f) const
		{
		if (dimx != f.dimx) return(false);
		if (dimy != f.dimy) return(false);
		return(true);
		}


	// Checks if any point has a bad_data_double value
	long Check_missing_data(double bad_data_double) const
		{
		long count=0,ixl;

		for (ixl=0;ixl<size();ixl++)
			if (data[ixl]==bad_data_double) count++;

		return (count);
		}

	// performs convolution with radius R. It does not support global domains.
	// MISSING data and OUT_OF_DOMAIN areas are treated as areas vith precipitation value=0mm
	CField convolution_fast_using_summed_fields(int R) const
		{
		CField fsum=*this;

		// replace missing values with zeros
		fsum.Replace_values(BAD_DATA_FLOAT,0);

		// calculate the summed fields
		// sum by rows first - we can skip the left row since it allready contains the valeus from f1
		#pragma omp parallel for
		for (int iy=0; iy < dimy; iy++)
		for (int ix=1; ix < dimx; ix++)
			fsum.set(ix,iy,fsum.get(ix-1,iy) + fsum.get(ix,iy));

		// make table of sqrt(R*R-ix*ix) values so that it will be faster
		vector <int> temp(R+1,0);
		for (int ix=0;ix<=R;ix++)
			temp[ix]=(int)floor(sqrt((double)R*(double)R-(double)ix*(double)ix));

		CField fconv=*this;
		#pragma omp parallel for
		for (int iy=0; iy < dimy; iy++)
		for (int ix=0; ix < dimx; ix++)
			{
			double sum=0;
			double count=0;
			for (int dy=-R;dy<=R;dy++)
				{
				count+=2*temp[abs(dy)]+1;
				// inside domain
				if (check_coordianes(ix, iy+dy))
					{
					int x1=ix-temp[abs(dy)]-1;
					int x2=ix+temp[abs(dy)];

					if (x2 > dimx-1)
						sum+=fsum.get(dimx-1,iy+dy);
					else
						sum+=fsum.get(x2,iy+dy);

					if (x1 < 0 )
						sum-=0;
					else
						sum-=fsum.get(x1,iy+dy);
					}
				}
			fconv.set(ix,iy,sum/count);
			}

		return (fconv);
		}

	// performs convolution with radius R. It does not support global domains.
	// MISSING data and OUT_OF_DOMAIN areas are treated as areas vith precipitation value=0mm
	CField convolution(int R) const
		{
		int ix,iy,dx,dy;
		vector <int> temp(R+1,0);
		int count;
		double t,xxx;
		CField f;

		/*if (Check_missing_data(BAD_DATA_FLOAT) != 0)
			{
			cout << "ERROR: Convoloved field has missing data \n";
			exit(1);}*/

		f=create_a_copy_with_data();

		// make table of sqrt(R*R-ix*ix) values so that it will be faster
		for (ix=0;ix<=R;ix++)
			temp[ix]=(int)floor(sqrt((double)R*(double)R-(double)ix*(double)ix));

		/*for (ix=0;ix<=R;ix++)
			cout << temp[ix] << "\n";*/

		for (iy=0;iy<dimy;iy++)
			{
			//cout << iy << endl;
			for (ix=0;ix<dimx;ix++)
				{
				//missing value check
				if (get(ix,iy)==BAD_DATA_FLOAT)
					f.set(ix,iy,BAD_DATA_FLOAT);
				else
					{
					t=0;
					count=0;
					for (dy=-R;dy<=R;dy++)
					for (dx=-temp[abs(dy)];dx<=temp[abs(dy)];dx++)
						{
						// out of domain check
						if (ix+dx<0 || ix+dx>=dimx || iy+dy<0 || iy+dy>=dimy )
							count++;
						else
							{
							xxx=get(ix+dx,iy+dy);
							// missing value check
							if (xxx!=BAD_DATA_FLOAT)
								{t=t+xxx;}
							count++;
							}
						}
					f.set(ix,iy,t/(double)count);
					}
				}
			}

		return (f);
		}



	// performs time convolution with radius Rs and Rt. It does not support global domains.
	// MISSING data and OUT_OF_DOMAIN areas are treated as areas vith precipitation value=0mm
	CField simple_convolution_in_time_also(int Rs, int Rt, const vector <CField> &fields) const
		{
		int ix,iy,dx,dy;
		vector <int> temp(Rs+1,0);
		int count;
		double t,xxx;
		CField f;

		if (Rt != 1)
			error(AT,FU, " Rt != 1. Only that is supported so far!");

		if (fields.size() != 3)
			error(AT,FU, " fields.size() != 3. Only that is supported so far!");

		for (ix=0;ix<(int)fields.size();ix++)
			if (!fields[ix].check_if_initilized())
				error(AT,FU, " !fields[ix].check_if_initilized()");

		for (ix=1;ix<(int)fields.size();ix++)
			if (!fields[ix].compare_have_they_the_same_dimensions_and_coordinates(fields[ix-1]))
				error(AT,FU, " !fields[ix].compare_have_they_the_same_dimensions_and_coordinates(fields[ix-1])");

		if (*this != fields[Rt])
			error(AT,FU, " *this != fields[Rt]!");

		f=create_a_copy_with_data();

		// make table of sqrt(R*R-ix*ix) values so that it will be faster
		for (ix=0;ix<=Rs;ix++)
			temp[ix]=(int)floor(sqrt((double)Rs*(double)Rs-(double)ix*(double)ix));

		for (ix=0;ix<=Rs;ix++)
			cout << temp[ix] << "\n";
		exit(1);

		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			//missing value check
			if (get(ix,iy)==BAD_DATA_FLOAT)
				f.set(ix,iy,BAD_DATA_FLOAT);
			else
				{
				t=0;
				count=0;
				// ------------------------------
				// current time
				// ------------------------------
				for (dy=-Rs;dy<=Rs;dy++)
				for (dx=-temp[abs(dy)];dx<=temp[abs(dy)];dx++)
					{
					// out of domain check
					if (ix+dx<0 || ix+dx>=dimx || iy+dy<0 || iy+dy>=dimy )
						count++;
					else
						{
						xxx=get(ix+dx,iy+dy);
						// missing value check
						if (xxx!=BAD_DATA_FLOAT)
							{t=t+xxx;}
						count++;
						}
					}

				// ------------------------------
				// next time step
				// ------------------------------
				xxx=fields[Rt+1].get(ix,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				count++;

				// surronding pixels
				xxx=fields[Rt+1].get_but_check_coordianes(ix+1,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt+1].get_but_check_coordianes(ix-1,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt+1].get_but_check_coordianes(ix,iy+1);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt+1].get_but_check_coordianes(ix,iy-1);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				count+=4;

				// ------------------------------
				// previous time step
				// ------------------------------
				xxx=fields[Rt-1].get(ix,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				count++;

				// surronding pixels
				xxx=fields[Rt-1].get_but_check_coordianes(ix+1,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt-1].get_but_check_coordianes(ix-1,iy);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt-1].get_but_check_coordianes(ix,iy+1);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				xxx=fields[Rt-1].get_but_check_coordianes(ix,iy-1);
				if (xxx!=BAD_DATA_FLOAT)
					t=t+xxx;
				count+=4;

				f.set(ix,iy,t/(double)count);
				}
			}

		return (f);
		}



	// performs convolution with convolution_radius - in a non equdistand grid. The grid axis still have to be at parallel. It does not support global domains.
	// MISSING data and OUT_OF_DOMAIN areas are treated as areas vith precipitation value=0mm
	CField convolution_in_non_equidistand_grid(double convolution_radius) const
		{
		int ix,iy,dx,dy,iz;
		int count;
		double t,xxx;
		CField f;
		int Rx,Ry;
		int r_temp;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->convolution_in_non_equidistand_grid ERROR: lat_data.size()==0 || lon_data.size()==0 ""\n";exit(1);}

		f=create_a_copy_with_data();

		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			//missing value check
			if (get(ix,iy)==BAD_DATA_FLOAT)
				f.set(ix,iy,BAD_DATA_FLOAT);
			else
				{
				// determining Rx
				iz=ix;
				if (ix < dimx/2)
					while (fabs(lon_data[ix] -lon_data[iz+1]) <= convolution_radius)
						{iz++;}
				if (ix >= dimx/2)
					while (fabs(lon_data[ix] -lon_data[iz-1]) <= convolution_radius)
						{iz--;}
				Rx=abs(ix-iz);

				// determining Ry
				iz=iy;
				if (iy < dimy/2)
					while (fabs(lat_data[iy] -lat_data[iz+1]) <= convolution_radius)
						{iz++;}
				if (iy >= dimy/2)
					while (fabs(lat_data[iy] -lat_data[iz-1]) <= convolution_radius)
						{iz--;}
				Ry=abs(iy-iz);

				//cout << ix << " " << iy << " " << Rx << " " << Ry << endl;

				t=0;
				count=0;
				for (dy=-Ry;dy<=Ry;dy++)
					{

					if (Rx==0 || Ry == 0) r_temp=0;
					else r_temp=(int)floor((double)Rx* sqrt(1-(double)dy*(double)dy/((double)Ry*(double)Ry)));

					for (dx=-r_temp;dx<=r_temp;dx++)
						{
						// out of domain check
						if (ix+dx<0 || ix+dx>=dimx || iy+dy<0 || iy+dy>=dimy )
							count++;
						else
							{
							xxx=get(ix+dx,iy+dy);
							// missing value check
							if (xxx!=BAD_DATA_FLOAT)
								{t=t+xxx;}
							count++;
							}
						}
					}
				f.set(ix,iy,t/(double)count);
				}
			}

		return (f);
		}
	CField convolution_in_non_equidistand_grid2(double convolution_radius) const
		{
		int ix,iy,dx,dy,iz;
		int count;
		double t,xxx;
		CField f;
		int Rx,Ry;
		int r_temp;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->convolution_in_non_equidistand_grid ERROR: lat_data.size()==0 || lon_data.size()==0 ""\n";exit(1);}

		f=create_a_copy_with_data();

		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			//missing value check
			if (get(ix,iy)==BAD_DATA_FLOAT)
				f.set(ix,iy,BAD_DATA_FLOAT);
			else
				{
				// determining Rx
				iz=ix;
				if (ix < dimx/2)
					while (fabs(lon_data[ix] -lon_data[iz+1]) <= convolution_radius)
						{iz++;}
				if (ix >= dimx/2)
					while (fabs(lon_data[ix] -lon_data[iz-1]) <= convolution_radius)
						{iz--;}
				Rx=abs(ix-iz);

				// determining Ry
				iz=iy;
				if (iy < dimy/2)
					while (fabs(lat_data[iy] -lat_data[iz+1]) <= convolution_radius)
						{iz++;}
				if (iy >= dimy/2)
					while (fabs(lat_data[iy] -lat_data[iz-1]) <= convolution_radius)
						{iz--;}
				Ry=abs(iy-iz);

				//cout << ix << " " << iy << " " << Rx << " " << Ry << endl;

				t=0;
				count=0;
				for (dy=-Ry;dy<=Ry;dy++)
					{

					if (Rx==0 || Ry == 0) r_temp=0;
					else r_temp=(int)floor((double)Rx* sqrt(1-(double)dy*(double)dy/((double)Ry*(double)Ry)));

					for (dx=-r_temp;dx<=r_temp;dx++)
						{
						count++;
						xxx=get_but_check_coordianes(ix+dx,iy+dy);
						if (xxx!=BAD_DATA_FLOAT)
							t=t+xxx;
						}
					}
				f.set(ix,iy,t/(double)count);
				}
			}

		return (f);
		}


	// performs convolution with convolution_radius - in a non equdistand grid. The grid axis still have to be at parallel. It assumes the domain is global in x direction.
	// MISSING data and OUT_OF_DOMAIN areas are treated as areas vith precipitation value=0mm
	CField convolution_in_non_equidistand_grid_with_global_domain_support(double convolution_radius, bool domain_global_in_x) const
		{
		int ix,iy,dx,dy,iz;
		int count;
		double t,xxx;
		CField f;
		//int Rx,Ry;
		int r_temp;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->convolution_in_non_equidistand_grid ERROR: lat_data.size()==0 || lon_data.size()==0 ""\n";exit(1);}

		f=create_a_copy_with_data();

		// perform nonglobal convolution
		//f=convolution_in_non_equidistand_grid(convolution_radius);

		// -----------------------
		// tabulate Ry,Rx
		// -----------------------
		vector <int> Ry;
		vector <int> Rx;
		for (iy=0;iy<dimy;iy++)
			{
			iz=iy;
			if (iy < dimy/2)
				while (fabs(lat_data[iy] -lat_data[iz+1]) <= convolution_radius)
					{iz++;}
			if (iy >= dimy/2)
				while (fabs(lat_data[iy] -lat_data[iz-1]) <= convolution_radius)
					{iz--;}
			Ry.push_back(abs(iy-iz));
			}
		for (ix=0;ix<dimx;ix++)
			{
			// determining Rx
			iz=ix;
			if (ix < dimx/2)
				while (fabs(lon_data[ix] -lon_data[iz+1]) <= convolution_radius)
					{iz++;}
			if (ix >= dimx/2)
				while (fabs(lon_data[ix] -lon_data[iz-1]) <= convolution_radius)
					{iz--;}
			Rx.push_back(abs(ix-iz));
			}


		// -----------------------
		// perform convolution
		// -----------------------
		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			//missing value check
			if (get(ix,iy)==BAD_DATA_FLOAT)
				f.set(ix,iy,BAD_DATA_FLOAT);
			else
				{
				//cout << ix << " " << iy << " " << Rx[ix] << " " << Ry[iy] << endl;

				t=0;
				count=0;
				for (dy=-Ry[iy];dy<=Ry[iy];dy++)
					{
					if (Rx[ix]==0 || Ry[iy] == 0) r_temp=0;
					else r_temp=(int)floor((double)Rx[ix]* sqrt(1-(double)dy*(double)dy/((double)Ry[iy]*(double)Ry[iy])));

					for (dx=-r_temp;dx<=r_temp;dx++)
						{
						count++;
						xxx=get_but_check_coordianes(ix+dx,iy+dy);
						if (xxx!=BAD_DATA_FLOAT)
							t=t+xxx;
						}
					}
				f.set(ix,iy,t/(double)count);
				}
			}

		return (f);
		}


	// sets the values < threshold to -1. It does not change missing values
	CField Threshold_cut(double threshold, CField &f_before_convolution) const
		{
		CField f;
		long ixl;


		f=create_a_copy_with_data();
		f.Set_all_values(BAD_DATA_FLOAT);

		// changes threshold if threshold is 0 - since there are no valid values below zero
		if (threshold==0) threshold=0.0001;

//		cout << threshold << endl;

		for (ixl=0;ixl<size();ixl++)
			{
			if (data[ixl]==BAD_DATA_FLOAT) f.data[ixl]=BAD_DATA_FLOAT;
			else if (data[ixl]<threshold) f.data[ixl]=-1;
			else f.data[ixl]=f_before_convolution.data[ixl];
			}

		return (f);
		}

	// Identifies and Enumerates objects - a new much faster object identification algorithm - developed july 2010
	CField Identify_Objects(CGroup_of_Identified_Objects_2D &Group_of_Identified_Objects_2D, bool is_domain_periodic_in_x) const
		{
		CField f_obj;

		if (lat_data.size()==0 || lon_data.size()==0)
			error(AT,FU, " lat_data.size()==0 || lon_data.size()==0");


		long il;

		f_obj=create_a_copy_with_data();
		//f_obj.allocatememory();

		//f_obj.Set_all_values(-1);
		for (il=0; il < f_obj.size(); il++)
			{
			if (data[il] == BAD_DATA_FLOAT) f_obj.data[il]=BAD_DATA_FLOAT;
			else f_obj.data[il]=-1;
			}

		int ix,iy;
		double count=0;
		double top,left,right;
		//long ixl;

		// initialize vector memory for first 1000 potential objects
		CIdentified_Object_2D t0,t2,t1;
		Group_of_Identified_Objects_2D.freememory();
		Group_of_Identified_Objects_2D.Identified_Object_2D.assign(1000,t0);

		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			// non-missing and inside object check
			if (get(ix,iy)>-1)
				{
				//cout << "aaa";
				top=f_obj.get_but_check_coordianes(ix,iy-1);
				left=f_obj.get_but_check_coordianes(ix-1,iy);

				// both have different objects
				if (top > -1 && left > -1 && top != left)
					{
					f_obj.set(ix,iy,top);

					// replace values in f_obj field
					t2 = Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left);

					//cout  << ix  << " " << iy << endl;
					//cout << left << " " <<  top << " " <<  t2.id << " " << t2.xmin  << " " <<t2.xmax<< " " <<t2.ymin<< " " <<t2. ymax << endl;


					f_obj.Replace_values_only_in_certain_area(left,top,t2.xmin,t2.xmax,t2.ymin,t2. ymax);

					// calculate the new object area
					t1=Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)top);
					if (t2.xmax > t1.xmax) t1.xmax=t2.xmax;
					if (t2.xmin < t1.xmin) t1.xmin=t2.xmin;
					if (t2.ymax > t1.ymax) t1.ymax=t2.ymax;
					if (t2.ymin < t1.ymin) t1.ymin=t2.ymin;
					Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)top)=t1;

					// erase the left object from vector
					Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left).id=-1;
					}

				// both have same object
				else if (top > -1 && left > -1 && top==left)
					{
					f_obj.set(ix,iy,top);
					}

				// only top has object
				else if (top > -1 && left < 0)
					{
					f_obj.set(ix,iy,top);
					Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)top).ymax=iy;
					}

				// only left has object
				else if (top < 0 && left > -1)
					{
					f_obj.set(ix,iy,left);

					// check if left object does not wind around top to larger xmax than ix. If not, then reset the xmax to ix
					if (Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left).xmax < ix)
						Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left).xmax=ix;
					}

				// new object
				else
					{
					f_obj.set(ix,iy,count);

					// resize vector by 1000 if more objects are found than are defined in vector
					if (count > Group_of_Identified_Objects_2D.Identified_Object_2D.size() - 2)
						Group_of_Identified_Objects_2D.Identified_Object_2D.resize(Group_of_Identified_Objects_2D.Identified_Object_2D.size() + 1000,t0);

					// add new object
					Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)count).set_object_borders(ix,ix,iy,iy);
					Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)count).id=count;

					count++;
					}
				}
			//else  f_obj.set(ix,iy,get(ix,iy));
			}


		// Modify if domain is periodic - join objects on left-right borders
		if (is_domain_periodic_in_x)
			for (iy=0; iy < f_obj.dimy; iy++)
				{
				left=f_obj.get(0,iy);
				right=f_obj.get(f_obj.dimx-1,iy);
				if (left != BAD_DATA_FLOAT && right != BAD_DATA_FLOAT && left > -1 && right > -1)
					if (left != right)
						{
						cout << "Global join: " << iy << "\t" << left << "\t" << right << "\t" << endl;

						// get object info
						t1=Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left);
						t2=Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)right);

						// Replace values in f_obj
						f_obj.Replace_values_only_in_certain_area(left,right,t1.xmin,t1.xmax,t1.ymin,t1.ymax);

						// calculate the new object area
						t2.xmin=0;
						t2.xmax=f_obj.dimx-1;
						if (t1.ymax > t2.ymax) t2.ymax=t1.ymax;
						if (t1.ymin < t2.ymin) t2.ymin=t1.ymin;
						Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)right)=t2;

						// erase the left object from vector
						Group_of_Identified_Objects_2D.Identified_Object_2D.at((long)left).id=-1;
						}
				}


/*		double counter2;
		for (ix=0; ix < Group_of_Identified_Objects_2D.Identified_Object_2D.size(); ix++)
			{
			counter2=0;
			for (il=0; il < f_obj.size(); il++)
				if (f_obj.data[il]==ix) counter2++;
			if (counter2 > 0)
			 	cout << ix <<  "\t" << counter2 << endl;
			}
*/
		Group_of_Identified_Objects_2D.compress();


/*
		for (ix=0; ix < Group_of_Identified_Objects_2D.Identified_Object_2D.size(); ix++)
			//if (Group_of_Identified_Objects_2D.Identified_Object_2D[ix].id != -1)
				{
				count=0;
				for (ixl=0; ixl < f_obj.size(); ixl++)
					if (f_obj.data[ixl]==Group_of_Identified_Objects_2D.Identified_Object_2D[ix].id) count++;
			 	t2=Group_of_Identified_Objects_2D.Identified_Object_2D[ix];
			 	cout << Group_of_Identified_Objects_2D.Identified_Object_2D[ix].id <<  "\t" << count << "\t" << t2.xmin  << " " <<t2.xmax<< " " <<t2.ymin<< " " <<t2. ymax << endl;
				}
		*/
		return(f_obj);
		}

	// old is the original_value_which_is_to_be_replaced, new_value in the new value
	void floodfill(int x, int y, double old, double new_value)
		{
		if (get(x,y)==old)
			{
			set(x,y,new_value);
			if (x < dimx - 1) floodfill(x+1, y  , old, new_value);
			if (x > 0       ) floodfill(x-1, y  , old, new_value);
			if (y < dimy - 1) floodfill(x  , y+1, old, new_value);
			if (y > 0       ) floodfill(x  , y-1, old, new_value);
			}
		}

	// old is the original_value_which_is_to_be_replaced, new_value in the new value
	void floodfill_periodic_sphere(int x, int y, double old, double new_value)
		{
		if (get_periodic_sphere(x,y)==old)
			{
			set_periodic_sphere(x,y,new_value);
			floodfill_periodic_sphere(x+1, y  , old, new_value);
			floodfill_periodic_sphere(x-1, y  , old, new_value);
			floodfill_periodic_sphere(x  , y+1, old, new_value);
			floodfill_periodic_sphere(x  , y-1, old, new_value);
			}
		}


	// old is the original_value_which_is_to_be_replaced, new_value in the new value
	void floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(int x, int y, double old, double new_value, vector <double> &touching_values_list, double dont_log_values_lower_than_this, vector <long> &list_of_affected_points)
		{
		double val=get(x,y);
		if (val==old)
			{
			set(x,y,new_value);
			list_of_affected_points.push_back(XY2DataLoc(x,y));
			if (x < dimx - 1) floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(x+1, y  , old, new_value, touching_values_list, dont_log_values_lower_than_this, list_of_affected_points);
			if (x > 0       ) floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(x-1, y  , old, new_value, touching_values_list, dont_log_values_lower_than_this, list_of_affected_points);
			if (y < dimy - 1) floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(x  , y+1, old, new_value, touching_values_list, dont_log_values_lower_than_this, list_of_affected_points);
			if (y > 0       ) floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(x  , y-1, old, new_value, touching_values_list, dont_log_values_lower_than_this, list_of_affected_points);
			}
		// else log the value - dont log if this value is dont_log_this_value
		else if (val >= dont_log_values_lower_than_this && val != new_value)
			{
			// search if this value was allready logged
			if (!binary_search(touching_values_list.begin(), touching_values_list.end(),val))
				{
				//cout << val << endl;
				// add the value to list
				touching_values_list.push_back(val);
				// sort - so that we can later use the very fast binary_search
				sort(touching_values_list.begin(), touching_values_list.end());
				}
			}
		}



	void Identify_2D_objects_using_floodfill_and_cascading_threshold(vector <double> &thrs_list, double lowest_object_id)
		{
		long il,it;
		double val;

		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		if (thrs_list.size() < 1)
			error(AT,FU, "thrs_list.size() < 1");

		// check if lowest_object_id is high enough
		if (thrs_list.size() +1 > lowest_object_id)
			error(AT,FU, "thrs_list.size() +1 > lowest_object_id");

		// check if all values fell between -1 and thrs_list.size()
		for (il=0; il < size(); il++)
			if (data[il] < -1 || data[il] > thrs_list.size() - 1)
				error(AT,FU, "data[il] < -1 || data[il] > thrs_list.size()");

		// test write fields
		//ostringstream s1;
		//s1.str("");
		//s1 << "aaa";
		//WriteField(s1.str()+".nc",f);
		//WriteField_PNG(s1.str()+".png",  f);

		// do the cascading object identification
		double counter=lowest_object_id;
		vector <double> touching_values_list;
		vector <long> list_of_affected_points;
		vector <long> list_of_affected_points_new;
		for (it=0; it < (long)thrs_list.size(); it++)
			{
			val=it;
			for (int x=0; x < dimx; x++)
			for (int y=0; y < dimy; y++)
				if (get(x,y) == val)
					{
					// do the flood fill
					touching_values_list.clear();
					list_of_affected_points.clear();
					floodfill_and_log_the_touching_values_and_make_a_list_of_all_affected_points(x,y, val, counter,touching_values_list, lowest_object_id, list_of_affected_points);

					//cout << val << "\t" << x << ":" << y << "\t" << list_of_affected_points.size() << "\t" << touching_values_list.size() << ":" ;
					//for (long io=0; io < (long)touching_values_list.size(); io++)
					//	cout << touching_values_list[io] << ",";
					//cout << endl;

					// if there is only one touching value - set to id of the touching value
					if (touching_values_list.size() == 1)
						for (il=0; il < (long)list_of_affected_points.size(); il++)
							data[list_of_affected_points[il]]=touching_values_list.back();

					// if there is more than one touching object - set the ids to -2 so that the growing algorithm can be applied later
					if (touching_values_list.size() > 1)
						for (il=0; il < (long)list_of_affected_points.size(); il++)
							data[list_of_affected_points[il]]=-2;

					// increase counter only if it is a new object
					if (touching_values_list.size() == 0)
						counter++;
					}

			Grow_object_area_set_by_allow_grow_index(-2, lowest_object_id, 1E100);

			// test write fields
			//ostringstream s1;
			//s1.str("");
			//s1 << "aaa" << it;
			//WriteField(s1.str()+".nc",f);
			//WriteField_PNG(s1.str()+".png",  f);
			}

		//return(f);
		}


	void Calculate_object_attributes(vector <CIdentified_Object_2D> &objlist, double lower_obj_index_limit) const
		{
		long il;

		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		// check if any 1 is present in the dataset
		for (il=0; il < size(); il++)
			if (data[il] < -1 )
				error(AT,FU, "data[il] < -1");

		// get max_object_id
		double max_object_id=get_maximum();

		// assign memory to object list
		CIdentified_Object_2D temp_o;
		temp_o.xmin=30000;
		temp_o.xmax=-1;
		temp_o.ymin=30000;
		temp_o.ymax=-1;

		CIdentified_Object_2D *op;
		objlist.assign((long)max_object_id+1,temp_o);

		// Area size and center
		for (int x=0; x < dimx; x++)
		for (int y=0; y < dimy; y++)
			if (get(x,y) >= lower_obj_index_limit)
				{
				op=&objlist[(long)get(x,y)];
				if (x < op->xmin) op->xmin=x;
				if (x > op->xmax) op->xmax=x;
				if (y < op->ymin) op->ymin=y;
				if (y > op->ymax) op->ymax=y;
				op->id=get(x,y);
				}
		}

	// check that all the id-s found in CGroup_of_Identified_Objects_2D are found are found in f_obj and vice versa
	void check_consistency_of_f_obj_Field_with_Group_of_Identified_Objects_2D(CGroup_of_Identified_Objects_2D &Group_of_Identified_Objects_2D) const
		{
		bool passed=true;

		long il, ik;

		long index;
		CIdentified_Object_2D t1;

		bool found;

		for (il=0; il < size(); il++)
			if (data[il] > -1)
				if (!Group_of_Identified_Objects_2D.check_if_object_is_defined(data[il], t1, index))
					{
					cout << "Could not find object: " << data[il] << " in Group_of_Identified_Objects_2D!" << endl;
					passed=false;
					}

		for (ik=0; ik < (long)Group_of_Identified_Objects_2D.Identified_Object_2D.size(); ik++)
			if (Group_of_Identified_Objects_2D.Identified_Object_2D[ik].id != -1)
				{
				found=false;
				for (il=0; il < size(); il++)
					if (data[il] == Group_of_Identified_Objects_2D.Identified_Object_2D[ik].id )
						found=true;

				if (found==false)
					{
					cout << "Could not find object: " << Group_of_Identified_Objects_2D.Identified_Object_2D[ik].id << " in f_obj!" << endl;
					passed=false;
					}
				}

		if (passed==false)
			error(AT,FU, " a problem!");

		}

	bool Check_if_any_neigboring_point_has_a_value_in_interval(int x, int y, double lower_limit, double upper_limit, double &v) const
		{
		if (x < dimx - 1) {v=get(x+1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		if (x > 0       ) {v=get(x-1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		if (y < dimy - 1) {v=get(x  , y+1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		if (y > 0       ) {v=get(x  , y-1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		return(false);
		}

	bool Check_if_any_neigboring_point_has_a_value_in_interval_periodic_x_y(int x, int y, double lower_limit, double upper_limit, double &v) const
		{
		v=get_periodic_x_y(x+1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_x_y(x-1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_x_y(x  , y+1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_x_y(x  , y-1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		return(false);
		}

	bool Check_if_any_neigboring_point_has_a_value_in_interval_periodic_x(int x, int y, double lower_limit, double upper_limit, double &v) const
		{
		v=get_periodic_x_y(x+1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_x_y(x-1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		if (y < dimy - 1) {v=get(x  , y+1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		if (y > 0       ) {v=get(x  , y-1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);}
		return(false);
		}

	bool Check_if_any_neigboring_point_has_a_value_in_interval_periodic_sphere(int x, int y, double lower_limit, double upper_limit, double &v) const
		{
		v=get_periodic_sphere(x+1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_sphere(x-1, y  ); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_sphere(x  , y+1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		v=get_periodic_sphere(x  , y-1); if (is_inside_interval(lower_limit, upper_limit, v)) return(true);
		return(false);
		}

	// in field f the allowed growth area is defined with allow_grow_index
	void Grow_object_area_set_by_allow_grow_index(double allow_grow_index, double lower_obj_index_limit, double upper_obj_index_limit)
		{
		long il;
		int x,y;
		double temp_obj;

		vector <Point> old_points;
		vector <Point> new_points;
		Point p;

		// identify starting points
		for (y=0;y<dimy;y++)
		for (x=0;x<dimx;x++)
			if (get(x,y) == allow_grow_index)
				if (Check_if_any_neigboring_point_has_a_value_in_interval(x,y,lower_obj_index_limit,upper_obj_index_limit,temp_obj))
					{
					p.set(x,y);
					p.lat=temp_obj;
					old_points.push_back(p.deep_copy());
					//cout << p.output() << endl;
					}


		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;
			new_points.clear();
			for (il=0; il < (long)old_points.size(); il++)
				{
				//cout << il << endl;
				x=old_points[il].x;
				y=old_points[il].y;
				temp_obj=old_points[il].lat;


				//check if the value was allready set by one of the previous objects - so we dont repeat ourselves
				if (get(x,y) == allow_grow_index)
					{
					// set the value of old_point
					set(x,y,temp_obj);

					// check neigboring points
					if (x < dimx - 1) if (get(x+1, y  ) == allow_grow_index) {p.set(x+1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (x > 0       ) if (get(x-1, y  ) == allow_grow_index) {p.set(x-1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (y < dimy - 1) if (get(x  , y+1) == allow_grow_index) {p.set(x  , y+1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (y > 0       ) if (get(x  , y-1) == allow_grow_index) {p.set(x  , y-1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					}
				}
			old_points=new_points;
			}
		}

	// in field f the allowed growth area is defined with allow_grow_index - periodic x and y
	void Grow_object_area_set_by_allow_grow_index_periodic_x_y(double allow_grow_index, double lower_obj_index_limit, double upper_obj_index_limit)
		{
		long il;
		int x,y;
		double temp_obj;

		vector <Point> old_points;
		vector <Point> new_points;
		Point p;

		// identify starting points
		for (y=0;y<dimy;y++)
		for (x=0;x<dimx;x++)
			if (get(x,y) == allow_grow_index)
				if (Check_if_any_neigboring_point_has_a_value_in_interval_periodic_x_y(x,y,lower_obj_index_limit,upper_obj_index_limit,temp_obj))
					{
					p.set(x,y);
					p.lat=temp_obj;
					old_points.push_back(p.deep_copy());
					//cout << p.output() << endl;
					}


		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;
			new_points.clear();
			for (il=0; il < (long)old_points.size(); il++)
				{
				//cout << il << endl;
				x=old_points[il].x;
				y=old_points[il].y;
				temp_obj=old_points[il].lat;


				//check if the value was allready set by one of the previous objects - so we dont repeat ourselves
				if (get(x,y) == allow_grow_index)
					{
					// set the value of old_point
					set(x,y,temp_obj);

					// check neigboring points
					if (get_periodic_x_y(x+1, y  ) == allow_grow_index) {p.set(x+1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_x_y(x-1, y  ) == allow_grow_index) {p.set(x-1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_x_y(x  , y+1) == allow_grow_index) {p.set(x  , y+1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_x_y(x  , y-1) == allow_grow_index) {p.set(x  , y-1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					}
				}
			old_points=new_points;
			}
		}


	// in field f the allowed growth area is defined with allow_grow_index - periodic x
	void Grow_object_area_set_by_allow_grow_index_periodic_x(double allow_grow_index, double lower_obj_index_limit, double upper_obj_index_limit)
		{
		long il;
		int x,y;
		double temp_obj;

		vector <Point> old_points;
		vector <Point> new_points;
		Point p;

		// identify starting points
		for (y=0;y<dimy;y++)
		for (x=0;x<dimx;x++)
			if (get(x,y) == allow_grow_index)
				if (Check_if_any_neigboring_point_has_a_value_in_interval_periodic_x(x,y,lower_obj_index_limit,upper_obj_index_limit,temp_obj))
					{
					p.set(x,y);
					p.lat=temp_obj;
					old_points.push_back(p.deep_copy());
					//cout << p.output() << endl;
					}


		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;
			new_points.clear();
			for (il=0; il < (long)old_points.size(); il++)
				{
				//cout << il << endl;
				x=old_points[il].x;
				y=old_points[il].y;
				temp_obj=old_points[il].lat;


				//check if the value was allready set by one of the previous objects - so we dont repeat ourselves
				if (get(x,y) == allow_grow_index)
					{
					// set the value of old_point
					set(x,y,temp_obj);

					// check neigboring points
					if (get_periodic_x_y(x+1, y  ) == allow_grow_index) {p.set(x+1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_x_y(x-1, y  ) == allow_grow_index) {p.set(x-1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}

					if (y < dimy - 1) if (get(x  , y+1) == allow_grow_index) {p.set(x  , y+1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (y > 0       ) if (get(x  , y-1) == allow_grow_index) {p.set(x  , y-1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					}
				}
			old_points=new_points;
			}
		}

	// in field f the allowed growth area is defined with allow_grow_index - periodic x and y
	void Grow_object_area_set_by_allow_grow_index_periodic_sphere(double allow_grow_index, double lower_obj_index_limit, double upper_obj_index_limit)
		{
		long il;
		int x,y;
		double temp_obj;

		vector <Point> old_points;
		vector <Point> new_points;
		Point p;

		// identify starting points
		for (y=0;y<dimy;y++)
		for (x=0;x<dimx;x++)
			if (get_periodic_sphere(x,y) == allow_grow_index)
				if (Check_if_any_neigboring_point_has_a_value_in_interval_periodic_sphere(x,y,lower_obj_index_limit,upper_obj_index_limit,temp_obj))
					{
					p.set(x,y);
					p.lat=temp_obj;
					old_points.push_back(p.deep_copy());
					//cout << p.output() << endl;
					}


		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;
			new_points.clear();
			for (il=0; il < (long)old_points.size(); il++)
				{
				//cout << il << endl;
				x=old_points[il].x;
				y=old_points[il].y;
				temp_obj=old_points[il].lat;


				//check if the value was allready set by one of the previous objects - so we dont repeat ourselves
				if (get_periodic_sphere(x,y) == allow_grow_index)
					{
					// set the value of old_point
					set_periodic_sphere(x,y,temp_obj);

					// check neigboring points
					if (get_periodic_sphere(x+1, y  ) == allow_grow_index) {p.set(x+1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_sphere(x-1, y  ) == allow_grow_index) {p.set(x-1, y  ); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_sphere(x  , y+1) == allow_grow_index) {p.set(x  , y+1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					if (get_periodic_sphere(x  , y-1) == allow_grow_index) {p.set(x  , y-1); p.lat=temp_obj; new_points.push_back(p.deep_copy());}
					}
				}

			// reindex new points so that they fall into the original domain
			for (il=0; il < (long)new_points.size(); il++)
				{
				if (new_points[il].x < 0)
					new_points[il].x=dimx-(-(long)new_points[il].x)%dimx;
				else
					new_points[il].x=(long)new_points[il].x%dimx;

				if (new_points[il].y == -1)
					new_points[il].y=0;
				if (new_points[il].y == dimy)
					new_points[il].y=dimy-1;
				}

			old_points=new_points;
			}
		}


 // Grows the area of objects to include touching areas of less-than-threshold preciptiation
	// The algorithm grows objects sequentially for 1 grid point a sequence
	CField Grow_object_area(CField &f_crop, double grow_threshold) const
		{
		CField f_grow;
		f_grow=create_a_copy_with_data();
		long il;
		int ix,iy;

		double temp_obj;

		vector <Point> old_points;
		vector <Point> new_points;
		Point p;

		//double x, y, lat, lon;


		// identify starting points
		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			// object not defined but not missing and some preciptiation does occur
			if (get(ix,iy) == -1 && f_crop.get(ix,iy) > grow_threshold)
					{
					// check neighborung grid points for objects
					temp_obj=-1;
					if (check_coordianes(ix-1,iy)) {if (get(ix-1,iy) > -1) temp_obj=get(ix-1,iy);}
					if (check_coordianes(ix+1,iy)) {if (get(ix+1,iy) > -1) temp_obj=get(ix+1,iy);}
					if (check_coordianes(ix,iy-1)) {if (get(ix,iy-1) > -1) temp_obj=get(ix,iy-1);}
					if (check_coordianes(ix,iy+1)) {if (get(ix,iy+1) > -1) temp_obj=get(ix,iy+1);}

					// add to the old_points
					if (temp_obj > -1)
						old_points.push_back(p.create((double)ix,(double)iy,temp_obj,0));
					}
			}

		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;

			new_points.clear();
			for (il=0; il < (long)old_points.size(); il++)
				{
				ix=(int)old_points[il].x;
				iy=(int)old_points[il].y;
				temp_obj=old_points[il].lat;

				if (!check_coordianes(ix,iy))
					error(AT,FU, " a problem!");


				//check if the value was allready set by one of the previous objects
				if (f_grow.get(ix,iy) == -1)
					{
					// set the value of old_point
					f_grow.set(ix,iy,temp_obj);

					// check neigboring points
					if (check_coordianes(ix-1,iy)) if (f_grow.get(ix-1,iy) == -1 && f_crop.get(ix-1,iy) > grow_threshold) new_points.push_back(p.create((double)ix-1,(double)iy,temp_obj,0));
					if (check_coordianes(ix+1,iy)) if (f_grow.get(ix+1,iy) == -1 && f_crop.get(ix+1,iy) > grow_threshold) new_points.push_back(p.create((double)ix+1,(double)iy,temp_obj,0));
					if (check_coordianes(ix,iy-1)) if (f_grow.get(ix,iy-1) == -1 && f_crop.get(ix,iy-1) > grow_threshold) new_points.push_back(p.create((double)ix,(double)iy-1,temp_obj,0));
					if (check_coordianes(ix,iy+1)) if (f_grow.get(ix,iy+1) == -1 && f_crop.get(ix,iy+1) > grow_threshold) new_points.push_back(p.create((double)ix,(double)iy+1,temp_obj,0));
					}
				}
			old_points.clear();
			old_points=new_points;
			}

		//for (il=0; il < old_points.size(); il++)
		//	cout << il << "\t" << old_points[il].x  << "\t" << old_points[il].y  << "\t" << old_points[il].lat << endl;

		return(f_grow);

		}


	// Grows the area of objects to include touching areas of less-than-threshold preciptiation
	// The algorithm grows objects sequentially for 1 grid point a sequence
	// Grows until distance from the original object is not greater than some threshold distance - parameter max_distance
	CField Grow_object_area_until_distance(CField &f_crop, double grow_threshold, double max_distance) const
		{
		CField f_grow;
		f_grow=create_a_copy_with_data();

		if (max_distance < 0 || max_distance > 99999 )
			error(AT,FU, " max_distance < 0 || max_distance > 99999");

		// unlimited growth if max_distance = 99999
		if (max_distance == 99999)
			{
			return(Grow_object_area(f_crop, grow_threshold));
			}


		long il,ip;
		int ix,iy,x2,y2;
		double min_distance, temp_distance;

		double temp_obj;

		// -----------------------------
		// Get the list of object border points
		// -----------------------------
		// list of border point for each object - first index is the obj_id, for each obj_id there is a vector list of border offsets
		vector < vector<long> > object_border_points(get_max_value()+1,vector<long>(0));
		// identify objects border
		for (il=0;il<size();il++)
			if (data[il] > -1)
				{
				DataLoc2XY(il, ix,iy);
				// check for border
				if (get_but_check_coordianes(ix-1,iy) != data[il]
					|| get_but_check_coordianes(ix+1,iy) != data[il]
					|| get_but_check_coordianes(ix,iy-1) != data[il]
					|| get_but_check_coordianes(ix,iy+1) != data[il])
					object_border_points[(long)data[il]].push_back(il);
				}



		vector <Point> old_points;
		vector <Point> temp_points;
		vector <Point> new_points;
		Point p;


		// -----------------------------
		// identify starting points
		// -----------------------------
		for (iy=0;iy<dimy;iy++)
		for (ix=0;ix<dimx;ix++)
			{
			// object not defined but not missing and some preciptiation does occur
			if (get(ix,iy) == -1 && f_crop.get(ix,iy) > grow_threshold)
					{
					// check neighborung grid points for objects
					temp_obj=-1;
					if (check_coordianes(ix-1,iy)) {if (get(ix-1,iy) > -1) temp_obj=get(ix-1,iy);}
					if (check_coordianes(ix+1,iy)) {if (get(ix+1,iy) > -1) temp_obj=get(ix+1,iy);}
					if (check_coordianes(ix,iy-1)) {if (get(ix,iy-1) > -1) temp_obj=get(ix,iy-1);}
					if (check_coordianes(ix,iy+1)) {if (get(ix,iy+1) > -1) temp_obj=get(ix,iy+1);}

					// add to the old_points
					if (temp_obj > -1)
						old_points.push_back(p.create((double)ix,(double)iy,temp_obj,0));
					}
			}


		// -----------------------------
		// iterate points
		// -----------------------------
		while (old_points.size() != 0)
			{
			//cout << old_points.size() << endl;

			// remove old_points that are more than a distance away from their non-grown object
			temp_points.clear();
			for (ip=0; ip < (long)old_points.size(); ip++)
				{
				ix=(int)old_points[ip].x;
				iy=(int)old_points[ip].y;
				temp_obj=old_points[ip].lat;

				//cout << ip << "\t" << ix  << "\t" << iy  << "\t" << temp_obj << endl;

				// find minimal distance to non-grown object -> check the list of all border points for this object
				min_distance=1E200;
				// loop over all border points for this object
				for (il=0; il < (long)object_border_points[(long)temp_obj].size(); il++)
					{
					DataLoc2XY(object_border_points[(long)temp_obj][il], x2,y2);
					temp_distance=calculate_distance_on_earth_in_km_between_two_grid_points(ix, iy, x2, y2);
					//cout << ip << "\t" << temp_distance << endl;
					if ( temp_distance < min_distance ) min_distance=temp_distance;
					}
				if (min_distance <= max_distance)
					temp_points.push_back(old_points[ip]);
				//cout << ip << "\t" << min_distance << endl;
				}

			//exit(1);
			//cout << old_points.size() <<  "\t" << temp_points.size() << endl;
			old_points=temp_points;

			new_points.clear();
			for (ip=0; ip < (long)old_points.size(); ip++)
				{
				ix=(int)old_points[ip].x;
				iy=(int)old_points[ip].y;
				temp_obj=old_points[ip].lat;

				if (!check_coordianes(ix,iy))
					error(AT,FU, " a problem!");


				//check if the value was allready set by one of the previous objects
				if (f_grow.get(ix,iy) == -1)
					{
					// set the value of old_point
					f_grow.set(ix,iy,temp_obj);

					// check neigboring points
					if (check_coordianes(ix-1,iy)) if (f_grow.get(ix-1,iy) == -1 && f_crop.get(ix-1,iy) > grow_threshold) new_points.push_back(p.create((double)ix-1,(double)iy,temp_obj,0));
					if (check_coordianes(ix+1,iy)) if (f_grow.get(ix+1,iy) == -1 && f_crop.get(ix+1,iy) > grow_threshold) new_points.push_back(p.create((double)ix+1,(double)iy,temp_obj,0));
					if (check_coordianes(ix,iy-1)) if (f_grow.get(ix,iy-1) == -1 && f_crop.get(ix,iy-1) > grow_threshold) new_points.push_back(p.create((double)ix,(double)iy-1,temp_obj,0));
					if (check_coordianes(ix,iy+1)) if (f_grow.get(ix,iy+1) == -1 && f_crop.get(ix,iy+1) > grow_threshold) new_points.push_back(p.create((double)ix,(double)iy+1,temp_obj,0));
					}
				}
			old_points.clear();
			old_points=new_points;
			}

		//for (il=0; il < old_points.size(); il++)
		//	cout << il << "\t" << old_points[il].x  << "\t" << old_points[il].y  << "\t" << old_points[il].lat << endl;

		return(f_grow);

		}





	//void find_all_objects_around_location(int x,int y, double radius_in_deg, vector <double> &list)
	double find_closest_object_around_location(int x,int y, double radius_in_deg, double &dmin) const
		{
		int ix,iy,io;

		// this algorithm so far supports only equaly spaced grids
		for (ix=0;ix<dimx-2;ix++)
			if (lon_data[ix]-lon_data[ix+1] != lon_data[ix+1]-lon_data[ix+2])
				error(AT,FU, "this algorithm so far supports only equaly spaced grids");
		for (ix=0;ix<dimy-2;ix++)
			if (lat_data[ix]-lat_data[ix+1] != lat_data[ix+1]-lat_data[ix+2])
				error(AT,FU, "this algorithm so far supports only equaly spaced grids");

		int iR;
		iR=(int)round(radius_in_deg/fabs(lon_data[0]-lon_data[1]));

		double object_index;
		double d;


		vector <double> list;
		list.clear();
		bool found;

		dmin=1E99;
		object_index=-1;
		for (iy=y-iR;iy<= y+iR;iy++)
		for (ix=x-iR;ix<= x+iR;ix++)
		// inside domain
		if (check_coordianes(ix, iy))
		// inside circle
		if ( (double)(ix-x)*(double)(ix-x)+(double)(iy-y)*(double)(iy-y) <= (double)iR*(double)iR)
		// found object
		if (get(ix,iy) > -1)
			{
			d=(double)(ix-x)*(double)(ix-x)+(double)(iy-y)*(double)(iy-y);
			if (d < dmin)
				{
				dmin=d;
				object_index=get(ix,iy);
				}



			found=false;
			// check list
			for (io=0; io < (long)list.size(); io++)
				if (get(ix,iy) == list[io])
					found=true;
			if (found==false)
				list.push_back(get(ix,iy));
			}

		dmin=sqrt(dmin);


		/*cout << "************** " << object_index << "\t" << dmin << endl;
		for (ix=0;ix< list.size();ix++)
			cout << list[ix] << endl;
		cout << endl;
		*/
		return(object_index);

		//exit(1);


		}


	//void find_all_objects_around_location(int x,int y, double radius_in_deg, vector <double> &list)
	vector <double> find_objects_around_location(int x,int y, double radius_in_deg, vector <double> &distance) const
		{
		int ix,iy;
		long il,ip;
		double d;


		// this algorithm so far supports only equaly spaced grids
		for (ix=0;ix<dimx-2;ix++)
			if (lon_data[ix]-lon_data[ix+1] != lon_data[ix+1]-lon_data[ix+2])
				error(AT,FU, "this algorithm so far supports only equaly spaced grids");
		for (ix=0;ix<dimy-2;ix++)
			if (lat_data[ix]-lat_data[ix+1] != lat_data[ix+1]-lat_data[ix+2])
				error(AT,FU, "this algorithm so far supports only equaly spaced grids");

		// -----------------------------
		// Get the list of object border points
		// -----------------------------
		vector <long> object_border_points;
		for (il=0;il<size();il++)
			if (data[il] > -1)
				{
				DataLoc2XY(il, ix,iy);
				// check for border
				if (get_but_check_coordianes(ix-1,iy) != data[il]
					|| get_but_check_coordianes(ix+1,iy) != data[il]
					|| get_but_check_coordianes(ix,iy-1) != data[il]
					|| get_but_check_coordianes(ix,iy+1) != data[il])
					object_border_points.push_back(il);
				}

		vector <double> object_min_distance((long)get_max_value()+1,1E100);
		for (ip=0;ip< (long)object_border_points.size();ip++)
			{
			DataLoc2XY(object_border_points[ip], ix,iy);

			//cout << object_border_points[ip] <<" " <<ix << " " << iy << endl;
			d=(double)(ix-x)*(double)(ix-x)+(double)(iy-y)*(double)(iy-y);
			if (d < object_min_distance[data[object_border_points[ip]]])
				object_min_distance[data[object_border_points[ip]]]=d;
			}

			//exit(1);


		// a special case when location is allready inside objects - just checking the border would give an uncorrect too large result for distance
		if (get(x,y) > -1)
			object_min_distance[(long)get(x,y)] = 0;

		//for (ip=0;ip< (long)object_min_distance.size();ip++)
		//	if (object_min_distance[ip] < 1E100)
		//		cout << ip << "\t" << sqrt(object_min_distance[ip]) << endl;

		// make a list of existing objects only
		vector <double> object_list;
		vector <double> min_distance_to_object;
		for (ip=0;ip< (long)object_min_distance.size();ip++)
			if (object_min_distance[ip] < 1E100)
				{
				object_list.push_back(ip);
				min_distance_to_object.push_back(sqrt(object_min_distance[ip]));
				}
		/*// convert distances to degrees
		for (ip=0;ip< (long)object_list.size();ip++)
			min_distance_to_object[ip] = min_distance_to_object[ip]/4.0;*/

		// sort by object distance
		Clist_of_values_to_sort list_of_values_to_sort;
		list_of_values_to_sort.sort_two_vectors_of_doubles_by_first_index(min_distance_to_object,object_list);

		//for (ip=0;ip< (long)object_list.size();ip++)
		//	cout << ip << "\t" << object_list[ip] << "\t" << min_distance_to_object[ip] << endl;

		// remove objects further than the threshold
		for (ip=0;ip< (long)object_list.size();ip++)
			if (min_distance_to_object[ip]/4.0 > radius_in_deg)
				{
				min_distance_to_object.erase(min_distance_to_object.begin()+ip,min_distance_to_object.end());
				object_list.erase(object_list.begin()+ip,object_list.end());
				}
		//for (ip=0;ip< (long)object_list.size();ip++)
		//	cout << ip << "\t" << object_list[ip] << "\t" << min_distance_to_object[ip] << endl;

		distance=min_distance_to_object;
		return(object_list);
		}



	// Sets value to missing data if any of two fields has a missing data
	void Missing_data_overlay(CField &original_field)
		{
		long ixl;

		// check dimensions
		if (!check_compare_dimensions(original_field))
			error(AT,FU,"!check_compare_dimensions(original_field)");
		if (!original_field.check_if_initilized())
			error(AT,FU,"!original_field.check_if_initilized()");
		if (!check_if_initilized())
			error(AT,FU,"!check_if_initilized()");


		for (ixl=0;ixl<size();ixl++)
			if (data[ixl] == BAD_DATA_FLOAT || original_field.data[ixl] == BAD_DATA_FLOAT )
				{
				data[ixl] = BAD_DATA_FLOAT;
				original_field.data[ixl] = BAD_DATA_FLOAT;
				}
		}

	// creates a o copy scaled by scale, only generate lon_data and lat_data, does not alocatememory or copy data
	CField create_a_scaled_copy_without_data(double scale) const
		{
		CField f;
		f.dimx=(int)((double)dimx*scale);
		f.dimy=(int)((double)dimy*scale);
		f.unique_timestamp=unique_timestamp;

		int ix;

		if (lat_data.size()==0 || lon_data.size()==0)
			{cout << "ERROR: CField->create_a_scaled_copy_without_data ERROR: lat_data.size()==0 || lon_data.size()==0 ""\n";exit(1);}

	    f.allocatememory_lat_lon_array();

		double double_oldix;

	    // scale lon data
	    for (ix=0; ix < f.dimx; ix++)
	    	{
	      	double_oldix=(double)ix*(double)(dimx-1)/(double)(f.dimx-1);

			f.lon_data[ix] = lon_data[(int)floor(double_oldix)] +
				(lon_data[(int)ceil(double_oldix)] - lon_data[(int)floor(double_oldix)])*
				(double_oldix-floor(double_oldix));

			}


	    for (ix=0; ix < f.dimy; ix++)
	    	{
	      	double_oldix=(double)ix*(double)(dimy-1)/(double)(f.dimy-1);

			f.lat_data[ix] = lat_data[(int)floor(double_oldix)] +
				(lat_data[(int)ceil(double_oldix)] - lat_data[(int)floor(double_oldix)])*
				(double_oldix-floor(double_oldix));

			}

		return(f);
		}


};
