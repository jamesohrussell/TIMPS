

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int cn_PnPoly_vec( double x, double y,  const vector <double> &x_pg,  const vector <double> &y_pg)
	{
    int    cn = 0;    // the crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i< (long)x_pg.size()-1; i++) {    // edge from V[i] to V[i+1]
       if (((y_pg[i] <= y) && (y_pg[i+1] > y))    // an upward crossing
        || ((y_pg[i] > y) && (y_pg[i+1] <= y))) { // a downward crossing
            // compute the actual edge-ray intersect x-coordinate
            float vt = (float)(y - y_pg[i]) / (y_pg[i+1] - y_pg[i]);
            if (x < x_pg[i] + vt * (x_pg[i+1] - x_pg[i])) // x < intersect
                ++cn;   // a valid crossing of y=y right of x
        }
    }
    return (cn&1);    // 0 if even (out), and 1 if odd (in)
	}

//------------------------------------------------------------------------
// class for polygon defined domain
class CPolygon_domain
	{
	public:
	string name;
	vector <double> lat_pg;
	vector <double> lon_pg;
	int land_sea_mask_switch;
	int color_seed;

	CPolygon_domain() {freememory();}
	~CPolygon_domain() {freememory();}

	void freememory()
		{
		name="";
		land_sea_mask_switch=-1;
		lon_pg.clear();
		lat_pg.clear();
		color_seed=-1;
		}

	int get_land_sea_mask_switch_from_string(string land_sea_mask_switch_str)  const
		{
		if (land_sea_mask_switch_str == "all") return(0);
		else if (land_sea_mask_switch_str == "sea") return(1);
		else if (land_sea_mask_switch_str == "land") return(2);
		else error(AT,FU, "get_land_sea_mask_switch_from_string");
		return(0);
		}

	string get_string_from_land_sea_mask_switch(int lsmw) const
		{
		if (land_sea_mask_switch == 0) return("all");
		else if (land_sea_mask_switch == 1) return("sea");
		else if (land_sea_mask_switch == 2) return("land");
		else error(AT,FU, "get_string_from_land_sea_mask_switch");
		return("error");
		}

	string output() const
		{
		long il;

		ostringstream s1;
		s1.str("");

		s1 << name << "\t" << get_string_from_land_sea_mask_switch(land_sea_mask_switch) << endl;
		for (il=0; il < (long)lat_pg.size(); il ++)
			s1 << il << "\t" << lon_pg[il]  << "\t" << lat_pg[il] << endl;
		return(s1.str());
		}

	bool check_if_initialized_and_consistent() const
		{
		if (lon_pg.size() < 4 || lon_pg.size() != lat_pg.size())
			return(false);
		return(true);
		}


	// this is needed to so that some grid poins do not fall into two neigboring polyogns
	// for example grid point lon=15.00000 would fall into the right and left polygon with
	// border at lon=15 deg. If 0.0001 deg is added to the lon value the grid point would
	// only fall into the right polygon
	void add_00001_deg_to_lat_and_lon_so_that_polygon_algorithm_works_better()
		{
		if (!check_if_initialized_and_consistent())
			error(AT,FU, "!check_if_initialized_and_consistent()");

		long il;
		for (il=0; il < (long)lon_pg.size(); il ++)
			{
			lon_pg[il]+=0.0001;
			lat_pg[il]+=0.0001;
			}
		}


	void initilize_rectangular_domain(string nameX, double lon_min, double lon_max, double lat_min, double lat_max, string land_sea_mask_switch_str, int color_seedX)
		{
		if (lon_min < 0) lon_min+=360;
		if (lon_max < 0) lon_max+=360;

		freememory();
		name=nameX;
		land_sea_mask_switch=get_land_sea_mask_switch_from_string(land_sea_mask_switch_str);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_min);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_max);
		lon_pg.push_back(lon_max);
		lat_pg.push_back(lat_max);
		lon_pg.push_back(lon_max);
		lat_pg.push_back(lat_min);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_min);
		add_00001_deg_to_lat_and_lon_so_that_polygon_algorithm_works_better();
		color_seed=color_seedX;
		}

	void initilize_rectangular_domain_dont_correct_negative_longitudes(string nameX, double lon_min, double lon_max, double lat_min, double lat_max, string land_sea_mask_switch_str, int color_seedX)
		{
		freememory();
		name=nameX;
		land_sea_mask_switch=get_land_sea_mask_switch_from_string(land_sea_mask_switch_str);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_min);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_max);
		lon_pg.push_back(lon_max);
		lat_pg.push_back(lat_max);
		lon_pg.push_back(lon_max);
		lat_pg.push_back(lat_min);
		lon_pg.push_back(lon_min);
		lat_pg.push_back(lat_min);
		add_00001_deg_to_lat_and_lon_so_that_polygon_algorithm_works_better();
		color_seed=color_seedX;
		}


	void initilize_domain_with_polygon(string nameX,  const vector <double> &lon_vec,  const vector <double> &lat_vec,  string land_sea_mask_switch_str, int color_seedX)
		{
		freememory();
		name=nameX;
		land_sea_mask_switch=get_land_sea_mask_switch_from_string(land_sea_mask_switch_str);
		if (lon_vec.size() < 4 || lon_vec.size() != lat_vec.size())
			error(AT,FU, "lon_vec.size() < 4 || lon_vec.size() != lat_vec.size()");

		long il;
		for (il=0; il < (long)lat_vec.size(); il ++)
			{

			if (lon_vec[il] < 0)
				lon_pg.push_back(lon_vec[il]+360);
			else
				lon_pg.push_back(lon_vec[il]);
			lat_pg.push_back(lat_vec[il]);
			}
		add_00001_deg_to_lat_and_lon_so_that_polygon_algorithm_works_better();
		color_seed=color_seedX;
		}

	void initilize_domain_with_polygon_dont_correct_negative_longitudes(string nameX,  const vector <double> &lon_vec,  const vector <double> &lat_vec,  string land_sea_mask_switch_str, int color_seedX)
		{
		freememory();
		name=nameX;
		land_sea_mask_switch=get_land_sea_mask_switch_from_string(land_sea_mask_switch_str);
		if (lon_vec.size() < 4 || lon_vec.size() != lat_vec.size())
			error(AT,FU, "lon_vec.size() < 4 || lon_vec.size() != lat_vec.size()");

		long il;
		for (il=0; il < (long)lat_vec.size(); il ++)
			{
			lon_pg.push_back(lon_vec[il]);
			lat_pg.push_back(lat_vec[il]);
			}
		add_00001_deg_to_lat_and_lon_so_that_polygon_algorithm_works_better();
		color_seed=color_seedX;
		}


	bool is_inside_domain_f_landfrac_on_same_grid(int ix, int iy,  const CField &f_landfrac) const
		{
		double landfrac;
		bool landfrac_test_ok=false;
		double lon,lat;

		landfrac=f_landfrac.get(ix,iy);
		if (land_sea_mask_switch==0)
			landfrac_test_ok=true;
		else if (land_sea_mask_switch==1)
			{
			if (landfrac <= 0.5)
				landfrac_test_ok=true;
			}
		else if (land_sea_mask_switch==2)
			{
			if (landfrac > 0.5)
				landfrac_test_ok=true;
			}

		// break if landfrac_test failed
		if (landfrac_test_ok==false) return (false);

		f_landfrac.XY2LatLon(ix,iy,lon,lat);
		// otherwise test also if inside poligon
		if(cn_PnPoly_vec( lon, lat, lon_pg, lat_pg) == 1)
			return(true);

		return(false);
		}

	bool is_inside_domain_f_landfrac_different_grid(double lon, double lat,  const CField &f_landfrac) const
		{
		double landfrac;
		bool landfrac_test_ok=false;

		landfrac=f_landfrac.get_by_lon_and_lat(lon, lat);
		if (land_sea_mask_switch==0)
			landfrac_test_ok=true;
		else if (land_sea_mask_switch==1)
			{
			if (landfrac <= 0.5)
				landfrac_test_ok=true;
			}
		else if (land_sea_mask_switch==2)
			{
			if (landfrac > 0.5)
				landfrac_test_ok=true;
			}

		// break if landfrac_test failed
		if (landfrac_test_ok==false) return (false);

		// otherwise test also if inside poligon
		if(cn_PnPoly_vec( lon, lat, lon_pg, lat_pg) == 1)
			return(true);

		return(false);

		}

	};


//------------------------------------------------------------------------
// class for longitudal average of domain using a land mask
class CDomain
	{
	public:
	string name;
	double lon_min;
	double lon_max;
	double lat_min;
	double lat_max;
	string land_sea_mask_swithch;

	vector <double> lat_out;
	vector <double> value_out;
	vector <double> counter;

	CDomain()
		{
		freememory();
		}

	~CDomain()
		{
		freememory();
		}

	void freememory()
		{
		name="";
		lon_min=BAD_DATA_FLOAT;
		lat_max=BAD_DATA_FLOAT;
		lon_max=BAD_DATA_FLOAT;
		lat_min=BAD_DATA_FLOAT;
		land_sea_mask_swithch="";
		lat_out.clear();
		value_out.clear();
		counter.clear();
		}

	string output()
		{
		if (!check_if_initialized())
			error(AT,FU, "!check_if_initilized()");

		ostringstream s1;
		s1.str("");
		long iy;
		for (iy=0; iy < (long)lat_out.size(); iy++)
			s1 << lat_out[iy] << "\t" << value_out[iy] << "\t" << counter[iy] << endl;
		return(s1.str());
		}


	bool check_if_initialized() const
		{
		if (lon_min == BAD_DATA_FLOAT || lon_max == BAD_DATA_FLOAT || lat_max == BAD_DATA_FLOAT || lat_min == BAD_DATA_FLOAT)
		 return(false);

		return(true);
		}


	CDomain initialize(string nameX, double lon_minX, double lon_maxX, double lat_minX, double lat_maxX, string land_sea_mask_swithchX)
		{
		CDomain cd;
		cd.name=nameX;

		// if domain specified with - longitudes
		if (lon_minX < 0) lon_minX=360+lon_minX;
		if (lon_maxX < 0) lon_maxX=360+lon_maxX;

		cd.lon_min=lon_minX;
		cd.lon_max=lon_maxX;
		cd.lat_min=lat_minX;
		cd.lat_max=lat_maxX;
		cd.land_sea_mask_swithch=land_sea_mask_swithchX;
		return(cd);
		}


	void calculate_longitude_average( const CField &f,  const CField &f_land_sea_mask)
		{
		int ix,iy;
		int x,y;
		double lat,lon;
		//double ave;

		lat_out.clear();
		value_out.clear();
		counter.clear();

		for (iy=0; iy < f.dimy; iy++)
			{
			lat_out.push_back(f.lat_data[iy]);
			value_out.push_back(0);
			counter.push_back(0);
			}

		int ok_switch;
		for (iy=0; iy < f.dimy; iy++)
			{
			//ave=0;
			counter[iy]=0;
			for (ix=0; ix < f.dimx; ix++)
				{
				f.XY2LatLon(ix, iy, lon, lat);

				ok_switch=0;
				// domain crosses longitue 0
				if (lon_min > lon_max)
					if (lon >= lon_min || lon <= lon_max)
						ok_switch=1;

				// domain does not cross longitue 0
				if  (lon_min < lon_max)
					if (lon >= lon_min && lon <= lon_max)
						ok_switch=1;

				if (ok_switch==1)
					{

					// check only sea
					if (land_sea_mask_swithch=="sea")
						{
						f_land_sea_mask.LatLon2XY_new(lon, lat, x, y);
						if (f_land_sea_mask.get_but_check_coordianes(x, y) == 0)
							// check missing data
							if (f.get_but_check_coordianes(ix,iy) != BAD_DATA_FLOAT)
								{
								value_out[iy]+=f.get_but_check_coordianes(ix,iy);
								counter[iy]++;
								}
						}
					// check only land
					if (land_sea_mask_swithch=="land")
						{
						f_land_sea_mask.LatLon2XY_new(lon, lat, x, y);
						if (f_land_sea_mask.get_but_check_coordianes(x, y) == 1)
							// check missing data
							if (f.get_but_check_coordianes(ix,iy) != BAD_DATA_FLOAT)
							{
							value_out[iy]+=f.get_but_check_coordianes(ix,iy);
							counter[iy]++;
							}
						}
					// all
					if (land_sea_mask_swithch=="all")
						{
						// check missing data
						if (f.get_but_check_coordianes(ix,iy) != BAD_DATA_FLOAT)
							{
							value_out[iy]+=f.get_but_check_coordianes(ix,iy);
							counter[iy]++;
							}
						}
					}
				}
			if (counter[iy]!=0)
				value_out[iy]=value_out[iy]/counter[iy];
			else value_out[iy]=BAD_DATA_FLOAT;
			}

		vector <double> l;
		vector <double> v;
		vector <double> c;

		l=lat_out;
		v=value_out;
		c=counter;

		lat_out.clear();
		value_out.clear();
		counter.clear();


		// delete lat out of domain
		for (iy=0; iy < f.dimy; iy++)
			{
			if (lat_max >= l.at(iy)  && lat_min <= l.at(iy) )
				{
				lat_out.push_back(l.at(iy));
				value_out.push_back(v.at(iy));
				counter.push_back(c.at(iy));
				}
			}


		//cout << name <<  endl;
		for (iy=0; iy < (long)lat_out.size(); iy++)
			{
			//cout << iy << "\t" << lat_out[iy] << "\t"  <<  value_out[iy]  << "\t"  << counter[iy] << endl;
			}
		}


	void calculate_longitude_average_and_write_to_file( const CField &f,  const CField &f_land_sea_mask, string filename)
		{
		ofstream myfile;

//		calculate_longitude_average(f, f_land_sea_mask, lat_out, value_out, counter);
//
//		myfile.open (filename.c_str());
//
//		for (ix=0;ix<lat_out.size();ix++)
//			{
//			myfile << lat_out << "\t";
//			myfile << value_out << "\t";
//			myfile << counter << "\t";
//			myfile << endl;
//			}
//		myfile.close();
		}

	};




