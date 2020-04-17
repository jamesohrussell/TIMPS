

void draw_a_cross_to_GDimage(gdImagePtr image, int x,int y, int radius, int color)
	{
	gdImageLine( image , x-radius , y-radius , x+radius ,y+radius , color );
	gdImageLine( image , x-radius , y+radius , x+radius ,y-radius , color );
	}


class CRGBcolor {
	public:
	double r,g,b;

	// constructor
	CRGBcolor() {freememory();}
	~CRGBcolor() {freememory();}

	CRGBcolor(double r_,double g_, double b_)
		{
		set( r_, g_,  b_);
		}

	CRGBcolor create(double r_,double g_, double b_)
		{
		CRGBcolor out;
		out.set( r_, g_,  b_);
		return(out);
		}

	CRGBcolor create_by_int(int r_,int g_, int b_)
		{
		CRGBcolor out;
		out.set( (double)r_/255.0, (double)g_/255.0,  (double)b_/255.0);
		return(out);
		}

	CRGBcolor create_greyscale_by_int(int x_)
		{
		CRGBcolor out;
		out.set( (double)x_/255.0, (double)x_/255.0,  (double)x_/255.0);
		return(out);
		}


	// frees memory
	void freememory()
		{
		r=b=g=-1;
		}

	bool compare_are_they_the_same(const CRGBcolor &f) const
		{
		if (r != f.r) return(false);
		if (g != f.g) return(false);
		if (b != f.b) return(false);
		return(true);
		}

	// overload operators == nad !=
	bool operator== (const CRGBcolor& f) const
		{return(compare_are_they_the_same(f));}
	bool operator!= (const CRGBcolor& f) const
		{return(!compare_are_they_the_same(f));}

	CRGBcolor& operator= (const CRGBcolor &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			r=rhs.r;
			g=rhs.g;
			b=rhs.b;
			}
		return *this;    // Return ref for multiple assignment
		}

	void set_to_grayscale(double level)
		{
		r=g=b=level;
		}

	void set(double r_,double g_, double b_)
		{
		r=r_;
		g=g_;
		b=b_;
		}

	// distance is between 0 (color 1) and 1 (color 2)
	void set_RGB_color_betwwen_two_colors(CRGBcolor const &c1, CRGBcolor const &c2, double distance)
		{
		r=c1.r*(1.0-distance)+c2.r*distance;
		g=c1.g*(1.0-distance)+c2.g*distance;
		b=c1.b*(1.0-distance)+c2.b*distance;
		}

	// seed defines the color
	void set_random_RGB_color_by_seed(double seed)
		{
		long idum;
		idum=(long)(-seed);

		r=ran2(&idum);
		g=ran2(&idum);
		b=ran2(&idum);
		}

	// seed defines the color, min and max brightness should be between 0 and 1
	void set_random_RGB_color_by_seed_and_limit_brightness (double seed, double min_brightness, double max_brightness)
		{
		//long idum;

		//idum=(long)(-seed);

		set_random_RGB_color_by_seed(seed);
		// eliminate darkest colors - since they are too similar to the black background
		while (r+g+b < min_brightness*3 || r+g+b > max_brightness*3)
			set_random_RGB_color_by_seed(seed);
		}

	void set_RGB_color_from_colortable_continuous_colors(double value, vector <double> &levels, vector <CRGBcolor> &colors)
		{
		double distance;
		vector<double>::iterator low;

		if (levels.size() != colors.size())
			error(AT,FU, "levels.size() != colors.size()");

		if (value <= levels.front())
			*this=colors.front();
		else if (value >= levels.back())
			*this=colors.back();
		else
			{
			low=lower_bound (levels.begin(), levels.end(), value);
			distance=(*low-value)/(*low-*(low-1));
			set_RGB_color_betwwen_two_colors(colors[low- levels.begin()], colors[low- levels.begin()-1],  distance);
			}

		}

	void set_RGB_color_from_colortable_discrete_colors(double value, vector <double> &levels, vector <CRGBcolor> &colors)
		{
		if (levels.size() != colors.size() - 1)
			error(AT,FU, "levels.size() != colors.size() - 1");

		vector<double>::iterator low;

		if (value <= levels.front())
			*this=colors.front();
		else if (value > levels.back())
			*this=colors.back();
		else
			{
			low=lower_bound (levels.begin(), levels.end(), value);
			*this=colors[low- levels.begin()];
			}
		}


	void set_RGB_color_from_colortable_discrete_colors_nonmatching_number_of_colors_and_levels(double value, const vector <double> &levels, const vector <CRGBcolor> &colors)
		{
		vector<double>::iterator low;

		//cout << levels.size() << "\t" << colors.size() << endl;
		ERRORIF(levels.size()+1 > colors.size());

		long index=lower_bound (levels.begin(), levels.end(), value)- levels.begin();
		long color_index=floor((double)(index)*(double)(colors.size()-1)/((double)levels.size()));

		//cout << value << "\t" << index << "\t" << color_index << endl;

		*this=colors[color_index];
		}

	void set_RGB_color_from_colortable_using_two_levels_and_linear_scale(double value, const vector <double> &levels, const vector <CRGBcolor> &colors)
		{
		vector<double>::iterator low;

		//cout << levels.size() << "\t" << colors.size() << endl;
		ERRORIF(levels.size() != 2);
		ERRORIF(colors.size() < 1);

		double x=(value-levels[0])/(levels[1]-levels[0]);
		if (x < 0)
			*this=colors.front();
		else if (x >= 1)
			*this=colors.back();
		else
			*this=colors[floor( x * (double)colors.size())];
		}


	void read_colortable_from_NCL_colortable_file(string fname, vector <CRGBcolor> &colors)
		{
		CRGBcolor rgb1;

		colors.clear();
		vector <string> file_line = get_each_line_as_string_from_txt_file(fname);

		for (long il = 0; il < (long)file_line.size(); il++)
			{
			string line=file_line[il];

			// remove unnecesarry whitespace
			Remove_repeating_whitespace_characters_from_string(line);

			// remove all charaters after #
			std::size_t found = line.find("#");
			if (found!=std::string::npos)
				line.erase(found,string::npos);

			if (line.size() > 0)
				{
				// remove space if it is the first character
				if (line[0] == ' ')
					line.erase(0,1);

				// remove space if it is the last character
				if (line[line.size()-1] == ' ')
					line.erase(line.size()-1,string::npos);

				vector<string> fields=StringExplode(line, " ");
				if (fields.size() == 3)
				if (is_string_an_Integer_test(fields[0]) && is_string_an_Integer_test(fields[1]) && is_string_an_Integer_test(fields[2]))
					{
					colors.push_back(rgb1.create_by_int( atoi((fields[0]).c_str()), atoi((fields[1]).c_str()) , atoi((fields[2]).c_str()) ));
					}
				}
			}
		}

	void darken_or_lighten_the_color_by_multiplying(double factor)
		{
		r*=factor;
		g*=factor;
		b*=factor;
		if (r>1) r=1;
		if (g>1) g=1;
		if (b>1) b=1;
		}


	string output() const
		{
		ostringstream s1;
		s1.str("");
		s1 << r << "," << g << "," << b;
		return(s1.str());
		}

	};



void Create_marked_coast_CField_from_land_sea_fraction_same_grid_for_f_lsf_and_f_in_and_out( const CField &f_lsf,  CField &f_in_and_out)
	{
	int ix, iy;
	if (!f_lsf.check_if_initilized())
		error(AT,FU, "!f_lsm.check_if_initilized()");
	if (!f_in_and_out.check_if_initilized())
		error(AT,FU, "!f_in_and_out.check_if_initilized()");
	if (f_lsf.dimx != f_in_and_out.dimx || f_lsf.dimy != f_in_and_out.dimy)
		error(AT,FU, "f_lsf.dimx != f_in_and_out.dimx || f_lsf.dimy != f_in_and_out.dimy");
	for (ix=0; ix < f_lsf.dimx; ix++)
		if (f_lsf.lon_data[ix] != f_in_and_out.lon_data[ix])
			error(AT,FU, "f_lsf.lon_data[ix] != f_in_and_out.lon_data[ix]");
	for (ix=0; ix < f_lsf.dimy; ix++)
		if (f_lsf.lat_data[ix] != f_in_and_out.lat_data[ix])
			error(AT,FU, "f_lsf.lat_data[ix] != f_in_and_out.lat_data[ix]");

	f_in_and_out.Set_all_values(0);

	// identify coasts in lsf
	for (iy=1;iy<f_lsf.dimy;iy++)
	for (ix=1;ix<f_lsf.dimx;ix++)
		{
		// check if land sea changes in neighboring points in x and y directions
		if (f_lsf.get(ix,iy-1) > 0.5 && f_lsf.get(ix,iy) <= 0.5)
			f_in_and_out.set(ix,iy,1);
		if (f_lsf.get(ix,iy-1) < 0.5 && f_lsf.get(ix,iy) >= 0.5)
			f_in_and_out.set(ix,iy,1);
		if (f_lsf.get(ix-1,iy) > 0.5 && f_lsf.get(ix,iy) <= 0.5)
			f_in_and_out.set(ix,iy,1);
		if (f_lsf.get(ix-1,iy) < 0.5 && f_lsf.get(ix,iy) >= 0.5)
			f_in_and_out.set(ix,iy,1);
		}
	}


void Create_marked_coast_CField_from_land_sea_fraction( const CField &f_lsf,  CField &f_in_and_out)
	{

	if (!f_lsf.check_if_initilized())
		error(AT,FU, "!f_lsm.check_if_initilized()");
	if (!f_in_and_out.check_if_initilized())
		error(AT,FU, "!f_in_and_out.check_if_initilized()");

	int ix, iy;

	CField f_temp;
	f_temp=f_lsf.create_a_copy_with_data();
	f_temp.Set_all_values(0);

	// identify coasts in lsf
	for (iy=1;iy<f_lsf.dimy;iy++)
	for (ix=1;ix<f_lsf.dimx;ix++)
		{
		// check if land sea changes in neighboring points in x and y directions
		if (f_lsf.get(ix,iy-1) > 0.5 && f_lsf.get(ix,iy) <= 0.5)
			f_temp.set(ix,iy,1);
		if (f_lsf.get(ix,iy-1) < 0.5 && f_lsf.get(ix,iy) >= 0.5)
			f_temp.set(ix,iy,1);
		if (f_lsf.get(ix-1,iy) > 0.5 && f_lsf.get(ix,iy) <= 0.5)
			f_temp.set(ix,iy,1);
		if (f_lsf.get(ix-1,iy) < 0.5 && f_lsf.get(ix,iy) >= 0.5)
			f_temp.set(ix,iy,1);
		}

	f_in_and_out.Set_all_values(0);
	// map coasts into f_in_and_out
	for (iy=0;iy<f_lsf.dimy;iy++)
		for (ix=0;ix<f_lsf.dimx;ix++)
			// check that the ix and iy fall inside the f_in_and_out domain
			if (f_lsf.lon_data[ix] >= f_in_and_out.lon_data[0] && f_lsf.lon_data[ix] <= f_in_and_out.lon_data[f_in_and_out.dimx-1] && f_lsf.lat_data[iy] >= f_in_and_out.lat_data[0] && f_lsf.lat_data[iy] <= f_in_and_out.lat_data[f_in_and_out.dimy-1])
				if (f_temp.get(ix,iy) == 1)
					f_in_and_out.set_by_lon_and_lat(f_lsf.lon_data[ix], f_lsf.lat_data[iy], 1);

	f_temp.freememory();

	//WriteField("lsf.nc",f_temp);
	//WriteField("f_in_and_out.nc",f_in_and_out);
	}

/*
// intended for PNG output of f_obj CFields. Each object number generates different random number
void WriteField_f_obj_to_PNG_with_landseamask(string fname, CField f_obj, CField f_lsm)
	{
	int ix, iy;

	double t1, t2, t3;
	double r,g,b;

	pngwriter image(f_obj.dimx, f_obj.dimy, 1.0, fname.c_str());

	long idum=-1;

		for (iy=0;iy<f.dimy;iy++)
		for (ix=0;ix<f.dimx;ix++)
			{
			//missing value check
			t1=f.get(ix,iy);

			// missing value
			if (t1==BAD_DATA_FLOAT) image.plot(ix+1, iy+1, 1.0, 1.0, 1.0);

			// no object
			else if (t1==-1) image.plot(ix+1, iy+1, 0.0, 0.0, 0.0);

			// object found
			else
				{
				r=g=b=0;
				idum=(long)-t1-1;

				// cout << ix << "\t" << iy << "\t" << idum << endl;

				// eliminate darkest colors - since they are too similar to the black background
				while (r+g+b < 0.5)
					{
					r=ran2(&idum);
					g=ran2(&idum);
					b=ran2(&idum);
					}
				// cout << ix << "\t" << iy << "\t" << idum << "\t" << r << "\t" << g << "\t" << b << endl;


				image.plot(ix+1, iy+1, r, g, b);
				}
			}
	image.close();
	}
*/

// intended for PNG output of f_obj CFields. Each object number generates different random number
#ifndef DISABLE_PNGWRITER
void WriteField_PNG(string fname,  const CField &f )
	{
	int ix, iy;

	double t1;
	double r,g,b;

	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());

	long idum=-1;

		for (iy=0;iy<f.dimy;iy++)
		for (ix=0;ix<f.dimx;ix++)
			{
			//missing value check
			t1=f.get(ix,iy);

			// missing value
			if (t1==BAD_DATA_FLOAT) image.plot(ix+1, iy+1, 1.0, 1.0, 1.0);

			// no object
			else if (t1==-1) image.plot(ix+1, iy+1, 0.0, 0.0, 0.0);

			// object found
			else
				{
				r=g=b=0;
				idum=(long)-t1-1;

				// cout << ix << "\t" << iy << "\t" << idum << endl;

				// eliminate darkest colors - since they are too similar to the black background
				while (r+g+b < 0.5)
					{
					r=ran2(&idum);
					g=ran2(&idum);
					b=ran2(&idum);
					}
				// cout << ix << "\t" << iy << "\t" << idum << "\t" << r << "\t" << g << "\t" << b << endl;


				image.plot(ix+1, iy+1, r, g, b);
				}
			}
	image.close();
	}
#endif

// intended for PNG output of f_obj CFields. Each object number generates different random number
#ifndef DISABLE_PNGWRITER
void WriteField_PNG_v2(string fname,  const CField &f )
	{
	int ix, iy;

	double t1;
	double r,g,b;

	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());

	long idum=-1;

		for (iy=0;iy<f.dimy;iy++)
		for (ix=0;ix<f.dimx;ix++)
			{
			//missing value check
			t1=f.get(ix,iy);

			// missing value
			if (t1==BAD_DATA_FLOAT) image.plot(ix+1, iy+1, 0.5, 0.5, 0.5);

			// no object
			else if (t1==-1) image.plot(ix+1, iy+1, 1.0, 1.0, 1.0);

			// object found
			else
				{
				r=g=b=0;
				idum=(long)-t1-1;

				// cout << ix << "\t" << iy << "\t" << idum << endl;

				// eliminate darkest colors - since they are too similar to the black background
				while (r+g+b < 1 || r+g+b > 2.5)
					{
					r=ran2(&idum);
					g=ran2(&idum);
					b=ran2(&idum);
					}
				// cout << ix << "\t" << iy << "\t" << idum << "\t" << r << "\t" << g << "\t" << b << endl;


				image.plot(ix+1, iy+1, r, g, b);
				}
			}
	image.close();
	}
#endif





// intended for PNG output of f_obj CFields. Each object number generates different random number
#ifndef DISABLE_PNGWRITER
void WriteField_PNG_blend(string fname,  const CField &f_obj,  const CField &f_crop,  const CField &f_coasts, double f_crop_max_value, double opacity, double brightnes_of_low_preciptation,  const vector <CField> &vector_f_optional_yellow_overlay)
	{
	int ix, iy, ip;

	double t1, t2, t3;
	double r,g,b;
	double temp_opacity=opacity;

	pngwriter image(f_obj.dimx+2, f_obj.dimy+2, 0, fname.c_str());

	long idum=-1;


	// open coastal data area
	//CField f_coast=ReadField("/home/gskok/work/B_13_Polygon_point_overlay/polygon_world_countryborders.dat.Pacific1.nc");

	for (iy=0;iy<f_obj.dimy;iy++)
	for (ix=0;ix<f_obj.dimx;ix++)
		{
		//------------------------------
		// f_obj value
		t1=f_obj.get(ix,iy);

		// missing value
		if (t1==BAD_DATA_FLOAT) image.plot(ix+2, iy+2, 1.0, 1.0, 1.0);
			// no object
		else if (t1==-1) image.plot(ix+2, iy+2, 0.0, 0.0, 0.0);

		// object found
		else
			{
			r=g=b=0;
			idum=(long)-t1-1;

			// cout << ix << "\t" << iy << "\t" << idum << endl;
				// eliminate darkest colors - since they are too similar to the black background
			while (r+g+b < 0.5)
				{
				r=ran2(&idum);
				g=ran2(&idum);
				b=ran2(&idum);
				}
			// cout << ix << "\t" << iy << "\t" << idum << "\t" << r << "\t" << g << "\t" << b << endl;


			image.plot(ix+2, iy+2, r, g, b);
			}

		// ----------------------------
		// original f_crop field
		// f_crop value
		t2=f_crop.get(ix,iy);
		temp_opacity=opacity;

		// check for max value
		if (t2> f_crop_max_value)
			t2=f_crop_max_value;

		t3=(1-brightnes_of_low_preciptation)/f_crop_max_value*t2 +brightnes_of_low_preciptation;

		// t2=0
 		if (t2==0)
			temp_opacity=0;

		// missing value
 		if (t2< 0)
			temp_opacity=0;

		image.plot_blend(ix+2, iy+2, temp_opacity, t3, t3, t3);


		// ----------------------------
		// coastal data
		t2=f_coasts.get(ix,iy);
		if (t2 == 1)
		image.plot_blend(ix+2, iy+2, 0.5, 1.0, 1.0, 1.0);

		// ----------------------------
		// vector_f_optional_yellow_overlay
		for (ip=0; ip < (int)vector_f_optional_yellow_overlay.size(); ip++)
			if (vector_f_optional_yellow_overlay[ip].get(ix,iy) == 1)
				image.plot_blend(ix+2, iy+2, 1.0, 1.0, 0.0, 0.0);

		}

	// frame
	image.square_blend( 1, 1, f_obj.dimx+2, f_obj.dimy+2, 0.5, 1.0, 1.0, 1.0);

	image.close();

	//f_coast.freememory();

	}
#endif

// intended for PNG output of f_obj CFields. Each object number generates different random number
#ifndef DISABLE_PNGWRITER
void WriteField_PNG_contourfield_grayscale(string fname,  const CField &f, double minlevel, double maxlevel)
	{
	int ix, iy;
	double gray;

	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());


	for (iy=0;iy<f.dimy;iy++)
	for (ix=0;ix<f.dimx;ix++)
		{
		if (f.get(ix,iy) >= maxlevel)
			gray=0;
		else if (f.get(ix,iy) < minlevel)
			gray=1;
		else
			gray=(maxlevel-f.get(ix,iy))/(maxlevel-minlevel);

		image.plot(ix+1, iy+1, gray, gray, gray);
		}
	image.close();
	}
#endif

// intended for PNG output of f_obj CFields. Each object number generates different random number
#ifndef DISABLE_PNGWRITER
void WriteField_PNG_contourfield_with_levels_and_colorbar(string fname,  const CField &f, vector <double> &levels, vector <CRGBcolor> &colors )
	{
	int ix, iy;
	CRGBcolor rgb1;

	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());

	for (ix=0; ix < f.dimx; ix++)
	for (iy=0; iy < f.dimy; iy++)
		if (f.get(ix,iy) >= 0)
			{
			rgb1.set_RGB_color_from_colortable_discrete_colors(f.get(ix,iy), levels, colors);
			image.plot(ix+1, iy+1, rgb1.r, rgb1.g, rgb1.b);
			}
	image.close();
	}
#endif


#ifndef DISABLE_PNGWRITER
void draw_pngimage_of_binary_field(string fname, const CField &f )
	{
	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());
	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		// only 0 and 1 values are allowed
		ERRORIF(f.get(ix,iy) != 0 && f.get(ix,iy) != 1);
		if (f.get(ix,iy) > 0.5)
			image.plot(ix+1, iy+1, 0., 0., 0.);
		else
			image.plot(ix+1, iy+1, 1., 1.0, 1.0);
		}
	image.close();
	}
#endif

#ifndef DISABLE_PNGWRITER
void draw_pngimage_of_binary_field_reverse_color(string fname, const CField &f )
	{
	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());
	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		// only 0 and 1 values are allowed
		ERRORIF(f.get(ix,iy) != 0 && f.get(ix,iy) != 1);
		if (f.get(ix,iy) > 0.5)
			image.plot(ix+1, iy+1, 1., 1.0, 1.0);
		else
			image.plot(ix+1, iy+1, 0., 0., 0.);
		}
	image.close();
	}
#endif


#ifndef DISABLE_PNGWRITER
void draw_grayscale_pngimage_between_0_and_1(string fname, const CField &f )
	{
	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());
	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		// only values between 0 and 1 values are allowed
		ERRORIF(f.get(ix,iy) < 0 || f.get(ix,iy) > 1);
		double color=f.get(ix,iy);
		image.plot(ix+1, iy+1, color, color, color);
		}
	image.close();
	}
#endif


#ifndef DISABLE_PNGWRITER
void draw_grayscale_pngimage_between_minlevel_and_maxlevel(string fname, const CField &f,double minlevel, double maxlevel )
	{
	double dlev=fabs(maxlevel-minlevel);
	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());
	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		double color;
		if (maxlevel > minlevel)
			color=(f.get(ix,iy)-minlevel)/dlev;
		else
			color=1-(f.get(ix,iy)-maxlevel)/dlev;
		image.plot(ix+1, iy+1, color, color, color);
		}
	image.close();
	}
#endif

#ifndef DISABLE_PNGWRITER
void draw_grayscale_pngimage_between_minlevel_and_maxlevel2_invert(string fname, const CField &f,double minlevel, double maxlevel )
	{
	double dlev=maxlevel-minlevel;
	ERRORIF(dlev < 0);
	pngwriter image(f.dimx, f.dimy, 1.0, fname.c_str());
	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		double color;
		if (f.get(ix,iy) <= minlevel)
			color=0;
		else if (f.get(ix,iy) >= maxlevel)
			color=1;
		else
			color=(f.get(ix,iy)-minlevel)/dlev;

		// invert colors
		image.plot(ix+1, iy+1, 1-color, 1-color, 1-color);
		}
	image.close();
	}
#endif


#ifndef DISABLE_PNGWRITER
CField read_grayscale_pngimage_between_0_and_1_into_CField(string fname)
	{
	pngwriter image(1,1,0,"aaa.png");
	image.readfromfile(fname.c_str());

	CField f;
	f.allocatememory_and_set_dimx_dimy(image.getwidth(), image.getheight());
	f.allocatememory_lat_lon_array();

	for (long ix=0; ix < f.dimx; ix++)
	for (long iy=0; iy < f.dimy; iy++)
		{
		double r=image.dread(ix+1, iy+1,1);
		double g=image.dread(ix+1, iy+1,2);
		double b=image.dread(ix+1, iy+1,3);
		f.set(ix,iy,(r+g+b)/3);
		}
	// convert white areas to totally white
	f.Set_all_values_in_value_range(1, 0.98, 1.5);
	image.close();

	return(f);
	}
#endif

#ifndef DISABLE_PNGWRITER
pngwriter verticaly_merge_two_figures(pngwriter &im1, pngwriter &im2)
	{
	long dimx1=im1.getwidth();
	long dimy1=im1.getheight();
	long dimx2=im2.getwidth();
	long dimy2=im2.getheight();

	pngwriter im(max(dimx1,dimx2),dimy1+dimy2,1.,"");

	for (long ix=1; ix <= dimx1; ix++)
	for (long iy=1; iy <= dimy1; iy++)
		{
		double r=im1.dread(ix, iy,1);
		double g=im1.dread(ix, iy,2);
		double b=im1.dread(ix, iy,3);
		im.plot(ix,iy+dimy2,r,g,b);
		}

	for (long ix=1; ix <= dimx2; ix++)
	for (long iy=1; iy <= dimy2; iy++)
		{
		double r=im2.dread(ix, iy,1);
		double g=im2.dread(ix, iy,2);
		double b=im2.dread(ix, iy,3);
		im.plot(ix,iy,r,g,b);
		}
	return(im);
	}
#endif

