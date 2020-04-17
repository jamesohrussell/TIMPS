// ----------------------------------------------------------------------------------
// A class to represent properties of a 3D Object slice at a certain time.
// ----------------------------------------------------------------------------------
class CFrameObject
	{
	public:
	double object_id;
	long timestamp;
	int xmin;
	int xmax;
	int ymin;
	int ymax;

	double centerx,centery;
	double area;
	double area_cos;
	double pieces;
	double vx;
	double vy;

	CFrameObject() {freememory();}
	~CFrameObject() {freememory();}

	void freememory()
		{
		object_id=-1;
		timestamp=-1;
		centerx=centery=-1;
		area=area_cos=-1;
		pieces=-1;
		vx=0;
		vy=0;
		xmin=xmax=ymin=ymax=-1;
		}

	CFrameObject deep_copy() const
		{
		CFrameObject out;
		out=*this;
		return(out);
		}

	bool compare_are_they_the_same(const CFrameObject &f) const
		{
		if (object_id != f.object_id) return(false);
		if (timestamp != f.timestamp) return(false);
		if (xmin != f.xmin) return(false);
		if (xmax != f.xmax) return(false);
		if (ymin != f.ymin) return(false);
		if (ymax != f.ymax) return(false);

		if (centerx != f.centerx) return(false);
		if (centery != f.centery) return(false);
		if (area != f.area) return(false);
		if (area_cos != f.area_cos) return(false);
		if (pieces != f.pieces) return(false);
		if (vx != f.vx) return(false);
		if (vy != f.vy) return(false);
		return(true);
		}


	// overload operators == nad !=
	bool operator== (const CFrameObject& f) const
		{return(compare_are_they_the_same(f));}
	bool operator!= (const CFrameObject& f) const
		{return(!compare_are_they_the_same(f));}



	// calculates basic 2D object atributes such as center and area
	void Set_Object_Attributes(double id, long ts,  const CField &f,	int Xxmin, int Xxmax, int Xymin, int Xymax)
		{
		object_id=id;
		timestamp=ts;
		xmin=Xxmin;
		xmax=Xxmax;
		ymin=Xymin;
		ymax=Xymax;

		int ix, iy;
		double temp;

		centerx=centery=area=area_cos=0;

		area=0;
		area_cos=0;

		// Area size and center
		for (iy=ymin;iy<=ymax;iy++)
		for (ix=xmin;ix<=xmax;ix++)
			{
			temp=f.get(ix,iy);
			// object found
			if (temp == object_id)
				{
				area++;
				centerx+=(double)ix;
				centery+=(double)iy;

				area_cos+=f.area_of_grid_square_in_km2(ix,iy);
				}
			}

		if (area == 0)
				{
				cout << id << endl;
				error(AT,FU, "area == 0  !");
				}

		centerx=centerx/area;
		centery=centery/area;

		/*
		// modify the calculation of centerx if object spans more than dimx/2 - in this case it most likely crosses the periodic border - but the algorithm does not assume that. The algorithm then identifies all intervals of x where object does not exist, and points the x_break at the one which is the largest -> max_x_break
		int x_break;
		int first_x_break;
		bool found;
		int interval_length;
		int max_interval_length;
		int	max_x_break;
		long il;
		if (xmax - xmin > f.dimx/2)
			{
			cout << object_id << endl;

			centerx=0;
			// find where object is broken in x direction
			first_x_break=-1;
			for (ix=xmin;ix<=xmax;ix++)
				{
				found=false;
				for (iy=ymin;iy<=ymax;iy++)
					if (f.get(ix,iy) == object_id)
						{
						found=true;
						break;
						}
				if (found == false)
					{
					first_x_break=ix;
					break;
					}
				}


			if (first_x_break == -1)
				error(AT,FU, "first_x_break == -1: object intinterupted in x direction");

			interval_length=0;
			max_interval_length=0;
			for (il=0;il < f.dimx; il++)
				{
				found=false;
				ix=(first_x_break + il) % f.dimx;
				for (iy=ymin;iy<=ymax;iy++)
					if (f.get(ix,iy) == object_id)
						{
						found=true;
						break;
						}
				if (found == false)
					{
					if (interval_length == 0) x_break=ix;
					interval_length++;
					}

				else
					{
					if (interval_length > max_interval_length)
						{
						max_interval_length=interval_length;
						max_x_break=x_break;
						}
					interval_length=0;
					}
				}

			cout << first_x_break << " " << max_x_break << " " << max_interval_length << endl;

			// this is more of a sanity test - in fact such large objects could in theory exist but are very unlikely and would probably mean an error in code somewhere
			// max_interval_length is a maximal uninteruped swath of domain in x direction that does not contain this object
			if (max_interval_length < f.dimx/2)
				error(AT,FU, "max_interval_length < f.dimx/2");


			for (iy=ymin;iy<=ymax;iy++)
			for (ix=xmin;ix<=xmax;ix++)
				{
				temp=f.get(ix,iy);
				// object found
				if (temp == object_id)
					{
					if (ix < max_x_break) centerx+=(double)ix;
					else if (ix > max_x_break) centerx+=(double)(ix-f.dimx);
					}
				}
			centerx=centerx/area;
			cout << centerx << endl;
			if (centerx < 0)
				centerx=(double)f.dimx + centerx;
			}
		*/
		}


	};

// A class to represent properties of a 3D Object. Basically it is a container of CFrameObject-s with some properties such as lifespan, ....
class C3DObject
	{
	public:
	double object_id;
	long timestamp_start;
	long timestamp_end;

	vector <CFrameObject> FrameObject;

	C3DObject() {freememory();}
	~C3DObject() {freememory();}

	void freememory()
		{
		object_id=-1;
		timestamp_start=timestamp_end=-1;
		FrameObject.clear();
		}


	bool compare_are_they_the_same(const C3DObject &f) const
		{
		if (object_id != f.object_id) return(false);
		if (timestamp_start != f.timestamp_start) return(false);
		if (timestamp_end != f.timestamp_end) return(false);
		if (FrameObject != f.FrameObject) return(false);
		return(true);
		}

	// overload operators == nad !=
	bool operator== (const C3DObject& f) const
		{return(compare_are_they_the_same(f));}
	bool operator!= (const C3DObject& f) const
		{return(!compare_are_they_the_same(f));}



	long lifespan() const
		{
		return(FrameObject.size());
		}

	C3DObject deep_copy() const
		{
		C3DObject out;
		out=*this;
		return(out);
		}

	long get_object_subindex_from_timestamp(long timestamp) const
		{

		//cout << timestamp << " " << timestamp_start << " " << FrameObject.back().timestamp << endl;
		if (FrameObject.front().timestamp != timestamp_start )
			error(AT,FU, "FrameObject.front().timestamp != timestamp_start");


		if (timestamp < timestamp_start || timestamp > FrameObject.back().timestamp)
			error(AT,FU, "timestamp < timestamp_start || timestamp > FrameObject.back().timestamp");

		//cout << timestamp << "\t" << timestamp_start << "\t" <<  FrameObject.back().timestamp << endl;
		if (timestamp==timestamp_start) return(0);

		long time_stp=	(FrameObject.back().timestamp - timestamp_start)/(lifespan() - 1);

		return( (timestamp - timestamp_start)/(time_stp));
		}

	// A method that adds an additional frame to an 3D object. If this objects is new this will initialize the object. id specifies the id of object and can be only set at initialzation. timestamp is the timestampt of this frame. f is the current CField frame.
	// Velocity is only calculated if this is not an initialization. The velocity at frame_index=4, represents the displacent of centre from frame_index=3 to frame_index=4;
	void add_object_frame(double id, long timestamp,  const CField &f, int xmin, int xmax, int ymin, int ymax)
		{
		// set object id
		if (object_id == -1)
			object_id=id;

		else if (object_id != id)
			error(AT,FU, "object_id cannot be changed once assigned!");

		// if necesarry set timestamp_start
		if (timestamp_start == -1 || timestamp_start > timestamp) timestamp_start=timestamp;

		// if necesarry set timestamp_end
		if (timestamp_end == -1 || timestamp_end < timestamp) timestamp_end=timestamp;


		// get the current frame
		CFrameObject current_FrameObject;
		current_FrameObject.Set_Object_Attributes(object_id, timestamp, f, xmin, xmax, ymin, ymax);

		// add the new CFrameObject array to the end of vector
		FrameObject.push_back(current_FrameObject.deep_copy());

		current_FrameObject.freememory();

		// calculate velocity
		if (FrameObject.size() > 1)
			{
			FrameObject.back().vx=FrameObject.back().centerx - FrameObject.at(lifespan()-2).centerx;
			FrameObject.back().vy=FrameObject.back().centery - FrameObject.at(lifespan()-2).centery;
			//cout <<  FrameObject.back().centerx << "\t" << FrameObject.at(lifespan()-2).centerx << endl;
			}

		// make correction for domain periodic in x - if object crosses the x boundary
		if (FrameObject.back().vx > f.dimx/2) FrameObject.back().vx=-f.dimx+FrameObject.back().vx;
		if (FrameObject.back().vx < -f.dimx/2) FrameObject.back().vx=FrameObject.back().vx+f.dimx;
		}

	// returns area_cos of object when the object was the largest
	double return_max_area_size() const
		{
		unsigned int ix;
		double max_area=0;

		for (ix = 0; ix < FrameObject.size(); ix++)
			if (FrameObject[ix].area_cos > max_area )
				max_area=FrameObject[ix].area_cos;

		return(max_area);
		}


	};



// A class to represent properties of a 3D Object. Basically it is a container of CFrameObject-s with some properties such as lifespan, ....
class C3DObjectGroup
	{
	public:
	vector <C3DObject> Object3D;

	C3DObjectGroup() {freememory();}
	~C3DObjectGroup() {freememory();}

	void freememory()
		{
		Object3D.clear();
		}

	long number_of_initialized_objects() const
		{
		return(Object3D.size());
		}

	C3DObjectGroup deep_copy() const
		{
		C3DObjectGroup out;
		out=*this;
		return(out);
		}


	void add_frame_to_object(double id, long timestamp,  const CField &f, int xmin, int xmax, int ymin, int ymax)
		{
		C3DObject temp;

		// check if vector size is long enough and enlarge if necesarry
		if ((long)id + 1 > (long)Object3D.size() )
			Object3D.resize((long)id + 10000,temp);

		Object3D.at((long)id).add_object_frame(id, timestamp, f, xmin, xmax, ymin, ymax);
		}

	void Write_To_Binary_File(string filename) const
		{
		long il,ily;
		long lfspn;
		ofstream myFile (filename.c_str(), ios::out | ios::binary);

		// count initialized objects
		long object_counter=0;
		for (il=0; il < (long)Object3D.size(); il ++)
			if (Object3D[il].lifespan() > 0) object_counter++;

		// write object_counter to file
		myFile.write ((char *) &object_counter, sizeof(long));

		// write initialized Object3D-s to file
		for (il=0; il < (long)Object3D.size(); il ++)
			if (Object3D[il].lifespan() > 0)
				{
				myFile.write ((char *) &Object3D[il].object_id, sizeof(double));
				myFile.write ((char *) &Object3D[il].timestamp_start, sizeof(long));
				myFile.write ((char *) &Object3D[il].timestamp_start, sizeof(long));

				// write lifespan data
				lfspn=Object3D[il].lifespan();
				myFile.write ((char *) &lfspn, sizeof(long));

				for (ily=0; ily < Object3D[il].lifespan(); ily++)
					myFile.write ((char *) &Object3D[il].FrameObject[ily], sizeof(CFrameObject));
				}

		myFile.close();
		}

	void Read_From_Binary_File(string filename)
		{
		long lfspn;
		long il,ily;

		if (!FileExists(filename))
			error(AT,FU, filename + " file does not exist!");

		ifstream myFile (filename.c_str(), ios::out | ios::binary);

		long object_counter=0;
		myFile.read ((char *) &object_counter, sizeof(long));

		// read only object_counter number of Object3D-s
		C3DObject temp_3DObject;
		CFrameObject temp_FrameObject;
		for (il=0; il < object_counter; il ++)
			{
			temp_3DObject.freememory();
			myFile.read ((char *) &temp_3DObject.object_id, sizeof(double));
			myFile.read ((char *) &temp_3DObject.timestamp_start, sizeof(long));
			myFile.read ((char *) &temp_3DObject.timestamp_start, sizeof(long));
			Object3D.push_back(temp_3DObject.deep_copy());

			myFile.read ((char *) &lfspn, sizeof(long));

			//cout << " " <<il  << " " << Object3D.back().object_id  << " " << lfspn  << " " << endl;

			for (ily=0; ily < lfspn; ily++)
				{
				myFile.read ((char *) &temp_FrameObject, sizeof(CFrameObject));
				Object3D.at(il).FrameObject.push_back(temp_FrameObject.deep_copy());
				}
			}
		myFile.close();
		}

	C3DObjectGroup Reorganize_objects_so_that_index_equals_object_id()
		{
		long il;

		C3DObjectGroup out;

		// determine max_object_id
		double max_object_id=0;
		for (il=0; il < number_of_initialized_objects(); il ++)
			if (max_object_id < Object3D.at(il).object_id ) max_object_id=Object3D.at(il).object_id ;


		// initialize memory
		C3DObject temp_C3DObject;
		out.Object3D.assign((long)max_object_id+1,temp_C3DObject);

		// copy data
		for (il=0; il < number_of_initialized_objects(); il ++)
			{
			//cout << il << " " << Object3D.size() << " " << temp3DObjectGroup.Object3D.at(il).object_id << endl;
			out.Object3D.at((long)Object3D.at(il).object_id)=Object3D.at(il).deep_copy();
			}

		return(out);
		}

	// selects only objects that have at lest one center inside the defined area
	C3DObjectGroup create_a_new_object_subset_for_specific_lat_lon_area( const CField &f, double lon_min, double lon_max, double lat_min, double lat_max)  const
		{
		long ily,il;
		double lat,lon;

		C3DObjectGroup out;

		bool inside;
		CFrameObject temp_FrameObject;
		for (il=0; il < (long)Object3D.size(); il ++)
			if (Object3D[il].lifespan() > 0)
				{
				inside= false;

				for (ily=0; ily < Object3D[il].lifespan(); ily++)
					{
					temp_FrameObject=Object3D[il].FrameObject[ily];
					f.XY2LatLon((int)temp_FrameObject.centerx,(int)temp_FrameObject.centery, lon, lat);
					if (lon >= lon_min &&
						lon <= lon_max &&
						lat >= lat_min &&
						lat <= lat_max)
							inside=true;
					}

				// add object if inside
				if (inside)
					out.Object3D.push_back(Object3D[il].deep_copy());
				}

		return(out);
		}


	void Display_ascii_data_for_object_with_id(long id) const
		{
		long il,ixl;

		for (il=0; il < (long)Object3D.size() ; il ++)
			if (Object3D[il].object_id == id)
				{
				cout << il << "\t" << Object3D[il].object_id << "\t" << Object3D[il].lifespan() << endl;
				for (ixl=0; ixl < Object3D[il].lifespan(); ixl ++)
					{
					cout << ixl << "\t" <<
					Object3D[il].FrameObject[ixl].centerx << "\t" <<
					Object3D[il].FrameObject[ixl].centery << "\t" <<
					Object3D[il].FrameObject[ixl].area << "\t" <<
					Object3D[il].FrameObject[ixl].area_cos << "\t" <<
					output_date_string_to_hour_from_timestamp(Object3D[il].FrameObject[ixl].timestamp) << "\t" <<
					Object3D[il].FrameObject[ixl].vx << "\t" <<
					Object3D[il].FrameObject[ixl].vy << "\t" <<
					Object3D[il].FrameObject[ixl].xmin << "\t" <<
					Object3D[il].FrameObject[ixl].xmax << "\t" <<
					Object3D[il].FrameObject[ixl].ymin << "\t" <<
					Object3D[il].FrameObject[ixl].ymax << endl;
					}
				}
		}

	void Display_ascii_data_for_Nth_object(long n) const
		{
		if (n >= (long)Object3D.size())
			error(AT,FU, "n >= (long)Object3D.size()");

		Display_ascii_data_for_object_with_id(Object3D[n].object_id);
		}


	void Display_ascii_data_for_N_objects(long n) const
		{
		long il;


		//for (il=0; il < n ; il ++)
			//if (Object3D[il].lifespan > 0)
				//cout << il << "\t" << Object3D.at(il).object_id << "\t" << Object3D.at(il).lifespan() << endl;


		long counter=0;
		il=0;
		while (counter < n && il < number_of_initialized_objects() )
			{
			if (Object3D[il].lifespan() > 0)
				{
				Display_ascii_data_for_object_with_id((long)Object3D.at(il).object_id);
				counter++;
				}
			il++;
			}
		}


	#ifndef DISABLE_PNGWRITER
	void draw_trajectories(string png_name_template,  const CField &f_tempate) const
		{
		double r,g,b;
		double direction_r,direction_g,direction_b;

		CMy_time_format tt, ttemp;

		ostringstream s1,s2,s3,s4,s5,s6,s7,dataset;

		int object_lifespan;
		long ixl, iyl;
		long il;
		int x0,y0,x1,y1;
		double lon0,lat0,lon1, lat1;
		int ix,ik;

		pngwriter image;
		pngwriter image_direction;


		double max_detected_lifespan=-1;
		for (il=0; il < number_of_initialized_objects(); il ++)
		if (Object3D[il].lifespan() > max_detected_lifespan) max_detected_lifespan=Object3D[il].lifespan();

		CField f,f_coast;
		double scale;

		//*******************************************++
		// Drawing images


		// create high-resolution Cfield from tempate
		scale=1400/(double)f_tempate.dimx;
		f=f_tempate.create_a_scaled_copy_without_data(scale);
		f.allocatememory();
		f.Set_all_values(0);

		// coastal data
		Group_of_Polygons gp=ReadPolygonFromFile("/home/gskok/work/B_13_Polygon_point_overlay/polygon_world_countryborders.dat");
		f_coast=Point_Overlay_Polygon(gp,f,2);


		// write PNG images
		s1.str("");
		s1 << png_name_template << "_lifespan.png";
		image = pngwriter(f.dimx+2, f.dimy+2, 1.0, s1.str().c_str());
		s1.str("");
		s1 << png_name_template << "_direction.png";
		image_direction = pngwriter(f.dimx+2, f.dimy+2, 1.0, s1.str().c_str());



		// draw smaller object first loop
		for (object_lifespan=2; object_lifespan <= max_detected_lifespan; object_lifespan ++)
			{
			// draw lines
			for (ixl=0; ixl < number_of_initialized_objects(); ixl ++)
				{
				if (Object3D[ixl].lifespan() == object_lifespan)
					{
					//cout << object_lifespan << " " << ixl << endl;

					for (iyl=1; iyl < Object3D[ixl].lifespan(); iyl++)
						{
						ttemp.timestamp=Object3D[ixl].FrameObject[iyl].timestamp;
						ttemp.timestamp_to_date();

						// to color objects according to lifespan
						if (object_lifespan < 10)
							{
							r=0;
							g=0;
							b=1;
							}
						else if (object_lifespan >= 10 && object_lifespan < 20)
							{
							r=1;
							g=176.0/256.0;
							b=30.0/256.0;
							}
						else if (object_lifespan >= 20 && object_lifespan < 30)
							{
							r=0;
							g=158.0/256;
							b=0;
							}
						else
							{
							r=1;
							g=0;
							b=0;
							}


						if (Object3D[ixl].FrameObject[iyl].vx > 0)
							{
							direction_r=255.0/256.0;
							direction_g=176.0/256.0;
							direction_b=30.0/256.0;
							}
						else
							{
							direction_r=0.0/256.0;
							direction_g=158.0/256.0;
							direction_b=0.0;
							}

						f_tempate.XY2LatLon((int)(round(Object3D[ixl].FrameObject[iyl-1].centerx)),
											(int)(round(Object3D[ixl].FrameObject[iyl-1].centery)),
											lon0, lat0);
						//cout << X3DObjectGroup_inside_subdomain.Object3D[ixl].FrameObject[iyl-1].centerx*scale << " " << X3DObjectGroup_inside_subdomain.Object3D[ixl].FrameObject[iyl-1].centery*scale << " " << lon0 << " " << lat0 << endl;
						f.LatLon2XY(lon0, lat0, x0, y0);

						f_tempate.XY2LatLon((int)(round(Object3D[ixl].FrameObject[iyl].centerx)),
											(int)(round(Object3D[ixl].FrameObject[iyl].centery)),
											lon1, lat1);
						f.LatLon2XY(lon1, lat1, x1, y1);

						 // years
						image.line_blend(x0+2, y0+2, x1+2, y1+2, 1.0, r, g, b);
						image_direction.line_blend(x0+2, y0+2, x1+2, y1+2, 1.0, direction_r, direction_g, direction_b);


						}
					}
				}
			}


		// coasts
		for (ik=0;ik<f.dimy;ik++)
		for (ix=0;ix<f.dimx;ix++)
			if (f_coast.get(ix,ik)==1)
				{
				image.plot_blend(ix+2, ik+2, 1.0, 0, 0, 0);
				image_direction.plot_blend(ix+2, ik+2, 1.0, 0, 0, 0);
				}


		// frame
		image.square_blend( 1, 1, f.dimx+2, f.dimy+2, 0.5, 0, 0, 0);
		image_direction.square_blend( 1, 1, f.dimx+2, f.dimy+2, 0.5, 0, 0, 0);

		// write image
		image.close();
		image_direction.close();

		f.freememory();
		f_coast.freememory();
		}
	#endif

	};

