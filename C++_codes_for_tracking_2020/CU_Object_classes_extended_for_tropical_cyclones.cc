// ----------------------------------------------------------------------------------
// A class to represent properties of a 3D Object slice at a certain time.
// ----------------------------------------------------------------------------------
 class CFrameObject_extended
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
	double u;
	double v;
	double lat;
	double lon;
	double min_ps_value_land;
	double min_ps_value_ocean;
	double max_wind_value_land;
	double max_wind_value_ocean;
	double precipitation_volume;


	CFrameObject_extended() {freememory();}
	~CFrameObject_extended() {freememory();}

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
		u=BAD_DATA_FLOAT;
		v=BAD_DATA_FLOAT;
		lat=BAD_DATA_FLOAT;
		lon=BAD_DATA_FLOAT;
		min_ps_value_land=BAD_DATA_FLOAT;
		min_ps_value_ocean=BAD_DATA_FLOAT;
		max_wind_value_land=BAD_DATA_FLOAT;
		max_wind_value_ocean=BAD_DATA_FLOAT;
		precipitation_volume=BAD_DATA_FLOAT;
		}

	CFrameObject_extended deep_copy() const
		{
		CFrameObject_extended out;
		out=*this;
		return(out);
		}

	};

// A class to represent properties of a 3D Object. Basically it is a container of CFrameObject_extended-s with some properties such as lifespan, ....
class C3DObject_extended
	{
	public:
	double object_id;
	long timestamp_start;
	long timestamp_end;

	vector <CFrameObject_extended> FrameObject;

	C3DObject_extended() {freememory();}
	~C3DObject_extended() {freememory();}

	void freememory()
		{
		object_id=-1;
		timestamp_start=timestamp_end=-1;
		FrameObject.clear();
		}

	C3DObject_extended deep_copy() const
		{
		C3DObject_extended out;
		out=*this;
		return(out);
		}

	long lifespan() const
		{
		return(FrameObject.size());
		}

	bool did_object_exist_at_timestamp(long timestamp) const
		{
		if (timestamp < timestamp_start || timestamp > FrameObject.back().timestamp)
			return(false);
		return(true);
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

	// returns area_cos of object when the object was the largest
	double return_max_area_size() const
		{
		double max_area=0;

		for (unsigned long ix = 0; ix < FrameObject.size(); ix++)
			if (FrameObject[ix].area_cos > max_area )
				max_area=FrameObject[ix].area_cos;

		return(max_area);
		}

	// returns minimal min_ps_value_ocean of all frames
	double return_min_min_ps_value_ocean() const
		{
		double min=10000000;

		for (unsigned long ix  = 0; ix < FrameObject.size(); ix++)
			if (FrameObject[ix].min_ps_value_ocean < min )
				min=FrameObject[ix].min_ps_value_ocean;

		return(min);
		}

	// returns max max_wind_value_ocean of all frames
	double return_max_max_wind_value_ocean() const
		{
		double max=-1;

		for (unsigned long ix  = 0; ix < FrameObject.size(); ix++)
			if (FrameObject[ix].max_wind_value_ocean > max )
				max=FrameObject[ix].max_wind_value_ocean;

		return(max);
		}

	// returns max max_wind_value_ocean of all frames
	double return_total_precipitation_volume() const
		{
		double prec=0;

		for (unsigned long ix  = 0; ix < FrameObject.size(); ix++)
			prec+=FrameObject[ix].precipitation_volume;

		return(prec);
		}


	};



// A class to represent properties of a 3D Object. Basically it is a container of CFrameObject_extended-s with some properties such as lifespan, ....
class C3DObjectGroup_extended
	{
	public:
	vector <C3DObject_extended> Object3D;
	vector <long> index;

	C3DObjectGroup_extended() {freememory();}
	~C3DObjectGroup_extended() {freememory();}

	void freememory()
		{
		Object3D.clear();
		index.clear();
		}

	unsigned long number_of_initialized_objects() const
		{
		return(Object3D.size());
		}


	// this is index for faster search of data
	void set_index()
		{
		long il;
		long max_index = -1;

		for (il=0; il < (long)Object3D.size(); il++)
			if (max_index < Object3D[il].object_id ) max_index = (long)Object3D[il].object_id;

		index.assign(max_index+1,-1);
		for (il=0; il < (long)Object3D.size(); il++)
			index[(long)Object3D[il].object_id]=il;
		}


	long get_index(long il) const
		{
		if (index[il] < 0)
			error(AT,FU, "index[il] < 0");
		return(index[il]);
		}

	C3DObjectGroup_extended deep_copy() const
		{
		C3DObjectGroup_extended out;
		out=*this;
		return(out);
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
					myFile.write ((char *) &Object3D[il].FrameObject[ily], sizeof(CFrameObject_extended));
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
		C3DObject_extended temp_3DObject;
		CFrameObject_extended temp_FrameObject;
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
				myFile.read ((char *) &temp_FrameObject, sizeof(CFrameObject_extended));
				Object3D.at(il).FrameObject.push_back(temp_FrameObject.deep_copy());
				}
			}
		myFile.close();
		}



	void Display_ascii_data_for_object_with_id(long id) const
		{
		long il,ixl;

		if (index.size() == 0)
			error(AT,FU, "index.size() == 0");
		if (index[id] == -1)
			error(AT,FU, "index[id] == -1");

		il=index[id];


		cout << il << "\t" << Object3D[il].object_id << "\t" << Object3D[il].lifespan() << endl;
		for (ixl=0; ixl < Object3D[il].lifespan(); ixl ++)
			{
			cout << ixl << "\t" <<
			Object3D[il].FrameObject[ixl].centerx << "\t" <<
			Object3D[il].FrameObject[ixl].centery << "\t" <<
			Object3D[il].FrameObject[ixl].area << "\t" <<
			Object3D[il].FrameObject[ixl].area_cos << "\t" <<
			Object3D[il].FrameObject[ixl].timestamp << "\t" <<

			Object3D[il].FrameObject[ixl].vx << "\t" <<
			Object3D[il].FrameObject[ixl].vy << "\t" <<
			Object3D[il].FrameObject[ixl].xmin << "\t" <<
			Object3D[il].FrameObject[ixl].xmax << "\t" <<
			Object3D[il].FrameObject[ixl].ymin << "\t" <<
			Object3D[il].FrameObject[ixl].ymax << "\t" <<

			Object3D[il].FrameObject[ixl].u << "\t" <<
			Object3D[il].FrameObject[ixl].v << "\t" <<
			Object3D[il].FrameObject[ixl].lat << "\t" <<
			Object3D[il].FrameObject[ixl].lon << "\t" <<
			Object3D[il].FrameObject[ixl].min_ps_value_land << "\t" <<
			Object3D[il].FrameObject[ixl].min_ps_value_ocean << "\t" <<
			Object3D[il].FrameObject[ixl].max_wind_value_land << "\t" <<
			Object3D[il].FrameObject[ixl].max_wind_value_ocean << "\t" <<
			Object3D[il].FrameObject[ixl].precipitation_volume << "\t" <<
			endl;
			}
		}

	void Display_ascii_data_for_N_objects(long n) const
		{
		long il;

		//for (il=0; il < n ; il ++)
			//if (Object3D[il].lifespan > 0)
				//cout << il << "\t" << Object3D.at(il).object_id << "\t" << Object3D.at(il).lifespan() << endl;

		long counter=0;
		il=0;
		while (counter < n && il < (long)number_of_initialized_objects() )
			{
			if (Object3D[il].lifespan() > 0)
				{
				Display_ascii_data_for_object_with_id((long)Object3D.at(il).object_id);
				counter++;
				}
			il++;
			}
		}

	};

