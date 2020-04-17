

void read_tropical_cyclone_measured_tracks_from_txt_files(string in_dir, vector<vector<double> > &lon3, vector<vector<double> > &lat3, vector<vector<double> > &time3, vector<string> &name3, long cutoff_timestamp, long cutoff_end_timestamp)
	{
	int ix;
	long il;

	CMy_time_format tt, ttemp;

	vector <string> line;
	vector <string> value_string;
	vector <double> value_number;

	vector<vector<double> > lon;
	vector<vector<double> > lat;
	vector<vector<double> > time;
	vector<string> name;

	vector<vector<double> > lon2;
	vector<vector<double> > lat2;
	vector<vector<double> > time2;

	// read lon
	line=get_each_line_as_string_from_txt_file(in_dir + "tropical_cyclones_tracks_lon.txt");
	for (il=0; il < (long)line.size(); il++)
		{
		value_string=StringExplode(line[il],",");
		value_number.assign(value_string.size(),BAD_DATA_FLOAT);
		for (ix=0; ix < (long)value_string.size(); ix++)
			{
			value_number[ix]=atof(value_string[ix].c_str())/100;

			// remap lon to positive values
			if (value_number[ix] < 0) value_number[ix]+=360;
			}


		lon.push_back(value_number);
		}

	// read lat
	line=get_each_line_as_string_from_txt_file(in_dir + "tropical_cyclones_tracks_lat.txt");
	for (il=0; il < (long)line.size(); il++)
		{
		value_string=StringExplode(line[il],",");
		value_number.assign(value_string.size(),BAD_DATA_FLOAT);
		for (ix=0; ix < (long)value_string.size(); ix++)
			value_number[ix]=atof(value_string[ix].c_str())/100;
		lat.push_back(value_number);
		}

	// read time
	line=get_each_line_as_string_from_txt_file(in_dir + "tropical_cyclones_tracks_time.txt");
	for (il=0; il < (long)line.size(); il++)
		{
		value_string=StringExplode(line[il],",");
		value_number.assign(value_string.size(),BAD_DATA_FLOAT);
		for (ix=0; ix < (long)value_string.size(); ix++)
			{
			tt.set_by_modified_julian_day(atof(value_string[ix].c_str()));
			value_number[ix]=tt.timestamp;
			}
		time.push_back(value_number);
		}

	// read name
	line=get_each_line_as_string_from_txt_file(in_dir + "tropical_cyclones_tracks_name.txt");
	for (il=0; il < (long)line.size(); il++)
		name.push_back(line[il]);

	// make a check
	if (lon.size()==0 || lon.size() != lat.size() || lon.size() != time.size() || name.size() != time.size())
		error(AT,FU, "lon.size()==0 || lon.size() != lat.size() || lon.size() != time.size()");


	// since tc data is 6h, interpolate 3h in the middle of the 6h interval
	for (il=0; il < (long)lon.size(); il++)
		{
		lon2.push_back(vector <double> (lon[il].size()*2-1,0));
		lat2.push_back(vector <double> (lat[il].size()*2-1,0));
		time2.push_back(vector <double> (time[il].size()*2-1,0));

		for (ix=0; ix < (long)lon[il].size(); ix++)
			{
			// remap 6h values
			lon2[il][ix*2]=lon[il][ix];
			lat2[il][ix*2]=lat[il][ix];
			time2[il][ix*2]=time[il][ix];
			// interpolate intermidiate values
			if (ix > 0)
				{
				lon2[il][ix*2-1]=(lon[il][ix]+lon[il][ix-1])/2;
				lat2[il][ix*2-1]=(lat[il][ix]+lat[il][ix-1])/2;
				time2[il][ix*2-1]=(time[il][ix]+time[il][ix-1])/2;
				}
			}

		}


	lon3.clear();
	lat3.clear();
	time3.clear();
	name3.clear();

	// select only cyclones that existed after tbegin and before tend
	for (il=0; il < (long)time2.size(); il++)
		if ( time2[il].back() >= cutoff_timestamp && time2[il].front() <= cutoff_end_timestamp)
			{
			time3.push_back(time2[il]);
			lon3.push_back(lon2[il]);
			lat3.push_back(lat2[il]);
			name3.push_back(name[il]);
			}

	/*for (ix=0; ix < time3.size(); ix++)
		{
		cout << ix << "\t" << name3[ix] << "\t" << time3[ix].size() << endl;
		if (name3[ix]=="FERNANDA")
			for (il=0; il < time3[ix].size(); il++)
				{
				tt.timestamp=time3[ix][il];
				tt.timestamp_to_date();
				cout << tt.output_date_string() << " " ;
				}
		cout << endl;
		}

	*/
	//exit(1);
	/*cout << endl;
	for (ix=0; ix < lon2[1].size(); ix++)
		cout << time2[1][ix] << "\t";
	cout << endl;
	exit(1);
*/
	}

class CObject_tc
	{
	public:
	double id;
	bool is_defined;
	bool is_tc;
	vector <double> tc_index_list;
	double lifespan_of_associated_tc;
	double u;
	double ps;
	double how_many_times_idetified_as_tc;
	double total_obj_prec;
	double lifespan;

	CObject_tc() {freememory();}
	~CObject_tc() {freememory();}

	void freememory()
		{
		id=-1;
		is_defined=false;
		is_tc=false;
		tc_index_list.clear();
		lifespan_of_associated_tc=0;
		u=BAD_DATA_FLOAT;
		ps=BAD_DATA_FLOAT;
		how_many_times_idetified_as_tc=0;
		total_obj_prec=0;
		lifespan=0;
		}


	string output()
		{
		ostringstream s1;
		s1.str("");
		s1
		<< (long)id << "\t"
		<< is_defined << "\t"
		<< is_tc << "\t"
		<< tc_index_list.size() << "\t"
		<< lifespan_of_associated_tc << "\t";

		for (long il=0; il < (long)tc_index_list.size(); il++)
			s1 << tc_index_list[il] << ",";
		s1 << "\t";

		s1 << how_many_times_idetified_as_tc << "\t"
		<< total_obj_prec << "\t"
		<<  lifespan << "\t"
		<< u << "\t"
		<< ps << "\t";
		return(s1.str());
		}


	};


class CTropical_cyclone_frame
	{
	public:
	double lat;
	double lon;
	double prec;
	long timestamp;
	vector <double> obj_id;

	CTropical_cyclone_frame() {freememory();}
	~CTropical_cyclone_frame() {freememory();}

	void freememory()
		{
		lat=-1;
		lon=-1;
		prec=0;
		timestamp=-1;
		obj_id.clear();
		}

	template <typename T>
	void archive_object(T & myFile)
		{
		Archive(myFile, lat);
		Archive(myFile, lon);
		Archive(myFile, prec);
		Archive(myFile, timestamp);
		Archive(myFile, obj_id);
		}

	// overload operators == nad !=
	bool operator== (const CTropical_cyclone_frame& f) const
		{
		if (lat==f.lat && lon==f.lon && prec==f.prec && timestamp==f.timestamp && obj_id==f.obj_id) return(true);
		return(false);
		}


	CTropical_cyclone_frame deep_copy()
		{
		CTropical_cyclone_frame out;
		out=*this;
		out.obj_id=obj_id;
		return(out);
		}


	string output() const
		{
		ostringstream s1;
		CMy_time_format tt;
		s1.str("");

		tt.timestamp=timestamp;
		tt.timestamp_to_date();
		s1 << tt.output_date_string() << "\t" << lat << "\t" << lon << "\t" << prec << "\t" << obj_id.size() << ":";
		for (long il=0; il < (long)obj_id.size(); il++)
			s1 << obj_id[il] << ",";
		return(s1.str());
		}

	};

class CTropical_cyclone
	{
	public:
	vector <CTropical_cyclone_frame> frame;
	string name;

	CTropical_cyclone() {freememory();}
	~CTropical_cyclone() {freememory();}

	void freememory()
		{
		frame.clear();
		name="";
		}

	string output() const
		{
		ostringstream s1;
		s1.str("");

		int ix;
		s1 << name << endl;
		for (ix=0; ix < (long)frame.size(); ix++)
			s1 << ix << "\t" << frame[ix].output()  <<endl;
		return(s1.str());
		}

	CTropical_cyclone deep_copy()
		{
		CTropical_cyclone out;
		out=*this;
		out.frame=frame;
		return(out);
		}

	// overload operators == nad !=
	bool operator== (const CTropical_cyclone& f) const
		{
		if (name==f.name && frame==f.frame) return(true);
		return(false);
		}

	template <typename T>
	void archive_object(T & myFile)
		{
		Archive(myFile, name);
		Archive(myFile, frame);
		}

	bool did_tc_exist_at_timestamp(long timestamp) const
		{
		if (timestamp >= frame.front().timestamp && timestamp <= frame.back().timestamp)
			return(true);
		return(false);
		}

	int get_tc_subindex_from_timestamp(long timestamp) const
		{
		if (!did_tc_exist_at_timestamp(timestamp))
			error(AT,FU, "!did_tc_exist_at_timestamp(timestamp)");
		// avoid devision by 0 for size() == 1 tcs
		if (frame.size() < 2) return (0);
		return ((timestamp - frame.front().timestamp)/((frame.back().timestamp - frame.front().timestamp)/(frame.size() - 1)));
		}


	double get_precipitation() const
		{
		double out=0;
		int ix;
		for (ix=0; ix < (long)frame.size(); ix++)
			out+=frame[ix].prec;
		return(out);
		}
	};

class CTropical_cyclone_dataset
	{
	public:
	vector <CTropical_cyclone> tc;

	CTropical_cyclone_dataset() {freememory();}
	~CTropical_cyclone_dataset() {freememory();}
	void freememory()
		{
		tc.clear();
		}

	unsigned long number_of_tcs()
		{return(tc.size());}

	// overload operators == nad !=
	bool operator== (const CTropical_cyclone_dataset& f) const
		{
		if (tc==f.tc) return(true);
		return(false);
		}

	void get_tc_list_from_timestamp(long timestamp,  vector <long> &tc_index_list) const
		{
		long il;
		tc_index_list.clear();
		for (il=0; il < (long)tc.size(); il++)
			if (tc[il].did_tc_exist_at_timestamp(timestamp))
				tc_index_list.push_back(il);
		}


	void get_object_list_from_timestamp(long timestamp, vector <long> &obj_index_list) const
		{
		long il,it;
		long tc_subindex;
		vector <long> tc_index_list;
		get_tc_list_from_timestamp(timestamp, tc_index_list);

		obj_index_list.clear();
		for (il=0; il < (long)tc_index_list.size(); il++)
			{
			tc_subindex=tc[tc_index_list[il]].get_tc_subindex_from_timestamp(timestamp);
				for (it=0; it < (long)tc[tc_index_list[il]].frame[tc_subindex].obj_id.size(); it++)
					// comment below to enable or disable usingn only the closest object
					//if (it==0)
						obj_index_list.push_back((long)tc[tc_index_list[il]].frame[tc_subindex].obj_id[it]);
			}
		}

	// remove all frames that did not exist in time interval
	void remove_all_frames_that_did_not_exist_in_time_interval(long timestamp_start, long timestamp_end)
		{
		CTropical_cyclone_dataset tcd;
		CTropical_cyclone tc_temp;

		long il,ix;
		for (il=0; il < (long)tc.size(); il++)
			{
			tc_temp.freememory();
			tc_temp.name=tc[il].name;

			// populate the tc_temp with frames
			for (ix=0; ix < (long)tc[il].frame.size(); ix++)
				if (tc[il].frame[ix].timestamp >= timestamp_start && tc[il].frame[ix].timestamp <= timestamp_end)
					tc_temp.frame.push_back(tc[il].frame[ix].deep_copy());

			// add to new database if tc_temp is not empty
			if (tc_temp.frame.size() > 0)
				tcd.tc.push_back(tc_temp);
			}
		(*this).freememory();
		*this=tcd;
		}



	// initilize from original txt files
	void initilize_from_original_txt_files(string in_dir, long cutoff_timestamp, long cutoff_end_timestamp )
		{
		long il;
		int ix;
		vector<vector<double> > lon, lat, time;
		vector<string> name;

		read_tropical_cyclone_measured_tracks_from_txt_files(in_dir, lon, lat, time, name, cutoff_timestamp, cutoff_end_timestamp);

		// first consistency check
		if (lon.size()==0 || lon.size() != lat.size() || lon.size() != time.size() || lon.size() != name.size())
			error(AT,FU, "lon.size()==0 || lon.size() != lat.size() || lon.size() != time.size() || lon.size() != name.size()");
		// second consistency check
		for (il=0; il < (long)lon.size(); il++)
			if (lon[il].size()==0 || lon[il].size() != lat[il].size() || lon[il].size() != time[il].size())
				error(AT,FU, "lon[il].size()==0 || lon[il].size() != lat[il].size() || lon[il].size() != time[il].size()");

		//populate the vector
		CTropical_cyclone_frame	temp_tcf;
		CTropical_cyclone temp_tc;

		tc.assign(lon.size(),temp_tc);
		for (il=0; il < (long)lon.size(); il++)
			{
			tc[il].name=name[il];
			for (ix=0; ix < (long)lon[il].size(); ix++)
				{
				temp_tcf.lat=lat[il][ix];
				temp_tcf.lon=lon[il][ix];
				temp_tcf.timestamp=(long)time[il][ix];
				tc[il].frame.push_back(temp_tcf.deep_copy());
				}
			}

		}

	template <typename T>
	void archive_object(T & myFile)
		{
		Archive(myFile, tc);
		}

	void Write_to_binary_file(string filename)
		{
		ofstream myFile (filename.c_str(), ios::out | ios::binary);
		archive_object(myFile);
		myFile.close();
		}

	void Read_from_binary_file(string filename)
		{
		ifstream myFile (filename.c_str(), ios::out | ios::binary);
		archive_object(myFile);
		myFile.close();
		}
	};


/*
// --------------------------------------------------------------------
// Cost function for matching TCs (IBRtRACs) with objects using distances
// --------------------------------------------------------------------
double cost_function_for_matching_TCs_with_objects_using_distances(CTropical_cyclone_dataset &tcd, C3DObjectGroup_extended &X3DObjectGroup_extended_selection, double distance_threshold_limit)
	{
	vector <double> min_distance;
	C3DObject_extended *obj;
	CTropical_cyclone *tc;
	bool found;
	double prec_volume;
	long timestamp;
	double temp1;
	double distance;
	long subindex;
	double dlon,dlat;
	long io,iof,it,itf;
	long counter1,counter2,counter3,counter4;
	counter1=counter2=counter3=counter4=0;

	// ----------------------------------------------------------------------------
	// add all minimal distances from the objects to closest IBTRACS points
	// ----------------------------------------------------------------------------
	for (io=0; io < (long)X3DObjectGroup_extended_selection.number_of_initialized_objects(); io++)
		{
		obj=&X3DObjectGroup_extended_selection.Object3D[io];
		prec_volume+=obj->return_total_precipitation_volume();

		for (iof=0; iof < obj->lifespan(); iof++)
			{
			timestamp=obj->FrameObject[iof].timestamp;
			temp1=distance_threshold_limit;
			found=false;
			// loop over all TCs
			for (it=0; it < (long)tcd.tc.size(); it++)
				{
				tc=&tcd.tc[it];
				// check if TC existed in this timestep
				if (tc->did_tc_exist_at_timestamp(timestamp))
					{
					found=true;
					subindex=tc->get_tc_subindex_from_timestamp(timestamp);
					dlat=obj->FrameObject[iof].lat- tc->frame[subindex].lat;
					dlon=obj->FrameObject[iof].lon- tc->frame[subindex].lon;

					dlon=fabs(dlon);
					// correct dlon if more than 180 deg
					if (dlon > 180) dlon=360-dlon;

					// calculate distance
					distance=sqrt(dlat*dlat+dlon*dlon);

					// check if new minimal distance
					if (distance < temp1) temp1=distance;

					//cout << found << "\t" << temp1 << "\t" << distance << endl;
					}
				}

			min_distance.push_back(temp1);

			// add minimal distance if IBTRaCS was present in the timestep
			if (found)
				{
				if (temp1 == distance_threshold_limit)
					counter4++;
				else
					counter1++;
				}
			// a punishment if IBTRaCS was not present in the timestep
			else
				counter2++;

			}
		}

	// ----------------------------------------------------------------------------
	// add all minimal distances from the IBTRACS points to closest objects
	// ----------------------------------------------------------------------------
	for (it=0; it < (long)tcd.tc.size(); it++)
		{
		tc=&tcd.tc[it];
		for (itf=0; itf < (long)tc->frame.size(); itf++)
			{
			timestamp=tc->frame[itf].timestamp;
			temp1=distance_threshold_limit;
			found=false;

			// find the closest object
			for (io=0; io < (long)X3DObjectGroup_extended_selection.number_of_initialized_objects(); io++)
				{
				obj=&X3DObjectGroup_extended_selection.Object3D[io];
				if (obj->did_object_exist_at_timestamp(timestamp))
					{
					found=true;
					subindex=obj->get_object_subindex_from_timestamp(timestamp);

					dlat=obj->FrameObject[subindex].lat- tc->frame[itf].lat;
					dlon=obj->FrameObject[subindex].lon- tc->frame[itf].lon;

					dlon=fabs(dlon);
					// correct dlon if more than 180 deg
					if (dlon > 180) dlon=360-dlon;

					// calculate distance
					distance=sqrt(dlat*dlat+dlon*dlon);

					// check if new minimal distance
					if (distance < temp1) temp1=distance;
					}
				}

			//tt.set_by_timestamp(timestamp);
			//cout << it << " " << temp1 << " " << tt.output_date_string_to_hour() << endl;

			min_distance.push_back(temp1);

			// add minimal distance if object was present in the timestep
			if (found)
				{
				if (temp1 == distance_threshold_limit)
					counter4++;
				else
					counter1++;
				}
			// a punishment if object was not present in the timestep
			else
				counter3++;

			}
		}

	//cout << "/" << round_to_digits(average_from_vector(min_distance),2) << "("<<X3DObjectGroup_extended_selection.number_of_initialized_objects()<<"/"<<min_distance.size()<<"/"<<counter1<<"/"<<counter4<<"/"<<counter2<<"/"<<counter3 << ")\t";


	return(average_from_vector(min_distance));
	}

*/


