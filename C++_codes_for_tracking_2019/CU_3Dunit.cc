



// this is old way used of identifying the 3D objects used until 2010-11-02 - it treats objects as full 3D objects - both into the future and into the past - any overlap means the same object. In order this to work you need to keep a large number of old f_obj in memory since merger of really old objects can occur -  dimz needs to be sufficiently large. A modified version of this code was written around 2010-11-02 - which merges objects only forward in time.

// This object is needed for identifying of 3D objects
class Cgroup_of_CFields
	{
	public:
	int dimz, dimx, dimy;
	double object_count_index;
	vector <CField> f_obj;
	vector <CGroup_of_Identified_Objects_2D> Group_of_Identified_Objects_2D;


	// constructor
	Cgroup_of_CFields()
		{
		freememory();
		}

	// destructor
	~Cgroup_of_CFields()
		{
		freememory();
		}

	// stevilo time stepov v spominu
	int current_number_of_timesteps_in_memory()
		{
		return(f_obj.size());
		}


	// freememory
	void freememory()
		{
		long il;
		object_count_index=-1;
		dimz=dimx=dimy-1;
		for (il=0; il < (long)f_obj.size(); il++)
			f_obj[il].freememory();
		f_obj.clear();
		for (il=0; il < (long)Group_of_Identified_Objects_2D.size(); il++)
			Group_of_Identified_Objects_2D[il].freememory();
		Group_of_Identified_Objects_2D.clear();
		}


	void initialize(CField f, long Xdimz)
		{

		if (f_obj.size() != 0)
			error(AT,FU, " f_obj.size() != 0");

		freememory();
		dimx=f.dimx;
		dimy=f.dimy;
		dimz=Xdimz;
		object_count_index=0;
		}


	void add_time_step(CField f_thrs, bool is_domain_periodic_in_x)
		{
		long il;
		int ix;
		CField temp_f_obj;
		long index,index2;
		ostringstream s1;
		double previous,current;

		bool object_found_at_previus_timestep;
		CIdentified_Object_2D t1,t2;

		// we need an offset to start id-ing 3D objects. This is needed because there should be no inedx collision (same index) between any 3D object and any 2D object returned by Identify_Objects routine.
		double start_id_offset_for_3d_objects=100000;

		// check dimensions
		if (dimx !=f_thrs.dimx || dimy !=f_thrs.dimy)
			error(AT,FU, " dimx !=f.dimx || dimy !=f.dimy");

		CGroup_of_Identified_Objects_2D gt1;

		// Perform 2D identification of objects in current time step from the thresholded field
		temp_f_obj=f_thrs.Identify_Objects(gt1,is_domain_periodic_in_x);


		//temp_f_obj.check_consistency_of_f_obj_Field_with_Group_of_Identified_Objects_2D(gt1);
		//WriteField("aaa1.nc",temp_f_obj);

		//gt1.output_ascii();

		// check if largest id number of object identified is larger than start_id_offset_for_3d_objects
		if (gt1.get_max_id_in_vector() > start_id_offset_for_3d_objects - 1)
			error(AT,FU, " gt1.get_max_id_in_vector() > start_id_offset_for_3d_objects - 1");

		// check the number_of_all_available_timesteps
		int number_of_all_available_timesteps = f_obj.size();

		// if this is the first time step reset object_count_index
		if (number_of_all_available_timesteps == 0)
			object_count_index=start_id_offset_for_3d_objects;

		// check all grid points for identified objects
		for (il=0; il < temp_f_obj.size(); il++)
			{
			if (temp_f_obj.data[il] > -1)
				{
				// check if there is a object on the same place in previous time step
				if (number_of_all_available_timesteps == 0) object_found_at_previus_timestep=false;
				else if (f_obj[0].data[il] > -1)
					object_found_at_previus_timestep=true;
				else
					object_found_at_previus_timestep=false;

				// -------------------------------------------------------------------
				// if object found in prevous time step theh join objects
				// -------------------------------------------------------------------
				if (object_found_at_previus_timestep == true)
					{
					previous = f_obj[0].data[il];
					current = temp_f_obj.data[il];

					// check if object at previous time step if different than in current time step
					if (previous != current)
						{
						// get object info
						if (!gt1.check_if_object_is_defined(current, t1, index))
							{
							s1.str("");
							s1 << current << " " << il<< " " << temp_f_obj.DataLoc2X(il)  << " " << temp_f_obj.DataLoc2Y(il) << " object not found!";
							gt1.output_ascii();
							error(AT,FU, s1.str());
							}

						// you have to go back through time to check for objects with current object id
						if ( current >= start_id_offset_for_3d_objects)
							for (ix=0; ix < number_of_all_available_timesteps; ix++)
								{
								//cout << "aaaa" << endl;

								if (!Group_of_Identified_Objects_2D[ix].check_if_object_is_defined(current, t2, index2)) break;

								f_obj[ix].Replace_values_only_in_certain_area(current,previous,t2.xmin,t2.xmax,t2.ymin,t2.ymax);

								Group_of_Identified_Objects_2D[ix].join_two_Identified_Object_2D_objects(current,previous);
								}

						//finaly reindex the current object
						temp_f_obj.Replace_values_only_in_certain_area(current,previous,t1.xmin,t1.xmax,t1.ymin,t1. ymax);
						gt1.join_two_Identified_Object_2D_objects(current,previous);
						}
					}



				// -------------------------------------------------------------------
				// if object was not found in prevous time step - create new object if necesarry
				// -------------------------------------------------------------------
				else if (object_found_at_previus_timestep == false)
					{
					current = temp_f_obj.data[il];


					// check if this is the first time we see this object - id smaller than start_id_offset_for_3d_objects
					if (current < start_id_offset_for_3d_objects)
						{
						//cout << current << endl;

						//cout << gt1.check_if_object_is_defined(temp_f_obj.data[il], t1, index) << "\tbbb\t" << il << endl;
						// get object info
						if (!gt1.check_if_object_is_defined(current, t1, index) )
							{
							s1.str("");
							s1 << current << " " << il<< " " << temp_f_obj.DataLoc2X(il)  << " " << temp_f_obj.DataLoc2Y(il) << " object not foundX";

							//gt1.output_ascii();
							error(AT,FU, s1.str());
							}

						temp_f_obj.Replace_values_only_in_certain_area(current,object_count_index,t1.xmin,t1.xmax,t1.ymin,t1.ymax);

						//reindex object as a new 3D object
						if (!gt1.change_id_of_object(current,object_count_index))
							{
							s1.str("");
							s1 << current << " " << object_count_index << " " << temp_f_obj.DataLoc2X(il)  << " " << temp_f_obj.DataLoc2Y(il) << " object not found";
							//gt1.output_ascii();
							error(AT,FU, s1.str());
							}

						//!!! Need to enlarge the xmin,yim,ymax,ymin areas;


						object_count_index++;
						}

					}


				}

			}

		// compress gt1
		gt1.compress();

		//gt1.output_ascii();


		//cout << f_obj.size() << "\t" << Group_of_Identified_Objects_2D.size() << endl;

		// delete last timestep
		if (number_of_all_available_timesteps == dimz)
			{
			f_obj[dimz-1].freememory();
			f_obj.erase(f_obj.end());
			Group_of_Identified_Objects_2D[dimz-1].freememory();
			Group_of_Identified_Objects_2D.erase(Group_of_Identified_Objects_2D.end());
			}

		// add the currente timestep
		f_obj.insert(f_obj.begin(), temp_f_obj);
		Group_of_Identified_Objects_2D.insert(Group_of_Identified_Objects_2D.begin(), gt1.deep_copy());

		//cout << f_obj.size() << "\t" << Group_of_Identified_Objects_2D.size() << endl;

		gt1.freememory();
		}





	};




