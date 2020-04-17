

//struct Point{double x, y;};
class Point
		{
		public:
		double x, y, lat, lon;

		Point() {x=y=lat=lon=BAD_DATA_FLOAT;}

		Point deep_copy()
			{
			Point p;
			p.x=x;
			p.y=y;
			p.lat=lat;
			p.lon=lon;
			return(p);
			}

		void set (double _x, double _y)
			{
			x=_x;
			y=_y;
			}


		Point create(double Xx, double Xy, double Xlat, double Xlon)
			{
			Point p;
			p.x=Xx;
			p.y=Xy;
			p.lat=Xlat;
			p.lon=Xlon;
			return(p);
			}

		Point& operator= (const Point &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			x=rhs.x;
			y=rhs.y;
			lat=rhs.lat;
			lon=rhs.lon;
			}
		return *this;    // Return ref for multiple assignment
		}


		};

//struct Point{double x, y;};
class Point_with_attributes
		{
		public:
		double x, y, lat, lon;
		vector <double> att;

		Point_with_attributes() {x=y=lat=lon=BAD_DATA_FLOAT; att.clear();}

		Point_with_attributes deep_copy()
			{
			Point_with_attributes p;
			p.x=x;
			p.y=y;
			p.lat=lat;
			p.lon=lon;
			p.att=att;
			return(p);
			}

		void set (double _x, double _y)
			{
			x=_x;
			y=_y;
			}


		Point_with_attributes create(double Xx, double Xy, double Xlat, double Xlon)
			{
			Point_with_attributes p;
			p.x=Xx;
			p.y=Xy;
			p.lat=Xlat;
			p.lon=Xlon;
			return(p);
			}

		Point_with_attributes& operator= (const Point_with_attributes &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			x=rhs.x;
			y=rhs.y;
			lat=rhs.lat;
			lon=rhs.lon;
			att=rhs.att;
			}
		return *this;    // Return ref for multiple assignment
		}


		};


class Cobject2D_overlap
	{
	public:
	double id;
	double obj_size;
	double overlap_size;

	Cobject2D_overlap() {freememory();}
	~Cobject2D_overlap() {freememory();}

	void freememory()
		{
		id=-1;
		obj_size=0;
		overlap_size=0;
		}

	Cobject2D_overlap& operator= (const Cobject2D_overlap &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			//if (!rhs.check_if_initilized())
			//	error(AT,FU, "!rhs.check_if_initilized()");
			id=rhs.id;
			obj_size=rhs.obj_size;
			overlap_size=rhs.overlap_size;
			}
		return *this;    // Return ref for multiple assignment
		}

	// less_than_equal_according_to_id
	bool operator<(const Cobject2D_overlap &rhs) const
		{
		if (id < rhs.id) return(true);
		return(false);
		}

	string output()
		{
		ostringstream s1;
		s1.str("");
		s1 << id << ":" << overlap_size ;
		return(s1.str());
		}


	};

class Cobject2D_overlap_list
	{
	public:
	vector <Cobject2D_overlap> object2D_overlap;

	Cobject2D_overlap_list() {freememory();}
	~Cobject2D_overlap_list() {freememory();}

	void freememory()
		{
		object2D_overlap.clear();
		}

	Cobject2D_overlap_list& operator= (const Cobject2D_overlap_list &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			//if (!rhs.check_if_initilized())
			//	error(AT,FU, "!rhs.check_if_initilized()");
			object2D_overlap=rhs.object2D_overlap;
			}
		return *this;    // Return ref for multiple assignment
		}

	void add_overlap(double id, double obj_size, Cobject2D_overlap &a)
		{
 		bool found=false;

 		vector<Cobject2D_overlap>::iterator low;

		// we need to use a object to do the lower_bound search
		a.id=id;
		// check if allready here
 		if (object2D_overlap.size() > 0)
 			{
	 		low=lower_bound_which_never_returns_an_out_of_bounds_result(object2D_overlap.begin(), object2D_overlap.end(),a);
			if (low->id==id)
				found=true;
			}

 		// new object
		if (found==false)
 			{
			// add
			Cobject2D_overlap b;
			b.id=id;
			b.obj_size=obj_size;
			b.overlap_size=1;
			object2D_overlap.push_back(b);
			// sort by id - so that we can later use the very fast binary_search
			sort(object2D_overlap.begin(), object2D_overlap.end());
			}

 		// old object
		else
			low->overlap_size++;
		}

	bool get_object_with_biggest_overlap_size(Cobject2D_overlap &obj)
		{
		double max=0;
		for (long il=0; il < (long)object2D_overlap.size(); il++)
			if (object2D_overlap[il].overlap_size > max)
				{
				obj=object2D_overlap[il];
				max=object2D_overlap[il].overlap_size;
				}
		if (max > 0) return(true);
		return(false);
		}

	string output()
		{
		ostringstream s1;
		s1.str("");
		s1 << object2D_overlap.size() << ":\t";
		for (long il=0; il < (long)object2D_overlap.size(); il++ )
			s1 << object2D_overlap[il].output() << ", ";
		//s1 << endl;
		return(s1.str());
		}


	};





// This object is needed by identifying of 2D and 3D objects
class CIdentified_Object_2D
	{
	public:
	double id;
	int xmin;
	int xmax;
	int ymin;
	int ymax;

	CIdentified_Object_2D()
		{
		freememory();
		}

	void freememory()
		{
		id=-1;
		}

	CIdentified_Object_2D& operator= (const CIdentified_Object_2D &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			//if (!rhs.check_if_initilized())
			//	error(AT,FU, "!rhs.check_if_initilized()");
			id=rhs.id;
			xmin=rhs.xmin;
			xmax=rhs.xmax;
			ymin=rhs.ymin;
			ymax=rhs.ymax;
			}
		return *this;    // Return ref for multiple assignment
		}


	CIdentified_Object_2D deep_copy() const
		{
		CIdentified_Object_2D out;
		out=*this;
		return(out);
		}

	void set_object_borders(int Cxmin, int Cxmax, int Cymin, int Cymax)
		{
		xmin=Cxmin;
		xmax=Cxmax;
		ymin=Cymin;
		ymax=Cymax;
		}

	string output_ascii_string()
		{
		ostringstream s1;
		s1.str("");
		s1 << id << "\t" << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax ;
		return(s1.str());
		}


	};

// This object is needed by identifying of 2D and 3D objects
class CGroup_of_Identified_Objects_2D
	{
	public:
	vector <CIdentified_Object_2D> Identified_Object_2D;
	//double offset;

	CGroup_of_Identified_Objects_2D()
		{
		freememory();
		}

	~CGroup_of_Identified_Objects_2D()
		{
		freememory();
		}

	void freememory()
		{
		Identified_Object_2D.clear();
		//offset=-1;
		}

	CGroup_of_Identified_Objects_2D deep_copy()
		{
		CGroup_of_Identified_Objects_2D out;
		out.Identified_Object_2D = Identified_Object_2D;
		//out.offset = offset;
		return(out);
		}

	CGroup_of_Identified_Objects_2D& operator= (const CGroup_of_Identified_Objects_2D &rhs)
		{
		if (this != &rhs) // make sure not same object
			{
			freememory();                     // Delete old name's memory.
			//if (!rhs.check_if_initilized())
			//	error(AT,FU, "!rhs.check_if_initilized()");
			Identified_Object_2D = rhs.Identified_Object_2D;
			}
		return *this;    // Return ref for multiple assignment
		}


	bool check_if_object_is_defined(double id, CIdentified_Object_2D &matched_Identified_Object_2D, long &index)
		{
		long il;
		//cout << "Looking for" << "\t" << id  << endl;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id == id)
				{
				//cout << "aaaaY" << endl;
				matched_Identified_Object_2D=Identified_Object_2D[il];
				index=il;
				return(true);
				}
		//cout << "aaaaX" << endl;
		return(false);
		}

	bool check_if_object_is_defined_index_only(double id)
		{
		long il;
		//cout << "Looking for" << "\t" << id  << endl;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id == id)
				{
				//cout << "aaaaY" << endl;
				return(true);
				}
		//cout << "aaaaX" << endl;
		return(false);
		}



	bool change_id_of_object(double oldid, double newid)
		{
		long il;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id == oldid)
				{
				Identified_Object_2D[il].id=newid;
				//cout << "Changed id ofA " << oldid << "->" << newid << endl;
				return(true);
				}
		cout << "Changed id ofA2 " << oldid << "->" << newid << endl;
		return(false);
		}

	void join_two_Identified_Object_2D_objects(double oldid, double newid)
		{
		long index_slave,index_master;
		CIdentified_Object_2D slave, master;

		if (!check_if_object_is_defined(oldid, slave, index_slave))
			error(AT,FU, "oldid not found!");

		// check weather the master is found in the timestep
		if (check_if_object_is_defined(newid, master, index_master))
			{
			// if master is found expand the masters border so include the slave
			if (slave.xmax > master.xmax) master.xmax=slave.xmax;
			if (slave.xmin < master.xmin) master.xmin=slave.xmin;
			if (slave.ymax > master.ymax) master.ymax=slave.ymax;
			if (slave.ymin < master.ymin) master.ymin=slave.ymin;

			// put master in vector
			Identified_Object_2D[index_master].freememory();
			Identified_Object_2D[index_master] = master.deep_copy();

			// remove slave from vector
			Identified_Object_2D[index_slave].freememory();
			Identified_Object_2D[index_slave].id = -1;

			//cout << "Joined " << oldid << "->" << newid << endl;
			}

		else
			{
			// only change the id of slave to master
			Identified_Object_2D[index_slave].id = newid;
			//cout << "Changed id ofB " << oldid << "->" << newid << endl;
			}

		}

	double get_max_object_id()
		{
		long il;
		// count defined objects
		double max_id=-1;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id > max_id)
				max_id=Identified_Object_2D[il].id;
		if (max_id == -1)
			error(AT,FU, "max_id == -1!");
		return(max_id);
		}

	long count_number_of_objects()
		{
		long il;
		// count defined objects
		long count=0;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id != -1)
				count++;
		return(count);
		}

	long size()
		{
		return(Identified_Object_2D.size());
		}



/*		bool check_if_object_is_defined_non_compressed(long id, CIdentified_Object_2D &matched_Identified_Object_2D)
		{
		if ( id >  Identified_Object_2D.size()) return false;
		else if (Identified_Object_2D[id].id == -1) return false;
		else true;
		}*/

	//	remove nondefined objects from the vector
	void compress()
		{
		long il;
		// count defined objects
		long count=0;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id != -1)
				count++;

		CIdentified_Object_2D temp;
		vector <CIdentified_Object_2D> temp_vector;
		temp_vector.assign(count, temp);

		count=0;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id != -1)
				{
				temp_vector.at(count)=Identified_Object_2D.at(il).deep_copy();
				count++;
				}

		Identified_Object_2D=temp_vector;
		temp_vector.clear();
		}


	//	uncompresses objects so that obj index equals object id
	void uncompress()
		{
		long il;

		// check first that the vector is really compressed
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id == -1)
				error(AT,FU, "Identified_Object_2D[il].id == -1 !");

		// find max id
		double max_id=-1;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id > max_id)
				max_id=Identified_Object_2D[il].id;

		CIdentified_Object_2D temp;
		vector <CIdentified_Object_2D> temp_vector;
		temp_vector.assign((long)(max_id+1), temp);

		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			temp_vector.at((long)Identified_Object_2D[il].id) = Identified_Object_2D.at(il).deep_copy();

		Identified_Object_2D=temp_vector;
		temp_vector.clear();
		}


	//	get max id of objec in the vector
	double get_max_id_in_vector()
		{
		long il;
		// count defined objects
		double max_id=-1;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id > max_id)
				max_id = Identified_Object_2D[il].id;

		return(max_id);
		}


	void output_ascii()
		{
		long il;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			cout << il << "\t" << Identified_Object_2D[il].output_ascii_string() << endl;
		}

	void output_ascii_for_id(double id)
		{
		long il;
		for (il=0; il < (long)Identified_Object_2D.size(); il++)
			if (Identified_Object_2D[il].id == id)
			cout << il << "\t" << Identified_Object_2D[il].output_ascii_string() << endl;
		}



	};





