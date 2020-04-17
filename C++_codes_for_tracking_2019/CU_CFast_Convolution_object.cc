

class CFast_Convolution_grid_point
	{
	public:
	vector <long> long_add_points;
	vector <long> long_remove_points;
	long previous_id;
	long count;
	int Rx;
	int Ry;

	CFast_Convolution_grid_point()
		{
		freememory();
		}

	~CFast_Convolution_grid_point()
		{
		freememory();
		}

	void freememory()
		{
		long_add_points.clear();
		long_remove_points.clear();
		previous_id=-1;
		Rx=-1;
		Ry=-1;
		count=-1;
		}

	CFast_Convolution_grid_point deep_copy() const
		{
		CFast_Convolution_grid_point out;
		out.long_add_points=long_add_points;
		out.long_remove_points=long_remove_points;
		out.previous_id=previous_id;
		out.Rx=Rx;
		out.Ry=Ry;
		out.count=count;
		return(out);
		}

	void display_all_points_txt() const
		{
		cout << "previous_id:" << previous_id << endl;
		cout << "Rx:" << Rx << "\tRy:" << Ry << endl;
		cout << "count: " << count << endl;
		if (long_add_points.size() > 0) cout << "long_add_points:" << endl;
		for (long il=0; il < (long)long_add_points.size(); il++)
			cout << long_add_points[il] << endl;
		if (long_remove_points.size() > 0) cout << "long_remove_points:" << endl;
		for (long il=0; il < (long)long_remove_points.size(); il++)
			cout << long_remove_points[il] << endl;
		}


	};


class CFast_Convolution_object
	{
	public:
	vector <CFast_Convolution_grid_point> fcgp;

	~CFast_Convolution_object()
		{
		freememory();
		}

	void freememory()
		{
		fcgp.clear();
		}

	void create_grid_point_database(const CField &f, double convolution_radius, bool is_domain_periodic_in_x)
		{
		// free memory
		freememory();

		int ix,iy,dx,dy,iz;
		int previous_ix,previous_iy;
		int r_temp;
		long temp_ix;
		long il,ol,hl;
		long previous_id;
		vector <long> previous_vector;
		vector <long> current_vector;

		CFast_Convolution_grid_point empty_Fast_Convolution_grid_point;

		if (!f.check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		vector <bool> point_inventory(f.size(),false);

		for (il=0;il<f.size();il++)
			{
			// determine ix and iy
			f.DataLoc2XY(il, ix, iy);

			// add new point
			fcgp.push_back(empty_Fast_Convolution_grid_point.deep_copy());

			// determining Rx
			iz=ix;
			if (ix < f.dimx/2)
				while (fabs(f.lon_data[ix] -f.lon_data[iz+1]) <= convolution_radius)
					{iz++;}
			if (ix >= f.dimx/2)
				while (fabs(f.lon_data[ix] -f.lon_data[iz-1]) <= convolution_radius)
					{iz--;}
			fcgp[il].Rx=abs(ix-iz);

			// determining Ry
			iz=iy;
			if (iy < f.dimy/2)
				while (fabs(f.lat_data[iy] -f.lat_data[iz+1]) <= convolution_radius)
					{iz++;}
			if (iy >= f.dimy/2)
				while (fabs(f.lat_data[iy] -f.lat_data[iz-1]) <= convolution_radius)
					{iz--;}
			fcgp[il].Ry=abs(iy-iz);

			// define which previous point to use
			// first line - point on the left
			previous_id=il-1;
			// next lines - point on the top
			if (il >= f.dimx)
				previous_id=il-f.dimx;

			fcgp[il].previous_id=previous_id;

			// clear vectors
			previous_vector.clear();
			current_vector.clear();

			// determine previous_ix and previous_iy
			f.DataLoc2XY(previous_id, previous_ix, previous_iy);

			// populate previous vector
			if (il != 0)
				for (dy=-fcgp[previous_id].Ry;dy<=fcgp[previous_id].Ry;dy++)
					{
					if (fcgp[previous_id].Rx==0 || fcgp[previous_id].Ry == 0) r_temp=0;
					else r_temp=(int)floor((double)fcgp[previous_id].Rx* sqrt(1-(double)dy*(double)dy/((double)fcgp[previous_id].Ry*(double)fcgp[previous_id].Ry)));

					for (dx=-r_temp;dx<=r_temp;dx++)
						{
						// non periodic domain
						if ( !is_domain_periodic_in_x )
							{
							if (f.check_coordianes(previous_ix+dx,previous_iy+dy))
								previous_vector.push_back(f.XY2DataLoc(previous_ix+dx,previous_iy+dy));
							}

						// periodic domain in x
						else
							if (previous_iy+dy >= 0 && previous_iy+dy < f.dimy)
								{
								temp_ix=previous_ix+dx;
								if (temp_ix < 0) temp_ix=f.dimx+temp_ix;
								else if (temp_ix > f.dimx-1) temp_ix= temp_ix - f.dimx;
								previous_vector.push_back(f.XY2DataLoc(temp_ix,previous_iy+dy));
								}

						}
					}


			// populate current vector
			fcgp[il].count=0;
			for (dy=-fcgp[il].Ry;dy<=fcgp[il].Ry;dy++)
				{
				if (fcgp[il].Rx==0 || fcgp[il].Ry == 0) r_temp=0;
				else r_temp=(int)floor((double)fcgp[il].Rx* sqrt(1-(double)dy*(double)dy/((double)fcgp[il].Ry*(double)fcgp[il].Ry)));

				for (dx=-r_temp;dx<=r_temp;dx++)
					{
					fcgp[il].count++;
					// non periodic domain
					if ( !is_domain_periodic_in_x )
						{
						if (f.check_coordianes(ix+dx,iy+dy))
							current_vector.push_back(f.XY2DataLoc(ix+dx,iy+dy));
						}

						// periodic domain in x
						else
							if (iy+dy >= 0 && iy+dy < f.dimy)
								{
								temp_ix=ix+dx;
								if (temp_ix < 0) temp_ix=f.dimx+temp_ix;
								else if (temp_ix > f.dimx-1) temp_ix= temp_ix - f.dimx;
								current_vector.push_back(f.XY2DataLoc(temp_ix,iy+dy));
								}

					}
				}




			// --------------------------------------------------
			// determine add/remove points
			// --------------------------------------------------

			// determine remove points - the points that are present in previous grid_point but not in the cuurent grid_point
			// populate point_inventory with current vector
			for (hl=0; hl < (long)current_vector.size(); hl++)
				point_inventory[current_vector[hl]]=true;
			for (ol=0; ol < (long)previous_vector.size(); ol++)
				if (point_inventory[previous_vector[ol]] == false)
					fcgp[il].long_remove_points.push_back(previous_vector[ol]);
			// depopulate point_inventory
			for (hl=0; hl < (long)current_vector.size(); hl++)
				point_inventory[current_vector[hl]]=false;

			// determine add points - the points that are present in cuurent grid_point but not in the previous grid_point
			for (hl=0; hl < (long)previous_vector.size(); hl++)
				point_inventory[previous_vector[hl]]=true;
			for (ol=0; ol < (long)current_vector.size(); ol++)
				if (point_inventory[current_vector[ol]] == false)
					fcgp[il].long_add_points.push_back(current_vector[ol]);
			// depopulate point_inventory
			for (hl=0; hl < (long)previous_vector.size(); hl++)
				point_inventory[previous_vector[hl]]=false;





			/*cout << "---------------------------------" << endl;
			cout << "point: " << il << endl;
			cout << "previous: " << endl;
			for (ol=0;ol<previous_vector.size();ol++)
				cout << previous_vector[ol] << endl;
			cout << "current: " << endl;
			for (ol=0;ol<current_vector.size();ol++)
				cout << current_vector[ol] << endl;

 			fcgp[il].display_all_points_txt();

			if (il==1441) exit(1);

			*/
			if (il % ((long)f.size()/100) == 0 )
				cout << (100*il) / f.size() << " % convolution initialization complete" << endl;
			}



		}

	CField perform_convolution( const CField &f)  const
		{
		CField f_conv, f_copy;
		long il,ol;


		if (!f.check_if_initilized())
			error(AT,FU, "!check_if_initilized()");


		// replace all missig values with 0 mm/3m
		f_copy=f.create_a_copy_with_data();
		f_copy.Replace_values(BAD_DATA_FLOAT,0);

		// create f_conv and set al values as 0 mm/3h
		f_conv=f.create_a_copy_with_data();
		f_conv.Set_all_values(0);

		double remove_points_contribution;
		double add_points_contribution;

		for (il=0;il<f.size();il++)
			{
			// start from previous value
			if (il == 0) f_conv.data[il]=0;
			else f_conv.data[il]=f_conv.data[fcgp[il].previous_id]*(double)fcgp[fcgp[il].previous_id].count/(double)fcgp[il].count;

			// removed points
			remove_points_contribution=0;
			for (ol=0; ol < (long)fcgp[il].long_remove_points.size(); ol++ )
				if (fcgp[il].long_remove_points[ol] != BAD_DATA_FLOAT)
					remove_points_contribution+=f_copy.data[fcgp[il].long_remove_points[ol]];
			if ((long)fcgp[il].long_remove_points.size() > 0)


			// added points
			add_points_contribution=0;
			for (ol=0; ol < (long)fcgp[il].long_add_points.size(); ol++ )
				if (fcgp[il].long_add_points[ol] != BAD_DATA_FLOAT)
					add_points_contribution+=f_copy.data[fcgp[il].long_add_points[ol]];
			if (fcgp[il].long_add_points.size() > 0)

			// final value
			f_conv.data[il]+=add_points_contribution/(double)fcgp[il].count - remove_points_contribution/(double)fcgp[il].count;

			// correct final value if very close to 0mm but rounded a but wrong
			if (f_conv.data[il] < 0.000001 &&  f_conv.data[il] > -0.000001)
				f_conv.data[il]=0;


			/*cout << "Point: " << il << endl;
			fcgp[il].display_all_points_txt();
			cout << "----------------------------------" << endl;


			if (il == 1441) exit(1);*/
			}

		// restore missing values
		for (il=0;il<f.size();il++)
			if (f.data[il]==BAD_DATA_FLOAT)
				f_conv.data[il]=BAD_DATA_FLOAT;

		// freememory
		f_copy.freememory();

		return(f_conv);
		}
	};

