

// ******************************************************
// ******************************************************
// Classes and function used for testing if point is inside polygon
//


class Closed_Polygon
	{
	public:
	long number_of_points;
	Point* V;

	Closed_Polygon()
		{
		V=NULL;
		}

	void allocatememory()
		{
		V = new Point [number_of_points+1];
		}

	void freememory()
		{
		delete[] V;
		V=NULL;
		}


	};

class Group_of_Polygons
	{
	public:
	long number_of_polygons;
	Closed_Polygon* cp;

	Group_of_Polygons()
		{
		cp=NULL;
		}

	void allocatememory()
		{
		cp = new Closed_Polygon [number_of_polygons];
		}

	void freememory()
		{
		delete[] cp;
		cp=NULL;
		}


	};

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int cn_PnPoly(  const Point &P,  const Point* V, int n )
	{
    int    cn = 0;    // the crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i] to V[i+1]
       if (((V[i].y <= P.y) && (V[i+1].y > P.y))    // an upward crossing
        || ((V[i].y > P.y) && (V[i+1].y <= P.y))) { // a downward crossing
            // compute the actual edge-ray intersect x-coordinate
            float vt = (float)(P.y - V[i].y) / (V[i+1].y - V[i].y);
            if (P.x < V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn&1);    // 0 if even (out), and 1 if odd (in)
	}


// the same as cn_PnPoly, but for a group of closed polygons
int cn_GroupPoly(  const Point &P,  const Group_of_Polygons &gp)
	{
    int cn=0;
    int ix;

    for (ix=0; ix < gp.number_of_polygons; ix++)
		if ( cn_PnPoly( P, gp.cp[ix].V, gp.cp[ix].number_of_points-1 ) == 1 ) cn=1;

    return (cn);    // 0 if even (out), and 1 if odd (in)
	}


Group_of_Polygons ReadPolygonFromFile(string fname)
	{
	Group_of_Polygons gp;
	FILE * pFile;
	pFile = fopen (fname.c_str(),"r");
	if (!pFile) {
			cout << "Unable to open file";
			exit(1); // terminate with error
		}

	int ix;
	int objectcount=0;
	int point_number=0;
	float x,y;
	int object_old=-9999, object_new;
	long tocka;

	// counts all polygons
	while( fscanf( pFile,"%d,%ld,%f,%f", &object_new, &tocka, &x, &y) !=EOF )
		{
		if (object_new != object_old) objectcount++;
		object_old=object_new;
		}

	// Group_of_Polygons class and allocate memory
	gp.number_of_polygons=objectcount;
	gp.allocatememory();

	// set number of vertexes of each Closed_Polygon to zero
	for (ix=0; ix < objectcount; ix++)
		gp.cp[ix].number_of_points=0;

	// counts number of vertexes for each polygon
	rewind (pFile);
	objectcount=0;
	object_old=-9999;
	while( fscanf( pFile,"%d,%ld,%f,%f", &object_new, &tocka, &x, &y) !=EOF )
		{
		if (object_new != object_old) objectcount++;
		object_old=object_new;
		gp.cp[objectcount-1].number_of_points++;
		}

	// adds another vertex since the last vertex has to be identical to the first
	for (ix=0; ix < objectcount; ix++)
		gp.cp[ix].number_of_points++;

	// reservers memory for polygon vertexes
	for (ix=0; ix < objectcount; ix++)
		gp.cp[ix].allocatememory();

	// reads the x,y data
	rewind (pFile);
	objectcount=0;
	object_old=-9999;
	while( fscanf( pFile,"%d,%ld,%f,%f", &object_new, &tocka, &x, &y) !=EOF )
		{
		if (object_new != object_old) {objectcount++;point_number=0;}
		object_old=object_new;

		// !!!!!  EDIT HOW LONGITUDE IS TREATED!!!
		gp.cp[objectcount-1].V[point_number].x=x;
//		if (x<0) gp.cp[objectcount-1].V[point_number].x=360+x;
//		else gp.cp[objectcount-1].V[point_number].x=x;

		gp.cp[objectcount-1].V[point_number].y=y;
		point_number++;
		}

	// sets the values for additional vertex since the last vertex has to be identical to the first
	for (ix=0; ix < objectcount; ix++)
		{
		gp.cp[ix].V[gp.cp[ix].number_of_points-1].x=gp.cp[ix].V[0].x;
		gp.cp[ix].V[gp.cp[ix].number_of_points-1].y=gp.cp[ix].V[0].y;
		}

/*
	// Modifies longitude of vertexes of each object so they are on the same plane together.
	// For example: if object spans across lon=0, the point_in_polygon routine will have problems, since some vertexes will be defined as lat="close to 360" and some lat="close to zero" - the object will WRONGLY seem to span accros all the world
	// The loop remaps each vertex on the same plane - either close to lat=0 or close to lat=360 - the plane is defined by the first vertex
	for (ix =0 ; ix < gp.number_of_polygons; ix++)
		for (iy =0 ; iy < gp.cp[ix].number_of_points; iy++)
		{
		// if they span more than 180 degress
		if (fabs(gp.cp[ix].V[0].x - gp.cp[ix].V[iy].x) > 180)
			{
			if (gp.cp[ix].V[0].x > 180) gp.cp[ix].V[iy].x = gp.cp[ix].V[iy].x + 360;
			else gp.cp[ix].V[iy].x = gp.cp[ix].V[iy].x - 360;
			}
		}
*/

	fclose (pFile);
	return(gp);
	}
/*
	// Mask Polygon - Creates Cfield where all points inside Group_of_Polygons are decleared missing values and 0 elsewhere
	// Can be used in conjunction with "add" to mask out the selected areas as missing values
	CField Create_Mask_Field_from_Polygon(Group_of_Polygons gp, CField ftemp)
		{
		CField f;
		long ixl;
		float lon,lat;
		Point P1,P2,P3;
		int ix,iy;

		f=ftemp.create_a_copy();
		f.allocatememory();

	for (iy=0;iy<dimy;iy++)
		{
		cout << "line: "<< iy << "\n";
		for (ix=0;ix<dimx;ix++)
			{
			XY2LatLon(ix, iy, lon, lat);
			ixl=XY2DataLoc(ix, iy);

			P1.y=P2.y=P3.y=lat;
			// checks three possible lon - so no objects are missed
			P1.x=lon;
			P2.x=lon+360;
			P3.x=lon-360;
			if (cn_GroupPoly( P1, gp) != 0) f.data[ixl] = BAD_DATA_FLOAT;
			else if (cn_GroupPoly( P2, gp) != 0) f.data[ixl] = BAD_DATA_FLOAT;
			else if (cn_GroupPoly( P3, gp) != 0) f.data[ixl] = BAD_DATA_FLOAT;
			else f.data[ixl]= 0;
			}
		}
	return (f);
	}
*/

	// Point Overlay Polygon - Creates Cfield where all points coresponding to lat, lon in the polzgon datafile are marked as 1. Everzthing else is marked as 0
	// Can be used to draw borders in PNG files
	// CField template is given as a source of start_lat and start_lon parameteres...
	CField Point_Overlay_Polygon( const Group_of_Polygons &gp,  const CField &ftemp, int line_width)
		{
		CField f;
		Point p, p2;
		long ipolygon, ipoint;


		if (!ftemp.check_if_initilized())
			{cout << "ERROR: Point_Overlay_Polygon ERROR: !ftemp.check_if_initilized()""\n";exit(1);}

		f=ftemp.create_a_copy_with_data();
		f.Set_all_values(0);



//		cout << p.lat << " " << p.lon << endl;

		// go over all points
		for (ipolygon=0;ipolygon < gp.number_of_polygons; ipolygon++)
			for (ipoint=0;ipoint < gp.cp[ipolygon].number_of_points; ipoint++)
				{
				for (p.lon=gp.cp[ipolygon].V[ipoint].x - 360; p.lon <= gp.cp[ipolygon].V[ipoint].x + 360 + 0.01; p.lon+=360)
					{
					p.lat=gp.cp[ipolygon].V[ipoint].y;

					//cout << "qqq:  " << ipolygon << " " << ipoint << " " << p.lat << " " << p.lon  << endl;


					//cout << p.lat << "\t" << f.lat_start << "\t" << << p.lon << endl;


					// draw a point
					if (p.lon >= ftemp.lon_data[0] && p.lon <= ftemp.lon_data[ftemp.dimx-1]
					 && p.lat >= ftemp.lat_data[0] && p.lat <= ftemp.lat_data[ftemp.dimy-1] )
					 	{
						//cout << p.lat << " " << p.lon << endl;
						f.Point_LatLon2XY(p);
						//cout << "ccc" << endl;

						// draw a line
						if (ipoint!=0)
							{
							p2.lon=gp.cp[ipolygon].V[ipoint-1].x - gp.cp[ipolygon].V[ipoint].x + p.lon;
							p2.lat=gp.cp[ipolygon].V[ipoint-1].y;

							if (p2.lon >= ftemp.lon_data[0] && p2.lon <= ftemp.lon_data[ftemp.dimx-1]
							 && p2.lat >= ftemp.lat_data[0] && p2.lat <= ftemp.lat_data[ftemp.dimy-1] )
								{
								//cout << p.lat << " " << p.lon  << " " << p.x << " " << p.y << endl;
								f.Point_LatLon2XY(p2);
								if (p2.x >= 0 && p2.x <= (float)f.dimx - 1 && p2.y >= 0 && p2.y <= (float)f.dimy - 1 )
									{
									f.lineImproved((int)p.x, (int)p.y, (int)p2.x, (int)p2.y, 1, line_width);
									}

								}

							}


						}



					}
				}

		return (f);
		}



