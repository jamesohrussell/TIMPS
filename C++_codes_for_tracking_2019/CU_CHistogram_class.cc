// ******************************************************
// ******************************************************
// Class CField. Stores one 2-D field with dimension atributes and raw data.
//
class CHistogram {
	public:
	vector <double> levels;
	vector <double> sum;
	string name;

	// constructor
	CHistogram() {freememory();}
	~CHistogram() {freememory();}

	void freememory()
		{
		levels.clear();
		sum.clear();
		name="";
		}

	CHistogram deep_copy()  const
		{
		CHistogram out;
		out.levels=levels;
		out.sum=sum;
		return(out);
		}


	void create_equidistant_levels(double start,double end, double interval_size, string Xname)
		{
		double l;

		freememory();

		name=Xname;

		l=start;
		while (l <= end )
			{
			levels.push_back(l);
			l+=interval_size;
			//cout << l << endl;
			}

		// initialize sum
		sum.assign(levels.size(),0);
		}

	void create_levels_from_vector(	vector <double> &levelsX)
		{
		levels=levelsX;
		// initialize sum
		sum.assign(levels.size(),0);
		}


	void add_item(double level_value, double sum_value)
		{
		long il;


		if (level_value < levels.front() ) sum.front()+=sum_value;
		else if (level_value >= levels.back() ) sum.back()+=sum_value;
		else
			for (il=0; il < (long)levels.size() - 1 ; il++)
				if (level_value >= levels[il] && level_value < levels[il+1])
					{
					//cout << level_value << "\t" << il << endl;
					sum[il]+=sum_value;
					}

		//cout << level_value << "\t" << il << endl;
		//cout << levels.front() << "\t" << levels.back() << endl;
		}

	void add_item_with_sum_one( double level_value)
		{
		add_item(level_value,1);
		}

	void add_vector( const vector <double> &level_value,  const vector <double> &sum_value)
		{
		long il;

		if (level_value.size() != sum_value.size() )
			error(AT,FU, "level_value.size() != sum_value.size()");

		for (il=0; il < (long)sum_value.size() ; il++)
			add_item(level_value[il],sum_value[il]);
		}

	void add_vector_with_sum_one( const vector <double> &level_value)
		{
		long il;

		for (il=0; il < (long)level_value.size() ; il++)
			add_item(level_value[il],1);
		}

	void display() const
		{
		long il;
		cout << "Histogram: " << name << endl;
		for (il=0; il < (long)levels.size() ; il++)
			cout << levels[il] << "\t" << sum[il] << endl;
		}

	string output() const
		{
		ostringstream s1;
		long il;
		s1 << "Histogram: " << name << endl;
		for (il=0; il < (long)levels.size() ; il++)
			s1 << levels[il] << "\t" << sum[il] << endl;
		return(s1.str());
		}


	double cumulative_probability(long indeks)
		{
		if (indeks > (long) sum.size() - 1 )
			error(AT,FU, "indeks > (long) sum.size() - 1");

		double sumcum=0;
		for (long il=0; il <= indeks ; il++)
			sumcum+=sum[il];
		return(sumcum/sum_vector(sum));
		}


	};
