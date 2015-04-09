/*
 
 LICENSE AGREEMENT
 
 This software is released under the MIT license
 
 Copyright (C) 2012
 
 James Rosindell (Imperial College London, University of Leeds)
 Luke Harmon (University of Idaho)
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 Users of the software are kindly requested to cite the related paper ("A unified model of species immigration, extinction and abundance on islands" Rosindell & Harmon, Journal of Biogeography - publication date 2012 or 2013)
 
 Thank you for your interest in our software
 
 */

/*
 
 INSTUCTIONS FOR USE
 
 This software simulates dynamics on an island under ecological neutral theory
 there is a choice of metacommunity distrubtion - uniform or read in from a file - two example files are included with the code.
 This version of the code is for the *protracted immigration* model - the basic immigration and clustered immigration models are simulated with code available separately.
 Instructions on how to edit and run the code are found in the 'main' routine below.
 Use will require a basic knowledge of C++
 
 You will need a random number generator.
 The attached file uses a standard
 choice included in C++ but has instructions on
 how to substitute this for your own preference.
 You are very strongly advised to do this because a large
 proportion of the computations executed by this code will
 be calls for random numbers and both speed and
 reliability of the generator are thus extremely important.
 We recommend numerical recipes 'ran2' routine which you
 will find available on line
 http://www.nr.com/oldverswitcher.html
 We use this routine for all our simulations and recommend it
 very highly but we do not attempt to distribute this code
 as it is written by others
 
 */



/************************************************************
 INCLUDES
 ************************************************************/

// standard inludes
# include <stdio.h>
# include <fstream>
# include <vector>
# include <iostream>
# include <string>
# include <math.h>
# include <time.h>
# include <ctime>

using namespace std;

/************************************************************
 RANDOM NUMBER GENERATOR WRAPPER OBJECT
 ************************************************************/

# include "RandomWrap.cpp"

class NRrand {
	
private:
	// an object that contains the random number generator of choice
	RandomWrap X;
	// the last result (for normal deviates)
	double lastresult;
	// when doing normal deviates and values are in pairs
	// true when a new pair is needed, false when lastresult can be used
	bool normflag;
	
public:
	
	NRrand()
	{
		normflag = true;
	}
	
	void set_seed(long seed)
	{
		X.set_seed(seed);
	}
	
	double d01()
	{
		return X.d01();
	}  
	
	long i0(long max)
	// integer between 0 and max inclusive
	{
		return (long(floor(d01()*(max+1))));
	}
	
	double norm()
	// normal deviate
	{
		if (normflag)
		{
			double r2 = 2;
			double xx;
			double yy;
			while (r2>1)
			{
				xx=2.0*d01()-1.0;
				yy=2.0*d01()-1.0;
				r2=(xx*xx)+(yy*yy);
			}
			double fac=sqrt(-2.0*log(r2)/r2);
			lastresult = xx*fac;
			double result = yy*fac;
			normflag = false;
			return result;
		}
		else
		{
			normflag = true;
			return lastresult;
		}
	}
	
	double fattail(double z)
	// 2D fat tailed deviate
	{
		double result;
		result = pow((pow(d01(),(1.0/(1.0-z)))-1.0),0.5);
		return result;
	}
	
	double direction()
	// direction (for 2D fat tailed deviate)
	{
		double xx = 1.0 , yy = 1.0;
		while (xx*xx+yy*yy>1.0)
		{
			xx = d01();
			yy = d01();
		}
		return pow((xx/yy),2.0);
	}
	// to reconstruct distribution, use x = fattail/squrt(1+direction) , y = fattail/squrt(1+(direction^-1))
	
	bool event(double probin)
	{
		if (probin < 0.000001)
		{
			if (X.d01() <= 0.000001)
			{
				return (event(probin * 1000000.0));
			}
			else 
			{
				return false;
			}
		}
		else 
		{
			if (probin > 0.999999)
			{
				return (!(event(1.0-probin)));
			}
			else
			{
				return (X.d01() <= probin);
			}
		}
	}
	
};

class metacommunity
{
	
private:
	
	vector<double> probdist;
	long maxab;
	
public:
	
	metacommunity() 
	{
		reset();
	}
	
	void reset()
	{
		probdist.clear();
		probdist.push_back(0);
		maxab = -1;
	}
	
	long numspec()
	{
		return probdist.size()-1;
		// the [0] entry is reserved to hold a 0;
	}
	
	double get_prob(long i)
	{
		if ((i >= probdist.size())||(i <= 0))
		{
			return 0;
		}
		else 
		{
			return(probdist[i]-probdist[i-1]);
		}
	}
	
	long get_maxab()
	{
		return maxab;
	}
	
	void calc_max()
	{
		double maxprob = 0;
		for (long i = 0 ; i < probdist.size() ; i ++)
		{
			if( get_prob(i) > maxprob)
			{
				maxprob = get_prob(i);
				maxab = i;
			}
		}
	}
	
	long randspec(double randin)
	{
		// returns a random species based on a random input in (0,1)
		// we are looking in the list for the index of the number which is greater than randin, whilst the number with an index one less is less than randin
		// use interval bisection for this.
		long min = 0;
		long max = probdist.size() -1;
		long med = long(floor((double(min)+double(max))/2.0));
		while ((max-min)>1)
		{
			if (probdist[med] >= randin)
			{
				max = med;
			}
			else 
			{
				min = med;
			}
			med = long(floor((double(min)+double(max))/2.0));
		}
		return max;
	}
	
	void readinDLS()
	{
		// this routine reads in a metacommunity file.
		vector<long> abundances;
		abundances.clear();
		ifstream in;
		in.open("metacommunity_DLS3.txt");
		while (!in.eof())
		{
			long templong;
			in >> templong;
			if (!in.eof())
			{
				abundances.push_back(templong);
			}
		}
		
		in.close();
		double totind = 0.0;
		for (long i = 0 ; i < abundances.size() ; i ++)
		{
			totind += double(abundances[i]);
		}
		
		double gekko = 0.0;
		
		double totsofar = 0.0;
		for (long i = 0 ; i < abundances.size() ; i ++)
		{
			totsofar = totsofar+(double(abundances[i])/totind);
			probdist.push_back(totsofar);
		}
		calc_max();
	}
	
	void readinLS()
	{
		// this routine reads in a metacommunity file.
		vector<long> abundances;
		abundances.clear();
		ifstream in;
		in.open("metacommunity_LS4.txt");
		while (!in.eof())
		{
			long templong;
			in >> templong;
			if (!in.eof())
			{
				abundances.push_back(templong);
			}
		}
		in.close();
		double totind = 0.0;
		for (long i = 0 ; i < abundances.size() ; i ++)
		{
			totind += double(abundances[i]);
		}
		double totsofar = 0.0;
		for (long i = 0 ; i < abundances.size() ; i ++)
		{
			totsofar = totsofar+(double(abundances[i])/totind);
			probdist.push_back(totsofar);
		}
		calc_max();
	}
	
	void setup(short flag)
	{
		reset();
		if (flag == 1)
		{
			readinDLS();
			// difference log series
		}
		if (flag == 2)
		{
			readinLS();
			// log series
		}
		if (flag == 3)
		{
			// uniform distribution: 1000 species each with abundance 1000,000
			for (int i = 1 ; i <= 1000 ; i ++)
			{
				probdist.push_back(double(i)*0.001);
			}
			maxab = 1;
		}
	}
	
};


class dynamicsim
{
	
private:
	
	// variables of the simulation
	double immigration;
	// immigration rate
	long local;
	// local community size
	short meta_type;
	// metacommunity type	
	double immtau;
	// the number of generations a species is requied to be present for before it is counted
	
	// simulation parameters
	vector<int> data;
	// the simulation data - individuals
	vector<bool> original;
	// true for each inidividual if it is a descendent of origial island population
	long numtrue;
	// the number of values that are 'true' in the above object
	// this is a tool to check to see when equilibrium has been reached.
	int numcycles;
	// the number of complete turnovers of all individuals
	// not generations but total extinction of all lineages and replacement from mainland
	
	vector<long> abundance;
	vector<long> abundance_inoct;
	// stores a list of the abundances of all the species in the system (meta.getsize) - the second vector stores in octaves
	vector<double> numgenspec;
	// stores the number of generations that this particular species has been on the island.
	// generations are meaured in terms of births and deaths of that species - not overall births and deaths.
	// we will only record the abundanecs of species fully established on the island 
	//(others are ignored for abundances aswell as for immigration and extinction rates)
	vector<bool> fullspec;
	// is true when the species is considered to have been present long enough to be classed as really there
		
	long richness_tot;
	// the total richness in the system
	
	metacommunity meta;
	
	vector<double> extrate_av;
	vector<double> immrate_av;
	vector<double> numread;
	// this stores the average extinction and immigration rates as a function of richness
	vector< vector<double> > abund_vs_rich;
	// this stores the average numbers in each abundance class rates as a function of richness
	
	double gen;
	long numsteps;
	// the number of generations that have passed.
	
	// RANDOM NUMBER GENERATOR
	long seed;
	NRrand NR;
	bool seeded;
	// variables todo with the random number generator.
	
	bool debug;
	// debug variable
	
	
public:
	
	dynamicsim() 
	{ 
		reset(); 
	}
	
	void reset() 
	{
		data.clear();
		original.clear();
		numtrue = -1;
		numcycles = 0;
		abundance.clear();
		abundance_inoct.clear();
		numgenspec.clear();
		fullspec.clear();
		richness_tot = 0;
		
		extrate_av.clear();
		immrate_av.clear();
		numread.clear();
		abund_vs_rich.clear();
		
		numsteps = 0;
		gen = 0.0;
		seeded = false;
		debug = false;		
	}
	
	void set_seed(long seedin)
	{
		if (!seeded)
		{
			seed = seedin;
			seeded = true;
			NR.set_seed(seed);
		}
	}
	
	void debugon()
	{
		debug = true;
	}
	
	int octave_sort(long ab_in)
	{
		int result;
		if(ab_in <= 0)
		{
			result = 0;
		}
		else 
		{
			long min = 1;
			long max = 2;
			result = 1;
			while(!((ab_in < max)&&(ab_in >= min)))
			{
				min = min*2;
				max = max*2;
				result ++;
			}
		}
		return result;
	}
	
	void resetup(	double immigrationin , long localin , short meta_typein , bool highrich, double immtauin)
	{
		reset();
		
		immigration = immigrationin;
		local = localin;
		meta_type = meta_typein;
		immtau = immtauin;
		
		meta.setup(meta_type);
		
		vector<double> tempclass;
		tempclass.clear();
		int tempsize;
		tempsize = floor(log(double(localin))/log(2.0));
		for (long i = 0 ; i < tempsize+5 ; i ++)
		{
			tempclass.push_back(0.0);
			abundance_inoct.push_back(0.0);
		}
		for (long i = 0 ; i < meta.numspec()+1 ; i ++)
		{
			extrate_av.push_back(0.0);
			immrate_av.push_back(0.0);
			numread.push_back(0.0);
			abund_vs_rich.push_back(tempclass);
		}
		
		gen = 0;
		numcycles = 0;
				
		if (highrich)
		{
			for (int i = 0 ; i < meta.numspec()+1 ; i ++)
			{
				abundance.push_back(0);
				numgenspec.push_back(0);
				fullspec.push_back(false);
			}
			for (int i = 0 ; i < local ; i ++)
			{
				data.push_back(meta.randspec(NR.d01()));
				original.push_back(true);
				abundance[data[data.size()-1]] ++;
			}
			richness_tot = 0;
			for (int i = 0 ; i < meta.numspec()+1 ; i ++)
			{
				if (abundance[i] > 0)
				{
					richness_tot ++;
					abundance_inoct[octave_sort(abundance[i])] ++;
					numgenspec[i] = (immtauin*2.0);
					fullspec[i] = true;
					// initially all species are assumed to have been present for sufficiently long
				}
			}
			numtrue = local;
		}
		else 
		{
			for (int i = 0 ; i < local ; i ++)
			{
				data.push_back(meta.get_maxab());
				original.push_back(true);
			}
			for (int i = 0 ; i < meta.numspec()+1 ; i ++)
			{
				abundance.push_back(0);
				numgenspec.push_back(0);
				fullspec.push_back(false);
			}
			abundance[meta.get_maxab()] = local;
			numgenspec[meta.get_maxab()] = (immtauin*2.0); 
			fullspec[meta.get_maxab()] = true;
			abundance_inoct[octave_sort(abundance[meta.get_maxab()])] ++;
			numtrue = local;
			richness_tot = 1;
			
		}
		
	}
	
	
	void step()
	{
		// increment number of steps
		numsteps ++;
		gen += 1.0/local;
						
		// now choose individual to die
		long chosen;
		chosen = NR.i0(local-1);

		// record the important information 
		long richness_tot_temp;
		richness_tot_temp = richness_tot; // store the richness at the start
		numread[richness_tot_temp] += 1.0;  // increment the number of readings
		for (int i = 1 ; i < abundance_inoct.size() ; i ++)
		{	
			abund_vs_rich[richness_tot_temp][i] += abundance_inoct[i]; 
		}
		
		// update for death
		long spec_infocus = data[chosen];
		numgenspec[spec_infocus] += (0.5/(abundance[spec_infocus]));
		if (((!fullspec[spec_infocus])&&(numgenspec[spec_infocus] > immtau))&&(abundance[spec_infocus] >1))
		{
			fullspec[spec_infocus] = true;
			immrate_av[richness_tot_temp] += 1.0;
			richness_tot ++;
			abundance[spec_infocus] --;
			abundance_inoct[octave_sort(abundance[spec_infocus])] ++;
		}
		else 
		{
			if (fullspec[spec_infocus])
			{
				abundance_inoct[octave_sort(abundance[spec_infocus])] --;
				abundance[spec_infocus] --;
				abundance_inoct[octave_sort(abundance[spec_infocus])] ++;
            }
			else 
			{
				abundance[spec_infocus] --;
			}
		}
		
		if (original[chosen])
		{
			numtrue --;
		}
			
		
		
		if (NR.d01() < immigration)
		{
			// immigration
			long toimm;
			toimm = meta.randspec(NR.d01());
			data[chosen] = toimm;
			if (abundance[toimm] == 0)
			{
				numgenspec[toimm] += 0.5;
			}
			else
			{
				numgenspec[toimm] += (0.5/(abundance[toimm]));
			}
								
			if ((!fullspec[toimm])&&(numgenspec[toimm] > immtau))
			{
				fullspec[toimm] = true;
				immrate_av[richness_tot_temp] += 1.0;
				richness_tot ++;
				abundance[toimm] ++;
				abundance_inoct[octave_sort(abundance[toimm])] ++;
			}
			else 
			{
				if (fullspec[toimm])
				{
					abundance_inoct[octave_sort(abundance[toimm])] --;
					abundance[toimm] ++;
					abundance_inoct[octave_sort(abundance[toimm])] ++;
				}
				else 
				{
					abundance[toimm] ++;
				}
			}
			
			original[chosen] = false;
						
		}
		else
		{
			// birth
			long torep;
			torep = chosen;
			while(torep == chosen)
			{
				torep = NR.i0(local-1);
			}
			data[chosen] = data[torep];
			numgenspec[data[torep]] += (0.5/(abundance[data[torep]]));
			// checker to add here										   
			original[chosen] = original[torep];
			
			if (original[chosen])
			{
				numtrue ++;
			}
													   
			if ((!fullspec[data[torep]])&&(numgenspec[data[torep]] > immtau))
			{
				fullspec[data[torep]] = true;
				immrate_av[richness_tot_temp] += 1.0;
				richness_tot ++;
				abundance[data[torep]] ++;
				abundance_inoct[octave_sort(abundance[data[torep]])] ++;
			}
			else 
			{
				if (fullspec[data[torep]])
				{
					abundance_inoct[octave_sort(abundance[data[torep]])] --;
					abundance[data[torep]] ++;
					abundance_inoct[octave_sort(abundance[data[torep]])] ++;
				}
				else 
				{
					abundance[data[torep]] ++;
				}
			}
			
			
			
		}
		
		if (abundance[spec_infocus] == 0)
		{
			if (fullspec[spec_infocus])
			{
				richness_tot --;
				extrate_av[richness_tot_temp] += 1.0;
			}
			numgenspec[spec_infocus] = 0;
			fullspec[spec_infocus] = false;
		} 
		
		
		// manage cycle count
		if (numtrue == 0)
		{
			numcycles ++;
			for (int i = 0 ; i < original.size() ; i ++)
			{
				original[i] = true;
				numtrue = local;
			}
		}
	}
	
	
	void sim(double immigrationin , long localin , short meta_typein , double timeallowed , char* file_out , char* file_log , double immtauin)
	{
		ofstream out;
		out.open( file_log );
		out << " simulating immigrationin = " << immigrationin << " , localin = " << localin << " , meta_typein = " << meta_typein << " , immtau = " << immtauin << " , timeallowed = " << timeallowed << " , \n";
		cout << " simulating immigrationin = " << immigrationin << " , localin = " << localin << " , meta_typein = " << meta_typein << " , immtau = " << immtauin << " , timeallowed = " << timeallowed << " , \n";
		out.close();
		time_t start , end;
		time(&start);
		
		// vars for M&W dynaics
		
		double gen2e = 0;
		double gen2eread = 0;
		double gen2e2 = 0;
		
		vector<double> numreadX_pre_full;
		vector<double> numreadX_pre_fullB;
		vector<double> numreadX_post_full;
		
		vector<double> immrateX_pre_full;
		vector<double> immrateX_pre_fullB;
		vector<double> immrateX_post_full;
		
		vector<double> extrateX_pre_full;
		vector<double> extrateX_pre_fullB;
		vector<double> extrateX_post_full;
		
		vector<double> immrateX_pre_full2;
		vector<double> immrateX_pre_fullB2;
		vector<double> immrateX_post_full2;
		
		vector<double> extrateX_pre_full2;
		vector<double> extrateX_pre_fullB2;
		vector<double> extrateX_post_full2;
		
		vector< vector<double> > abvsri_pre;
		vector< vector<double> > abvsri_preB;
		vector< vector<double> > abvsri_post;
		
		
		numreadX_pre_full.clear();
		numreadX_pre_fullB.clear();
		numreadX_post_full.clear();
		
		immrateX_pre_full.clear();
		immrateX_pre_fullB.clear();
		immrateX_post_full.clear();
		
		extrateX_pre_full.clear();
		extrateX_pre_fullB.clear();
		extrateX_post_full.clear();
		
		immrateX_pre_full2.clear();
		immrateX_pre_fullB2.clear();
		immrateX_post_full2.clear();
		
		extrateX_pre_full2.clear();
		extrateX_pre_fullB2.clear();
		extrateX_post_full2.clear();
		
		abvsri_pre.clear();
		abvsri_preB.clear();
		abvsri_post.clear();
		
		resetup(immigrationin , localin , meta_typein , false , immtauin);
		// do not move this line down or you get a seg error
		
		vector<double> tempclass;
		tempclass.clear();
		int tempsize;
		tempsize = floor(log(double(localin))/log(2.0));
		for (long i = 0 ; i < tempsize+5 ; i ++)
		{
			tempclass.push_back(0.0);
		}
		
		for (long i = 0 ; i < meta.numspec()+1 ; i ++)
		{
			numreadX_pre_full.push_back(0.0);
			numreadX_pre_fullB.push_back(0.0);
			numreadX_post_full.push_back(0.0);
			
			immrateX_pre_full.push_back(0.0);
			immrateX_pre_fullB.push_back(0.0);
			immrateX_post_full.push_back(0.0);
			
			extrateX_pre_full.push_back(0.0);
			extrateX_pre_fullB.push_back(0.0);
			extrateX_post_full.push_back(0.0);
			
			immrateX_pre_full2.push_back(0.0);
			immrateX_pre_fullB2.push_back(0.0);
			immrateX_post_full2.push_back(0.0);
			
			extrateX_pre_full2.push_back(0.0);
			extrateX_pre_fullB2.push_back(0.0);
			extrateX_post_full2.push_back(0.0);
			
			abvsri_pre.push_back(tempclass);
			abvsri_preB.push_back(tempclass);
			abvsri_post.push_back(tempclass);
		}
		
		bool first = true;
		
		time(&end);
		while (difftime(end,start) < timeallowed)
		{
			out.open( file_log );
			out << " simulating immigrationin = " << immigrationin << " , localin = " << localin << " , meta_typein = " << meta_typein << " , immtau = " << immtauin  << " , timeallowed = " << timeallowed << " , \n";
			out << "run for " << difftime(end,start) << " of " << timeallowed << " total\n";
			cout << " simulating immigrationin = " << immigrationin << " , localin = " << localin << " , meta_typein = " << meta_typein << " , immtau = " << immtauin  << " , timeallowed = " << timeallowed << " , \n";
			cout << "run for " << difftime(end,start) << " of " << timeallowed << " total\n";
			out.close();
			resetup(immigrationin , localin , meta_typein , false , immtauin);
			while((numcycles == 0)&&(difftime(end,start) < timeallowed))
			{
				step();
				time(&end);
			}
			
			for (int i = 0 ; i < numread.size() ; i ++)
			{
				numreadX_pre_full[i]+=(numread[i]);
				immrateX_pre_full[i]+=(immrate_av[i]);
				extrateX_pre_full[i]+=(extrate_av[i]);
				immrateX_pre_full2[i]+=((immrate_av[i])*(immrate_av[i]));
				extrateX_pre_full2[i]+=((extrate_av[i])*(extrate_av[i]));
				numread[i]=0;
				immrate_av[i]=0;
				extrate_av[i]=0;
				
				for (int j = 0 ; j < (abund_vs_rich[i]).size() ; j ++)
				{
					// record abundance here
					abvsri_pre[i][j] += abund_vs_rich[i][j];
					abund_vs_rich[i][j] = 0;
				}
			}
			gen2eread += 1.0;
			gen2e += gen;
			gen2e2 += (gen*gen);
			while((numcycles < 2)&&(difftime(end,start) < timeallowed))
			{
				step();
				time(&end);
            }
			for (int i = 0 ; i < numread.size() ; i ++)
			{
				numreadX_post_full[i]+=(numread[i]);
				immrateX_post_full[i]+=(immrate_av[i]);
				extrateX_post_full[i]+=(extrate_av[i]);
				immrateX_post_full2[i]+=((immrate_av[i])*(immrate_av[i]));
				extrateX_post_full2[i]+=((extrate_av[i])*(extrate_av[i]));
				numread[i]=0;
				immrate_av[i]=0;
				extrate_av[i]=0;
				
				for (int j = 0 ; j < (abund_vs_rich[i]).size() ; j ++)
				{
					// record abundance stuff here
					abvsri_post[i][j] += abund_vs_rich[i][j];
					abund_vs_rich[i][j] = 0;
				}
			}
			resetup(immigrationin , localin , meta_typein , true , immtauin);
			while((numcycles == 0)&&(difftime(end,start) < timeallowed))
			{
				step();
				time(&end);
			}
			for (int i = 0 ; i < numread.size() ; i ++)
			{
				numreadX_pre_fullB[i]+=(numread[i]);
				immrateX_pre_fullB[i]+=(immrate_av[i]);
				extrateX_pre_fullB[i]+=(extrate_av[i]);
				immrateX_pre_fullB2[i]+=((immrate_av[i])*(immrate_av[i]));
				extrateX_pre_fullB2[i]+=((extrate_av[i])*(extrate_av[i]));
				numread[i]=0;
				immrate_av[i]=0;
				extrate_av[i]=0;
				
				for (int j = 0 ; j < (abund_vs_rich[i]).size() ; j ++)
				{
					// record abundance stuff here
					abvsri_preB[i][j] += abund_vs_rich[i][j];
					abund_vs_rich[i][j] = 0;
				}
			}
			gen2eread += 1.0;
			gen2e += gen;
			gen2e2 += (gen*gen);
			while((numcycles < 2)&&(difftime(end,start) < timeallowed))
			{
				step();
				time(&end);
			}
			for (int i = 0 ; i < numread.size() ; i ++)
			{
				numreadX_post_full[i]+=(numread[i]);
				immrateX_post_full[i]+=(immrate_av[i]);
				extrateX_post_full[i]+=(extrate_av[i]);
				immrateX_post_full2[i]+=((immrate_av[i])*(immrate_av[i]));
				extrateX_post_full2[i]+=((extrate_av[i])*(extrate_av[i]));
				numread[i]=0;
				immrate_av[i]=0;
				extrate_av[i]=0;
				
				for (int j = 0 ; j < (abund_vs_rich[i]).size() ; j ++)
				{
					// record abundance stuff here
					abvsri_post[i][j] += abund_vs_rich[i][j];
					abund_vs_rich[i][j] = 0;
				}
				
			}
			time(&end);
		}
		
		out.open( file_out , ofstream:: app);
		out << " , , immigration= , " << immigrationin << " , local= , " << localin << " , metatype= , " << meta_typein << " , immtau= , " << immtauin ;
		out << " , readings= , " << gen2eread << " , genM= , " << gen2e/gen2eread << " , genV= , " << ((gen2e2/gen2eread)-((gen2e/gen2eread)*(gen2e/gen2eread))) << " , START_DATA , \n";
		for (long i = 0 ; i < numreadX_pre_full.size() ; i ++)
		{
			out << meta_typein << " , "; // each line starts with the meta type to make it easier for gnuplot later
			out << numreadX_pre_full[i] <<  " , ";
			out << immrateX_pre_full[i] <<  " , ";
			out << extrateX_pre_full[i] <<  " , ";
			out << immrateX_pre_full2[i] <<  " , ";
			out << extrateX_pre_full2[i] <<  " , ";
			
			out << numreadX_post_full[i] <<  " , ";
			out << immrateX_post_full[i] <<  " , ";
			out << extrateX_post_full[i] <<  " , ";
			out << immrateX_post_full2[i] <<  " , ";
			out << extrateX_post_full2[i] <<  " , ";
			
			out << numreadX_pre_fullB[i] <<  " , ";
			out << immrateX_pre_fullB[i] <<  " , ";
			out << extrateX_pre_fullB[i] <<  " , ";
			out << immrateX_pre_fullB2[i] <<  " , ";
			out << extrateX_pre_fullB2[i] <<  " , ";
			
			out << i << " , ";
			if (numreadX_pre_full[i] > 0)
			{
				out << immrateX_pre_full[i] / numreadX_pre_full[i] <<  " , ";
				out << extrateX_pre_full[i] / numreadX_pre_full[i] <<  " , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << abvsri_pre[i][j] / numreadX_pre_full[i] <<  " , ";
				}
			}
			else 
			{
				out << " ,  , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << " , ";
				}
			}
			if (numreadX_post_full[i] > 0)
			{
				out << immrateX_post_full[i] / numreadX_post_full[i] << " , ";
				out << extrateX_post_full[i] / numreadX_post_full[i] << " , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << abvsri_post[i][j] / numreadX_post_full[i] <<  " , ";
				}
			}
			else 
			{
				out << " ,  , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << " , ";
				}
			}
			if (numreadX_pre_fullB[i] > 0)
			{
				out << immrateX_pre_fullB[i] / numreadX_pre_fullB[i] <<  " , ";
				out << extrateX_pre_fullB[i] / numreadX_pre_fullB[i] <<  " , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << abvsri_preB[i][j] / numreadX_pre_fullB[i] <<  " , ";
				}
			}
			else 
			{
				out << " ,  , ";
				for (int j = 1 ; j < (abvsri_pre[i]).size() ; j ++)
				{
					out << " , ";
				}
			}
			
			if (numreadX_pre_full[i] > 0)
			{
				out << pow(((immrateX_pre_full2[i] / numreadX_pre_full[i])-((immrateX_pre_full[i] / numreadX_pre_full[i])*(immrateX_pre_full[i] / numreadX_pre_full[i]))),0.5) <<  " , ";
				out << pow(((extrateX_pre_full2[i] / numreadX_pre_full[i])-((extrateX_pre_full[i] / numreadX_pre_full[i])*(extrateX_pre_full[i] / numreadX_pre_full[i]))),0.5) <<  " , ";
			}
			else 
			{
				out << " ,  , ";
			}
			if (numreadX_post_full[i] > 0)
			{
				out << pow(((immrateX_post_full2[i] / numreadX_post_full[i])-((immrateX_post_full[i] / numreadX_post_full[i])*(immrateX_post_full[i] / numreadX_post_full[i]))),0.5) << " , ";
				out << pow(((extrateX_post_full2[i] / numreadX_post_full[i])-((extrateX_post_full[i] / numreadX_post_full[i])*(extrateX_post_full[i] / numreadX_post_full[i]))),0.5) << " , ";
			}
			else 
			{
				out << " ,  , ";
			}
			if (numreadX_pre_fullB[i] > 0)
			{
				out << pow(((immrateX_pre_fullB2[i] / numreadX_pre_fullB[i])-((immrateX_pre_fullB[i] / numreadX_pre_fullB[i])*(immrateX_pre_fullB[i] / numreadX_pre_fullB[i]))),0.5) <<  " , ";
				out << pow(((extrateX_pre_fullB2[i] / numreadX_pre_fullB[i])-((extrateX_pre_fullB[i] / numreadX_pre_fullB[i])*(extrateX_pre_fullB[i] / numreadX_pre_fullB[i]))),0.5) <<  " , ";
			}
			else 
			{
				out << " ,  , ";
			}
			
			
			out << " \n";
			
		}
		out.close();
		
	}
};

int charconvertor(char charin)
{
	switch (charin) 
	{
		case '0': return (0);
		case '1': return (1);
		case '2': return (2);
		case '3': return (3);
		case '4': return (4);
		case '5': return (5);
		case '6': return (6);
		case '7': return (7);
		case '8': return (8);
		case '9': return (9);
		default: return (-1);
	}
}

int jobconvertor(char* argin)
{
	
	int maxind = 0;
	while (charconvertor(argin[maxind]) != -1) 
	{
		maxind ++;
	}
	int jobtoret = 0;
	int pow10 = 1;
	for (int i = maxind-1 ; i >=0 ; i --)
	{
		jobtoret += (pow10*charconvertor(argin[i]));
		pow10 = pow10*10;
	}
	return jobtoret;
}

int main(int argc) 
{
	dynamicsim test; // set seed here
	test.set_seed( 300 );
	double time_budget = 60*60/60; // 8 hours simulation time
	
    // these lines create the file names for output
    char file_out_name[50];
	char file_log_name[50];
	sprintf (file_out_name, "MWCI_out_%i.txt", 300);
	sprintf (file_log_name, "MWCI_log_%i.txt", 300);
    
    // run the simulation - use the following form
    
    /*
     test.sim( immigration rate - since this is clustered immigration you might choose to divide this by the cluster size,
     local community size,
     metacommunity type 1 = Difference log series 2 = log series 3 = uniform,
     simulation run time (seconds),
     file_out_name,
     file_log_name,
     number of generations needed for immigrant species to become 'good immigrants');
     */
    
	test.sim( 0.001 , 10000 , 1 , time_budget , file_out_name , file_log_name , 10 );
	
	
	
}