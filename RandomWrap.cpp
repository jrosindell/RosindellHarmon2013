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
 
 This file uses a standard
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


# include <stdlib.h>

using namespace std;

class RandomWrap
{
    
private:
    
	bool seeded;
    // declare variables needed for your random number generator here.
    
public:
    
	RandomWrap()
	{
		seeded = false;
	}
    
	void set_seed(long seed)
	{
		if (!seeded)
		{
            srand ( seed ); // edit this line to set the seed as 'seed' in your random number generator of choice
			seeded = true;
		}
	}
    
	double d01()
	{
        return(double(rand())/double(RAND_MAX)); // edit this line to return a random number between 0 and 1 according to a uniform distribution using your random number generator
	}  
    
};

