/* 
   msff - apply a frequency filter to output from ms

   Copyright (C) 2002 Kevin Thornton

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include <vector>
// OS X
/*
#ifndef _UNISTD_H_ 
#include "getopt.h"
#endif
*/
#include <getopt.h>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/PolyTableFunctions.hpp>

using namespace std;
using namespace Sequence;
enum filtertype{MINOR,DERIVED};
struct msffargs
{
  filtertype filter;
  double freq;
};

void parseargs(int argc, char *argv[],msffargs *args);
void usage(void);

int main(int argc, char *argv[]) {
	msffargs args;
  	parseargs(argc,argv,&args);

  	SimParams p;
  	cin >> p;
  	SimData d(p.totsam());
  	unsigned nruns = p.runs();
 	unsigned totsam = p.totsam();
	vector<string> filtered_haps; // we will use this vector to store new haplotypes
	unsigned dercounts;
	
  	std::ios_base::sync_with_stdio(true);

  	int rv;
  	while( (rv=d.fromfile(stdin)) != EOF)
	{
		// counts derived alleles at first site ONLY
		dercounts = 0;
		for (unsigned j = 0 ; j < d.size() ; ++j)
	  	{
	    		dercounts += d[j][0] == '1' ? 1 : 0;
	  	}

		//depending on what we filter, print out individuals that have major allele
		switch (args.filter)
	  	{
	  		//if we are interested in minor allele, do:
	  		case MINOR:
				double minor_freq;
				minor_freq = 1.0-double(dercounts)/double(totsam) < double(dercounts)/double(totsam) ? 1.0-double(dercounts)/double(totsam) : double(dercounts)/double(totsam);
				//if minor allele at site 0 is > specified frequency
	    			if( 1.0-minor_freq > args.freq ) 
	      			{
					//identify which allele is minor 0 or 1; if tie assign 0 as major and 1 as minor
					char major_allele = dercounts > totsam-dercounts ? '1' : '0';

					//iterate over individuals (rows)	
					for(unsigned i = 0 ; i < totsam ; ++i)
		  			{
						//assign haplotype to vector of keepers if starts with major_allele
		 				if( d[i][0] == major_allele )
						{	
							filtered_haps.push_back(d[i]);	
						}
		  			}
	      			}
	    		break;
	  
	  		//if we are interested in derived allele, do:
	  		case DERIVED:
				//if derived allele at site 0 is > specified frequency
	    			if(double(dercounts)/double(totsam)>args.freq )
	      			{
					//iterate over individuals (rows)	
					for(unsigned i = 0 ; i < totsam ; ++i)
		  			{
						//assign haplotype to vector of keepers if starts with 1
		 				if( d[i][0] == '1' )
						{	
							filtered_haps.push_back(d[i]);	
						}
		  			}
	      			}
	    		break;	
	  	}

		//declare new SimData filtered_data and fill with the haplotypes that you want to keep
		SimData filtered_data;
		filtered_data.assign( &*d.pbegin(),d.numsites(),&filtered_haps[0],filtered_haps.size() );

		//Remove invariant sites
		RemoveInvariantColumns(&filtered_data);

		// print header of ms sim, with positions etc.
		fprintf(stdout,"//\nsegsites: %d\npositions: ",filtered_data.numsites());
		for(unsigned j = 0 ; j < filtered_data.numsites(); ++j) 
		{
			fprintf(stdout,"%lf ",filtered_data.position(j)); 
		}
		fprintf(stdout,"\n");
	
	  	for(unsigned i = 0 ; i < filtered_data.size() ; ++i)
  		{		
			cout << filtered_data[i] << endl;
			// fprintf(stdout,"%c",filtered_data[i][j]); // too lazy to write for loop for C style
  		}
  	}
	filtered_haps.clear();
}

void parseargs(int argc, char *argv[],msffargs *args)
{
  if(argc==1){
    usage();
    exit(1);
  }

  //assign some defaults
  args->filter=MINOR;
  args->freq = 999.0;
  int c;

  while ((c = getopt (argc, argv, "m:d:h")) != -1)
    {
      switch (c)
	{
	case 'm':
	  args->filter=MINOR;
	  args->freq = atof(optarg);
	  break;
	case 'd':
	  args->filter = DERIVED;
	  args->freq = atof(optarg);
	  break;
	case 'h':
	  usage();
	  exit(0);
	  break;
	case '?':
	  cerr<<"error: unknown argument "<< argv[optind-1]<<endl;
	  usage();
	  exit(0);
	  break;
	}
    }
  if (args->freq >= 1.0)
    {
      usage();
      exit(1);
    }
}

void usage(void)
{
  cerr<<"usage:\n";
  cerr<<"To filter based on minor allele frequency:\n";
  cerr<<"msff -m freq\n";
  cerr<<"To filter based on derived allele frequency,\n";
  cerr<<"use -d instead of -m"<<endl;
}
