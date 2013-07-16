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
  SimData d2(p.totsam());	
  unsigned nruns = p.runs();

  cout << p << '\n';

  unsigned totsam = p.totsam();

//  unsigned *indexes = static_cast<unsigned *>(malloc(MAXSITES*sizeof(unsigned)));
//  unsigned nindexes,dercounts;
  unsigned dercounts;
	
  std::ios_base::sync_with_stdio(true);

  int rv;
  while( (rv=d.fromfile(stdin)) != EOF)
    {
	
	// print header of ms sim, with positions etc.
	fprintf(stdout,"//\nsegsites: %d\npositions: ",d.numsites());
	for(unsigned j = 0 ; j < d.numsites() ; ++j) //JRI
	  {
	    fprintf(stdout,"%lf ",d.position(j)); //JRI
	  }
	fprintf(stdout,"\n");

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
	    if (1.0-double(dercounts)/double(totsam)>args.freq&&
		double(dercounts)/double(totsam)>args.freq  )
	      {
		//indexes[nindexes]=i; //JRI
		for(unsigned i = 0 ; i < totsam ; ++i) 
	  	{
	    		for(unsigned j = 0 ; j < d.numsites() ; ++j) //JRI
	      		{
				fprintf(stdout,"%c",d[i][j]); //JRI
	      		}
	    		fprintf(stdout,"\n");
	  	}	
	      }
	    break;
	  
	  //if we are interested in derived allele, do:
	  case DERIVED:
	    if(double(dercounts)/double(totsam)>args.freq )
	      {
		//indexes[pnindexes]=i; //JRI
		unsigned newdudes=0;
		for(unsigned i = 0 ; i < totsam ; ++i)
		  {
		 	if( d[i][0] == '1' ){	
			    //iterate over individuals
			    for(unsigned j = 0 ; j < d.numsites() ; ++j) //JRI
			      {
				//print each SNP j for individual i	
				d2[newdudes][j]=d[i][j];
			      }
				//cout << endl;
			    newdudes++;
			}
		  }
	      }
	    break;
	  }

cout << d2[0][13] << endl << d2[1][13];

//  RemoveInvariantColumns(&d2);
  for(unsigned i = 0 ; i < totsam ; ++i)
  {
	for(unsigned j = 0 ; j < d.numsites() ; ++j) //
	{
		//cout << i << " " << j << endl;	
//		cout << i << " " << j << " " << d2[i][j] << endl;
		//fprintf(stdout,"%c",d2[i][j]); //JRI
	}
	fprintf(stdout,"\n");

  }


  }
  //free(indexes);

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
