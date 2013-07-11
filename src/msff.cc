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

  cout << p << '\n';

  unsigned totsam = p.totsam();

  unsigned MAXSITES = 100;
  unsigned *indexes = static_cast<unsigned *>(malloc(MAXSITES*sizeof(unsigned)));
  unsigned nindexes,dercounts;

  std::ios_base::sync_with_stdio(true);

  int rv;
  while( (rv=d.fromfile(stdin)) != EOF)
    {
    if (d.numsites() > MAXSITES)
      {
	MAXSITES = d.numsites()+1;
	indexes = static_cast<unsigned *>(realloc(indexes,MAXSITES*sizeof(unsigned)));
      }
    nindexes=0;
	
    for (unsigned i = 0 ; i < d.numsites() ; ++i)
//      for (unsigned i = 0; i < 1; ++i)
      {
	dercounts = 0;
	for (unsigned j = 0 ; j < d.size() ; ++j)
	  {
	    dercounts += d[j][i] == '1' ? 1 : 0;
	  }
	switch (args.filter)
	  {
	  case MINOR:
	    if (1.0-double(dercounts)/double(totsam)>args.freq&&
		double(dercounts)/double(totsam)>args.freq || i>0 )
	      {
		indexes[nindexes]=i; //JRI
//		cerr << "WTF " << i << " " << nindexes << endl; 
		if( i == 0 ){ nindexes++; } //JRI
	      }
	    break;
	  case DERIVED:
	    if(double(dercounts)/double(totsam)>args.freq || i>0 )
	      {
		indexes[nindexes]=i; //JRI
		if( i == 0 ){ nindexes++; } //JRI
	      }
	    break;
	  }
      }

    //print out the new, filtered, gametes
    //used C-style I/O b/c it can be much faster,
    //and speed is important here
    if (nindexes > 0)
      {
	fprintf(stdout,"//\nsegsites: %d\npositions: ",nindexes);
//	for(unsigned j = 0 ; j < nindexes ; ++j)
	for(unsigned j = 0 ; j < d.numsites() ; ++j) //JRI
	
	  {
	 //   fprintf(stdout,"%lf ",d.position(indexes[j]));
	    fprintf(stdout,"%lf ",d.position(j)); //JRI
	  }
	fprintf(stdout,"\n");
	for(unsigned i = 0 ; i < totsam ; ++i)
	  {
//	    for(unsigned j = 0 ; j < nindexes ; ++j)
	    for(unsigned j = 0 ; j < d.numsites() ; ++j) //JRI
	      {
//		fprintf(stdout,"%c",d[i][indexes[j]]);
		fprintf(stdout,"%c",d[i][j]); //JRI
	      }
	    fprintf(stdout,"\n");
	  }
      }
    else
      {
	fprintf(stdout,"//\nsegsites: 0\n\n");
      }
  }
  free(indexes);
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
