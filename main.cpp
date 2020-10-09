#include <random>
#include <iostream>

using namespace std;

const int sample_size=90;
const int njack=90;
const int nboot=1000;
const double err=0.1;
const int npoints=10;
vector<mt19937_64> gen(npoints);

int iexp;

//! return an integer random number in the range [min,max)
double get_gauss(const double ave,const double sig,mt19937_64& g)
{
  return normal_distribution<>(ave,sig)(g);
}

//! return an integer random number in the range [min,max)
double get_int(const int min,const int max,mt19937_64& g)
{
  return uniform_int_distribution<>(min,max-1)(g);
}

struct data
{
  double x;
  
  vector<double> jack;
  
  vector<double> boot;
  
  double jack_ave() const
  {
    double ave=0;
    
    for(auto& j : jack)
      ave+=j;
    
    ave/=njack;
    
    return ave;
  }
  
  double jack_err() const
  {
    double ave=0,ave2=0;
    
    for(auto& j : jack)
      {
	ave+=j;
	ave2+=j*j;
      }
    
    ave/=njack;
    ave2/=njack;
    
    ave2-=ave*ave;
    
    return sqrt(ave2)*sqrt(njack-1);
  }
  
  double boot_ave() const
  {
    double ave=0;
    
    for(auto& b : boot)
	ave+=b;
    
    ave/=nboot;
    
    return ave;
  }
  
  double boot_err() const
  {
    double ave=0,ave2=0;
    
    for(auto& b : boot)
      {
	ave+=b;
	ave2+=b*b;
      }
    
    ave/=nboot;
    ave2/=nboot;
    
    ave2-=ave*ave;
    
    return sqrt(ave2);
  }
  
  data(const double& x,mt19937_64& g) : x(x),jack(njack),boot(nboot)
  {
    vector<double> r(sample_size);
    
    for(auto& ri : r)
      ri=get_gauss(x,err,g);
    
    if(iexp==0)
      for(auto& b : boot)
      {
	b=0;
	for(int i=0;i<sample_size;i++)
	  b+=r[get_int(0,sample_size,g)];
	b/=sample_size;
      }
    
    const int clust_size=sample_size/njack;
    for(int j=0;j<njack;j++)
      {
	jack[j]=0;
	for(int i=j*clust_size;i<(j+1)*clust_size;i++)
	  jack[j]+=r[i];
	jack[j]/=clust_size;
      }
    
    double s=0;
    for(auto& j : jack)
      s+=j;
    
    for(auto& j : jack)
      j=(s-j)/(njack-1);
  }
};

int main(int narg,char** arg)
{
  mt19937_64 glb_gen;
  
  if(narg>1)
    glb_gen.seed(atoi(arg[1]));
  
  const int nexp=1000;
  
  for(iexp=0;iexp<nexp;iexp++)
    {
      vector<data> n;
      
      for(int ipoint=0;ipoint<npoints;ipoint++)
	{
	  mt19937_64& g=gen[ipoint];
	  g.seed(get_int(0,1<<31,glb_gen));
	  //gen.seed(iatoi(arg[1]));
	  n.emplace_back(get_gauss(0,1.0,g),g);
	}
      
      double ca;
      double ce;
      
      if(iexp==0)
	{
	  ca=0;
	  ce=0;
	  
	  for(int ijack=0;ijack<njack;ijack++)
	    {
	      double c=0;
	      for(int i=0;i<npoints;i++)
		{
		  const double x=(n[i].jack[ijack]-n[i].x)/n[i].jack_err();
		  c+=x*x;
		}
	      
	      c/=npoints;
	      
	      cout<<"jch2: "<<c<<endl;
	      
	      ca+=c;
	      ce+=c*c;
	    }
	  
	  ca/=njack;
	  ce/=njack;
	  ce-=ca*ca;
	  cout<<"Jack ch2: "<<ca<<" +- "<<sqrt(ce/njack*(njack-1))<<endl;
	  
	  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
	}
      
      if(iexp==0)
	{
	  ca=0;
	  ce=0;
	  
	  for(int iboot=0;iboot<nboot;iboot++)
	    {
	      double c=0;
	      for(int i=0;i<npoints;i++)
		{
		  const double x=(n[i].boot[iboot]-n[i].x)/n[i].boot_err();
		  c+=x*x;
		}
	      
	      c/=npoints;
	      
	      cout<<"bch2: "<<c<<endl;
	      
	      ca+=c;
	      ce+=c*c;
	    }
	  
	  ca/=nboot;
	  ce/=nboot;
	  ce-=ca*ca;
	  
	  cout<<"Boot ch2: "<<ca<<" +- "<<sqrt(ce)<<endl;
	  
	  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
	}
      
      double c=0;
      for(int i=0;i<npoints;i++)
	{
	  const double x=(n[i].jack_ave()-n[i].x)/n[i].jack_err();
	  c+=x*x;
	}
      
      c/=npoints;
      
      cout<<"ch2: "<<c<<endl;
      
      // for(int i=0;i<npoints;i++)
      // 	cout<<n[i].x<<"\t"
      // 	    <<" "<<(n[i].boot_ave()-n[i].x)/n[i].boot_err()<<"\t"
      // 	    <<" "<<(n[i].jack_ave()-n[i].x)/n[i].jack_err()<<"\t"
      // 	    <<" "<<n[i].boot_err()/n[i].jack_err()<<endl;
    }
  
  return 0;
}
