#ifndef ADAPT_TRAPEZOID_H
#define ADAPT_TRAPEZOID_H
#include <list>
#include <utility>
#include <cassert>
//#include <iostream>

//using namespace std;

template<typename T1,typename T2,typename T3>
class triple
{
 public:
  T1 first;
  T2 second;
  T3 third;
  triple(T1 x1,T2 x2,T3 x3)
    :first(x1),second(x2),third(x3)
    {
    }
};


template<typename T1,typename T2,typename T3>
triple<T1,T2,T3> make_triple(T1 x1,T2 x2,T3 x3)
{
  return triple<T1,T2,T3>(x1,x2,x3);
}


template <typename T>
T adapt_trapezoid(T (*fun)(T,void*),T x1,T x2,T err_limit,void* param)
{
  //  const err_limit=.001;
  typedef triple<T,T,bool> interval;
  /*T for interval type,
    bool for state trur for still to be updated,
    false for do not need to be updated
   */
  std::list<interval> interval_list;
  T current_sum=((fun(x1,param)+fun(x2,param))/2.*(x2-x1));
  interval_list.push_back(make_triple(x1,current_sum,true));
  interval_list.push_back(make_triple(x2,(T)0.,true));
  bool int_state=1;
  int n_intervals=1;
  while(int_state)
    {
      //std::cout<<n_intervals<<std::endl;
      int_state=0;
      typename std::list<interval>::iterator i1=interval_list.begin();
      typename std::list<interval>::iterator i2=interval_list.begin();
      i2++;
      for(;i2!=interval_list.end();i1++,i2=i1,i2++)
	{
	  //cout<<i1->first<<"\t"<<i2->first<<endl;
	  //assert(i2->first>i1->first);
	  if(i1->third)
	    {
	      interval new_interval((i1->first+i2->first)/2,0,true);
	      
	      T sum1,sum2;
	      sum1=(fun(new_interval.first,param)+fun(i1->first,param))/2*(new_interval.first-i1->first);
	      sum2=(fun(new_interval.first,param)+fun(i2->first,param))/2*(i2->first-new_interval.first);
	      new_interval.second=sum2;
	      T err;
	      err=i1->second-sum1-sum2;
	      err=err<0?-err:err;

	      if(err>err_limit/n_intervals)
		{
		  i1->second=sum1;
		  interval_list.insert(i2,new_interval);
		  n_intervals++;
		  if(n_intervals>10e6)
		    {
		      
		      break;
		    }
		}
	      else
		{
		  i1->third=false;
		}
	      int_state=1;
	    }
	}
      
    }
  T result=0;
  for(typename std::list<interval>::iterator i=interval_list.begin();i!=interval_list.end();i++)
    {
      result+=i->second;
    }
  return result;
}

#endif
//end of the file
