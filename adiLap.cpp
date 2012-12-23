/* Computational Physics: Apply the ADI(alternating direction implicit scheme) method to the Laplace problem with M = N = 5. Suppose the boundary value conditions as, for example, u5,j = 0 , u1,j = 1 for j = 1, ..., 5, and those values of points along the boundaries j = 1 and j = 5 changing linearly */




#include<iostream>
#include <iomanip> 
using namespace std;
void main()
{
	double A[9][9],T[9],v[9],w[9],vN[9],wN[9],b1w[9],b2w[9],b1v[9],b2v[9],g[8],h[8];
	double bv[9]={-0.25,-0.5,-1.75,0,0,-1,-0.25,-0.5,-1.75};//matrix bv initialization
	double bw[9]={-0.25,0,-0.25,-0.5,0,-0.5,-1.75,-1,-1.75};//matrix bw initialization
	int i,j,t,t1;
	for(j=0;j<=8;j++)                                       //matrix A initialization
	{
		for(i=0;i<=8;i++)
		{
			if(i==j)
				A[j][i]=-2;
			else if((i+j!=11)&&(i+j!=5)&&((i-j==1)||(i-j==-1)))
				A[j][i]=1;
			else
				A[j][i]=0;
		}
	}   
    double x[8]={1,1,0,1,1,0,1,1};
	cout<<"请输入9个初始条件:"<<endl;
	for(i=0;i<9;i++)
		cin>>v[i];                                             //matrix v initialization
    w[0]=v[0];							       //ralation of matrix v and matrix w
	w[1]=v[3];
	w[2]=v[6];
    w[3]=v[1];
	w[4]=v[4];
	w[5]=v[7];
	w[6]=v[2];
	w[7]=v[5];
	w[8]=v[8];	
	cout<<"please input the number of iteration t:";
	cin>>t;
	cout<<endl;
	for(t1=0;t1<t;t1++)
	{  
		for(i=0;i<9;i++)
	   {  
			b1w[i]=0;
		    for(j=0;j<9;j++)
				b1w[i]=b1w[i]+A[j][i]*w[j];
	   }
	
        b1v[0]=b1w[0];b1v[1]=b1w[3];b1v[2]=b1w[6];
	    b1v[3]=b1w[1];b1v[4]=b1w[4];b1v[5]=b1w[7];
     	b1v[6]=b1w[2];b1v[7]=b1w[5];b1v[8]=b1w[8];        
		for (i=0;i<9;i++)
			T[i]=bv[i]-b1v[i]-0.6283*v[i];	
    	g[7]=0.3805;
    	h[7]=T[8]/-2.6283;
    	for (i=6;i>=0;i--)
		{
			g[i]=-x[i]/(-2.6283+x[i+1]*g[i+1]);
			h[i]=(T[i+1]-x[i+1]*h[i+1])/(-2.6283+x[i+1]*g[i+1]);
		}
		vN[0]=(T[0]-h[0])/(-2.6283+g[0]);
		for (i=1;i<9;i++)
			vN[i]=g[i-1]*vN[i-1]+h[i-1];
		wN[0]=vN[0];
		wN[3]=vN[1];
		wN[6]=vN[2];
		wN[1]=vN[3];
		wN[4]=vN[4];
		wN[7]=vN[5];
		wN[2]=vN[6];
		wN[5]=vN[7];
		wN[8]=vN[8];
		for (i=0;i<9;i++)
			for(i=0;i<9;i++)
		   {
				b2v[i]=0;
				for(j=0;j<9;j++)
					b2v[i]=b2v[i]+A[j][i]*vN[j];
			}
        
		b2w[0]=b2v[0];b2w[1]=b2v[3];b2w[2]=b2v[6];
	    b2w[3]=b2v[1];b2w[4]=b2v[4];b2w[5]=b2v[7];
	    b2w[6]=b2v[2];b2w[7]=b2v[5];b2w[8]=b2v[8];
        
		for (i=0;i<9;i++)
			T[i]=bw[i]-b2w[i]-0.6283*wN[i];
		
		g[7]=0.3805;
    	h[7]=T[8]/-2.6283;
    	for (i=6;i>=0;i--)
		{
		g[i]=-x[i]/(-2.6283+x[i+1]*g[i+1]);
		h[i]=(T[i+1]-x[i+1]*h[i+1])/(-2.6283+x[i+1]*g[i+1]);
	
		}
		w[0]=(T[0]-h[0])/(-2.6283+g[0]);
        for (i=1;i<9;i++)
			w[i]=g[i-1]*w[i-1]+h[i-1];
        
		v[0]=w[0];v[1]=w[3];v[2]=w[6];                          // ralation of matrix v and matrix w
		v[3]=w[1];v[4]=w[4];v[5]=w[7];
		v[6]=w[2];v[7]=w[5];v[8]=w[8];
	}
	cout<<"0"<<setw(10)<<"0.25"<<setw(10)<<"0.5"<<setw(10)<<"0.75"<<setw(10)<<"1"<<endl;
	cout<<"0"<<setw(10)<<v[0]<<setw(10)<<v[1]<<setw(10)<<v[2]<<setw(10)<<"1"<<endl;
	cout<<"0"<<setw(10)<<v[3]<<setw(10)<<v[4]<<setw(10)<<v[5]<<setw(10)<<"1"<<endl;
    cout<<"0"<<setw(10)<<v[6]<<setw(10)<<v[7]<<setw(10)<<v[8]<<setw(10)<<"1"<<endl;
   	cout<<"0"<<setw(10)<<"0.25"<<setw(10)<<"0.5"<<setw(10)<<"0.75"<<setw(10)<<"1"<<endl;
}