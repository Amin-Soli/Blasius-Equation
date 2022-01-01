#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>

using namespace std;

#ifdef _WIN32
#include <windows.h>
#include <sys/stat.h>
#include <dirent.h>

#define DIV 1048576 
#define WIDTH 7
#endif

#ifdef linux
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#endif

/********************************************************************************************************/
/*********************************** prototype of functions *********************************************/
/********************************************************************************************************/

double ** product_three_matrix(double **, double **, double **);
double * product_matrix_and_vector(double **, double *);
double * product_two_matrix_and_vector(double **, double **, double *);
double ** inverse_matrix (double **);
void keller_box_method_using_thomas_algorithm(double *,double *,double *,int ,double ,double *);

/********************************************************************************************************/
/********************************************************************************************************/


int main()
{
	/************* define variables **************/

    int i,n,count=0;
    double dx,k,max_eta ;
    double *f, *u, *v, *eta, *h, *copy_f, *copy_u, *copy_v;
    double error_f,error_u,error_v, max_Error, alfa, ErrorLimit;

    ////////////////// input data section given by user /////////////////

    dx = 0.25;  // the first grid size

    k = 1.15;  // ratio magnitude (best choose = 1.15)

    max_eta = 14;   // maximum of eta

    alfa = 1;

    ErrorLimit = 1.0e-6;

    /////////////////////////////////////////////////////////////////////

    /************* number of grids **************/

    n=log(1+ max_eta*(k-1)/dx)/log(k) + 1;

    /********************************************/

    f = new double[n];
    u = new double[n];
    v = new double[n];
    h = new double[n];
    eta = new double[n];
    copy_f = new double[n];
    copy_u = new double[n];
    copy_v = new double[n];

    /************* calculation of eta[i] **************/

    for(i=0;i<n;i++)
        eta[i]=dx*(pow(k,i)-1)/(k-1);

    /************* calculation of different grid size **************/

    for(i=1;i<n;i++)
    {
        if(i==1)
            h[i]=dx;

        else
            h[i]=k*h[i-1];
    }


    /************* initial guess for f,u,v **************/

    for(i=0;i<n;i++)
        f[i]=eta[i];

    for(i=0;i<n;i++)
        u[i]=1.0*i/(n-1);

    for(i=0;i<n;i++)
        v[i]=exp(-eta[i]*eta[i]/2.0);


    /************* solution using keller box method using tomas algorithm **************/

    do
    {

        count++;

        for(i=0;i<n;i++)
        {
            copy_f[i]=f[i];
            copy_u[i]=u[i];
            copy_v[i]=v[i];
        }

        keller_box_method_using_thomas_algorithm(f,u,v,n,alfa,h);

		max_Error = 0.0;

        for(i=0;i<n;i++)
        {
            error_f=fabs(f[i]-copy_f[i]);
            error_u=fabs(u[i]-copy_u[i]);
            error_v=fabs(v[i]-copy_v[i]);

            if(error_f > max_Error)
                max_Error = error_f;

            if(error_u > max_Error)
                max_Error = error_u;

            if(error_v > max_Error)
                max_Error = error_v;

        }


    }while(max_Error > ErrorLimit);

    /************************* write solutions in output *******************************/
    
    cout << "input data of the problem:" << endl << endl;
    cout << "the first grid size: " << dx << endl; 
    cout << "ratio magnitude : " << k << endl; 
    cout << "maximum of eta : " << max_eta << endl; 
    cout << "alpha : " << alfa << endl;
	cout <<"maximum error for solving : "<< ErrorLimit << endl << endl;
	
	#ifdef _WIN32

	if (mkdir("results") != -1)
        cout << "Directory of results was created." << endl;
        
	#endif
	
	#ifdef linux
	
	if (opendir("results"))
		system("rm -r results");

	if (mkdir("results", 0777) != -1)
        cout << "Directory of results was created." << endl;
        
	#endif

    ofstream file1;
    file1.open("results/f.txt");

    file1 << "eta" << '\t' << '\t' << '\t' << "f" << endl << endl;

    for(i=0;i<n;i++)
        file1 << std::fixed << std::setprecision(6) << eta[i] <<'\t' << '\t' << f[i] << endl;

    file1.close();

    ofstream file2;
    file2.open("results/u.txt");

    file2 << "eta" << '\t' << '\t' << '\t' << "u" << endl << endl;

    for(i=0;i<n;i++)
        file2 << std::fixed << std::setprecision(6) << eta[i] <<'\t' << '\t' << u[i] << endl;

    file2.close();

    ofstream file3;
    file3.open("results/v.txt");

    file3 << "eta" << '\t' << '\t' << '\t' << "v" << endl << endl;

    for(i=0;i<n;i++)
        file3 << std::fixed << std::setprecision(6) << eta[i] <<'\t' << '\t' << v[i] << endl;

    file3.close();

	cout << "result data was written." << endl;
	
	cout << "number of iterations to solve the problem: " << count << endl << endl;

    cout<<"end program."<<endl;

    cin.get();

}

/********************************************************************************************************/
/********************************************************************************************************/
/************************************ definition of functions *******************************************/
/********************************************************************************************************/
/********************************************************************************************************/


void keller_box_method_using_thomas_algorithm(double *f,double *u,double *v,int n,double alfa,double *h)
{
    int i,j,k, m = 3;

    double ***B, ***A, ***C, **R, **X;

    double **inverseB, **product_A_and_inverseB_and_C, *product_A_and_inverseB_and_R,

           *product_inverseB_and_R, *product_inverseB_and_C_and_X;

    B = new double **[n];
    A = new double **[n-1];
    C = new double **[n-1];
    R = new double *[n];
    X = new double *[n];


    for (i=0;i<n;i++)
    {
		 B[i] = new double *[m];
		 R[i] = new double [m];
		 X[i] = new double [m];
		 if (i!=n-1)
			{
				A[i] = new double *[m];
				C[i] = new double *[m];
			}
	}

	for (i=0;i<n;i++)
		for (j=0;j<m;j++)
		{
			B[i][j] = new double [m];
			if (i!=n-1)
				{
					A[i][j] = new double [m];
					C[i][j] = new double [m];
				}
		}


    /**************** set up R matrix ********************/

    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
        {
            if(j==0 && i!=0)
               R[i][j]=f[i-1]-f[i] + h[i]*(u[i]+u[i-1])/2.0;

            else if(j==1 && i!=0)
               R[i][j]=v[i-1]-v[i] - alfa*h[i]*(f[i]*v[i]+f[i-1]*v[i-1])/4.0;

            else if(j==2 && i!=n-1)
                R[i][j]=u[i] - u[i+1]+ h[i+1]*(v[i]+v[i+1])/2.0 ;

            else
                R[i][j]=0 ;

        }


   /**************** set up B matrix ********************/

   for(i=0;i<n;i++)
    {
       if(i==0)
        for(j=0;j<m;j++)
         for(k=0;k<m;k++)
         {
            if(j==k && j!=2)
               B[i][j][k]=1;

            else if(j==k && j==2)
               B[i][j][k]=-h[i+1]/2.0;

            else if(j==k+1 && j==2)
                B[i][j][k]=-1;

            else
                B[i][j][k]=0;
         }


       else if(i==n-1)
        for(j=0;j<m;j++)
         for(k=0;k<m;k++)
          {
            if(j==k && j==0)
               B[i][j][k]=1;

            else if(j==k+1 && j==1)
               B[i][j][k]=alfa*h[i]*v[i]/4.0;

            else if(j==k+1 && j==2)
                B[i][j][k]=1;

            else if(j==k-1 && j==0)
               B[i][j][k]=-1*h[i]/2.0;

            else if(j==k-1 && j==1)
                B[i][j][k]=1+alfa*h[i]*f[i]/4.0;

            else
                B[i][j][k]=0;
           }

      else
        for(j=0;j<m;j++)
         for(k=0;k<m;k++)
          {
            if(j==k && j==0)
               B[i][j][k]=1;

            else if(j==k && j==2)
               B[i][j][k]=-1*h[i+1]/2.0;

            else if(j==k+1 && j==1)
               B[i][j][k]=alfa*h[i]*v[i]/4.0;

            else if(j==k+1 && j==2)
                B[i][j][k]=-1;

            else if(j==k-1 && j==0)
               B[i][j][k]=-1*h[i]/2.0;

            else if(j==k-1 && j==1)
                B[i][j][k]=1+alfa*h[i]*f[i]/4.0;

            else
                B[i][j][k]=0;
         }

    }

    /**************** set up A matrix ********************/

    for(i=0;i<n-1;i++)
     for(j=0;j<m;j++)
      for(k=0;k<m;k++)
        {
            if(j==k && j==0)
               A[i][j][k]=-1;

            else if(j==k+1 && j==1)
               A[i][j][k]=alfa*h[i]*v[i-1]/4.0;

            else if(j==k-1 && j==0)
               A[i][j][k]=-1*h[i]/2.0;

            else if(j==k-1 && j==1)
                A[i][j][k]=-1+alfa*h[i]*f[i-1]/4.0;

            else
                A[i][j][k]=0;
        }

    /**************** set up C matrix ********************/

    for(i=0;i<n-1;i++)
     for(j=0;j<m;j++)
      for(k=0;k<m;k++)
        {
          if(j==k && j==2)
            C[i][j][k]=-1*h[i+1]/2.0;

          else if(j==k+1 && j==2)
             C[i][j][k]=1;

          else
            C[i][j][k]=0;
        }

    /**************** forward elimination ********************/

    for(i=1;i<n;i++)
    {
        inverseB = inverse_matrix(B[i-1]);
        product_A_and_inverseB_and_C = product_three_matrix(A[i-1] ,inverseB , C[i-1]);
        product_A_and_inverseB_and_R = product_two_matrix_and_vector(A[i-1] ,inverseB , R[i-1]);

        for(j=0;j<m;j++)
            for(k=0;k<m;k++)
              B[i][j][k] = B[i][j][k] - product_A_and_inverseB_and_C[j][k];


        for(j=0;j<m;j++)
              R[i][j] = R[i][j] - product_A_and_inverseB_and_R[j];

    }

    /**************** backward subsitution to calculate X ********************/

    for(i=n-1;i>=0;i--)
    {
        if(i==n-1)
        {
            inverseB = inverse_matrix(B[i]);
            X[i] = product_matrix_and_vector(inverseB, R[i]);
        }

        else
        {
            inverseB = inverse_matrix(B[i]);
            product_inverseB_and_R = product_matrix_and_vector(inverseB, R[i]);
            product_inverseB_and_C_and_X = product_two_matrix_and_vector(inverseB, C[i] ,X[i+1]);

            for(j=0;j<m;j++)
              X[i][j] = product_inverseB_and_R [j] - product_inverseB_and_C_and_X[j];

        }

    }

    /**************** calculate f, u, and v using X ********************/

    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
    {
        if(j==0)
          f[i]=f[i]+X[i][j];

        else if(j==1)
          u[i]=u[i]+X[i][j];

        else
          v[i]=v[i]+X[i][j];
    }

}

/********************************************************************************************************/
/********************************************************************************************************/

double ** product_three_matrix(double **A, double **inverseB, double **C)
{
	int i,j,k, m=3;
    double a[m][m];
    double **x;
    x = new double *[m];
    for (i=0;i<m;i++)
		x[i] = new double [m];

    for(i=0;i<m;i++)
     for(j=0;j<m;j++)
       a[i][j]=0;

    for(i=0;i<m;i++)
     for(j=0;j<m;j++)
       x[i][j]=0;

    for(i=0;i<m;i++)
     for(j=0;j<m;j++)
       for(k=0;k<m;k++)
           a[i][j]+=A[i][k]*inverseB[k][j];


    for(i=0;i<m;i++)
     for(j=0;j<m;j++)
       for(k=0;k<m;k++)
           x[i][j]+=a[i][k]*C[k][j];

	return x;
}

/********************************************************************************************************/
/********************************************************************************************************/

double * product_matrix_and_vector(double **inverseB, double *R)
{
	int i,j, m=3;
    double *x;
    x = new double [m];

    for(i=0;i<m;i++)
        x[i]=0;

    for(i=0;i<m;i++)
       for(j=0;j<m;j++)
           x[i]+=inverseB[i][j]* R[j];

    return x;

}

/********************************************************************************************************/
/********************************************************************************************************/

double * product_two_matrix_and_vector(double **y, double **z, double *w)
{
	int i,j, m=3;
    double a[3], *x;
    x = new double [3];

    for(i=0;i<m;i++)
       a[i]=0;

    for(i=0;i<m;i++)
       x[i]=0;

    for(i=0;i<m;i++)
       for(j=0;j<m;j++)
           a[i]+=z[i][j]* w[j];

    for(i=0;i<m;i++)
       for(j=0;j<m;j++)
           x[i]+=y[i][j]* a[j];

    return x;
}

/********************************************************************************************************/
/********************************************************************************************************/

double ** inverse_matrix (double **B)
{
    int i, j, k, flag, w, z, v, m=3;
    double a[3][2*m], t, b;
    double **x;
    x = new double *[m];
    for (i=0;i<m;i++)
		x[i] = new double[m];

    for(i=0;i<m;i++)
        for(j=0;j<m;j++)
        a[i][j]=B[i][j];

    for(i=0;i<m;i++)
     for (j=m;j<2*m;j++)
        {
          a[i][j]=0;
          if(i+m==j)
            a[i][j]=1;
        }

    for(i=0;i<m-1;i++)
    for(j=i+1;j<m;j++)
    {
        if(fabs(a[i][i])<fabs(a[j][i]))
            for(k=0;k<2*m;k++)
            {
             t=a[i][k];
             a[i][k]=a[j][k];
             a[j][k]=t;
            }
    }

    for(i=0;i<m-1;i++)
    {
      flag=1;
      for(w=i,z=1;(w<m-1 && flag==1);w++,z++)
      {
        flag=0;
        if(a[i][i]==0)
        {
            for(v=0;v<2*m;v++)
            {
            t=a[i][v];
            a[i][v]=a[i+z][v];
            a[i+z][v]=t;
            }
            flag=1;

        }

      }

      for(j=i+1;j<m;j++)
       {
        b=a[j][i]/a[i][i];
        for(k=i;k<2*m;k++)
            a[j][k]-=a[i][k]*b;
        }
    }

    for(i=m-1;i>0;i--)
    for(j=i-1;j>=0;j--)
    {
     b=a[j][i]/a[i][i];
     for(k=i;k<2*m;k++)
        a[j][k]-=a[i][k]*b;
    }

    for(i=0;i<m;i++)
     for(j=m;j<2*m;j++)
        a[i][j]/=a[i][i];

    for(i=0;i<m;i++)
     for(j=m,k=0;j<2*m;j++,k++)
        x[i][k]=a[i][j];

	return x;

}


