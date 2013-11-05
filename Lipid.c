
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

//------physical parameters-------------------

double L  = 40.0;        // length of domain 

double pi = 3.141592;    // pi is pi 

int    n = 150;          // subdivisions
double dx;               // 1/(n-1)
double delT  = 0.501;    // time step
double delT0 = 0.501;    // time step
int    NElements;        // n  - 1
int    NVariables;       // n*4*2
int    NPositions  = 16; //make sure this number agrees with the one used 
                         //in the definition of struct ele 
int MaxNeighbor = 24;
double r  = 2.5;         // radius
double R  = 3.0;         // height
double h0 = 1.5;         // lipid length at rest 

char command[100];
char loadfilename[100]; 
int LOAD = 0;            // 0 = no, 1 = yes
int CALLIBRATION = 0;    // 0 = no, 1 = yes...run arclength calibration...30 or so steps at delT = 0.01 works
double ihs0, its0, ots0, ohs0; //inner outer head tail arclength 

double Kt       = 1.0;    // tilt modulus
double Kb       = 1.0;    // bending modulus
double innerK0  = 0.05;   // spontaneous curvature 
double outerK0  = 0.15;   // spontaneous curvature 
double ST       = 0.0000; // surface tenstion 
double STinner  = 0.0000; // surface tenstion 
double STouter  = 0.0000; // surface tenstion 

double NewtTol = 1.0e-7; 
int NewtStep, MaxNewtStep = 10;
int TimeStep, MaxTimeStep = 500;
  
//--------------------------------------------

struct vec {double x;double y;};
struct elem { int e[16]; };
struct quint {double *D;double *U;double *UU;double *L;double *LL;};

int SurfEngDen (double rhol, double rhor, double phi, double *w);
int BendingEnergy(double ** Variables,\
						double *Variables0,\
						double *dS,\
						int * iEle, int * jEle, int *Ele,\
						double *E,\
						double *EGrad,\
						double * aHess, int * iHess, int * jHess);
int BendingEnergyDensity(double ihlx, double ihly, double itlx, double itly, \
								 double ihrx, double ihry, double itrx, double itry, \
								 double ohlx, double ohly, double otlx, double otly, \
								 double ohrx, double ohry, double otrx, double otry, \
								 double ihlx0, double ihly0, double itlx0, double itly0, \
								 double ohlx0, double ohly0, double otlx0, double otly0, \
								 double ihrx0, double ihry0, double itrx0, double itry0, \
								 double ohrx0, double ohry0, double otrx0, double otry0, \
								 double dSConstraint, 
								 double *E);
double Volume(struct vec hl, struct vec hr,\
				  struct vec tl, struct vec tr);
int QuintDiagSolve(double *D,\
						 double *U,  double *L,\
						 double *UU, double *LL,\
						 double *b,  int n); 
int Elements(int *iEle, int *jEle, double *aE,\
				 int *ia, int *ja, double *aa);
int Transpose(int m, int n, double *a, int *ia, int *ja, double **b, int **ib, int **jb);
int ShowMatrix(int m, int n, double *a, int *ia, int *ja, char dot);
int ATimesB(int m, int n, int p, double *a, int *ia, int *ja, double *b, int *ib, int *jb, double **c, int **ci, int **cj, int mem, int eval);
int ZeroOut(int n, double *kill, int *ia, int *ja, double *aa); 
int SSOR(double *aa, int *ja, int *ia, int n, double *b, double w);
int PCG(double *aa, int *ja, int *ia, int n, double *x, double *b, double *r, double *d, double *q, double tol);
int Load(int NVariables,\
			struct vec *innerh,\
			struct vec *innert,\
			struct vec *outerh,\
			struct vec *outert);
int GetArcLengths(struct vec * innerh, struct vec * innert,\
						struct vec * outerh, struct vec * outert);
int main(int argc, char **argv) { 

  
	if(argc < 9) {
	  
	  printf("\nargc = %d; please supply input: r R L delT maxiter innerK0 outerK0 h0\n\n",argc);
	  exit(1);
	}
  

	r           = atof(argv[1]); 
	R           = atof(argv[2]);
	L           = atof(argv[3]);
	delT        = atof(argv[4]);
	MaxTimeStep = atoi(argv[5]);
	innerK0     = atof(argv[6]);
	outerK0     = atof(argv[7]);
	h0          = atof(argv[8]);

	sprintf(loadfilename,"loadprofile.dat");

	if(strcmp(argv[9],"L")==0) {
		LOAD = 1;
	}

	char infofilename[100]; 
	char profilefilename[100]; 
  
	sprintf(infofilename,"%g-%g-%g-%g-%d-%g-%g-%g.info",\
			  r,R,L,delT,MaxTimeStep,innerK0,outerK0,h0);

	sprintf(profilefilename,"%g-%g-%g-%g-%d-%g-%g-%g.profile",\
			  r,R,L,delT,MaxTimeStep,innerK0,outerK0,h0);

	FILE * infofile = fopen(infofilename, "w");

	fprintf(infofile,"r = %g, R = %g, L = %g, delT = %g, maxiter = %d\n", r,R,L,delT,MaxTimeStep);
	fprintf(infofile,"Kb = %g, Kt = %g, innerK0 = %g, outerK0 = %g, h0 = %g\n", Kb,Kt,innerK0, outerK0, h0);

	fclose(infofile);

	//mesh parameter n 
	//mesh size dx 
	dx = 1.0/((double)n-1.0);
	int dim = 2;   

	/* 

     DEFINE THE INITIAL GEOMETRY
     OF THE FUSION PORE

	*/

	struct vec *innerh = malloc(sizeof(struct vec)*n); //positions of head group on inner leaflet
	struct vec *innert = malloc(sizeof(struct vec)*n); //positions of tail group on inner leaflet 
	struct vec *outerh = malloc(sizeof(struct vec)*n); //positions of head group on outer leaflet 
	struct vec *outert = malloc(sizeof(struct vec)*n); //positions of tail group on outer leaflet 

	//construct an initial toroidal shape 
	/*

   length apportioned to flat portion  = L - r - R
   length apportioned to round portion = pi*R/2
   a*n points in flat portion 
   b*n points in round portion 
   a + b = 1

   b = pi*R/(L - R - r + pi*R) gives a spacing at the terminal round portion
   equal to the spacing on the beginning flat portion  
 
	*/

	double phiend, *phi = malloc(sizeof(double)*n);
	phiend = pi - atan( R / ( L - (r + R) )); //end angle for phi
	double inc = (-R*cos(phiend)/sin(phiend) + R*pi/2.0)/((double)(n-1));
	int i,j,k;
	for(i = 0; inc*((double)i) <=  R*pi/2.0; ++i) phi[i] = inc*((double) i)/R;  
	for(j = i; j < n ;++j) phi[j] = atan(1.0 / (pi/2.0 - inc*((double) j)/R)) + pi;    

	double rho, dx, dy;

	for(i = 0; i < n; ++i) {

		if(phi[i] <= pi/2.0) {
			rho = R;
			dx = h0*cos(phi[i]); dy = h0*sin(phi[i]);
		}

		else {
			rho = R/sin(phi[i]);
			dx = 0.0; dy = h0;
		}

		//parametrize inner torus shape 
		innert[i].x = r + R - rho*cos(phi[i]);
		innert[i].y = rho*sin(phi[i]);
		innerh[i].x = innert[i].x - dx;
		innerh[i].y = innert[i].y + dy;

		//parametrize outer torus shape
		outert[i].x = r + R - rho*cos(phi[i]);
		outert[i].y = rho*sin(phi[i]);
		outerh[i].x = outert[i].x + dx;
		outerh[i].y = outert[i].y - dy;

	}


	double A = 0.8, B = 0.2;  

	double * dS = malloc(sizeof(double)*(n-1));
  
	double b = ((2.0*A+B)*pi*R/2.0)/(L - R - r + pi*R);
	double *sphi = malloc(sizeof(double)*n);

	double *st = malloc(sizeof(double)*n);

	for(i=0;i<b*n;++i) st[i] = ((double) i)/( (b*((double) n)-1.0) );
  
	for(i=0;i<b*n;++i) sphi[i] = (A*st[i]*st[i]+B*st[i])*pi/2.0;
	inc = (L-R-r)/((1.0-b)*((double)n));

	for(i=0;i<b*n;++i) {

		dx = h0*cos(sphi[i]); dy = h0*sin(sphi[i]);

		//parametrize inner torus shape 
		innert[i].x = r + R - R*cos(sphi[i]);
		innert[i].y = R*sin(sphi[i]);
		innerh[i].x = innert[i].x - dx;
		innerh[i].y = innert[i].y + dy;

		//parametrize outer torus shape
		outert[i].x = r + R - R*cos(sphi[i]);
		outert[i].y = R*sin(sphi[i]);
		outerh[i].x = outert[i].x + dx;
		outerh[i].y = outert[i].y - dy;

	}

	for(i=floor(b*n);i<n;++i) {

		//parametrize inner torus shape 
		innert[i].x = r + R + (0.5 + i-floor(b*n))*inc;
		innert[i].y = R;
		innerh[i].x = innert[i].x;
		innerh[i].y = innert[i].y + h0;

		//parametrize outer torus shape
		outert[i].x = r + R + (0.5 + i-floor(b*n))*inc;
		outert[i].y = R;
		outerh[i].x = outert[i].x;
		outerh[i].y = outert[i].y - h0;

	}

	for(i=0;i<n-1;++i) {
    
		dx = innert[i+1].x-innert[i].x;
		dy = innert[i+1].y-innert[i].y;

		dS[i] = sqrt(dx*dx + dy*dy);
    
	}
  

	// print out initial torus
	FILE * pfile = fopen("profile.dat", "w");
	for(i=0;i<n;++i) {
		fprintf(pfile,"%g %g %g %g %g %g %g %g\n",\
				  innerh[i].x,innerh[i].y,\
				  innert[i].x,innert[i].y,\
				  outerh[i].x,outerh[i].y,\
				  outert[i].x,outert[i].y);
	} 
	fclose(pfile);

	/*
    ASSIGN THE VARIABLES TO GLOBAL POINTERS    
	*/  

	NVariables = 4*2*n; //n pairs of lipids, 4 points per pair, 2 coordinates per point
	double ** Variables  = malloc(sizeof(double *)*NVariables);

	int step = 0;
	for(i=0;i<n;++i) {

		//go through each lipid group and assign
		//the pointer to each of the eight possible
		//variables 

		//     inner   outer
		//    H~~~~~T T~~~~~H
		//  0,1   2,3 4,5   5,7 
		//
		// scheme for order...pay attention

		Variables[step]=&(innerh[i].x); ++step; 
		Variables[step]=&(innerh[i].y); ++step;
		Variables[step]=&(innert[i].x); ++step;
		Variables[step]=&(innert[i].y); ++step;
		Variables[step]=&(outerh[i].x); ++step;
		Variables[step]=&(outerh[i].y); ++step;
		Variables[step]=&(outert[i].x); ++step;
		Variables[step]=&(outert[i].y); ++step;
	}

	//POSITION OF THE ESSENTIAL BOUNDARY VALUES
	//the essential boundary values are the positions of the lipids
	//at the parallel portion of the membrane
	// they appear in lipid pair 0
	// i.e. the coordinates 0:7
	// and in lipid pair n-1
	// i.e. the coordinates 8*(n-1):8(n-1)+8  

	double * Essential = malloc(sizeof(double)*NVariables);
	for(i=0;i<NVariables;++i) Essential[i] = 1.0;
	for(i=0;i<8; ++i) Essential[i] = 0.0;
	for(i=8*(n-1);i<NVariables;++i) Essential[i] = 0.0;

	//here we adjust for the fact that the outer and inner head
	//can change their x coordinate but not their y coordinate
	//which should be indentically 0
	// ihx ihy itx ity otx oty ohx ohy  
	//  1   0   1   0   0   0   1   0
	//  0   1   2   3   4   5   6   7

	Essential[0] = 1.0;
	Essential[2] = 1.0;
	Essential[6] = 1.0;

	/*                                  *
     
    CONSTRUCTING THE ADJACENCY MATRIX

	 *                                  */
  
	//go through the rhomboids and assign 
	//the index of the 16 variables in each element 
	//we've order the first 8 elements the same 
	//as in the assignment to variable, namely 
	//0    1    2    3    4    5    6    7
	//ih.x ih.y it.x it.y oh.x oh.y ot.x ot.y 

	NElements   = n - 1;  
	int e; //e for element   
	int * iEle = malloc(sizeof(int)*NElements);
	int * jEle = malloc(sizeof(int)*NElements*NPositions);
	int *  Ele = malloc(sizeof(int)*NElements*NPositions);

	iEle[0] = 0;

	for(e = 0; e < n-1; ++e) {
    
		iEle[e+1] = iEle[e] + 16;

		for(i = 0; i < 16;++i) {
			jEle[iEle[e]+i] = 8*e + i;
			Ele[iEle[e]+i]  = i;
		}

	}

	//Transpose of adjacency matrix 
	int *iEleT;
	int *jEleT;

	Transpose(NElements, NVariables, NULL, iEle, jEle, NULL, &iEleT, &jEleT);
	//ShowMatrix(NElements, NVariables, NULL, iEle, jEle, '.');
	//ShowMatrix(NVariables, NElements, NULL, iEleT, jEleT, '.');


	double * aHess = malloc(sizeof(double)*NVariables*MaxNeighbor);
	int    * iHess = malloc(sizeof(int)*NVariables*MaxNeighbor);
	int    * jHess = malloc(sizeof(int)*NVariables*MaxNeighbor);

	ATimesB(NVariables, NElements, NVariables, NULL, iEleT, jEleT, NULL, iEle, jEle, NULL, &iHess, &jHess, 1, 0);
	//ShowMatrix(NVariables, NVariables, NULL, iHess, jHess, '.');

	double Energy;
	double *EGrad = malloc(sizeof(double)*NVariables); 

	double *  Variables0 = malloc(sizeof(double)*NVariables); //variables for the previous time step
	double *  dVariables = malloc(sizeof(double)*NVariables); //variables for the previous time step
	double *  dd = malloc(sizeof(double)*NVariables); //variables for the previous time step
	double *  qq = malloc(sizeof(double)*NVariables); //variables for the previous time step
	double *  rr = malloc(sizeof(double)*NVariables); //variables for the previous time step

	//PCG variables 
	//here we don't need the structure because we're just going to point to the previous values and 
	//then copy them over 

	double residual;

	//USE DATA LOADED FROM loadprofile.dat INSTEAD!
	if(LOAD) { 

		Load(NVariables,\
			  innerh,\
			  innert,\
			  outerh,\
			  outert);

	}

	//start time stepping 
	double residual0;
	for(TimeStep = 0; TimeStep<MaxTimeStep;++TimeStep) {
		//initialize old variable 
		for(i=0;i<NVariables;++i) Variables0[i] = *(Variables[i]);

		NewtStep = 0;
		delT0 = delT;

		while(1) {

			if(NewtStep > 10 | residual*(NewtStep>2) > residual0) {

				delT *= 0.5;

				NewtStep = 0;

				for(i=0;i<NVariables;++i) *(Variables[i]) = Variables0[i];

			}

			residual0 = residual;
      
			GetArcLengths(innerh, innert,    \
							  outerh, outert);

			BendingEnergy(Variables,\
							  Variables0,\
							  dS,\
							  iEle, jEle, Ele,\
							  &Energy,\
							  EGrad,\
							  aHess, iHess, jHess); 

			ZeroOut(NVariables,Essential,iHess,jHess,aHess); //zero out essential boundary elements of hessian
			for(i=0; i < NVariables; ++i) EGrad[i] *= Essential[i]; //zero out essential boundary elements of gradient
      
			residual = 0.0;

			for(i=0;i<NVariables;++i) {
				dVariables[i] = 0.0;
				residual += EGrad[i]*EGrad[i];
			}

			residual = sqrt(residual);
      
			if(residual < NewtTol) break; //Newton iteration completed within tolerance
      
//			printf("Time Step = %g; Residual = %g\n", delT, residual);
			    
			//solve Ax = b
			//where A = Hess, x = Variables and b = EGrad

			PCG(aHess, jHess, iHess, NVariables,\
				 dVariables,\
				 EGrad,\
				 rr,\
				 dd,\
				 qq,\
				 1.0e-25);
      
			for(i=0;i<NVariables;++i) *(Variables[i]) -= dVariables[i]; 
      
			++NewtStep;

		}

		if(NewtStep < 4) {
      
			delT *= 2.0;
      
		}

		if(TimeStep % 25 ==0)
			printf("%d r: %g  Energy = %g\n", TimeStep, r,  Energy);

		sprintf(command,"echo %d %g >> %s", TimeStep, Energy, infofilename);  
    
		system(command);
    
		//save current profile to "profilefilename" 
		//place where profile is stored 
		FILE *profilefile;

		profilefile = fopen(profilefilename,"w");
		for(i=0;i<n;++i) {
			fprintf(profilefile,"%g %g %g %g %g %g %g %g\n",\
					  innerh[i].x,innerh[i].y,\
					  innert[i].x,innert[i].y,\
					  outerh[i].x,outerh[i].y,\
					  outert[i].x,outert[i].y);
		} 

		fclose(profilefile);

	}

	//write to profile.dat for shifting 

	FILE *profilefile;

	profilefile = fopen("profile.dat","w");
	for(i=0;i<n;++i) {
		fprintf(profilefile,"%g %g %g %g %g %g %g %g\n",\
				  innerh[i].x,innerh[i].y,\
				  innert[i].x,innert[i].y,\
				  outerh[i].x,outerh[i].y,\
				  outert[i].x,outert[i].y);
	} 
  
	fclose(profilefile);  
	printf("done\n");

  
	return(0);
}

int function(double var1,\
				 double var2,\
				 double var3,\
				 double *out) {

	*out = (var1 - var2)*(var1 - var2) + (var3 - var2)*(var3 - var2) + (var1 - var3)*(var1 - var3);

	return(0);
}


int Load(int NVariables,\
			struct vec *innerh,\
			struct vec *innert,\
			struct vec *outerh,\
			struct vec *outert) {  
  
	FILE * pfile = fopen(loadfilename, "r");
	int i;

	printf("Loading from %s\n",loadfilename);

	for(i=0;i<NVariables;++i) {
		fscanf(pfile,"%lg %lg %lg %lg %lg %lg %lg %lg\n",\
				 &(innerh[i].x),&(innerh[i].y),\
				 &(innert[i].x),&(innert[i].y),\
				 &(outerh[i].x),&(outerh[i].y),\
				 &(outert[i].x),&(outert[i].y));
	} 
	fclose(pfile);

	return(0);

}

int ZeroOut(int n, double *kill, int *ia, int *ja, double *aa) {
	//zero out the elements of the matrix in the index ste kill
	//and make the diagonal = 1 there

	int i,ii,j,jj,k;

	for(i=0;i<n;++i) {
		for(j=ia[i];j<ia[i+1];++j){
			jj = ja[j];
			//consider the i,jj element      
			aa[j] = kill[i]*kill[jj]*aa[j] + (kill[i]==0)*(kill[jj]==0)*(i==jj);
		}
	}
	return(0);
}

int ATimesB(int m, int n, int p, double *a, int *ia, int *ja, double *b, int *ib, int *jb, double **c, int **ci, int **cj, int mem, int eval) {

	//added row size n, so we can also use only upper left parts 
	//of matrices 

	int maxdim;
  
	maxdim = (m>=n)*m + (m<n)*n;

	if(mem) {
		*ci = malloc(sizeof(int)*(maxdim + 1)); 
		if(*ci == NULL) return(0); 
	}

	int i; 
	for(i = 0; i <= maxdim; ++i) (*ci)[i] = -1;

	int j, jj, k, kk, nnz;
	/*count total occurences of connections
	  ci keeps track of number of occurences per row*/
	nnz = 0;
	for(i = 0; i < m; ++i) {
		for(j = ia[i]; j < ia[i+1]; ++j) {
			jj = ja[j];

			if(jj<n) {
				for(k = ib[jj]; k < ib[jj+1]; ++k){
					kk = jb[k];
					if(kk < p) {
						//kk already seen? then ci[kk] = i 
						if((*ci)[kk] != i) {
							++nnz;
							(*ci)[kk] = i; //kk seen this turn, don't count again
						}
					}//end if
				}
			}//end if

		}
	}

	if(mem) {
		*cj = malloc(sizeof(int)*nnz);
		if(*cj == NULL) return(0);
	}
	if(mem&&eval) {
		*c = malloc(sizeof(double)*nnz);
		if(*c == NULL) return(0);
	}


	if(eval) {
		for(i = 0; i < nnz; ++i) (*c)[i] = 0.0;
	}

	int l, ll, event;
	(*ci)[0] = 0;
	//go through again and retreive products
	for(i = 0; i < m; ++i) {
		(*ci)[i+1] = (*ci)[i];
		for(j = ia[i]; j < ia[i+1]; ++j) {
			jj = ja[j];

			if(jj<n) {
				for(k = ib[jj]; k < ib[jj+1]; ++k){
					kk = jb[k];
					    
					if(kk<p) {
						//kk already seen? then add product
						event = 0;
						for(l = (*ci)[i]; l < (*ci)[i+1]; ++l) {
							ll = (*cj)[l];
							if(kk == ll) {
								event = 1;  
								if(eval) (*c)[l] = (*c)[l] + a[j]*b[k];
							}
						}               
						//kk not seen? then iterate ia and add to ja and add product
						if(event == 0) {
							(*cj)[(*ci)[i+1]] = kk;
							(*ci)[i+1] = (*ci)[i+1]+1; 
							if(eval) (*c)[l] = (*c)[l] + a[j]*b[k];
						}
					} //end if
				}
			}//end if
		}
	}

	return(0);

}

int Transpose(int m, int n, double *a, int *ia, int *ja, double **b, int **ib, int **jb) {

	int i, j, jj, k, nnz;
	double temp;
	nnz = ia[m];
   
	int mem = 1;
	int eval = 0;

	if(mem&&eval) *b = malloc(sizeof(double)*nnz); 
	if(mem) *ib = malloc(sizeof(int)*(n+1)); 
	if(mem) *jb = malloc(sizeof(int)*nnz); 

	for(i=0;i<n+1;++i) (*ib)[i] = 0;
  
	for(i=0;i<m;++i) {
		for(j=ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			(*ib)[jj+1] = (*ib)[jj+1] + 1;
		}
	}

	(*ib)[0] = 0;
	for(i=0;i<n;++i) {
		(*ib)[i+1] = (*ib)[i+1] + (*ib)[i];    
	}

	for(i=0;i<nnz;++i) (*jb)[i]=-1;

	for(i=0;i<m;++i) {
		for(j=ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			if(eval) temp = a[j];
			for(k=(*ib)[jj];k<(*ib)[jj+1];++k){

				/*printf("%d,%d,%d,%d, %d,%d\n",i,jj,(*ib)[jj],k,(*ib)[jj+1],(*jb)[k]);*/
				if((*jb)[k]==-1) {
					(*jb)[k] = i;
					if(eval) (*b)[k] = temp;
					break;
				}
				/*printf("%d,%d,%d,%d,%d,%d\n",i,jj,(*ib)[jj],k,(*ib)[jj+1],(*jb)[k]);*/
			}
		}
	}
	return(0);
}


int ShowMatrix(int m, int n, double *a, int *ia, int *ja, char dot) {

	int i, j, jj; 
	double values[n];
  
	if(dot=='M') {
    
		for(i=0;i<m;++i) {
			for(j=ia[i];j<ia[i+1];++j) {
				printf("%d %d %g;\n", i, ja[j], a[j]);
			}
		}

	}

	else if(dot=='s'||dot=='S') {
		printf("\n%d by %d with %d nonzeros\n",m,n,ia[m]);

		printf("ia = [");
		for(i=0;i<m+1;++i){
			printf("%d ",ia[i]);
		}
		printf("];\nja= [");
		for(i=0;i<m;++i){
			for(j=ia[i];j<ia[i+1];++j){
				printf("%d ",ja[j]);       
			}
			printf(" | ");
		}
		printf("];\na= [");
		for(i=0;i<m;++i){
			for(j=ia[i];j<ia[i+1];++j){
				if(dot =='S') printf("%02.2g ",a[j]);
			}
			printf(" | ");
		}
		printf("];\n");
	}
	else {
		for(i=0;i<m;++i) {
			for(j=0;j<n;++j) values[j]=0.0;
			for(j=ia[i];j<ia[i+1];++j){
				jj=ja[j];
				if(dot == '.') values[jj]=1.0;
				else values[jj]=a[j];
			}
			for(j=0;j<n;++j) {
				if(dot=='.') {
					if(values[j]==1.0) printf("n ");
					else printf(". ");
				}
				else {
					if(dot=='D') printf("%-10.5f ",values[j]);      
					if(dot=='d') printf("%-6.1f ",values[j]);      

				}
			}
			printf(";");    printf("\n");

		}
	}

	return(0);

}

int PutAij(int i, int j, double *a, int *ia, int *ja, double z) {

	int k, fail;
	fail = 1;
	for(k=ia[i];k<ia[i+1];++k) {
		if(j == ja[k]) {
			a[k] += z;    
			fail = 0;
		}
	}
	if(fail) {
		printf("column %d not found in row %d\n",j,i);
		exit(1);
	}
	return(0);
}


int Elements(int *iEle, int *jEle, double *aE,\
				 int *ia, int *ja, double *aa) {



	return(0);

}


int BendingEnergy(double ** Variables,\
						double *Variables0,\
						double *dS,\
						int * iEle, int * jEle, int *Ele,\
						double *E,\
						double *EGrad,\
						double * aHess, int * iHess, int * jHess) {

  
	int e,f,i,j,jj,k; //e,f for element

	*E = 0.0;

	//zero out the gradient and the Hessian 
	for(i=0;i<NVariables;++i) {
		EGrad[i] = 0.0;
		for(j=iHess[i];j<iHess[i+1];++j) {
			aHess[j] = 0.0;
		}
	}

	//NPositions is the number of nodes in an element
	//the elements are the rhomboidal regions spanning two
	//pairs of lipids
  
	int In0[NPositions], In[NPositions], Position[NPositions];
	double LocalVariables[NPositions];
	double LocalVariables0[NPositions];

	double h = 0.0001; //perturbation parameter 
	double z[9], z0;
	for(e = 0; e < NElements; ++e) {

		//find out indeces of vertices and their position
		//relative to the current element
		for(i = iEle[e];i < iEle[e+1];++i) {
			In0[i-iEle[e]] = jEle[i];
			Position[i-iEle[e]] = Ele[i];
		}
    
		//put the Indeces of the current element's
		//vertices in order 
		//put the values of Variables pertinent to this
		//element into a local variable LocalVariables
		for(i = 0; i < NPositions; ++i) {
			In[Position[i]] = In0[i]; 
			LocalVariables[Position[i]] = *(Variables[In0[i]]); 
			LocalVariables0[Position[i]] = Variables0[In0[i]]; 
		}

		//    printf("%d : ", e);
		//    for(i=0;i<NPositions;++i) printf("%d ", In[i]);
		//    printf("\n");
		//calculate the bending energy just by itself

		BendingEnergyDensity(
			LocalVariables[0],\
			LocalVariables[1],\
			LocalVariables[2],\
			LocalVariables[3],\
			LocalVariables[4],\
			LocalVariables[5],\
			LocalVariables[6],\
			LocalVariables[7],\
			LocalVariables[8],\
			LocalVariables[9],\
			LocalVariables[10],\
			LocalVariables[11],\
			LocalVariables[12],\
			LocalVariables[13],\
			LocalVariables[14],\
			LocalVariables[15],\
			LocalVariables0[0],\
			LocalVariables0[1],\
			LocalVariables0[2],\
			LocalVariables0[3],\
			LocalVariables0[4],\
			LocalVariables0[5],\
			LocalVariables0[6],\
			LocalVariables0[7],\
			LocalVariables0[8],\
			LocalVariables0[9],\
			LocalVariables0[10],\
			LocalVariables0[11],\
			LocalVariables0[12],\
			LocalVariables0[13],\
			LocalVariables0[14],\
			LocalVariables0[15],\
			dS[e],\
			&z0); 
		*E += z0;

		//go through indeces of current element 
		//pairwise 
		for(i=0;i<NPositions;++i) {     

			//now calculate the gradient w.r.t. the ith position 
			for(k=0; k<3; ++k) {    
				BendingEnergyDensity(LocalVariables[0] + h*(i==0)*(k-1),\
											LocalVariables[1] + h*(i==1)*(k-1),\
											LocalVariables[2] + h*(i==2)*(k-1),\
											LocalVariables[3] + h*(i==3)*(k-1),\
											LocalVariables[4] + h*(i==4)*(k-1),\
											LocalVariables[5] + h*(i==5)*(k-1),\
											LocalVariables[6] + h*(i==6)*(k-1),\
											LocalVariables[7] + h*(i==7)*(k-1),\
											LocalVariables[8] + h*(i==8)*(k-1),\
											LocalVariables[9] + h*(i==9)*(k-1),\
											LocalVariables[10] + h*(i==10)*(k-1),\
											LocalVariables[11] + h*(i==11)*(k-1),\
											LocalVariables[12] + h*(i==12)*(k-1),\
											LocalVariables[13] + h*(i==13)*(k-1),\
											LocalVariables[14] + h*(i==14)*(k-1),\
											LocalVariables[15] + h*(i==15)*(k-1),\
											LocalVariables0[0],\
											LocalVariables0[1],\
											LocalVariables0[2],\
											LocalVariables0[3],\
											LocalVariables0[4],\
											LocalVariables0[5],\
											LocalVariables0[6],\
											LocalVariables0[7],\
											LocalVariables0[8],\
											LocalVariables0[9],\
											LocalVariables0[10],\
											LocalVariables0[11],\
											LocalVariables0[12],\
											LocalVariables0[13],\
											LocalVariables0[14],\
											LocalVariables0[15],\
											dS[e],\
											&(z[k])); 
			}
			
			EGrad[In[i]] += (z[2]-z[0])/(2.0*h); 

			//stop here       
			for(j=0;j<NPositions;++j) {

				for(k=0;k<9;++k) z[k] = 0.0;
				if(i!=j) {//the finite difference formula for the cases i = j and i~=j are differenct
					for(k=0; k<9; ++k) {    

						// how the 9 possible changes are tabulated 
						/*
						               -h 0 h i 
											  -h  0 1 2
											     0  3 4 5 
												  h  6 7 8 = k   
												     j
						*/

						BendingEnergyDensity(LocalVariables[0] + h*((i==0)*(k%3-1) + (j==0)*(k/3-1)), \
													LocalVariables[1] + h*((i==1)*(k%3-1) + (j==1)*(k/3-1)), \
													LocalVariables[2] + h*((i==2)*(k%3-1) + (j==2)*(k/3-1)), \
													LocalVariables[3] + h*((i==3)*(k%3-1) + (j==3)*(k/3-1)), \
													LocalVariables[4] + h*((i==4)*(k%3-1) + (j==4)*(k/3-1)), \
													LocalVariables[5] + h*((i==5)*(k%3-1) + (j==5)*(k/3-1)), \
													LocalVariables[6] + h*((i==6)*(k%3-1) + (j==6)*(k/3-1)), \
													LocalVariables[7] + h*((i==7)*(k%3-1) + (j==7)*(k/3-1)), \
													LocalVariables[8] + h*((i==8)*(k%3-1) + (j==8)*(k/3-1)), \
													LocalVariables[9] + h*((i==9)*(k%3-1) + (j==9)*(k/3-1)), \
													LocalVariables[10] + h*((i==10)*(k%3-1) + (j==10)*(k/3-1)), \
													LocalVariables[11] + h*((i==11)*(k%3-1) + (j==11)*(k/3-1)), \
													LocalVariables[12] + h*((i==12)*(k%3-1) + (j==12)*(k/3-1)), \
													LocalVariables[13] + h*((i==13)*(k%3-1) + (j==13)*(k/3-1)), \
													LocalVariables[14] + h*((i==14)*(k%3-1) + (j==14)*(k/3-1)), \
													LocalVariables[15] + h*((i==15)*(k%3-1) + (j==15)*(k/3-1)), \
													LocalVariables0[0],\
													LocalVariables0[1],\
													LocalVariables0[2],\
													LocalVariables0[3],\
													LocalVariables0[4],\
													LocalVariables0[5],\
													LocalVariables0[6],\
													LocalVariables0[7],\
													LocalVariables0[8],\
													LocalVariables0[9],\
													LocalVariables0[10],\
													LocalVariables0[11],\
													LocalVariables0[12],\
													LocalVariables0[13],\
													LocalVariables0[14],\
													LocalVariables0[15],\
													dS[e],\
													&(z[k])); 
					}
				}

				if(i==j) { //see above for i~=j
					for(k=0; k<3; ++k) {    
						BendingEnergyDensity(LocalVariables[0] + h*(i==0)*(k-1),\
													LocalVariables[1] + h*(i==1)*(k-1),\
													LocalVariables[2] + h*(i==2)*(k-1),\
													LocalVariables[3] + h*(i==3)*(k-1),\
													LocalVariables[4] + h*(i==4)*(k-1),\
													LocalVariables[5] + h*(i==5)*(k-1),\
													LocalVariables[6] + h*(i==6)*(k-1),\
													LocalVariables[7] + h*(i==7)*(k-1),\
													LocalVariables[8] + h*(i==8)*(k-1),\
													LocalVariables[9] + h*(i==9)*(k-1),\
													LocalVariables[10] + h*(i==10)*(k-1), \
													LocalVariables[11] + h*(i==11)*(k-1), \
													LocalVariables[12] + h*(i==12)*(k-1), \
													LocalVariables[13] + h*(i==13)*(k-1), \
													LocalVariables[14] + h*(i==14)*(k-1), \
													LocalVariables[15] + h*(i==15)*(k-1), \
													LocalVariables0[0],       \
													LocalVariables0[1],       \
													LocalVariables0[2],       \
													LocalVariables0[3],       \
													LocalVariables0[4],       \
													LocalVariables0[5],       \
													LocalVariables0[6],       \
													LocalVariables0[7],       \
													LocalVariables0[8],       \
													LocalVariables0[9],       \
													LocalVariables0[10],       \
													LocalVariables0[11],       \
													LocalVariables0[12],       \
													LocalVariables0[13],       \
													LocalVariables0[14],       \
													LocalVariables0[15],       \
													dS[e],       \
													&(z[k])); 
					}
				}
				
				// how the 9 possible changes are tabulated 
				/* i~=j
					    -h 0 h i 
						  -h  0 1 2
						    0  3 4 5 
							 h  6 7 8 = k   
							   j

				*/

				/* i==j
					   
					   -h 0 h i
						    0 1 2

				*/
				z0 = (i==j)*(z[2] + z[0] - 2.0*z[1])/(h*h) \
					+  (i!=j)*(z[8]-z[6]-z[2]+z[0])/(4.0*h*h);

				PutAij(In[i], In[j], aHess, iHess, jHess, z0);

			}

		}
    
	}
  
	return(0);
  
}

int GetArcLengths(struct vec * innerh, struct vec * innert,\
						struct vec * outerh, struct vec * outert) { 
  
	double S0, P0, Q0, T0;
	struct vec S, P, Q, T;

	//calculate the average arclength sizes...
	// S0 along inner head layer
	// P0 along inner tail layer
	// Q0 along outer tail layer 
	// T0 along outer head layer 

	int i;

	S0 = 0.0; 
	P0 = 0.0; 
	Q0 = 0.0; 
	T0 = 0.0;

	for(i=0; i < n-1; ++i) {

		S.x  = innerh[i+1].x - innerh[i].x;                 S.y  = innerh[i+1].y - innerh[i].y;
		S0 += sqrt(S.x*S.x + S.y*S.y);
    
		P.x  = innert[i+1].x - innert[i].x;                 P.y  = innert[i+1].y - innert[i].y;
		P0 += sqrt(P.x*P.x + P.y*P.y);

		Q.x  = outert[i+1].x - outert[i].x;                 Q.y  = outert[i+1].y - outert[i].y;
		Q0 += sqrt(Q.x*Q.x + Q.y*Q.y);

		T.x  = outerh[i+1].x - outerh[i].x;                 T.y  = outerh[i+1].y - outerh[i].y;
		T0 += sqrt(T.x*T.x + T.y*T.y);    

	}
 
	ihs0 = S0/((double) n-1);
	its0 = P0/((double) n-1);
	ots0 = Q0/((double) n-1);
	ohs0 = T0/((double) n-1);

	if(CALLIBRATION==1) printf("running arclength calibration\n%g %g %g %g\n", S0, P0, Q0, T0);

	return(0);
}

int BendingEnergyDensity(double ihlx, double ihly, double itlx, double itly, \
								 double ohlx, double ohly, double otlx, double otly, \
								 double ihrx, double ihry, double itrx, double itry, \
								 double ohrx, double ohry, double otrx, double otry, \
								 double ihlx0, double ihly0, double itlx0, double itly0, \
								 double ohlx0, double ohly0, double otlx0, double otly0, \
								 double ihrx0, double ihry0, double itrx0, double itry0, \
								 double ohrx0, double ohry0, double otrx0, double otry0, \
								 double dSConstraint,\
								 double *E) {

	struct vec innerhl, innerhr, innertl, innertr; 
	struct vec outerhl, outerhr, outertl, outertr;
	struct vec P, Q, S, T, M, N, Dl, dl, Dr, dr, El, el, Er, er;
	double innerJ, outerJ; //inner and outer area element

	/*

    l = left, r = right
    D = inner director, E = outer director
    d = inner unit director, e = outer unit director
    M = inner normal, N = outer normal
    S = inner tangert, T = outer tangent
    J = inner area element, K = outer area element 

	*/


	//  printf("input to bde\n %g %g %g %g %g %g %g %g\n", ihlx, ihly, itlx, itly, ohlx, ohly, otlx, otly);

	innerhl.x = ihlx; innerhl.y = ihly; innertl.x = itlx; innertl.y=itly;
	innerhr.x = ihrx; innerhr.y = ihry; innertr.x = itrx; innertr.y=itry;
	outerhl.x = ohlx; outerhl.y = ohly; outertl.x = otlx; outertl.y=otly;
	outerhr.x = ohrx; outerhr.y = ohry; outertr.x = otrx; outertr.y=otry;

	Dl.x = innertl.x - innerhl.x;                 Dl.y = innertl.y - innerhl.y;
	dl.x = Dl.x/sqrt(Dl.x*Dl.x + Dl.y*Dl.y);      dl.y = Dl.y/sqrt(Dl.x*Dl.x + Dl.y*Dl.y);
	Dr.x = innertr.x - innerhr.x;                 Dr.y = innertr.y - innerhr.y;
	dr.x = Dr.x/sqrt(Dr.x*Dr.x + Dr.y*Dr.y);      dr.y = Dr.y/sqrt(Dr.x*Dr.x + Dr.y*Dr.y);

	El.x = outertl.x - outerhl.x;                 El.y = outertl.y - outerhl.y;
	el.x = El.x/sqrt(El.x*El.x + El.y*El.y);      el.y = El.y/sqrt(El.x*El.x + El.y*El.y);
	Er.x = outertr.x - outerhr.x;                 Er.y = outertr.y - outerhr.y;
	er.x = Er.x/sqrt(Er.x*Er.x + Er.y*Er.y);      er.y = Er.y/sqrt(Er.x*Er.x + Er.y*Er.y);

	S.x  = innerhr.x - innerhl.x;                 S.y  = innerhr.y - innerhl.y;
	double dS = sqrt(S.x*S.x + S.y*S.y);
	S.x = S.x/dS;                                 S.y  = S.y/dS;
	M.x = -S.y;                                   M.y  = S.x; //inner normals M

	T.x  = outerhr.x - outerhl.x;                 T.y  = outerhr.y - outerhl.y;
	double dT = sqrt(T.x*T.x + T.y*T.y);
	T.x = T.x/dT;                                 T.y  = T.y/dT;
	N.x = T.y;                                    N.y  = -T.x; //outer normals N

	//tail arc elements for calibration steps 
	P.x  = innertr.x - innertl.x;                 P.y  = innertr.y - innertl.y;
	double dP = sqrt(P.x*P.x + P.y*P.y);
	Q.x  = outertr.x - outertl.x;                 Q.y  = outertr.y - outertl.y;
	double dQ = sqrt(Q.x*Q.x + Q.y*Q.y);

	double R1 = innerhl.x, R2 = innerhr.x, h = innerhr.y - innerhl.y; 
	innerJ = pi*(R1 + R2)*sqrt((R2-R1)*(R2-R1)+h*h); //inner differential area
	R1 = outerhl.x, R2 = outerhr.x, h = outerhr.y - outerhl.y; 
	outerJ = pi*(R1 + R2)*sqrt((R2-R1)*(R2-R1)+h*h); //outer differential area

	double Splayinner = -( S.x*(dr.x - dl.x) + S.y*(dr.y-dl.y) )/dS \
		- (dl.x + dr.x)/(innerhl.x + innerhr.x);    //(-1) for the orientation 
	double Splayouter = ( T.x*(er.x - el.x) + T.y*(er.y-el.y) )/dT\
		+ (el.x + er.x)/(outerhl.x + outerhr.x);    
  
	//stretch of outer molecules
	double ihl = sqrt((ihlx-itlx)*(ihlx-itlx) + (ihly-itly)*(ihly-itly)); //left/right inner lipid length 
	double ihr = sqrt((ihrx-itrx)*(ihrx-itrx) + (ihry-itry)*(ihry-itry));

	double ohl = sqrt((ohlx-otlx)*(ohlx-otlx) + (ohly-otly)*(ohly-otly)); //left/right outer lipid length 
	double ohr = sqrt((ohrx-otrx)*(ohrx-otrx) + (ohry-otry)*(ohry-otry));
  

	//VERSION 1 extension not allowed
	//  double iStl = (sqrt(0.0001+(ihl-h0)*(ihl-h0)) + ihl-h0)/2.0; 
	//  double iStr = (sqrt(0.0001+(ihr-h0)*(ihr-h0)) + ihr-h0)/2.0; 

	//double oStl = (sqrt(0.0001+(ohl-h0)*(ohl-h0)) + ohl-h0)/2.0; 
	//  double oStr = (sqrt(0.0001+(ohr-h0)*(ohr-h0)) + ohr-h0)/2.0; 

	//VERSION 2 neither extension nor compression of lipid allowed

	double iStl = ihl-h0; 
	double iStr = ihr-h0; 

	double oStl = ohl-h0; 
	double oStr = ohr-h0; 

	double iStretch = iStl*iStl + iStr*iStr;  
	double oStretch = oStl*oStl + oStr*oStr;  

	//compression 

	double innerVol = Volume(innerhl, innerhr,\
									 innertl, innertr);

	double outerVol = Volume(outerhl, outerhr,\
									 outertl, outertr);

	double innerCom = (innerVol/innerJ - h0)*(innerVol/innerJ - h0)/(2.0*h0);
	double outerCom = (outerVol/outerJ - h0)*(outerVol/outerJ - h0)/(2.0*h0);

	//tilt of molecules 
  
	double idx = (ihlx + ihrx - itlx - itrx)/2.0;//average of directors 
	double idy = (ihly + ihry - itly - itry)/2.0;
	double ih = sqrt(idx*idx + idy*idy); //their length 

	double odx = (ohlx + ohrx - otlx - otrx)/2.0;
	double ody = (ohly + ohry - otly - otry)/2.0;
	double oh = sqrt(odx*odx + ody*ody);
  
	//tilt = |d/|d| - n|^2 = 2 - 2d.n/|d|
	double Tiltinner = 2.0 -2.0*(idx*M.x+idy*M.y)/ih; 
	double Tiltouter = 2.0 -2.0*(odx*N.x+ody*N.y)/oh; 

	//SADDLE SPLAY???
  
	double Extension = (innertl.x - outertl.x)*(innertl.x - outertl.x)\
		+ (innertl.y - outertl.y)*(innertl.y - outertl.y)\
		+ (innertr.x - outertr.x)*(innertr.x - outertr.x)\
		+ (innertr.y - outertr.y)*(innertr.y - outertr.y);
  
	//  printf("%g %g %g %g\n", innerJ, outerJ, Splayinner, Splayouter);

	if(CALLIBRATION==0) {  //business as usual 

		*E = (\
			(Kb*(Splayinner-innerK0)*(Splayinner-innerK0) - Kb*innerK0*innerK0)*innerJ \
			+ (Kb*(Splayouter-outerK0)*(Splayouter-outerK0) - Kb*outerK0*outerK0)*outerJ \
			+ (innerJ+outerJ)*Extension\
			+ Kt*innerJ*Tiltinner + Kt*outerJ*Tiltouter +\
			+ 3.0*Kb*innerJ*iStretch + 3.0*Kb*outerJ*oStretch\
			+ 0.0*Kb*innerJ*innerCom + 0.0*Kb*outerJ*outerCom\
			);
    
		//only penalize for deviations from the arclength differential
		//along the inner and outer tail portions

		*E += (dP - dSConstraint)*(dP - dSConstraint)\
			+  (dQ - dSConstraint)*(dQ - dSConstraint);

		//add time component
    
		*E += (\
			innerJ*((ihlx - ihlx0)*(ihlx - ihlx0) + (ihly - ihly0)*(ihly - ihly0) + \
					  (ihrx - ihrx0)*(ihrx - ihrx0) + (ihry - ihry0)*(ihry - ihry0) + \
					  (itlx - itlx0)*(itlx - itlx0) + (itly - itly0)*(itly - itly0) + \
					  (itrx - itrx0)*(itrx - itrx0) + (itry - itry0)*(itry - itry0)) + \
			outerJ*((ohlx - ohlx0)*(ohlx - ohlx0) + (ohly - ohly0)*(ohly - ohly0) + \
					  (ohrx - ohrx0)*(ohrx - ohrx0) + (ohry - ohry0)*(ohry - ohry0) + \
					  (otlx - otlx0)*(otlx - otlx0) + (otly - otly0)*(otly - otly0) + \
					  (otrx - otrx0)*(otrx - otrx0) + (otry - otry0)*(otry - otry0)) \
			)/(2.0*delT);
	}

	else {  //run the callibration steps 

		*E = (dS - ihs0)*(dS - ihs0)\
			+  (dP - its0)*(dP - its0)\
			+  (dQ - ots0)*(dQ - ots0)\
			+  (dT - ohs0)*(dT - ohs0);
    
		*E += (\
			((ihlx - ihlx0)*(ihlx - ihlx0) + (ihly - ihly0)*(ihly - ihly0) + \
			 (ihrx - ihrx0)*(ihrx - ihrx0) + (ihry - ihry0)*(ihry - ihry0) + \
			 (itlx - itlx0)*(itlx - itlx0) + (itly - itly0)*(itly - itly0) + \
			 (itrx - itrx0)*(itrx - itrx0) + (itry - itry0)*(itry - itry0)) + \
			((ohlx - ohlx0)*(ohlx - ohlx0) + (ohly - ohly0)*(ohly - ohly0) + \
			 (ohrx - ohrx0)*(ohrx - ohrx0) + (ohry - ohry0)*(ohry - ohry0) + \
			 (otlx - otlx0)*(otlx - otlx0) + (otly - otly0)*(otly - otly0) + \
			 (otrx - otrx0)*(otrx - otrx0) + (otry - otry0)*(otry - otry0)) \
			)/(2.0*delT);  
		*E *= dx;

	}

	return(0);
}


//calculate the volume of a rotated frustral wedge via 3 point trapezoidal rule
double Volume(struct vec hl, struct vec hr,\
				  struct vec tl, struct vec tr) {

	double V, J, K, L;

	double R1 = hl.x, R2 = hr.x, h = hr.y - hl.y; 
	J = pi*(R1 + R2)*sqrt((R2-R1)*(R2-R1)+h*h); //head differential area
	double T1 = tl.x, T2 = tr.x, k = tr.y - tl.y; 
	L = pi*(T1 + T2)*sqrt((T2-T1)*(T2-T1)+k*k); //tail differential area
	double S1 = (R1+T1)/2.0, S2 = (R2+T2)/2.0, j = (h+k)/2.0; 
	K = pi*(S1 + S2)*sqrt((S2-S1)*(S2-S1)+j*j); //middle differential area

	double D = sqrt((hl.x + hr.x - tl.x - tr.x)*(hl.x + hr.x - tl.x - tr.x)/4.0 + \
						 (hr.y + hl.y - tr.y - tl.y)*(hr.y + hl.y - tr.y - tl.y)/4.0);  //distance between midpoints
  
	return(D*(0.25*J + 0.5*K + 0.25*L));
}

/*
  !*************************************************************************************
  ! 
  !  PROGRAM: PCG = Preconditioned Conjugate Gradient (with Gauss Seidel Preconditioner)
  !  AUTHOR:  Rolf J. Ryham, Rice University 2009
  !  PURPOSE: Perform Conjugate Gradient with Symmetric Guass Seidel Preconditioning
  !
  !                      x <= inv(A)b              
  !
  !           where A is an n by n symmetric matrix.  Here the preconditioner M is 
  !
  !                      M = (D + L) inv(D) (D + U)              
  !
  !           where A = D + L + U  with diagonal D, upper triangular part U and lower 
  !           triangular part L.
  !
  !  DESCRIPTION: The matrix A is stored in compressed row storage
  !
  !           a   are the nonzero entries in A
  !           ja  are the column numbers of the nonzero entries
  !           ia  are the indices of a and ja corresponding to rows  
  !               ia(1) = 1 and ia(n+1) = nnz+1 where nnz are the number   
  !               of nonzeros
  !
  !  EXAMPLE OF COMPRESSED ROW STORAGE: 
  ! 
  !      | 1 0 3 |     a  = [1, 3, -1, 2]
  !   A =| 0 0 0 |     ja = [1, 3, 1, 2]
  !      |-1 2 0 |     ia = [1, 3, 3, 5] 
  !
  !  SUBROUTINES: 
  !          
  !       subroutine SGS(aa, ja, ia, n, b)   !Symmetric Guass Seidel
  !       subroutine SSOR(aa, ja, ia, n, b,w)   !Symmetric Succesive Over Relaxation
  !
  !***********************************************************************************
  !  
  !  INPUTS/OUTPUTS:
  !
  !  double aa             length nnz nonzero entries of A
  !  integer ja            length nnz column of entries of A
  !  integer ia            length n + 1 indices of nonzero entries in A rowwise 
  !  integer n             dimension of A
  !  double b              length n vector right hand side of the equation Ax = b
  !  double x              length n vector approximate solution of the equation Ax = b
  !  double r, d, q        length n vector search directions, r(1) = residual on exit
  !                        r(2) = number of iterations on exit
  !  
  !*********************************************************************************** */

int PCG(double *aa, int *ja, int *ia, int n, double *x, double *b, double *r, double *d, double *q, double tol) {

	double al, be, d0, dnew, dold;
	double temp, temp2;
	int i, j, k, l;

	int imax, restart;
	double e, w;
  
	w = 1.99;
	restart = (int)(sqrt((double)n));
	imax = n + 1;

	e = tol;
	i = 1;

	//      !r = b - Ax
	//      !d = r

	for(j=0;j<n;++j){
		temp = 0.0;
		for(k = ia[j];k < ia[j+1];++k){
			l = ja[k];
			temp = temp + aa[k]*x[l];
		}    
		r[j] = b[j] - temp;
		d[j] = r[j];
	}
	//      !d <= inv(M)r (d = r) 

	//SGS(aa, ja, ia, n, d);
	SSOR(aa, ja, ia, n, d, w);

	//      !d0 = r^T r
	dnew = 0.0;
  
	for(j = 0;j<n;++j){
		dnew = dnew + d[j]*r[j];
	}
	d0 = dnew + 1.0;

	if(dnew == 0.0) return(0);

	//      !PCG iteration
	while(i<imax) {  
		temp2 = 0.0;
		//      !q = Ad  
		for(j = 0;j<n;++j) {
			temp = 0.0;
			for(k = ia[j];k< ia[j+1];++k) {
				l = ja[k];
				temp = temp + aa[k]*d[l];
			}
			// !temp = q(j)
			//         !d = d
			q[j] = temp;
			temp2 = temp2 + d[j]*q[j] ; //!d^T q
		}

		//  !al = dnew/d^T q
		al = dnew/temp2;

		//    !x = x + al.d 
		for(j = 0;j<n;++j) {
			x[j] = x[j] + al*d[j];
		}

		if(i%restart==0) {
			//         !restart the residual
			//         !r = b - Ax
			for(j=0;j<n;++j) {
				temp = 0.0;
				for(k = ia[j];k< ia[j+1]; ++k) {
					l = ja[k];
					temp = temp + aa[k]*x[l];
				}
				r[j] = b[j] - temp;
			}

		}
		else { 
			//    !update the residual
			for(j=0;j<n;++j) {
				r[j] = r[j] - al*q[j];
			}
		}


		for(j = 0;j<n;++j) {
			q[j] = r[j];
		}

		//     !q <= inv(M)r (s = r) 

		//SGS(aa, ja, ia, n, q);       
		SSOR(aa, ja, ia, n, q, w);

		dold = dnew;      
		//      !d0 = r^T r
		dnew = 0.0;
		for(j = 0;j<n;++j){
			dnew = dnew + r[j]*q[j];
		}
		//printf("%d,%g<%g %d %g\n", i, dnew,d0*e, fabs(dnew) < fabs(d0*e),e);
		if(fabs(dnew) < fabs(d0*e)) { 
			break ;
		}



		be = dnew/dold;

		for(j = 0;j<n;++j){
			d[j] = q[j] + be*d[j];
		}

		i = i + 1;

		//printf("%d/%d,%g<%g %d\n", i, imax,dnew,d0*e, dnew < d0*e);
	}

	//!residual upon exit
	for(j = 0;j<n;++j){
		temp = 0.0;
		for(k = ia[j];k< ia[j+1];++k) {
			l = ja[k];
			temp = temp + aa[k]*x[l];
		}
		r[j] = b[j] - temp;         
	}

	temp = 0.0;
	for(j=0;j<n;++j) {
		temp = temp + r[j]*r[j];
	}
  
	//  !l2 residual norm upon exit
  
	r[0] = sqrt(temp);
	r[1] = (double)i;
  
	return(0);
  
}


/****************************************************************************
! 
!  PROGRAM: SSOR = Symmetric Succesive Over Relaxation
!  AUTHOR:  Rolf J. Ryham, Fordham University 2010
!  PURPOSE: Perform one iteration of the SSOR
!           operation 
!
!                      1/(2-w)(D/w + L) inv(D/w)(D/w + U) x = b => x              
!
!           where A = D + L + U is an n by n symmetric matrix with diagonal 
!           D, upper triangular part U and lower triangular part L.
!  
!  DESCRIPTION: The matrix A is stored in compressed row storage
!
!           a   are the nonzero entries in A
!           ja  are the column numbers of the nonzero entries
!           ia  are the indices of a and ja corresponding to rows  
!               ia(1) = 1 and ia(n+1) = nnz+1 where nnz are the number   
!               of nonzeros
!
!  EXAMPLE OF COMPRESSED ROW STORAGE: 
! 
!      | 1 0 3 |     a  = [1, 3, -1, 2]
!   A =| 0 0 0 |     ja = [1, 3, 1, 2]
!      |-1 2 0 |     ia = [1, 3, 3, 5] 
!
!****************************************************************************
!  
!  INPUTS/OUTPUTS:
!
!  double aa             length nnz
!  integer ja            length nnz
!  integer ia            length n + 1
!  integer n             dimension of A
!  double b              vector to be multiplies on in and the product on out         
!  double w              relaxation constant
!
!****************************************************************************/

int SSOR(double *aa, int *ja, int *ia, int n, double *b, double w) { //  Symmetric Guass Seidel
    
	int i, j, jj, k ;
	double temp, sum;

/**************************************
!
! PART I: solve for (D/w + L)*b = (2-w)b
!
! do i = 1,n
!    for j: A(i,j) != 0, j < i 
!        sum = sum + A(i,j)*b(j)
!    end
!    b(i) = w*((2-w)b(i) - sum)/A(i,i)
! end 
!
!**************************************/

	for(i = 0;i < n;++i) {

		for(j = ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			if(jj==i) {
				temp = aa[j]; //    !temp = A(i,i)
			}
		}

		//!find j: A(i,j).NE.0, j > i
		sum = 0.0;
		for(j = ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			if(jj<i) {   
				sum  += aa[j]*b[jj]; 
			}
		}
    
		b[i] = w*((2-w)*b[i] - sum)/temp;
    
	}

	/*      !**************************************
			    !
				   ! PART I: replace b with inv(I + w inv(D)U)*b
					  !
					    ! do i = n,1,-1
						   !    b(i) = b(i)/A(i,i)
							  !    for j: A(j,i) != 0, j < i 
							    !        b(j) = b(j) - A(j,i)*b(i)
								   !    end
									  ! end 
									    !
										 !************************************** */
  
	for(i = n-1; i>-1;--i) {    
		for(j = ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			if(jj==i) {
				temp = aa[j];//     !temp = A(i,i)
			}
		}
		// !find j: A(i,j).NE.0, j > i
		sum = 0.0;
		for(j = ia[i];j<ia[i+1];++j) {
			jj = ja[j];
			if(jj>i){
				sum += aa[j]*b[jj] ;
			}
		}

		b[i] = b[i] - w*sum/temp;

	}

	return(0);

}

