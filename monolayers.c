#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//__________Two Monolayers___________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//------parameters---------------------------
double l  = 20.0;       // length

double pi = 3.141592;   // pi is pi 

int    n  = 200;        // subdivisions
double r  = 2.0;        // radius
double R  = 2.0;        // height
double h0 = 1.0;        // distance from midplane to neutral surface
double alpha   = 2.0;    // base point to lumen
double Kt      = 1.0;    // tilt modulus
double Kb      = 1.0;    // bending modulus
double Jinner  = 0.00;   // spontaneous curvature 
double Jouter  = 0.00;    // spontaneous curvature 
double ST      = 0.0000; // surface tenstion 
double STinner = 0.0000; // surface tenstion 
double STouter = 0.0000; // surface tenstion 


//------convergence parameters---------------

int NDecayTest  = 10;
double DecayTol = 1.0e-3; 
double DecayC   = 1.0e-3 + 1.0; 
double EOld; 

//--------------------------------------------


//------------------------------------------
double phie;

struct quint {
  
  double *D;
  double *U;
  double *UU;
  double *L;
  double *LL;

};

int SurfEngDen (double rhol, double rhor, double phi, double *w);

int BendingEnergy (double *Phi, double *rho, double *theta, double *zeta, double  *W, \
  	   double *GradRho, double * GradTheta, double *GradZeta, double *dA, \
		   struct quint JacRho, struct quint JacTheta, struct quint JacZeta);

int BendingEngDen(double phil, double phir, double phim,	\
		  double rhol, double rhor, double rhom,	\
		  double thetal, double thetar,			\
		  double zetal,  double zetar,			\
		  double *w);

int QuintDiagSolve(double *D,				\
		   double *U,  double *L,		\
		   double *UU, double *LL,		\
		   double *b,  int n) ;

int main(int argc, char **argv) { 

  if(argc!=3) { 
    printf("\nargc = %d; Please provide input:  r  R\n\n", argc);
    exit(1);
  } 
     r = atof(argv[1]); 
     R = atof(argv[2]);

  printf("r: %g R: %g\n", r, R);  

  double delT = 0.1;
  int maxiter = 5;
  double errortol = 1.0e-7;

  //--------------------------------------------------------------------------------
  double *rhotorus= malloc(sizeof(double)*(n+1));
  double *rho0    = malloc(sizeof(double)*(n+1));
  double *rhoT    = malloc(sizeof(double)*(n+1));
  double *rho     = malloc(sizeof(double)*(n+1));
  double *delrho  = malloc(sizeof(double)*(n+1));
  double *dA      = malloc(sizeof(double)*(n+3));

  double *GradRho   = malloc(sizeof(double)*(n+1));
  
  struct quint JacRho; 

  JacRho.D   = malloc(sizeof(double)*(n+1));
  JacRho.U   = malloc(sizeof(double)*(n+1));
  JacRho.L   = malloc(sizeof(double)*(n+1));
  JacRho.UU  = malloc(sizeof(double)*(n+1));
  JacRho.LL  = malloc(sizeof(double)*(n+1));

  //--------------------------------------------------------------------------------
  // theta defined on edges, as opposed to rho 
  // which is defined on the points 

  double *theta0    = malloc(sizeof(double)*(n));
  double *theta     = malloc(sizeof(double)*(n));
  double *deltheta  = malloc(sizeof(double)*(n));

  double *GradTheta   = malloc(sizeof(double)*(n));
  
  struct quint JacTheta; 

  JacTheta.D   = malloc(sizeof(double)*(n));
  JacTheta.U   = malloc(sizeof(double)*(n));
  JacTheta.L   = malloc(sizeof(double)*(n));
  JacTheta.UU  = malloc(sizeof(double)*(n));
  JacTheta.LL  = malloc(sizeof(double)*(n));


  double *zeta0    = malloc(sizeof(double)*(n));
  double *zeta     = malloc(sizeof(double)*(n));
  double *delzeta  = malloc(sizeof(double)*(n));

  double *GradZeta   = malloc(sizeof(double)*(n));
  
  struct quint JacZeta; 

  JacZeta.D   = malloc(sizeof(double)*(n));
  JacZeta.U   = malloc(sizeof(double)*(n));
  JacZeta.L   = malloc(sizeof(double)*(n));
  JacZeta.UU  = malloc(sizeof(double)*(n));
  JacZeta.LL  = malloc(sizeof(double)*(n));

  //--------------------------------------------------------------------------------

  double *Phi = malloc(sizeof(double)*(n+1));
  double *PhiT = malloc(sizeof(double)*(n+1));
  
  double phi, phi0;
  double EMin, ETorus;
  char name[100];

  /*
  FILE *pfile; 
  pfile = fopen("profile.dat", "w");

  FILE *qfile; 
  qfile = fopen("energytemp.dat", "w");

  FILE *tfile; 
  tfile = fopen("profilee.dat", "w");
  */

  int i,j,k; 
  int I, J, swich;
  int iter = 0;
  double error;

  phie = pi - atan( R / ( l - (r + alpha) )); //end angle for phi
  
  double inc = (-R*cos(phie)/sin(phie) + R*pi/2.0)/((double)n);
  
  for(i = 0; inc*((double)i) <=  R*pi/2.0; ++i) Phi[i] = inc*((double) i)/R;
  
  for(j = i; j < n+1 ;++j) Phi[j] = atan(1.0 / (pi/2.0 - inc*((double) j)/R)) + pi;
    
  double PhiLast = Phi[n];
  for(k = 0; k < n+1; ++k) Phi[k] *= Phi[k]/PhiLast;


  FILE *sfile = fopen("profile.dat", "w");
  //initializing  torus
  for (i = 0; i < n+1; ++i) 
    {
      /*
	here, the phi is arranged such that 

	r+R - rho[i]cos(phi[i]) 
	
	is evenly spaced when rho[i] forms the profile of a line,
	namely 

	rho[i] = R/sin(phi(i))
      */

      phi  = Phi[i];

      /*
      //#####TORUS#####
      if (phi <= pi/2.0) 

	{
	  rhotorus[i] = R;
	}	

      else
      
	{
	  rhotorus[i] = R / sin(phi) ;
	}
      */

      //#####ELLIPSE#####
      if(phi <= pi/2.0)
	{
	  rhotorus[i] = 1 / ( sqrt(  cos(phi)*cos(phi)/( alpha*alpha ) + sin(phi)*sin(phi)/(R*R) ));
	}
      else
	{
	  rhotorus[i] = ( R / sin(phi) );
	}
      
      
      fprintf(sfile,"%g %g %g %g\n",phi,rhotorus[i],(r+R)-rhotorus[i]*cos(phi),rhotorus[i]*sin(phi));  
    }
  fclose(sfile);


  for(i = 0; i < n; ++i)  theta0[i] = 0.0;
  for(i = 0; i < n; ++i)  zeta0[i] = 0.0;
    
  for(i=0; i<n+1; ++i) rhoT[i] = 0.0;




  //torus initialized

  BendingEnergy (Phi, rhotorus, theta, zeta, &ETorus,	\
		 GradRho, GradTheta, GradZeta,	dA,	\
		 JacRho, JacTheta, JacZeta);     
  double Coef = 10.0;
  ETorus *= Coef;
  //  printf("energy torus = %16.16g \n", ETorus);
  printf("%16.16g\n", ETorus);

  /*

  char LoadThis[80]; 
  strcpy(LoadThis,"monolayers_NoTilt200/profile.r");
  strcat(LoadThis, argv[1]);
  strcat(LoadThis, "_R");
  strcat(LoadThis, argv[2]);
  strcat(LoadThis, ".dat");
  

  //usage of fscanf 
  FILE *pfile = fopen("monolayers_NoTilt200/profile.r1.5_R1.5.dat", "r"); //open with "r" to read 
    double dummy; 
  
  for(i=0;i<n+1;++i)
    {

      //the formatted outputs of fscanf go to a pointer 
      //since we don't need the last two columns of profilee.dat, 
      //we dump the values into a dummy variable 
      
      //rho + i is the pointer to rho[i]

      //we're printing the values scanned in to compare whether they are the same
      //as the ones that were printed by fprintf 

      fscanf(pfile, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg"
	     ,Phi+i, rho0+i, theta0+i,zeta0+i,			\
	     &dummy,&dummy, &dummy,&dummy, &dummy,&dummy,	\
	     &dummy,&dummy, &dummy,&dummy);
      
      //      printf("%g %g %g %g\n",Phi[i], rho0[i],theta0[i],zeta0[i]);

    }
  
  fclose(pfile);
  
  
  
  */

  //set initial data to simply the torus 
  for(i=0;i<n+1;++i) rho0[i] = rhotorus[i];

  /*  

  //bend up the tip 
  //Phi[n] = Phi[n-1] + (Phi[n] - Phi[n-1])*0.90;

  double Eng;
  BendingEnergy (Phi, rho0, theta0, zeta0, &Eng,	\
		 GradRho, GradTheta, GradZeta, dA,	\
		 JacRho, JacTheta, JacZeta);     


  printf("energy loaded = %16.16g \n", Eng);

  //---------END TEST-------------------------

  */

  while(iter < maxiter) 

    {

    error = errortol+1.0;

    for( i = 0; i < (n+1); ++i ) rho[i] = rho0[i];
    for( i = 0; i < (n); ++i ) theta[i] = theta0[i];
    for( i = 0; i < (n); ++i )  zeta[i] = zeta0[i];
   

    int NewtStep = 0;

    //----------- Rho Newton start -----------------------------------------------

    while(error > errortol) 
      {

	//      SurfaceEnergy (rho, &W, GradW, JacWD, JacWU, JacWL);

	BendingEnergy (Phi, rho, theta, zeta, &EMin,		\
		       GradRho, GradTheta, GradZeta, dA,	\
		       JacRho, JacTheta, JacZeta);
	
	//set rhs to equation	
	for(i=0;i<(n+1); ++i) 
	  {	    
	    delrho[i]   = -(rho[i]   - rho0[i]   + delT*GradRho[i]);
	  }
      
	//add identity to Jacobian
	
	for(i = 0; i < (n+1); ++i) 
	  {
	    
	    JacRho.D[i]  = delT*JacRho.D[i] + 1.0;
	    JacRho.U[i]  = delT*JacRho.U[i];
	    JacRho.L[i]  = delT*JacRho.L[i];
	    JacRho.UU[i] = delT*JacRho.UU[i];
	    JacRho.LL[i] = delT*JacRho.LL[i];
	    
	  }
	
	QuintDiagSolve(JacRho.D,			\
		       JacRho.U,  JacRho.L,		\
		       JacRho.UU, JacRho.LL,		\
		       delrho,  n+1);

	error = 0.0;
	
	for (i=0; i < (n+1); i++)
	  {
	    rho[i]   += delrho[i];	  
	    error    += delrho[i]*delrho[i];
	  }      
	
	error = sqrt(error); 

	++NewtStep;

	if(NewtStep > 12) { 
	  delT *= 0.7;
	  for( i = 0; i < (n+1); ++i ) rho[i] = rho0[i];
	  NewtStep -= 1;
	}
	
      }
    
    if(NewtStep < 6)  delT = 1.5*delT;
    if(delT > 0.5) delT = 0.5;

    //------------ Rho Newton done ----------------------------------------------
    if((iter%25)==0);
    printf( "%d %16.16g %g err = %g delt = %g %g %g\n", iter, EMin , ETorus, error, delT, DecayC/sqrt((double)iter), DecayTol);	       

    //----------- Theta  Newton start --------------------------------------------

    NewtStep = 0;
    double delS = 0.001;
    //never do this
    if(iter>0){
      for(l = 0; k < delT/delS; ++k) {	
	error = errortol + 1.0;
	while(error > errortol){
	  	
	// SurfaceEnergy (rho, &W, GradW, JacWD, JacWU, JacWL);

	BendingEnergy (Phi, rho, theta, zeta, &EMin,		\
		       GradRho, GradTheta, GradZeta, dA,	\
		       JacRho, JacTheta, JacZeta);
	
	//set rhs to equation
	
	for(i=0;i<(n); ++i) 
	  {	    
	    deltheta[i] = -(theta[i] - theta0[i] + delS*GradTheta[i]);	    
	    delzeta[i]  = -(zeta[i]  - zeta0[i]  + delS*GradZeta[i]);	    
	  }
      
	//add identity to Jacobian
	
	for(i = 0; i < (n); ++i) 
	  {
	    
	    JacTheta.D[i]  = delS*JacTheta.D[i] + 1.0;
	    JacTheta.U[i]  = delS*JacTheta.U[i];
	    JacTheta.L[i]  = delS*JacTheta.L[i];
	    JacTheta.UU[i] = delS*JacTheta.UU[i];
	    JacTheta.LL[i] = delS*JacTheta.LL[i];

	    JacZeta.D[i]  = delS*JacZeta.D[i] + 1.0;
	    JacZeta.U[i]  = delS*JacZeta.U[i];
	    JacZeta.L[i]  = delS*JacZeta.L[i];
	    JacZeta.UU[i] = delS*JacZeta.UU[i];
	    JacZeta.LL[i] = delS*JacZeta.LL[i];
	
	  }
	
	QuintDiagSolve(JacTheta.D,			\
		       JacTheta.U,  JacTheta.L,		\
		       JacTheta.UU, JacTheta.LL,	\
		       deltheta, n);

	QuintDiagSolve(JacZeta.D,			\
		       JacZeta.U,  JacZeta.L,		\
		       JacZeta.UU, JacZeta.LL,	\
		       delzeta, n);
      
	error = 0.0;
	
	for (i=0; i < (n); i++)
	  {
	    theta[i] += deltheta[i];	  
	    zeta[i]  += delzeta[i];	  
	    error    += deltheta[i] * deltheta[i];
	    error    += delzeta[i]  * delzeta[i];
	  }      
	
	error = sqrt(error); 
	//	printf( "subiter = %d, %g/%g, error = %g, delT = %g\n", k, EMin , ETorus, error, delS);
	++NewtStep;
	
	}
      }
    }
    
    //------------ Theta Newton done ----------------------------------------------

    ++iter; //printf("iter = %d maxiter = %d\n", iter, maxiter);

    for(i=0; i<n+1; i++) rho0[i] = rho[i];
    for(i=0; i<n; i++) theta0[i] = theta[i];
    for(i=0; i<n; i++) zeta0[i]  = zeta[i];
    

    if(iter > NDecayTest) {

      DecayC = pow((double)iter,1.5)*(EOld-EMin);

    }
    
    if(DecayC/sqrt((double)iter) < DecayTol) ;;//break ; //break from time stepping       

    EOld = EMin;  

    }//end of time iteration, minimum "found"
  iter = 0;

    
  //write torus energy
  char TorusEnergy[80];
  strcpy(TorusEnergy,"TorusEnergy.r");
  strcat(TorusEnergy, argv[1]);
  strcat(TorusEnergy, "._R");
  strcat(TorusEnergy, argv[2]);
  strcat(TorusEnergy, ".dat");
  FILE *tfile = fopen(TorusEnergy, "w");
  fprintf(tfile,"%g",ETorus);
  fclose(tfile);


  
  char str[80]; 
  strcpy(str,"profile.r");
  strcat(str, argv[1]);
  strcat(str, "._R");
  strcat(str, argv[2]);
  strcat(str, ".dat");
  FILE *pfile = fopen(str, "w");
  
  for (i = 0; i < n+1; ++i) 
    {

      phi = Phi[i];

      double TauX, TauY, NX, NY,nx,ny,nx2,ny2;
      j = i;
      if(i==n) j = n-1;
      double rhom = (rho0[j+1] + rho0[j])/2.0;
      double drho = (rho0[j+1] - rho0[j])/(Phi[j+1]-Phi[j]);
      TauX = (-drho*cos(phi) + rhom*sin(phi))/sqrt(drho*drho + rhom*rhom);
      TauY = ( drho*sin(phi) + rhom*cos(phi))/sqrt(drho*drho + rhom*rhom);
      NX   = -TauY; 
      NY   =  TauX;
 
      nx = NX*cos(theta[j]) - NY*sin(theta[j]);
      ny = NX*sin(theta[j]) + NY*cos(theta[j]);

      nx2 = -(NX*cos(zeta[j]) - NY*sin(zeta[j]));
      ny2 = -(NX*sin(zeta[j]) + NY*cos(zeta[j]));
      
      fprintf(pfile,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",	\
	      phi, rho0[i], theta0[j], zeta0[j],			\
	      (r+R) - rho0[i]*cos(phi), rho0[i]*sin(phi),		\
      	      NX,  NY,							\
      	      nx,  ny,							\
      	      nx2, ny2,							\
      	      (r+R) - rhotorus[i]*cos(phi), rhotorus[i]*sin(phi));
      
    }

   fclose(pfile);


  char str2[80];
  strcpy(str2, "energy.r");
  strcat(str2, argv[1]);
  strcat(str2, "._R");
  strcat(str2, argv[2]);
  strcat(str2, ".dat");
  
  FILE *qfile = fopen(str2, "w");
  fprintf(qfile, "%16.16g", EMin);
  fclose(qfile);
  
  
  char done[80];
  strcpy(done, "Done r: ");
  strcat(done, argv[1]);
  strcat(done, " R: ");
  strcat(done, argv[2]);
  printf("%s", done);
  
 

  return(0); 

}



int SurfaceEnergy (double *rho, double *W,		\
		   double *GradW,			\
		   double *JacWD, double *JacWU, double *JacWL)
{
  double dphi  = ( pi / 2.0 ) / ((double)n);
  double phi, phi0;
  double rhom, drho;
  double w0, w1,w2, w[9];
  double h = 1.0e-8, dS;
  
  int i,j,k,iter; 
  
  *W = 0.0;
  
  for(i=0;i<(n+1);++i) 
    {
    
      GradW[i] = 0.0;
      JacWD[i] = 0.0;   
      JacWL[i] = 0.0;
      JacWU[i] = 0.0;

    }


  for(i = 0; i< n; i++) // n intervals!!!!
    {
      phi = ((double)i + 0.5) * dphi;

      SurfEngDen (rho[i]+h, rho[i+1]+h, phi, w);
      SurfEngDen (rho[i]+h, rho[i+1],   phi, w+1);
      SurfEngDen (rho[i]+h, rho[i+1]-h, phi, w+2);
      SurfEngDen (rho[i],   rho[i+1]+h, phi, w+3);
      SurfEngDen (rho[i],   rho[i+1],   phi, w+4);
      SurfEngDen (rho[i],   rho[i+1]-h, phi, w+5);
      SurfEngDen (rho[i]-h, rho[i+1]+h, phi, w+6);
      SurfEngDen (rho[i]-h, rho[i+1],   phi, w+7);
      SurfEngDen (rho[i]-h, rho[i+1]-h, phi, w+8);
     
      *W += w[4];

      GradW[i]   += (w[1]-w[7])/(2.0*h);
      GradW[i+1] += (w[3]-w[5])/(2.0*h);

  
      JacWD[i]   += (w[1] + w[7] - 2.0*w[4])/(h*h);
      JacWD[i+1] += (w[3] + w[5] - 2.0*w[4])/(h*h);

      //W_i,i+1

      JacWU[i]   += (w[0] - w[6] - w[2] + w[8])/(4.0*h*h);
      JacWL[i]   += (w[0] - w[6] - w[2] + w[8])/(4.0*h*h);   
      
    }

  GradW[0] = 0.0;
  GradW[n] = 0.0;

  JacWU[0] = 0.0;
  JacWL[n-1] = 0.0;
  return(0);
}

int SurfEngDen (double rhol, double rhor, double phi, double *w)
{
  double rhom, drho, dphi;

  rhom = (rhol + rhor) / 2.0;
  drho = (rhor - rhol) / dphi; 

  *w = 2.0*pi*dphi*				\
    (r+R - rhom * cos(phi)) *			\
    sqrt(drho * drho + rhom * rhom);
      
  return(0);
}

int BendingEnergy (double *Phi, double *rho, double *theta, double *zeta, double *W, \
		   double *GradRho, double * GradTheta, double *GradZeta, double *dA, \
		   struct quint JacRho, struct quint JacTheta, struct quint JacZeta)
{
  double dphi  = ( pi / 2.0 ) / ((double)n);
  double phi, phi0;
  double phil, phim, phir;
  double rhom, rhol, rhor;
  double thetam, thetal, thetar;
  double zetam, zetal, zetar;
  double w0, w1,w2, w[3][3][3];
  double h = 1.0e-8, hpm[3];
  
  int i,j,k,l,iter; 
  
  *W = 0.0;
  
  hpm[0] = -h;
  hpm[1] =  0.0;
  hpm[2] =  h;

  for(i=0;i<(n+1);++i) {
    
    GradRho[i] = 0.0;
    JacRho.D[i] = 0.0;   
    JacRho.L[i] = 0.0;
    JacRho.U[i] = 0.0;
    JacRho.LL[i] = 0.0;
    JacRho.UU[i] = 0.0; 

  }
      
  for(i=0;i<(n);++i) {

      GradTheta[i] = 0.0;
      JacTheta.D[i] = 0.0;   
      JacTheta.L[i] = 0.0;
      JacTheta.U[i] = 0.0;
      JacTheta.LL[i] = 0.0;
      JacTheta.UU[i] = 0.0;

      GradZeta[i] = 0.0;
      JacZeta.D[i] = 0.0;   
      JacZeta.L[i] = 0.0;
      JacZeta.U[i] = 0.0;
      JacZeta.LL[i] = 0.0;
      JacZeta.UU[i] = 0.0;

    }

  //  printf("-----\n");
  for(i = 1; i < n; i++) //n+1 points!!!!!
    {

      phim = Phi[i];
      phil = Phi[i-1];
      phir = Phi[i+1];

      rhom = rho[i];      
      rhol = rho[i-1];      
      rhor = rho[i+1];

      thetal = theta[i-1];      
      thetar = theta[i];

      zetal = zeta[i-1];      
      zetar = zeta[i];

      double dphil = phim - phil;
      double dphir = phir - phim;
      double drhol = (rhom - rhol) / dphil;
      double drhor = (rhor - rhom) / dphir;      
      double drho  = (drhol + drhor) / 2.0;  

      dA[i] = 2.0*pi*					\
	( r + R - rhom * cos(phim)) *			\
	sqrt(drho * drho + rhom * rhom);
      
      ///- Form Hessian and Gradient w.r.t. rho-------------------------------------------------------------

     for(j=0;j<3;++j) {
	for(k=0;k<3;++k) {
	  for(l=0;l<3;++l) {
	    BendingEngDen(phil, phir, phim,				\
			  rhol+hpm[j], rhor+hpm[k], rhom + hpm[l],	\
			  thetal, thetar,				\
			  zetal,  zetar,				\
			  &(w[j][k][l]));
	  }
	}
      }

      *W += w[1][1][1];

      GradRho[i]   += (w[1][1][2] - w[1][1][0])/(2.0*h); // middle
      GradRho[i-1] += (w[2][1][1] - w[0][1][1])/(2.0*h); // left
      GradRho[i+1] += (w[1][2][1] - w[1][0][1])/(2.0*h); // right 

      JacRho.D[i]   += (w[1][1][2] +  w[1][1][0] - 2.0*w[1][1][1])/(h*h); // i,i
      JacRho.D[i-1] += (w[2][1][1] +  w[0][1][1] - 2.0*w[1][1][1])/(h*h); // i-1,i-1
      JacRho.D[i+1] += (w[1][2][1] +  w[1][0][1] - 2.0*w[1][1][1])/(h*h); // i+1,i+1

      JacRho.U[i-1]   += (w[2][1][2] -  w[0][1][2] - w[2][1][0] + w[0][1][0])/(4.0*h*h); // i-1,i; d/dp in p[i-1] with respect to changing p[i] 
      JacRho.U[i]     += (w[1][2][2] -  w[1][2][0] - w[1][0][2] + w[1][0][0])/(4.0*h*h); // i,i+1

      JacRho.UU[i-1]  += (w[2][2][1] -  w[0][2][1] - w[2][0][1] + w[0][0][1])/(4.0*h*h); // i-1,i+1

      JacRho.L[i-1]   += (w[2][1][2] -  w[0][1][2] - w[2][1][0] + w[0][1][0])/(4.0*h*h); // i-1,i
      JacRho.L[i]     += (w[1][2][2] -  w[1][2][0] - w[1][0][2] + w[1][0][0])/(4.0*h*h); // i,i+1

      JacRho.LL[i-1]  += (w[2][2][1] -  w[0][2][1] - w[2][0][1] + w[0][0][1])/(4.0*h*h); // i-1,i+1

      ///- Form Hessian and Gradient w.r.t. theta ------------------------------------------------------------

      for(j=0;j<3;++j) {
	for(k=0;k<3;++k) {
	  BendingEngDen(phil, phir, phim,				\
			rhol, rhor, rhom,				\
			thetal+hpm[j], thetar+hpm[k],			\
			zetal,  zetar,					\
			&(w[j][k][1]));
	}
      }

      GradTheta[i-1] += (w[2][1][1] - w[0][1][1])/(2.0*h); // left
      GradTheta[i]   += (w[1][2][1] - w[1][0][1])/(2.0*h); // right 
      
      JacTheta.D[i-1] += (w[2][1][1] +  w[0][1][1] - 2.0*w[1][1][1])/(h*h); // i-1,i-1
      JacTheta.D[i]   += (w[1][2][1] +  w[1][0][1] - 2.0*w[1][1][1])/(h*h); // i,i

      JacTheta.U[i-1]   += (w[2][2][1] -  w[0][2][1] - w[2][0][1] + w[0][0][1])/(4.0*h*h); // i-1,i
      JacTheta.L[i-1]   += (w[2][2][1] -  w[2][0][1] - w[0][2][1] + w[0][0][1])/(4.0*h*h); // i,i-1

      //in this set-up, JacTheta is actuarry tridiag (not sure what I was smoking 
      //not to see that before)

      for(j=0;j<3;++j) {
	for(k=0;k<3;++k) {
	  BendingEngDen(phil, phir, phim,				\
			rhol, rhor, rhom,				\
			thetal, thetar,					\
			zetal +hpm[j], zetar+hpm[k],			\
			&(w[j][k][1]));
	}
      }

      GradZeta[i-1] += (w[2][1][1] - w[0][1][1])/(2.0*h); // left
      GradZeta[i]   += (w[1][2][1] - w[1][0][1])/(2.0*h); // right 
      
      JacZeta.D[i-1] += (w[2][1][1] +  w[0][1][1] - 2.0*w[1][1][1])/(h*h); // i-1,i-1
      JacZeta.D[i]   += (w[1][2][1] +  w[1][0][1] - 2.0*w[1][1][1])/(h*h); // i,i

      JacZeta.U[i-1]   += (w[2][2][1] -  w[0][2][1] - w[2][0][1] + w[0][0][1])/(4.0*h*h); // i-1,i
      JacZeta.L[i-1]   += (w[2][2][1] -  w[2][0][1] - w[0][2][1] + w[0][0][1])/(4.0*h*h); // i,i-1

    }

  dA[0] = 1.0;
  dA[n] = 1.0;

  // i = 0 : 1 0 0 ... 0 | 0 
  //rho bc.s
  JacRho.D[0]  = 1.0;
  JacRho.U[0]  = 0.0;
  JacRho.UU[0] = 0.0;

  GradRho[0] = 0.0;

  // i = 1: 1 -1 0 ... 0 | = 0 
  

  JacRho.L[0]  = 1.0;
  JacRho.D[1]  = -1.0;
  JacRho.U[1]  = 0.0;
  JacRho.UU[1] = 0.0;

  GradRho[1] = 0.0;


  // i = n; 0...0 0 1 | 0
  //rho b.c.'s
  JacRho.LL[n-2]  = 0.0;
  JacRho.L[n-1]  = 0.0;
  JacRho.D[n]  = 1.0;

  GradRho[n] = 0.0;

  // i = n-1: 0..0 -1 1 | = 0 
  

  JacRho.LL[n-3]  = 0.0;
  JacRho.L[n-2]   = 0.0;
  JacRho.D[n-1]   = -1.0;
  JacRho.U[n-1]   = 1.0;

  GradRho[n-1] = 0.0;

  //theta bc.s
  JacTheta.D[0]  = 1.0;
  JacTheta.U[0]  = 0.0;

  GradTheta[0] = 0.0;

  JacTheta.L[n-2]  = 0.0;
  JacTheta.D[n-1]  = 1.0;

  GradTheta[n-1] = 0.0;



  //zeta bc.s
  JacZeta.D[0]  = 1.0;
  JacZeta.U[0]  = 0.0;

  GradZeta[0] = 0.0;

  JacZeta.L[n-2]  = 0.0;
  JacZeta.D[n-1]  = 1.0;

  GradZeta[n-1] = 0.0;

  return(0);

}

int BendingEngDen(double phil, double phir, double phim,	\
		  double rhol, double rhor, double rhom,	\
		  double thetal, double thetar,			\
		  double zetal,  double zetar,			\
		  double *w)
{
  double drhol, drhor, drho, J, J2, delJ, rr, rl;

  double dphil = phim - phil;
  double dphir = phir - phim;
  double dphih = (dphil+dphir)/2.0;

  drhol = (rhom - rhol) / dphil;
  drhor = (rhor - rhom) / dphir;

  drho  = (drhol + drhor) / 2.0;  


  //------coordinates of rho surface-------
  double xl,xm,xr,yl,ym,yr;
  double mlx,mly,mrx,mry;
  
  xl = (r + alpha - rhol*cos(phil)); yl = rhol*sin(phil);// coords of rhol
  xm = (r + alpha - rhom*cos(phim)); ym = rhom*sin(phim);// coords of rhom
  xr = (r + alpha - rhor*cos(phir)); yr = rhor*sin(phir);// coords of rhor

  mlx = (xm + xl)/2.0; mly = (ym + yl)/2.0;// midpoint of rhol/rhom
  mrx = (xr + xm)/2.0; mry = (yr + ym)/2.0;// midpoint of rhom/rhor

  //--------------------------------------


  /*
  J = 2.0*pi*					\
    ( r + R - rhom * cos(phim)) *			\
    sqrt(drho * drho + rhom * rhom);
  */

  double rhorm, rholm;
  rhorm = (rhom + rhor) / 2.0;
  rholm = (rhom + rhol) / 2.0;
  double phirm = phim + dphir/2.0;
  double philm = phim - dphil/2.0;
  
  rr = (r + alpha -rhorm*cos(phirm))*drhor/sqrt(drhor * drhor + rhorm * rhorm);
  rl = (r + alpha -rholm*cos(philm))*drhol/sqrt(drhol * drhol + rholm * rholm);
  
  double k = (drho * drho + rhom * rhom) / (rhom * rhom);
  
  delJ = 2.0*pi*(-cos(phim)*sqrt(drho * drho + rhom * rhom)  +		\
		 ((r + alpha -rhom*cos(phim))*rhom) /sqrt(drho * drho + rhom * rhom) - \
		 (rr-rl)/(phirm-philm)  );
  
  double H = 0.5*((delJ/J)*sqrt(k));

  double thetam = 0.5*(thetal + thetar);
  double zetam  = 0.5*(zetal  + zetar);

  double faux = (thetar-thetam)*(thetar-thetam)/(dphir*dphir)	\
    + (thetam-thetal)*(thetam-thetal)/(dphil*dphil);



  
  //Tilt energy 
  
  double Tilt = 2.0 - 2.0*cos(thetam);
  double Tilt2 = 2.0 - 2.0*cos(zetam) ;



  //Divergence Energy 
  // div n = n_phi * tau_phi / sqrt(rho'*rho' + rho*rho) + n_psi * tau_psi / (r+R - rho cos(phi))

  //*********************************************************************
  //n_psi * tau_psi = first coordinate of n in (r,z) space

  double TauX, TauY, NX, NY;

  TauX = (-drho*cos(phim) + rhom*sin(phim))/sqrt(drho*drho + rhom*rhom);
  TauY = ( drho*sin(phim) + rhom*cos(phim))/sqrt(drho*drho + rhom*rhom);

  NX   = -TauY; 
  NY   =  TauX;
  
  //2 being for the zeta's here 
  double npsi  = NX*cos(thetam) - NY*sin(thetam);
  double npsi2 = -(NX*cos(zetam) - NY*sin(zetam));

  //*********************************************************************
  //n_phi * tau_phi = ...it's complicated 
  
  double TaulX, TaulY, NlX, NlY;
  double TaurX, TaurY, NrX, NrY;
  double nrX, nrY, nlX, nlY;
  double nrX2, nrY2, nlX2, nlY2;


  TaulX = (-drhol*cos(philm) + rholm*sin(philm))/sqrt(drhol*drhol + rholm*rholm);
  TaulY = ( drhol*sin(philm) + rholm*cos(philm))/sqrt(drhol*drhol + rholm*rholm);

  TaurX = (-drhor*cos(phirm) + rhorm*sin(phirm))/sqrt(drhor*drhor + rhorm*rhorm);
  TaurY = ( drhor*sin(phirm) + rhorm*cos(phirm))/sqrt(drhor*drhor + rhorm*rhorm);

  NlX   = -TaulY;   NlY   =  TaulX;
  NrX   = -TaurY;   NrY   =  TaurX;

  nlX = NlX*cos(thetal) - NlY*sin(thetal);
  nlY = NlX*sin(thetal) + NlY*cos(thetal);

  nrX = NrX*cos(thetar) - NrY*sin(thetar);
  nrY = NrX*sin(thetar) + NrY*cos(thetar);


  //2 being for the zeta's here 
  
  nlX2 = -(NlX*cos(zetal) - NlY*sin(zetal));
  nlY2 = -(NlX*sin(zetal) + NlY*cos(zetal));

  nrX2 = -(NrX*cos(zetar) - NrY*sin(zetar));
  nrY2 = -(NrX*sin(zetar) + NrY*cos(zetar));


  //-----neutral surface------------
  double hlx, hly, hrx, hry;// coords for point on neutral surface
  double hlx2, hly2, hrx2, hry2;// coords for zeta side
  double rhohl,rhohr, rhoh, drhoh;  
  double rhohl2,rhohr2, rhoh2, drhoh2;  

  hlx = mlx + h0*cos(thetal)*NlX;//coords of left neutral surface 
  hly = mly + h0*cos(thetal)*NlY; 
  
  double a = hlx-(r + alpha);
  double b = hly; 
  rhohl = sqrt(a*a + b*b);//length of rho left on neutral surface

  hrx = mrx + h0*cos(thetar)*NrX;//coords of right neutral surface
  hry = mry + h0*cos(thetar)*NrY; 

  a = hrx-(r + alpha);
  b = hry; 
  rhohr = sqrt(a*a + b*b);//length of rho right on neutral surface

  rhoh = (rhohr + rhohl)/2.0;//length of rho[i] on upper/inner neutral surface


  hlx2 = mlx - h0*cos(zetal)*NlX;//coords of left neutral surface 
  hly2 = mly - h0*cos(zetal)*NlY; //zeta

  a = hlx2-(r+alpha);
  b = hly2; 
  rhohl2 = sqrt(a*a + b*b);//length of rho left on lower/outer neutral surface

  hrx2 = mrx - h0*cos(zetar)*NrX;//coords of right neutral surface
  hry2 = mry - h0*cos(zetar)*NrY; //zeta

  a = hrx2-(r+alpha);
  b = hry2; 
  rhohr2 = sqrt(a*a + b*b);//length of rho right on lower/outer neutral surface

  rhoh2 = (rhohr2 + rhohl2)/2.0;//length of rho[i] on lower/outer neutral surface


  drhoh  = (rhohr  - rhohl )/dphih;
  drhoh2 = (rhohr2 - rhohl2)/dphih;


  double tauhx, tauhy, norm;
  double tauhx2, tauhy2, norm2;

  tauhx = hrx - hlx;
  tauhy = hry - hly;

  norm = sqrt(tauhx*tauhx + tauhy*tauhy);

  tauhx = tauhx/norm;
  tauhy = tauhy/norm;


  tauhx2 = hrx2 - hlx2;
  tauhy2 = hry2 - hly2;

  norm2 = sqrt(tauhx2*tauhx2 + tauhy2*tauhy2);

  tauhx2 = tauhx2/norm2;
  tauhy2 = tauhy2/norm2;

  double Nhx, Nhy;
  double Nhx2, Nhy2;

  Nhx = -tauhy;
  Nhy = tauhx;

  Nhx2 = -tauhy2;
  Nhy2 = tauhx2;

  double delnX,delnY;
  double delnX2,delnY2;

  delnX = (nrX - nlX)/(phirm-philm);
  delnY = (nrY - nlY)/(phirm-philm);

  delnX2 = (nrX2 - nlX2)/(phirm-philm);
  delnY2 = (nrY2 - nlY2)/(phirm-philm);

  J = 2.0*pi*					\
    ( r + alpha - rhoh * cos(phim)) *			\
    sqrt(drhoh * drhoh + rhoh * rhoh);

  J2 = 2.0*pi*					\
    ( r + alpha - rhoh2 * cos(phim)) *			\
    sqrt(drhoh2 * drhoh2 + rhoh2 * rhoh2);

  /*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*
  double nphi  = (nrX - nlX)*TauX/(phirm-philm)   + (nrY - nlY)*TauY/(phirm-philm);
  double nphi2 = (nrX2 - nlX2)*TauX/(phirm-philm) + (nrY2 - nlY2)*TauY/(phirm-philm);

  double divn  = nphi/sqrt(drho*drho + rhom*rhom)  +  npsi/(r + R - rhom*cos(phim));
  double divn2 = nphi2/sqrt(drho*drho + rhom*rhom) +  npsi2/(r + R - rhom*cos(phim));
  *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*/
  

  //$$$$$$$$$$$$$ neutral surfaces $$$$$$$$$$$$$$$$$$$$$
  double nphi  = (delnX)*tauhx   + (delnY)*tauhy;
  double nphi2 = (delnX2)*tauhx2 + (delnY2)*tauhy2;

  double divn  = nphi/sqrt(drhoh*drhoh + rhoh*rhoh)  +  npsi/(r + alpha - rhoh*cos(phim));
  double divn2 = nphi2/sqrt(drhoh2*drhoh2 + rhoh2*rhoh2) +  npsi2/(r + alpha - rhoh2*cos(phim));
  
  double dphi = (phir - phil)/2.0;
  
  *w = (     ( 2.0*Kb*(divn -Jouter)*(divn -Jouter) + 0.5 * Kt*(Tilt))*J \
	     + ( 0.0*Kb*(divn2-Jinner)*(divn2-Jinner) + 0.5 * Kt*(Tilt2))*J2 \
	     + 0.00*faux						\
	     + 0.00*H*H							\
	     + STinner*J2 + STouter*J )					\
    *dphi;

  return(0);

}


//--RhoTransform----------------------------------------------------------------------------
    /*--------------------------------------------------------------------------------------
      deltaRh is amount that rho[n] inc/dec in transformation.
      0*deltaRh is added to rho[0] while 1*deltaRh is added to rho[n] 
      factor increases with ratio of current phiT (T for phi of transformation) to new end phi

      Rh added to rho[0] since point at which phi is taken shifted by Rh
      factor for Rh goes from 1 to 0, similar to other factor but cos()

      the graphs look good what the changes rh or Rh are small such as less than .1
      I took bigger changes such as rh=1 or Rh=1 to make sure the graphs ends were at the right placs
      but you get a little dippy or buldgey like curve around pi/2 area but i dont think we will be making
      changes of that size?? hopefull this is a "good" transformation
     */
int RhoTransform (double *Phi, double *PhiT, double *rho, double *rhoT, double rh, double Rh)
{
  double phien, phiend, alpha;
  double rhon, rhonn, deltaRh, deltarh;
  int i,j,k;  

  phiend = pi - atan( (R) / ( l - (r + R) )); //end angle 
  phien = pi - atan( (R + Rh) / ( l - (r + rh + R + Rh) )); //end angle with new R/r
  
  double inc = (-(R+Rh)*cos(phien)/sin(phien) + (R+Rh)*pi/2.0)/((double)n);
  for(i = 0; inc*((double)i) <= (R+Rh)*pi/2.0; ++i)
    PhiT[i] = inc*((double) i)/(R+Rh);
  for(j = i; j < n+1 ;++j) 
    PhiT[j] = atan(1.0 / (pi/2.0 - inc*((double) j)/(R+Rh) )) + pi;
    
  rhon = R / sin(phiend); //length of rho at phie
  rhonn = (R + Rh) / sin(phien); //new length of rho with Rh
  
  deltaRh = rhonn - rhon; //change in rho[n]
    
  for(i=0; i<n+1; i++){
    
    alpha = ((PhiT[i] / phien)*pi/2.0);
    
    rhoT[i] = rho[i] + deltaRh*(PhiT[i]/phien) + Rh*cos(alpha);
    
    //rhoT[i] = rho[i] + deltaRh*(PhiT[i]/phien) + Rh*( ((double)n+1.0 -(double)i) / ((double)n+1.0) );
    //rhoT[i] = rho[i] + deltaRh*( (double)i / (double)n ) + Rh*( ((double)n+1.0 -(double)i) / ((double)n+1.0) );
  
    //printf("%10.10g\t%g\t%g\t%g\n",PhiT[i], rho[i],rhoT[i], rhoT[i]-rho[i]);
    
  }

  

return 0;
}




/*
D = diag
U = upper
L = lower
x = unknown
b = rhs
n = dimen
*/

int TriDiagSolve(double *D, double *U, double *L, double *b, int n) 
{

  int i,j,k;

  double alpha;

  for(i=0;i<(n-1);++i) 
    {

      alpha   = -L[i]/D[i];
      D[i+1] += U[i]*alpha;    
      b[i+1] += b[i]*alpha;

    }


  b[n-1] = b[n-1]/D[n-1];

  for(i=(n-2);i>(-1);--i) 
    {

    b[i] = (b[i] - b[i+1]*U[i])/D[i];

    }

  return(0);      

}



int QuintDiagSolve(double *D,				\
		   double *U,  double *L,		\
		   double *UU, double *LL,		\
		   double *b,  int n) 
{

  int i,j,k;

  double alpha, beta;

  for(i=0;i<(n-2);++i) 
    {

      alpha   = - L[i]/D[i];
      beta    = -LL[i]/D[i];

      b[i+1] += b[i]*alpha;
      b[i+2] += b[i]*beta;

      D[i+1] += U[i]*alpha;    
      U[i+1] += UU[i]*alpha;    

      L[i+1] += U[i]*beta;    
      D[i+2] += UU[i]*beta;    

    }


  // the i = n - 2 case 
  
  i = n-2;

  alpha = -L[i]/D[i]; 

  b[i+1] += b[i]*alpha;
  D[i+1] += U[i]*alpha;



  b[n-1] = b[n-1]/D[n-1];
  b[n-2] = (b[n-2] - b[n-1]*U[n-2])/D[n-2];

  for(i=(n-3);i>(-1);--i) 
    {

    b[i] = (b[i] - b[i+1]*U[i] - b[i+2]*UU[i])/D[i];

    }

  return(0);      

}
