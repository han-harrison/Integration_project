 
#include <nag.h>
#include <stdio.h>
#include <nag_stdlib.h>
#include <math.h>
#include <nagd01.h>
#include <nagx01.h>
#include <stdlib.h>
#include <string.h>
#include <nags.h>
#include <fcntl.h>
#define NDIM 2

static double QU1( Integer ndim,  double var[]);
static double QU2( Integer ndim,  double var[]);
static double QU3( Integer ndim,  double var[]);
static double QU4( Integer ndim,  double var[]);
static double GAU(double x);

static double FUN(double x);
static double FUN_X(double x);
static double FUN_Y(double y);
static double LOR1(double t);
static double LOR2(double t);
/*static double bessel(double y);*/

const double pi=3.1415926;
const double qe=1.602177e-19*3e9;
double alpha, k, kx, ky, lam_e, om, xv, x0, theta, psi, beta, gam, dum, D;
double a, b, epsl, sigma_x,sigma_y,sigma_t;

typedef struct {double r, i;} complex;
typedef struct {complex x, y, z;} vector;

vector Vcreate(double xr, double xi, double yr, double yi, double zr, double zi);
vector Vadd(vector v1, vector v2);
vector Vmultiply(vector v, complex c);
vector Vconjugate(vector v);
complex Vdot(vector v1, vector v2);

complex Ccreate(double r, double i);
complex Cadd(complex a, complex b);
complex Cmultiply(complex a, complex b);
complex Cdivide(complex a, complex b);
complex Cconjugate(complex a);
vector G1,G2, G, nv, nvp, nvv;
complex G1x, G1y, G1z, G2x, G2y, G2z, Gx, Gy, Gz;
complex a1, a2, a1t, a2t, f5c;

int main()
{
  char outfile_name[64];
  FILE *out_file;
  char line[200], ans,ans1, ans2;
  const double pi=3.1415926;
  const double qe=1.602177e-19*3e9;
  extern double alpha, k, kx, ky, xv, x0, theta, psi, beta, gam, dum, D;
  extern double a,b,sigma_x,sigma_y, sigma_t, lam_e, om;
  int i, ii, ic, j, jj, l, m, n, liw, flag=0;
  int n0, int_thet_1, int_thet_2, jstart, ind, indend;
  Integer ndim=NDIM, minpts, maxpts;
 /* grating parameters below */
  int grat, Nt;
  double lgr, lgr_e, frac, alpha1, alpha2, hi, el, el1, el2, xx, wid;
 /* beam parameters below  */ 
  int pulse_shp, ngauss, philmt;
  double y0,y0_av,charlen,nemit, emit, x_rad, n_e, n_e_t, i_av_p ;
  double pulse_len, len50, epsl, sigma_tdum, i_p, x, x_e, apb;
  double sigma_t1, sigma_t2, sigma_t3, sigma_t4, sigma_t5, sigma_teff, ssigma_x;
  double mu1, mu2, mu3, mu4, mu5, g1, g2, g3, g4, g5, tot_int;
  double fwhm_x, fwhm_y, fwhm, fwhm_t1, fwhm_t2, fwhm_t3, fwhm_t4, fwhm_t5;
 /* optical system parameters below */ 
  double r, r2, domeg_m, theta_u, theta_d, d_theta, d_phi, thet_1, thet_2;
  double phi, phi0, phistep, lam ;
  double aven, aven_p1, aven_p2, aven_single, aven_single_p1, aven_single_p2;
 /* integration parameters below */
  double ll[2], ul[2];                           /* limits for the y,z integral*/ 
  double lol, upl, eps, acc, res1, res2, res3, res4 ;
  double epsabs, epsrel, abserr;
  double fmax, tmax, t, f, gtop, gbot;
  double u1, u2, u3, u4, u5;
  double ff, f1, f2, sf3, f_y_c, f_y_s, f_z_c, f_z_s, z, u, s_coh_long;
  double f5, f5_p1, f5_p2; 
  double s_coh_xy, s_coh_xy_p1, s_coh_xy_p2, incans, cohans, err;  
  double f4[179], en_inc[179];
  double x0_bin[50], inc_bin[50], coh_bin[50], inc_bin_p1[50], coh_bin_p1[50], inc_bin_p2[50], coh_bin_p2[50];
 /* general variables */
  double temp1, temp2, temp3, temp4, ratio, sum1, sum2, sum3, sum4, dumsum, percent;
  double sum_single, sum_single_p1, sum_single_p2 ;
  double sum_th, sum_th_p1, sum_th_p2, sum_th_single, sum_th_single_p1, sum_th_single_p2;
  double sum_tr, sum_tr_p1, sum_tr_p2, sum_tr_single, sum_tr_single_p1, sum_tr_single_p2;
  double sum, sum_p1, sum_p2; 
 /* output arrays listed below */
  double thet[179], en_coh[179], f_p[179];
  double s_coh_z[179][10];
  double s_inc_x[179][10], s_inc_x_p1[179][10], s_inc_x_p2[179][10];
  double s_coh_x[179][10], s_coh_x_p1[179][10], s_coh_x_p2[179][10];
  double s_coh_tr[179][10], s_coh_tr_p1[179][10], s_coh_tr_p2[179][10];
  double single_elec[179][15], single_elec_p1[179][15], single_elec_p2[179][15];
  double bin_single[179][15], bin_single_p1[179][15], bin_single_p2[179][15];
  double f_b[179][10], f_b_p1[179][10], f_b_p2[179][10];
  double bin[179][15], bin_p1[179][15], bin_p2[179][15];
  float outbin[171][8], outbin_single[171][8];  
  double dumm[3], hyperg[500];
  double bin_1[179][2],dum70[179][2];
  Nag_QuadProgress QP;
  Nag_QuadSubProgress QPSUB1,QPSUB2;
  Nag_User COMM1, COMM2, COMM3, COMMWC1, COMMWC2, COMMWC3, COMMWC4, COMMWC5, COMMWC6, COMMWC7, COMMWC8;
  Nag_TrigTransform WT_FUNC;
  static NagError fail;
  
  complex G1x, G1y, G1z, G2x, G2y, G2z, Gx, Gy, Gz;
  complex a1, a2, a1t, a2t, f5c;
  	
  printf("\n"); printf("\n");
  printf("  This code calculates the radiated energy, BOTH POLARIZATIONS, (order 1) from a grating\n");
  printf("                          !!!OF FINITE WIDTH!!!\n");
  printf("\n");  printf("\n");
  printf(" !!!!!!!!!!!!!! NB1. bunch length is now defined by its FWHM !!!!!!!!!!!!! \n");
  printf(" NB2. This version gives the total energy, in both polarisation directions \n");
  printf("\n");
  printf(" NB3. scale factor f1 DOES NOT include n_e, which has been moved to S_inc & S_coh \n");
  printf(" NB4. this version provides data for the single electron yields, diff. & total \n");
  printf(" NB5. This is the preferred version for the calculation of polarisation \n");
  printf(" NB6. All output is deposited automatically in 5 text files \n");
  printf(" NB7. Note the clarification regarding the meaning of the quantity x0 \n\n\n");
   
  m=1;                                                                /*default order 1*/
  flag=0;
  if (!flag) printf("Default order is 1. Do you want to change it? yes(y) or no(n) \n");
  line[0]= 'x';
  do { 
	if (line[0] != 'x') printf(" ");
	if (fgets(line, 200, stdin) == NULL) { printf("\nErr: fgets.\n"); exit(1); }
	}
	while (line[0] != 'y' && line[0] != 'n'); 
  ans = line[0];
  if (ans=='y')
   {printf("Enter order required; default is order 1 \n");
    scanf("%d",&m);
    if (m<=0) m=1;  }
  printf(" \n");

  lgr=40;            /* grating length*/
  printf(" Enter grating period in mm\n");
  scanf("%lf",&el);
  Nt= (int)lgr/el;                    /*Nt is the number of teeth for this grating*/
  printf(" A standard grating of 40mm length would have %5d", Nt); printf(" teeth\n");
  printf(" Enter grating width in mm; if <1, it will be set to the default value of 20mm \n");
  scanf("%lf",&wid);
  if (wid<=1) wid=20.0;


  printf(" Now select type of grating \n");
  printf(" Type:\n");
  printf("      -3 for strip grating\n");
  printf("      -1 for reverse blaze echelle grating \n");
  printf("       0 for idealized strip grating       \n");
  printf("       1 for normal blaze echelle grating  \n");
  printf("       2 for sawtooth grating              \n");
  printf(" Default is the sawtooth grating           \n");
  scanf("%d",&grat);
  printf("\n");
  if (grat==-3) 
    {printf(" This is a strip grating \n");
    alpha1=0;
    printf(" Enter height of vertical facet, in mm  \n");
    scanf("%lf",&hi);
    alpha2=0;
    printf(" Enter the fraction of the period taken by the strip \n");
    scanf("%lf",&frac);
    if (frac<0.0 || frac>1.0) frac=0.5; 
    el1=el*frac;
    el2=el-el1;}
  else if (grat==-1)   
    {printf("This is a reverse blaze echelle grating \n");
    frac=0;
    el1=0;
    alpha1=pi/2;
    el2=el;
    printf(" Enter second blaze angle of grating in degrees; MUST BE greater than 90 \n");
    scanf("%lf",&alpha2);
    alpha2=alpha2*pi/180;
    hi=fabs(el2*tan(alpha2));} 
  else if (grat==0)
    {printf(" This is an idealized  strip grating \n");
    alpha1=0;
    printf(" Enter the fraction of the period taken by the strip \n");
    scanf("%lf",&frac);
    if (frac<0.0 || frac>1.0) frac=0.5; 
    el1=el*frac;
    alpha2=0;} 
  else if (grat==1)
    {printf("This is an echelle-type grating \n");
    frac=1.0;
    el1=el;
    el2=0;
    alpha2= -pi/2;
    printf("Enter first blaze angle of grating, in degrees \n");
    scanf("%lf",&alpha1);
    alpha1=alpha1*pi/180;
    hi=fabs(el*tan(alpha1));}
  else
    {
    /*if (el==1.5) alpha1=30;
    else if (el==1.0) alpha1=35;SETUP SEPCIFICALLY FOR SLAC ANALYSIS
    else if(el==0.5) alpha1=40;
    else {printf("Enter first blaze angle of grating in degrees \n");*/
    printf("Enter first blaze angle of grating in degrees \n");
    scanf("%lf",&alpha1);
    alpha1=alpha1*pi/180;
    frac=pow(cos(alpha1),2);
    printf("Fraction of period taken by first blaze is= "); 
    printf("%-4.3f", frac, "\n");
    el1=el*frac;
    el2=el-el1;
    hi=fabs(el1*tan(alpha1));                  /*height of grating indentation*/    
    alpha2= -(0.5*pi- alpha1);                /* definition of alpha2 changed, 10/1/10 */ 
    /*alpha2=0.5*pi+ alpha1;*/}
		/*BEGIN DEBUG SECTION - WHY DOESN'T ECHELLE GRATING WORK?*/
  printf("\n");
  printf("frac = %f\n", frac);
  printf("el = %f\n", el);
  printf("el1 = %f\n", el1);
  printf("el2 = %f\n", el2);
  printf("hi = %f\n", hi);
  printf("alpha1 = %f\n", alpha1);
  printf("alpha2 = %f\n", alpha2);
		/*END DEBUG SECTION - WHY DOESN'T ECHELLE GRATING WORK?*/
 
  r=230;          /* distance of Winston cone entrance from beam axis, in mm */
  r2=21;         /* Winston cone entrance diameter, in mm */
                      /* Max. entry angle for Winston cone is 6.3 deg.*/
  domeg_m=pi*(r2/2)*(r2/2)/pow(r,2);                                      /*max. solid angle */
  printf("\n");
  printf(" With these dimensions ,we have: \n");
  printf("Max. solid angle at 90 is= "); printf("%-7.4f",domeg_m); printf(" sr \n");  

                              /*beam parameters  */

  printf(" Enter value of gamma; if <1 it will be set equal to 55773 \n");   /*default is SLAC */
  scanf("%lf",&gam);
  if (gam<=1) gam=55773;
  beta=sqrt(pow(gam,2)-1)/gam;
  printf("Type normalized emittance of the beam in mm.mrad \n");
  scanf("%le",&nemit);
  emit=1e-3*nemit/(beta*gam);                                  /*converted to mm.rad*/
  esc20:printf("Enter measured (or assumed) position (in mm) of beam centroid above the PEAK of the grating teeth\n");
  printf("NB. This is NOT the quantity x0, which is from the base of the teeth and which will be calculated by the code\n");
  scanf("%lf",&x_e);
  if (x_e<=0) x_e= 0.0;
  x0 = x_e + el* sin(alpha1)* cos(alpha1);
  printf("x0 for this grating is = %-5.3fmm\n", x0);
  printf("\n");
                              
  esc24:printf("Enter FWHM_X(in mm) of Gaussian beam at the centre of grating \n");
  scanf("%le",&fwhm_x);
  sigma_x= fwhm_x/(2*sqrt(2*log(2)));
  printf(" For this fwhm: \n");
  printf(" 0.50 of the beam is within a radius of "); printf("%-4.2e",0.6749*sigma_x); printf(" mm \n");
  printf(" 0.90 of the beam is within a radius of "); printf("%-4.2e",1.6450*sigma_x); printf(" mm \n");
  printf(" 0.99 of the beam is within a radius of "); printf("%-4.2e",2.5733*sigma_x); printf(" mm \n");
  y0=1.6450*sigma_x;                                  /*average value of radius over 1/2 grating*/
  charlen=pow(y0,2)/emit;
  sf3=0;   
  j=1;
  for (j=1; j<7; j++)                           /*first calculate beam profile over grating*/
    {z=(j-1)*lgr/10;
    x_rad=y0*sqrt(1+pow(z/charlen,2));                           /*varying beam half-width*/
    sf3=sf3+x_rad;}
  y0_av=sf3/6;                                  /*average value of radius over 1/2 grating*/
  sigma_x= y0_av/1.645;                         /*average value of sigma_x over 1/2 grating*/  
  printf("Average value of FWHM_X from centre to edge of grating is "); printf("%-4.2e",(2*sigma_x*sqrt(2*log(2)))); printf(" mm\n"); 
  printf("Enter FWHM_Y(in mm) of Gaussian beam at waist, at the centre of grating \n");
  scanf("%le",&fwhm_y);
  sigma_y= fwhm_y/(2*sqrt(2*log(2)));
  printf(" For this fwhm: \n");
  printf(" 0.50 of the beam is within a radius of "); printf("%-4.2e",0.6749*sigma_y); printf(" mm \n");
  printf(" 0.90 of the beam is within a radius of "); printf("%-4.2e",1.6450*sigma_y); printf(" mm \n");
  printf(" 0.99 of the beam is within a radius of "); printf("%-4.2e",2.5733*sigma_y); printf(" mm \n");

	esc1: printf( "\n");	
	res1 = 0.0;
  lol=x0-5*sigma_x;                                       /* GD modified 5/6/06*/
  if (lol<hi) lol=hi;
  upl=x0+5*sigma_x;
  epsabs=1.0e-12;
  epsrel=1.0e-6;
  liw=500;
 	d01sjc(GAU,lol,upl,epsabs,epsrel,liw,&res1,&abserr,&QP, &COMM1, &fail);
  res1=res1/(sigma_x*sqrt(2*pi));
  printf("For a fwhm of "); printf("%-4.2e",fwhm_x);
  printf("mm and a beam position of="); printf("%-4.3f",x_e); printf("mm above the peak of the teeth, \n");
  printf(" %-4.2f percent of the beam is above the grating \n\n",100*res1);
  NAG_FREE(QP.sub_int_beg_pts);
  NAG_FREE(QP.sub_int_end_pts);
  NAG_FREE(QP.sub_int_result);
  NAG_FREE(QP.sub_int_error);
    
                                            /* bunch parameters */
  
  printf("Enter total number of electrons in the bunch \n");
  scanf("%le",&n_e_t);
  printf("since only "); printf("%-4.1f",100*res1); printf(" percent is above the grating, I am setting N_e=");
  printf("%le",res1*n_e_t); printf("\n");
  n_e=n_e_t*res1;
  esc25: printf("Enter nominal bunch duration (FWHM), in ps \n");
  scanf("%lf",&pulse_len);
  i_av_p=n_e_t*1.602177*1e-7/pulse_len;
  printf("The average current in the pulse is= "); printf("%-4.3f",i_av_p); printf(" A \n");  
  
  printf("You now have the following  choice for profiles in z-direction \n");
  printf("Type one of the following: \n");
  printf(" 0 for a DC beam \n");
  printf(" 1 for a single Gaussian \n"); printf(" 2 for triangular \n"); printf(" 3 for double-exp. \n");
  printf(" 4 for parabolic \n"); printf(" 5 for cosine \n"); printf(" 6 NOT IN USE \n");
  printf(" 7 for Lorentz shape \n"); printf(" 8 for a multi-Gauss shape \n");
  printf(" which one do you want?? \n");
  scanf("%d",&pulse_shp);
  if (pulse_shp<0 || pulse_shp>8)
  { printf(" I do not recognize this \n");
   goto esc1;}
  if (pulse_shp==1)                     /*Gaussian*/
  {  
   printf("Type asymmetry factor epsl \n");
   scanf("%lf",&epsl);
   sigma_t= pulse_len/(sqrt(2*log(2))*(1+epsl));              /*length defined by FWHM 15/5/09*/
   printf(" the sigma of the leading edge is "); printf("%-4.3f",sigma_t); printf(" ps \n");
   printf(" and sigma of the trailing edge is "); 
   printf("%-4.3f",epsl*sigma_t); printf(" ps \n"); 
  }
  else if (pulse_shp==2)                      /*triangular, modified 11/05/09*/
  {
   printf("Type asymmetry factor epsl \n");
   scanf("%lf",&epsl);
   sigma_t= pulse_len/((1+epsl)*(1-0.5));                          /* now in terms of FWHM, 15/05/09 */
   printf(" For this bunch shape, the sigma_t is "); printf("%-4.3f",sigma_t); printf(" ps \n");
  }
  else if (pulse_shp==3)                     /*double exponential*/
  {  
   printf("Type asymmetry factor epsl \n");
   scanf("%lf",&epsl);
   sigma_t=pulse_len/((1+epsl)*log(1/0.5));             /*changed 15/05/09, now using FWHM */
   printf(" the sigma_t is "); printf("%-4.3f",sigma_t); printf(" ps \n");    
  }
  else if (pulse_shp==4)                        /*parabolic*/
  {
   printf("Type asymmetry factor epsl \n");
   scanf("%lf",&epsl);
   sigma_t=pulse_len/((1+epsl)*sqrt(1-0.5));       /*changed 15/05/09, using FWHM*/   
   printf(" For this bunch shape, sigma_t is "); printf("%-4.3f",sigma_t); printf(" ps \n");    
  }
  else if (pulse_shp==5)                     /*cosine*/
  {
   epsl=1;
   sigma_t=pulse_len*pi/(4*acos(0.5));                   /* changed 15/05/09, using FWHM*/   
   printf(" For this bunch shape, sigma_t is "); printf("%-4.3f",sigma_t); printf(" ps \n");
  }
  else if (pulse_shp==6)
  {
   goto esc1;
  }
  else if (pulse_shp==7)                    /*Lorentz*/
  {
   printf("Enter the asymmetry factor epsl; symmetric shape means eps=1 \n");
   scanf("%lf",&epsl);
   apb=pulse_len;                                 /* changed 15/5/09*/
   printf("a+b is = "); printf("%-4.3f",apb); printf(" ps \n");  
   a=apb/(1+epsl);       /* no need for sigma_t; shape is defined by a & b */
   b=a*epsl;
   printf(" a(leading edge)= "); printf("%-4.3f",a); printf("ps \n"); 
   printf(" b(trailing edge)= "); printf("%-4.3f",b); printf("ps \n");
  }
  else if (pulse_shp==8)        /*multiple Gaussian,  June 08, modified for FWHM 11/05/09*/
  {   
   esc81:fwhm_t1= fwhm_t2= fwhm_t3= fwhm_t4= fwhm_t5=0.0;
   sigma_t1= sigma_t2= sigma_t3= sigma_t4= sigma_t5=0.0;
   g1=g2=g3=g4=g5=0.0;
   mu1=mu2=mu3=mu4=mu5=0.0;
   gtop=gbot=len50=0.0;
   printf(" how many gaussians are you fitting? Min. is 2 and max. is 5 \n");
   scanf("%1d", &ngauss);
   res1=0; 
   printf(" for the 1st peak, enter fwhm(ps), amplitude & position(ps) \n");
   scanf("%lf %lf %lf", &fwhm_t1, &g1, &mu1);
   sigma_t1= fwhm_t1/(2*sqrt(2*log(2)));
   if (ngauss==2)
   {printf(" for the 2nd peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t2, &g2, &mu2);
    sigma_t2= fwhm_t2/(2*sqrt(2*log(2)));
    sigma_t3=sigma_t4=sigma_t5=0;
    g3=g4=g5=0;
    mu3=mu4=mu5=0;}
   else if (ngauss==3)
   {printf(" for the 2nd peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t2, &g2, &mu2);
    sigma_t2= fwhm_t2/(2*sqrt(2*log(2)));
    printf(" for the 3rd peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t3, &g3, &mu3);
    sigma_t3= fwhm_t3/(2*sqrt(2*log(2)));
    sigma_t4=sigma_t5=0;
    g4=g5=0;
    mu4=mu5=0;}
   else if (ngauss==4)
   {printf(" for the 2nd peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t2, &g2, &mu2);
    sigma_t2= fwhm_t2/(2*sqrt(2*log(2)));
    printf(" for the 3rd peak, enter fwhm(ps), amplitude & position(ps) \n");    
    scanf("%lf %lf %lf", &fwhm_t3, &g3, &mu3);
    sigma_t3= fwhm_t3/(2*sqrt(2*log(2)));
    printf(" for the 4th peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t4, &g4, &mu4);
    sigma_t4= fwhm_t4/(2*sqrt(2*log(2)));
    sigma_t5=0;
    g5=0;
    mu5=0;}
   else
   {printf(" for the 2nd peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t2, &g2, &mu2);
    sigma_t2= fwhm_t2/(2*sqrt(2*log(2)));
    printf(" for the 3rd peak, enter fwhm(ps), amplitude(ps) & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t3, &g3, &mu3);
    sigma_t3= fwhm_t3/(2*sqrt(2*log(2)));
    printf(" for the 4th peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t4, &g4, &mu4);
    sigma_t4= fwhm_t4/(2*sqrt(2*log(2)));
    printf(" for the 5th peak, enter fwhm(ps), amplitude & position(ps) \n");
    scanf("%lf %lf %lf", &fwhm_t5, &g5, &mu5);
    sigma_t5= fwhm_t5/(2*sqrt(2*log(2)));}

  /* check total value of this integral, i.e. normalisation factor */
   tot_int= sqrt(2*pi)*(g1*sigma_t1+g2*sigma_t2+g3*sigma_t3+g4*sigma_t4+g5*sigma_t5);
   printf(" the normalisation factor is  "); 
   printf("%-4.3e",tot_int); printf(" units \n");  
   sigma_teff= g1*sigma_t1+ g2*sigma_t2+g3*sigma_t3+ g4*sigma_t4+g5*sigma_t5; 
     /* find position & value of maximum of normalised distribution; total scan width=15ps */
   fmax=0;
   tmax=0;
   i=0;
   for (i==0; i<1500; i++)
   { 
   t= 1*sigma_teff- i*0.01;                  /* start well into the positive region*/  
   f=(g1*exp(-(t-mu1)*(t-mu1)/(2*sigma_t1*sigma_t1))+g2*exp(-(t-mu2)*(t-mu2)/(2*sigma_t2*sigma_t2))+
     g3*exp(-(t-mu3)*(t-mu3)/(2*sigma_t3*sigma_t3))+
     g4*exp(-(t-mu4)*(t-mu4)/(2*sigma_t4*sigma_t4))+
     g5*exp(-(t-mu5)*(t-mu5)/(2*sigma_t5*sigma_t5)))/tot_int;
   if (f>= fmax) {fmax=f; tmax=t;}
   }
   printf(" maximum for the (normalised) shape is at t="); 
   printf("%-4.3f", tmax); printf(" \n");
   printf(" and is equal to "); 
   printf("%-4.3f", fmax); printf(" \n");
   /* finished with this; now find points where value of function is about 0.5 */
   /* and use these as the limits for the evaluation of the normalised integral */
   /* start on positive side...*/
   i=0;
   for (i==0; i<1500; i++)
   { 
   t= tmax + i*0.01;   
   f=(g1*exp(-(t-mu1)*(t-mu1)/(2*sigma_t1*sigma_t1))+g2*exp(-(t-mu2)*(t-mu2)/(2*sigma_t2*sigma_t2))+
     g3*exp(-(t-mu3)*(t-mu3)/(2*sigma_t3*sigma_t3))+
     g4*exp(-(t-mu4)*(t-mu4)/(2*sigma_t4*sigma_t4))+
     g5*exp(-(t-mu5)*(t-mu5)/(2*sigma_t5*sigma_t5)))/tot_int;
   if (f<= 0.5*fmax) {gtop=t; goto esc82;}
   }
  esc82: printf("the upper 50 percent point of the distribution is "); printf("%-4.3f",gtop);
/* continue on negative...*/
   i=0;
   for (i==0; i<1500; i++)
   { 
   t= tmax - i*0.01;   
   f=(g1*exp(-(t-mu1)*(t-mu1)/(2*sigma_t1*sigma_t1))+g2*exp(-(t-mu2)*(t-mu2)/(2*sigma_t2*sigma_t2))+
     g3*exp(-(t-mu3)*(t-mu3)/(2*sigma_t3*sigma_t3))+
     g4*exp(-(t-mu4)*(t-mu4)/(2*sigma_t4*sigma_t4))+
     g5*exp(-(t-mu5)*(t-mu5)/(2*sigma_t5*sigma_t5)))/tot_int;
   if (f<= 0.5*fmax) {gbot=t; goto esc83;}
   }
   esc83: printf(" and the lower one "); printf("%-4.3f",gbot); printf("\n");
/* done*/        
   esc80: len50= gtop- gbot;    /* this is the assumed bunch length */
   printf(" Hence, the assumed bunch length (FWHM) is "); printf("%-4.3f",len50); printf(" ps \n");
   flag=0;
   if (!flag) printf(" Do you want to change anything? yes(y) or no(n) \n");
   line[0]= 'x';
   do { 
       if (line[0] != 'x') printf(" ");
       if (fgets(line, 200, stdin) == NULL) { printf("\nErr: fgets.\n"); exit(1); }
      }
   while (line[0] != 'y' && line[0] != 'n'); 
   ans = line[0];
   if (ans=='y')
   {printf("go back to input parameters \n");
   goto esc81;
   }
   printf(" \n");  
  }

                    /*estimate of grating factor R^2  and of the coherent & incoherent integrals*/
 printf("Alpha1 = %f\n", alpha1);
 printf("Alpha2 = %f\n", alpha2);
 printf("This code can consider up to 7 steps in the azimuthal direction phi \n");
 printf(" Define step size: min.=0.1deg, max.=10deg, default=1.0deg. \n");
 scanf("%lf", &phistep);
 if (phistep<0.1 || phistep>10) phistep=1.0;
 printf("OK, phi step is equal to= "); printf("%-3.1f", phistep); printf(" deg.\n");
    
/* set the y-limits for the y-z integrations*/
  ll[0]= 0.0; ul[0]= wid/2;                  /*y-integral; wid is the grating width */
  maxpts= 5000;
  eps=1e-4;                                /* acceptable error*/
  
  i=(int)(alpha1*180/pi);

	
  for ( i=(int)(alpha1*180/pi); i<171; i++)       /* angles lower than the blaze angle do not matter! */
  { 
    thet[i]=i;                                 /*scan angles from 10-170 deg.  */
    theta=thet[i]*pi/180;
    lam=el*(1/beta-cos(theta))/m;      /*lam is the wavelength of the emitted radiation*/
    k= 2*pi/lam;
    if (i==10||i==20||i==30|| i==40||i==50||i==60 || i==70 || i==80||i==90 || i==100 || i==110 || i==120
    ||i==130 || i==140 || i==150 || i==160 || i==170)    
    {printf("%-.1f  ", (theta*180/pi)); printf("\n");}
    ii=0;
    for (ii=0; ii<8; ii++)
   {
     phi0=phistep*ii*pi/180;                    /* azim. angle phi scanned up to 7xphistep deg.*/       
     kx=k*sin(theta)*cos(phi0);
     ky=k*sin(theta)*sin(phi0);
     lam_e=lam*beta*gam/(2*pi*sqrt(1+pow((beta*gam*sin(theta)*sin(phi0)),2)));  /*evanescent w/l */
     ff= 1/(pi*el);
/*  N.B. all distances from the base of the grating teeth.   */

                       /* calculate for a range of x0 values and store for integration */ 
     lol= x0-5*sigma_x;
     if (lol<=hi) lol=1.001*hi;      /* electrons must not hit the grating!!*/
     upl= x0+5*sigma_x;
     indend=20;                     /* number of steps in x, i.e. in 0.5sigma steps */
     ind=0;
	
     for (ind=0; ind<indend+1; ind++)
     { 	
        xv= lol+ ind*0.5*sigma_x; /*proceed in steps of 0.5sigma   first facet  */
        res1=res2=res3=res4=0.0;
        ll[1]=0.0; ul[1]=el1;  
          /* z limits for 1st facet */
        alpha=alpha1;
        psi=0;
        dum=0.0;      
                /* dum accounts for diff. expressions for impact param. for the two facets */
        D= k*(1/beta-cos(theta))-kx*tan(alpha);
	temp1=temp2=temp3=temp4=0.0;
        minpts=1;
        d01wcc(ndim, QU1, ll, ul, &minpts, maxpts, eps, &res1, &acc, &COMMWC1, &fail);  /* real part of the int. for the x&z components */ 
        minpts=1;
        	d01wcc(ndim, QU2, ll, ul, &minpts, maxpts, eps, &res2, &acc, &COMMWC2, &fail);   /* imag. part */
       G1x=Ccreate((ff*res1*tan(alpha1)), (ff*res2*tan(alpha1)));      /*changed 25/2/10*/
       G1z=Ccreate((ff*res1), (ff*res2));                              /*changed 25/2/10*/
       temp3=temp4=0.0;
       minpts=1;
       d01wcc(ndim, QU3, ll, ul, &minpts, maxpts, eps, &res3, &acc, &COMMWC3, &fail);  /*real part for the integr. of y component */
       minpts=1;
       d01wcc(ndim, QU4, ll, ul, &minpts, maxpts, eps, &res4, &acc, &COMMWC4, &fail);  /*imag. part of y-integral */
       G1y=Ccreate((ff*res3*tan(alpha1)), (-ff*res4*tan(alpha1)));             /*changed 25/2/10*/
	
/*                    second facet,including phase correction psi             */
			/*NO SECOND FACET IS NEEDED FOR ECHELLE GRATING*/

/* and create the complex vector G    */
      G= Vcreate( G1x.r, G1x.i, G1y.r, G1y.i, G1z.r, G1z.i);
      nv= Vcreate( (sin(theta)*cos(phi0)), 0, (sin(theta)*sin(phi0)), 0, cos(theta), 0);  /*dir. unit vector*/ 
      nvp= Vcreate( (cos(theta)*cos(phi0)), 0, (cos(theta)*sin(phi0)), 0, -sin(theta), 0); /*1st pol. unit vector*/
      nvv= Vcreate( -sin(phi0), 0, cos(phi0), 0, 0, 0);                                    /*2nd pol. unit vector*/ 
   ;
/*    first polarization e1, i.e vector nvp, in plane of n,z. Need G.nvp                                   */    
      a1=Vdot(G,nvp);
      a1t= Cmultiply(a1, Cconjugate(a1));        /*square of magn. in first polarisation*/ 
/*    second polarization e2, i.e. vector nvv, vertical to n and z                                         */
      a2=Vdot(G,nvv);
      a2t= Cmultiply(a2, Cconjugate(a2));        /*square of magn. in second polarisation*/ 
      f5c=Cadd(a1t,a2t);                /*total */
      f5= f5c.r;                       /* output is R^2, in both polarisations, at this specific value of x*/
      f5_p1= a1t.r;                           /* polar. p1  */
      f5_p2= a2t.r;                           /* .. and p2... */

      if (ind==10)
      {single_elec[i][ii]= f5;
       single_elec_p1[i][ii]= f5_p1;
       single_elec_p2[i][ii]= f5_p2;}            /*store the single-electron R^2 */  
      x0_bin[ind]= xv;
      inc_bin[ind]= f5*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
      inc_bin_p1[ind]= f5_p1*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
      inc_bin_p2[ind]= f5_p2*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
      coh_bin[ind]= sqrt(f5)*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
      coh_bin_p1[ind]= sqrt(f5_p1)*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
      coh_bin_p2[ind]= sqrt(f5_p2)*exp(-(xv-x0)*(xv-x0)/(2*sigma_x*sigma_x));
	
                                  
     }                              /*end of x0 scan  */ 

/* numerical integration over x0 follows, to find the coherent and incoherent integrals*/
/* the two arrays below contain the incoherent and coherent parts of the x-integrals */
/* arranged so that in col 0 data for phi=0, in 1 for phi=1(10), in 7 for phi=7(70)*/
    d01gac(indend, x0_bin, inc_bin, &incans, &err, &fail);
    s_inc_x[i][ii]= incans/(sqrt(2*pi)*sigma_x);    
    d01gac(indend, x0_bin, coh_bin, &cohans, &err, &fail);
    s_coh_x[i][ii]= cohans*cohans/(2*pi*sigma_x*sigma_x);
    
    d01gac(indend, x0_bin, inc_bin_p1, &incans, &err, &fail);
    s_inc_x_p1[i][ii]= incans/(sqrt(2*pi)*sigma_x);    
    d01gac(indend, x0_bin, coh_bin_p1, &cohans, &err, &fail);
    s_coh_x_p1[i][ii]= cohans*cohans/(2*pi*sigma_x*sigma_x);
    
    d01gac(indend, x0_bin, inc_bin_p2, &incans, &err, &fail);
    s_inc_x_p2[i][ii]= incans/(sqrt(2*pi)*sigma_x);    
    d01gac(indend, x0_bin, coh_bin_p2, &cohans, &err, &fail);
    s_coh_x_p2[i][ii]= cohans*cohans/(2*pi*sigma_x*sigma_x);
    
    }                             /* end of ii (phi) scan */
    
  }                               /* end of i (theta) scan */  
 
  
  printf(" The following is a list of  the values of R^2, total, P1, P2 averaged over x-dimension of bunch\n");
  printf(" NB!! The dimensionless quantity R^2 must be multiplied by the grating period in order to obtain\n");
  printf(" the true measure of the efficiencies of two different gratings.\n");
  printf(" See equation 8a of the Theory note for justification\n");
  printf(" Each column below is for a phi step of "); printf("%-3.1f",phistep); printf(" deg. \n");
  printf(" first column is the wavelength, in microns \n");
  printf(" last column is the average over the 3 values of phi on either side of zero, i.e.\n");
  printf(" between -%-3.1f & +%-3.1f deg. only, in this case \n",phistep*3, phistep*3);
  i=(int)(alpha1*180/pi);
  for ( i=(int)(alpha1*180/pi); i<171; i++) 
  {   printf("\n");
      lam= 1e3*el*(1-cos(i*pi/180))/m;
      printf("%6.1f  ", lam);          /* in microns */
      ii = 1;
      sum = sum_p1 = sum_p2 =0;
      for (ii = 1; ii<4; ii++)
      { sum = sum + s_coh_x[i][ii];
        sum_p1 = sum_p1 + s_coh_x_p1[i][ii];
	sum_p2 = sum_p2 + s_coh_x_p2[i][ii];}
        s_coh_x[i][8] = (s_coh_x[i][0] + 2*sum) /7;               /* average over phi, from -3 to +3 deg.*/
        s_coh_x_p1[i][8] = (s_coh_x_p1[i][0] + 2*sum_p1)/7;              /* average over phi, from -3 to +3 deg.*/
        s_coh_x_p2[i][8] = (s_coh_x_p2[i][0] + 2*sum_p2)/7;              /* average over phi, from -3 to +3 deg.*/
      for (ii = 0; ii<9; ii++)
      {  printf( "%-4.3e  ", s_coh_x[i][ii]);}
      printf(" \n"); printf("        ");
      for (ii = 0; ii<9; ii++)
      {  printf( "%-4.3e  ", s_coh_x_p1[i][ii]);}
      printf("\n"); printf("        ");
      for (ii = 0; ii<9; ii++)
      {  printf( "%-4.3e  ", s_coh_x_p2[i][ii]);}       
  }    
  printf("\n");
  printf(" Each column above is for a phi step of "); printf("%-3.1f",phistep); printf(" deg. \n");
  printf(" last column is the average over the 3 values of phi on either side of zero, i.e.\n");
  printf(" between -%-3.1f & +%-3.1f deg. only, in this case \n",phistep*3, phistep*3);
  printf(" Finished with the calculation of R^2; see above for interpretation \n");
  printf("\n"); 
	         /* calculation of coherence factors & of radiated energy */

  printf("\nCalculating energy per bunch, for this specific grating (4cm long)  \n");
  f1=2*pi*qe*qe*1e-7;                                           /* in Joules*/
  printf("The  scale factor f1="); printf("%le",f1); printf(" J.cm \n");
  f2=100/pow(el,2);                 /* geometry factor f2 in cm-2; grating period in mm */
  
  i=(int)(alpha1*180/pi);



  for ( i=(int)(alpha1*180/pi); i<171; i++)
  { thet[i]=i;                                 /*scan angles from 10-170 deg.  */
    theta=thet[i]*pi/180;
    lam=el*(1/beta-cos(theta))/m;
    k= 2*pi/lam;
    ii=0;
    for (ii=0; ii<8; ii++)
    {phi0=phistep*ii*pi/180;                   /* phi scanned in 7 steps, each of size=phistep*/       
     kx=k*sin(theta)*cos(phi0);
     ky=k*sin(theta)*sin(phi0);
     lam_e=lam*beta*gam/(2*pi*sqrt(1+pow((beta*gam*sin(theta)*sin(phi0)),2)));      
/*     find the angular distribution factor f4          */     
     f4[i]=pow(m,2)*pow((beta/(1-beta*cos(theta))),3);
     om=2*pi*0.3/lam;                          /*c=0.3 mm/ps, hence om is in ps^-1*/                           
                      /*calculate longitudinal coherence effects*/
     if (pulse_shp==0)
     {
       u=0;                         /*not necessary, but no coh. for DC beam*/
       f_z_c=0;
       f_z_s=0;
     }
     else if (pulse_shp==1)
     {u1= om*sigma_t;
      u2= om*sigma_t*epsl;
      temp1=0; temp2=0;
      f_z_c=(exp(-(u1*u1)/2)+ epsl*exp(-(u2*u2)/2))/(1+epsl);    /* yes, it is correct!*/
/* degenerate hypergeometric function below is taken from G&R page 480  */ 
      jj=1;
      hyperg[0]=1;
      for (jj==1; jj<201; jj++)
      {ratio= (u1*u1/2)*(2*jj-1)/(jj*(2*jj+1));
       hyperg[jj]= hyperg[jj-1]*ratio;
      }
      jj=1;
      sum=1;
      for (jj==1; jj<201; jj++)
      { dumsum= sum+ hyperg[jj];
      	sum=dumsum;
      }
      temp1= fabs(sum*u1*sigma_t*exp(-u1*u1/2));      
/*      now for the trailing part of the bunch */ 
      jj=1;
      hyperg[0]=1;
      for (jj==1; jj<201; jj++)
      { ratio= (u2*u2/2)*(2*jj-1)/(jj*(2*jj+1));
	 hyperg[jj]= hyperg[jj-1]*ratio;
      }
      jj=1;
      sum=1;
      for (jj==1; jj<201; jj++)
      { dumsum= sum+ hyperg[jj];
      	sum=dumsum;
      }
      temp2= fabs(sum*u2*sigma_t*epsl*exp(-u2*u2/2));
      f_z_s=2*(temp1-temp2)/(sqrt(2*pi)*sigma_t*(1+epsl));    /* including nor. factor of 2/sqrt(2pi)..  */	
      }        
      else if (pulse_shp==2)
      {
       u=om*sigma_t;
       f_z_c=2*(1-(epsl*cos(u)+cos(epsl*u))/(1+epsl))/(epsl*u*u);
       f_z_s=2*((sin(epsl*u)-epsl*sin(u))/(1+epsl))/(epsl*u*u);
      }
/* integrate to edge of pulse*/
      else if (pulse_shp==3)
      {
       u=om*sigma_t;
       f_z_c=(1+epsl*u*u)/(1+u*u)/(1+pow((epsl*u),2));
       f_z_s=(epsl*u-u)/(1+u*u)/(1+pow((epsl*u),2));
      }
/*integr. to infinity, so only sigma_t appears in result*/
      else if (pulse_shp==4)
      {
       u=om*sigma_t;
       f_z_c=3*(sin(u)+sin(epsl*u)/(epsl)-u*cos(u)-u*cos(epsl*u)/(epsl))/((1+epsl)*pow(u,3));
       f_z_s=3*(1/(epsl*epsl)-1+cos(u)-cos(epsl*u)/(epsl*epsl)+u*sin(u)-u*sin(epsl*u)/(epsl))/((1+epsl)*pow(u,3));
/* integrate to edge of pulse i.e. sigma_t*/
      }
      else if (pulse_shp==5)
      {
       u=om*sigma_t;
       f_z_c=pi*pi*cos(u)/(pi*pi-4*u*u);
       f_z_s=0;
/*must integrate to edge of pulse, i.e. sigma_t*/
      }
      else if (pulse_shp==6)
      {
       u=om*sigma_t;
       f_z_c=((1-u*u)*cos(u)+2*u*sin(u))/(pow((1+u*u),2));
       f_z_s=0;
/*must integrate from -inf. to edge of pulse, i.e. sigma_t*/
      }
      else if (pulse_shp==7)                      /*numerical integr. for Lorentz shape*/
      {
        /* a and b deal with asymmetric Lorentz shape. a is for leading edge*/
                            /* norm. factor is 2/pi/(a+b)  */
       f_z_c=(b*exp(-om*b)+a*exp(-om*a))/(apb);
       WT_FUNC=Nag_Sine;
       d01ssc(LOR1,0.0,om,WT_FUNC,50,500,1e-3,&res1,&abserr,&QPSUB1,&COMM2,&fail);
       d01ssc(LOR2,0.0,om,WT_FUNC,50,500,1e-3,&res2,&abserr,&QPSUB2,&COMM3,&fail);
       f_z_s=2*(b*b*res2-a*a*res1)/(pi*(apb));
       NAG_FREE(QPSUB1.interval_error);
       NAG_FREE(QPSUB1.interval_result);
       NAG_FREE(QPSUB1.interval_flag);
       NAG_FREE(QPSUB2.interval_error);
       NAG_FREE(QPSUB2.interval_result);
       NAG_FREE(QPSUB2.interval_flag);
      }
      else if (pulse_shp==8)                        /* multiple Gaussian, normalised, modified 15/3/08*/
      {sigma_teff = g1*sigma_t1+ g2*sigma_t2+g3*sigma_t3+ g4*sigma_t4+g5*sigma_t5;
       u1=om*sigma_t1; u2=om*sigma_t2;
       u3=om*sigma_t3; u4=om*sigma_t4; u5=om*sigma_t5;
       f_z_c=sqrt(2*pi)*(g1*sigma_t1*exp(-u1*u1/2)*cos(mu1*om)+
                         g2*sigma_t2*exp(-u2*u2/2)*cos(mu2*om)+
                         g3*sigma_t3*exp(-u3*u3/2)*cos(mu3*om)+
			 g4*sigma_t4*exp(-u4*u4/2)*cos(mu4*om)+
			 g5*sigma_t5*exp(-u5*u5/2)*cos(mu5*om))/tot_int;
       f_z_s=sqrt(2*pi)*(g1*sigma_t1*exp(-u1*u1/2)*sin(mu1*om)+
                         g2*sigma_t2*exp(-u2*u2/2)*sin(mu2*om)+
                         g3*sigma_t3*exp(-u3*u3/2)*sin(mu3*om)+
			 g4*sigma_t4*exp(-u4*u4/2)*sin(mu4*om)+
			 g5*sigma_t5*exp(-u5*u5/2)*sin(mu5*om))/tot_int;
      }		

      s_coh_long=pow(f_z_c,2)+pow(f_z_s,2);         /*finished with long. coherence factor*/
      s_coh_z[i][ii]= s_coh_long;                  /*store in separate array for future use*/
      f_y_c=exp(-pow((ky*sigma_y),2)/2);     /* NB integration here is to infinity; does it matter?? */
      f_y_s=0.0;
      s_coh_xy= (f_y_c*f_y_c+f_y_s*f_y_s)*s_coh_x[i][ii];
      s_coh_xy_p1= (f_y_c*f_y_c+f_y_s*f_y_s)*s_coh_x_p1[i][ii];
      s_coh_xy_p2= (f_y_c*f_y_c+f_y_s*f_y_s)*s_coh_x_p2[i][ii];
      s_coh_tr[i][ii]= s_coh_xy;             /*transverse coherence factor, in new array*/
      s_coh_tr_p1[i][ii]= s_coh_xy_p1;             /*transverse coherence factor, pol. 1,  in new array*/
      s_coh_tr_p2[i][ii]= s_coh_xy_p2;             /*transverse coherence factor, pol.2,  in new array*/
      f_b[i][ii]= n_e*s_inc_x[i][ii]+ n_e*n_e*s_coh_xy*s_coh_long;     /*overall effect of Ne electrons in a bunch*/
      f_b_p1[i][ii]= n_e*s_inc_x_p1[i][ii]+ n_e*n_e*s_coh_xy_p1*s_coh_long;     /*...in pol.1..*/
      f_b_p2[i][ii]= n_e*s_inc_x_p2[i][ii]+ n_e*n_e*s_coh_xy_p2*s_coh_long;     /*... and pol.2..*/
      bin[i][0]=i;
      bin[i][ii+1]= f1*f2*f4[i]*(n_e*s_inc_x[i][ii]+ n_e*n_e*s_coh_xy*s_coh_long);
      bin_p1[i][ii+1]= f1*f2*f4[i]*(n_e*s_inc_x_p1[i][ii]+ n_e*n_e*s_coh_xy_p1*s_coh_long);
      bin_p2[i][ii+1]= f1*f2*f4[i]*(n_e*s_inc_x_p2[i][ii]+ n_e*n_e*s_coh_xy_p2*s_coh_long);
      bin_single[i][0]= i;                      /* store for single electron yield */
      bin_single[i][ii+1]= f1*f2*f4[i]*single_elec[i][ii];
      bin_single_p1[i][ii+1]= f1*f2*f4[i]*single_elec_p1[i][ii];
      bin_single_p2[i][ii+1]= f1*f2*f4[i]*single_elec_p2[i][ii];
      /*NB!! in col.1 data for phi=0, in 2 for phi=1, in 8 for phi=7*/     
    }                                                   /* end of ii (phi) scan */                          
  }                                                    /* end of i (theta) scan */
  printf("\n"); printf("\n");
   
/*  printf("For Phi0=2 ONLY \n");
  printf(" Angle   Lambda   S_inc       S_coh_z     S_coh_tr   Bunch_factor   Energy \n");
  printf("          (mm)                                                      (J/sr) \n");
  i=(int)(alpha1*180/pi);
  for ( i=(int)(alpha1*180/pi); i<171; i++)
  { ii=2;
    thet[i]=i;                                 /*scan angles up to 170 deg.  */
/*    theta=thet[i]*pi/180;
    lam=el*(1/beta-cos(theta))/m;
    printf("%5.1f",thet[i]); printf("    ");
    printf("%4.3f",lam);
    printf("    "); printf("%-4.3e",s_inc_x[i][ii]); printf("   ");
    printf("%-4.3e",s_coh_z[i][ii]); printf("   ");   printf("%-4.3e",s_coh_tr[i][ii]); printf("   ");
    printf("%-4.3e",f_b[i][ii]); printf("   "); printf("%-4.3e",bin[i][ii+1]);
    printf("\n");
    printf("                  "); printf("%-4.3e",s_inc_x_p1[i][ii]); printf("   ");
    printf("%-4.3e",s_coh_z[i][ii]); printf("   ");   printf("%-4.3e",s_coh_tr_p1[i][ii]); printf("   ");
    printf("%-4.3e",f_b_p1[i][ii]); printf("   "); printf("%-4.3e",bin_p1[i][ii+1]);
    printf("\n");
    printf("                  "); printf("%-4.3e",s_inc_x_p2[i][ii]); printf("   ");
    printf("%-4.3e",s_coh_z[i][ii]); printf("   ");   printf("%-4.3e",s_coh_tr_p2[i][ii]); printf("   ");
    printf("%-4.3e",f_b_p2[i][ii]); printf("   "); printf("%-4.3e",bin_p2[i][ii+1]);
    printf("\n");
  }  
  printf(" Angle   Lambda   S_inc       S_coh_z     S_coh_tr   Bunch_factor   Energy \n");
  printf("          (mm)                                                      (J/sr) \n");
  printf("The above apply for Phi2=0 ONLY \n");
  printf("\n");  
  printf("the following is a list of the INCOHERENT & COHERENT parts of the x-integral; INCOHERENT is first");
  printf("\n");  
  i=(int)(alpha1*180/pi);
  for ( i=(int)(alpha1*180/pi); i<171; i++)
  {   printf("\n");
      printf("%3d  ", i);
      ii=0;
      for (ii=0; ii<8; ii++)
      { printf( "%.3e  ", s_inc_x[i][ii]); }
      printf("\n"); printf("     ");
      ii=0;
      for (ii=0; ii<8; ii++)
      { printf( "%.3e  ", s_coh_x[i][ii]); }         
  }      
  printf("\n");*/
  

  printf("      List of single-electron yields     \n");
  i=(int)(alpha1*180/pi);
  for ( i=(int)(alpha1*180/pi); i<171; i++)
  {   printf("\n");
      printf("%3d  ", i);
      ii=1;
      for (ii=1; ii<9; ii++)
      { printf( "%.3e  ", bin_single[i][ii]); } 
               
  }      
  printf("\n");
  
  i=0;
  for (i=0; i<151; i++)
  {dum70[i][0]=i+20;
   dum70[i][1]=bin[i+20][1];                   /*copy to DUM70 but only from theta=20 */   
  }
  l=151;                                                  /*start of sorting section, N=no. of rows*/
  n=1;                                                    /* no. of columns, i.e. 2 here*/
  ic=1;                                                   /* column index for sorting, i.e. 2nd here*/
  i=0;
  for (i=0; i<l; i++)
  { j=i+1; 
    for (j=i+1; j<l+1; j++)
     { if (dum70[i][ic]<dum70[j][ic])
       { jj=0;
	 for (jj=0; jj<n+1; jj++)
             { dumm[jj]=dum70[i][jj];
               dum70[i][jj]=dum70[j][jj];
               dum70[j][jj]=dumm[jj];
             }
        }
     }
  }
  printf("\n");
  printf(" Maximum for this case =%-4.3e",dum70[0][1]);
  printf(" J/sr/cm at theta =%6.2f\n",dum70[0][0]);

/*rough estimate of full width at half max. for the power distribution*/
  i=0;
  for (i=0; i<151; i++)
  {
    if (dum70[i][1]>=0.5*dum70[1][1]) 
    {
        bin_1[i][0]=dum70[i][0];
        bin_1[i][1]=dum70[i][1];
    }
    else
    {
        bin_1[i][0]=0;
        bin_1[i][1]=0;
    }
  }
  l=151;                                                  /*start of sorting section, N=no. of rows*/
  n=1;                                                    /* no. of columns, i.e. 2 here*/
  ic=0;                                                   /* column index for sorting, i.e. 1st here*/
  i=0;
  for (i=0; i<l; i++)
  { j=i+1;
    for (j=i+1; j<l+1; j++)
     { if (bin_1[i][ic]<bin_1[j][ic])
       { jj=0;
	 for (jj=0; jj<n+1; jj++)
         { dumm[jj]=bin_1[i][jj];
           bin_1[i][jj]=bin_1[j][jj];
           bin_1[j][jj]=dumm[jj];
         }
       }
     }
  }
  theta_u=bin_1[0][0];
  i=0;
  for (i=0; i<151; i++)
  {
    if (bin_1[i][1]<=0) goto esc2;
    theta_d=bin_1[i][0];
  }
  esc2:printf("\n");
  printf("Rough estimate of full-width, half maximum \n");
  printf(" Lower limit at= %6.2f\n",theta_d);
  printf(" Upper limit at= %6.2f\n",theta_u);
  printf(" \n");
                     /* rough estimate of total radiated energy, per single electron & per bunch */
  printf("I will now average over the azimuthal angle phi \n");
  printf("The min. number of steps in phi, INCLUDING phi=0, is 1; max is 8,");
  printf(" which would take you up to phi= %4.1f \n", (7*phistep));
  printf(" Enter number of steps;  default is 4, which would correspond to phi=3 (for steps of 1 deg.) \n");
  scanf("%1i",&philmt);         /*!! NB. different definition of philmt in this version!!*/
  if (philmt<1 || philmt>8) philmt=4;  
  printf(" I will average up to phi= %-3.1f ",((philmt-1)*phistep)); printf(" \n");
  i=(int)(alpha1*180/pi);
  for ( i=(int)(alpha1*180/pi); i<171; i++)       /* angles lower than the blaze angle do not matter! */
  {sum= sum_p1= sum_p2= sum_single= sum_single_p1= sum_single_p2= 0.0;
   ii=2;
   if (philmt==1) 
   {sum= bin[i][1]; sum_p1=bin_p1[i][1]; sum_p2=bin_p2[i][1]; sum_single= bin_single[i][1];
   sum_single_p1= bin_single_p1[i][1]; 
   sum_single_p2= bin_single_p2[i][1];}
   else
   for (ii=2; ii<philmt+1; ii++)
   {sum=sum+ 2*bin[i][ii];             /* start from col. index 2(phi=1xphistep) and double the sum to account for +-phi */
    sum_p1=sum_p1+ 2*bin_p1[i][ii];             /* ditto for pol. 1 */
    sum_p2=sum_p2+ 2*bin_p2[i][ii];             /* and for pol.2  */
    sum_single= sum_single+ 2*bin_single[i][ii];
    sum_single_p1= sum_single_p1+ 2*bin_single_p1[i][ii];
    sum_single_p2= sum_single_p2+ 2*bin_single_p2[i][ii];}
   if (philmt==1) 
   {bin[i][9]= bin[i][1]; 
    bin_p1[i][9]= bin_p1[i][1];
    bin_p2[i][9]= bin_p2[i][1];
    bin_single[i][9]= bin_single[i][1];
    bin_single_p1[i][9]= bin_single_p1[i][1];
    bin_single_p2[i][9]= bin_single_p2[i][1];}
   else 
   {bin[i][9]=(bin[i][1]+ sum)/(2*philmt-1);       /* col. 9 is the store for AVERAGE dE over all ACCEPTED phi's; in J/sr/cm */	 
    bin_p1[i][9]=(bin_p1[i][1]+ sum_p1)/(2*philmt-1);       /* ditto for pol. 1 */
    bin_p2[i][9]=(bin_p2[i][1]+ sum_p2)/(2*philmt-1);       /* and for pol. 2 */
    bin_single[i][9]= (bin_single[i][1]+ sum_single)/(2*philmt-1);
    bin_single_p1[i][9]= (bin_single_p1[i][1]+ sum_single_p1)/(2*philmt-1);
    bin_single_p2[i][9]= (bin_single_p2[i][1]+ sum_single_p2)/(2*philmt-1);}
  }
  printf("\n");

                              /* integrate over phi up to the max. value=7 */
  d_theta=1.745e-2;                                      /* d_theta step is always 1 degree */
  d_phi= phistep*1.745e-2;                                 /* d_phi expressed in rad */
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)             /* integrate dE over ALL phi & store in 10*/
  { sum= sum_p1= sum_p2= sum_single= sum_single_p1= sum_single_p2= 0.0; 
    ii=1;
    for (ii=2; ii<9; ii++)
    {sum = sum + sin(i*pi/180)*d_theta*d_phi*bin[i][ii];     
     sum_p1 = sum_p1 + sin(i*pi/180)*d_theta*d_phi*bin_p1[i][ii];           /* pol. 1 only */
     sum_p2 = sum_p2 + sin(i*pi/180)*d_theta*d_phi*bin_p2[i][ii];           /* and pol. 2*/
     sum_single = sum_single + sin(i*pi/180)*d_theta*d_phi*bin_single[i][ii];
     sum_single_p1 = sum_single_p1 + sin(i*pi/180)*d_theta*d_phi*bin_single_p1[i][ii];
     sum_single_p2 = sum_single_p2 + sin(i*pi/180)*d_theta*d_phi*bin_single_p2[i][ii];     
    }        
    bin[i][10]= sin(i*pi/180)*d_theta*d_phi*bin[i][1]+ 2*sum;
      /* col. 10 is a temp. store for phi-integration over ALL phi's, i.e. up to phi=7; in J/sr/cm */
      /* only used for the estimate of the total radiated energy, see below */

    bin_p1[i][10] = sin(i*pi/180)*d_theta*d_phi*bin_p1[i][1]+ 2*sum_p1;  /* pol. 1 only */
    bin_p2[i][10] = sin(i*pi/180)*d_theta*d_phi*bin_p2[i][1]+ 2*sum_p2;  /*  and pol. 2 */
    bin_single[i][10] = sin(i*pi/180)*d_theta*d_phi*bin_single[i][1]+ 2*sum_single;    
    bin_single_p1[i][10] = sin(i*pi/180)*d_theta*d_phi*bin_single_p1[i][1] + 2*sum_single_p1;
    bin_single_p2[i][10] = sin(i*pi/180)*d_theta*d_phi*bin_single_p2[i][1] + 2*sum_single_p2;
  }
                               /* now integrate over theta */
  i=(int)(alpha1*180/pi);
  sum1 = sum2 = sum3 = sum4 = 0.0;
  for (i=(int)(alpha1*180/pi); i<171; i++)           /* now truly in J! The length is in the f2 factor! */
  { sum1= sum1 + bin[i][10];
    sum3= sum3 + bin_p1[i][10];
    sum4= sum4 + bin_p2[i][10];
    sum2= sum2 + bin_single[i][10];
  }
  printf("\n"); 

  printf("                             For this grating                \n");
  printf("             Rough estimate of total radiated energy (single electron) over hemisphere = ");
  printf("%.3e", sum2*0.1*lgr);              /* no need for x2, since +-phi has been taken into account*/
  printf(" J"); printf("\n");
  printf("             Rough estimate of total radiated energy (bunch) over hemisphere = ");
  printf("%.3e", sum1*0.1*lgr);              /* no need for x2, since +-phi has been taken into account*/
  printf(" J"); printf("\n");
  percent=sum1*0.1*lgr/(n_e_t*1.602177e-19*(gam-1)*0.511*1e6);
  printf(" This is equal to "); printf("%.3e",percent); printf(" of the energy in the bunch \n");  
  printf(" Energy in polarization 1= "); printf("%.3e", sum3*0.1*lgr); printf(" J \n");
  printf(" Energy in polarization 2= "); printf("%.3e", sum4*0.1*lgr); printf(" J \n"); 
  printf("\n"); printf("\n");
  printf("This case had the following parameters:\n");
  printf("gamma= %-5.1f \n",gam);
  printf("beam height(mm)= %-6.3f \n",x0);
  printf("pulse length (FWHM,ps)= %-6.3f \n",pulse_len);
  printf("total number of electrons in the bunch= %-6.3e \n",n_e_t);
  printf("pulse shape code is= %1i \n",pulse_shp);
  if (pulse_shp==8) printf("no. of Gaussians= %1i \n", ngauss);
  else
  printf("asymmetry factor= %-3.1f \n",epsl);
  if (pulse_shp==7)
  {apb=pulse_len; a=apb/(1+epsl);
   printf("for Lorentz shape a+b= %-5.3f ",apb);
   printf(" ps and leading edge fraction = %-5.3f \n",a);}
  printf("grating period(mm)= %-5.3f \n",el);
  printf("emission order= %1i \n",m);
  printf("PHI step= %-4.2f \n",phistep);
  printf("limiting value of phi= %-3.1f ",(phistep*(philmt-1))); printf(" deg. \n");
  printf("\n");
    
  printf("Angle    Theta2    Theta1    energy/bunch        energy(p1)          energy(p2)    polarisation\n");
  printf("                                   J                J                   J\n");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  {
    thet[i]=i;
    theta=thet[i]*pi/180;
    lgr_e = fabs(21/sin(theta));           /* variable effective grating length, determined basically by 
                                                     the cone entrance diameter of 21mm*/
	/*lgr_e = lgr;*/					     
      
   /*thet_1 and thet_2 below are the max. possible angles for a given observation angle*/
  /* NB. !!!!definition changed from previous versions of the code!!!*/
  /* now determined by the Winston cone acceptance angle */
    thet_1 = theta + 6*pi/180;
    thet_2 = theta - 6*pi/180;
    if (thet_1<0)   thet_1=pi+thet_1;
    if (thet_2<0)   thet_2=pi+thet_2;
    int_thet_1=(int)180*thet_1/pi;
    int_thet_2=(int)180*thet_2/pi;
    if (180*thet_1/pi-int_thet_1>=0.5) int_thet_1++;    
    if (180*thet_2/pi-int_thet_2>=0.5) int_thet_2++;   
    outbin[i][0]=theta*180/pi; outbin_single[i][0]=theta*180/pi;
    outbin[i][1]=int_thet_2; outbin_single[i][1]=int_thet_2;
    outbin[i][2]=int_thet_1; outbin_single[i][2]=int_thet_1;  
/*  n0 is total no. of theta angles accepted, modified May 08; must be 13 values */
    n0=int_thet_1-int_thet_2;
    if (n0<=0) n0=-n0;       
    if (int_thet_1<=int_thet_2) jstart= int_thet_1;
    else
    jstart=int_thet_2;

/* now average over accepted values of theta the already averaged values over phi */
       /* remember, col.9 of bin holds the previous averaging over phi*/   
    sum_th= sum_th_p1= sum_th_p2= sum_th_single= sum_th_single_p1= sum_th_single_p2= 0.0;
    j = 1;  
    for (j = 1; j<n0+2; j++)
    { sum_th = sum_th + bin[jstart+j-1][9];
      sum_th_p1 = sum_th_p1 + bin_p1[jstart+j-1][9];
      sum_th_p2 = sum_th_p2 + bin_p2[jstart+j-1][9];
      sum_th_single = sum_th_single + bin_single[jstart+j-1][9];
      sum_th_single_p1 = sum_th_single_p1 + bin_single_p1[jstart+j-1][9];
      sum_th_single_p2 = sum_th_single_p2 + bin_single_p2[jstart+j-1][9];
    };

    aven = sum_th/(n0 + 1);        /*have averaged over accepted theta & phi angles; in J/sr/cm*/
    aven_p1 = sum_th_p1/(n0 + 1);                          /*ditto for pol.1*/
    aven_p2 = sum_th_p2/(n0 + 1);                          /*and for pol.2 */
    aven_single = sum_th_single/(n0 + 1);
    aven_single_p1 = sum_th_single_p1/(n0 + 1);
    aven_single_p2 = sum_th_single_p2/(n0 + 1); 

    sum_tr = aven * 0.1 * lgr_e * domeg_m;    /*variable grating length, Sept 12; 
                                     now in J*/
    sum_tr_single = aven_single * 0.1 * lgr_e * domeg_m;

    sum_tr_p1 = aven_p1 * 0.1 * lgr_e * domeg_m;
    sum_tr_p2 = aven_p2 * 0.1 * lgr_e * domeg_m;
    sum_tr_single = aven_single * 0.1 * lgr_e * domeg_m;
    sum_tr_single_p1 = aven_single_p1 * 0.1 * lgr_e * domeg_m;
    sum_tr_single_p2 = aven_single_p2 * 0.1 * lgr_e * domeg_m;
    outbin[i][3] = aven; 
    outbin_single[i][3] = aven_single;
    outbin[i][4] = sum_tr_p1; 
    outbin_single[i][4] = sum_tr_single_p1;   /* store single-electr. numbers but don't print */
    outbin[i][5] = sum_tr_p2;
    outbin_single[i][5] = sum_tr_single_p2;
    outbin[i][6] = (sum_tr_p1-sum_tr_p2)/(sum_tr_p1+sum_tr_p2);
    outbin_single[i][6] = (sum_tr_single_p1-sum_tr_single_p2)/(sum_tr_single_p1+sum_tr_single_p2);
    printf("%-5.1f",theta*180/pi); printf("   "); printf("%5d",int_thet_2); printf("      ");
    printf("%5d",int_thet_1); printf("       ");
    printf("%-6.3e    ",sum_tr); printf("     ");
    printf("%-6.3e",sum_tr_p1); printf("           "); printf("%-6.3e",sum_tr_p2); printf("    ");
    printf("%-6.3e\n", (sum_tr_p1-sum_tr_p2)/(sum_tr_p1+sum_tr_p2)); 
    
  };
  printf("Angle    Theta2    Theta1    energy/bunch        energy(p1)          energy(p2)    polarisation\n");
  printf("                                   J                J                   J\n\n");

                /* Start of output section; differential (total) energy first */

  out_file=fopen("dE.txt","w");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { fprintf(out_file,"\n");
    ii=0;
    for (ii=0; ii<10; ii++)
    {if (ii==0) fprintf(out_file, "%-5.1f  ", bin[i][ii]);
    else fprintf(out_file, " %.3e ", bin[i][ii]); }       
  }
  fprintf(out_file,"\n");
  fprintf(out_file,"gamma= %6.1f \n",gam);
  fprintf(out_file,"x0(mm)= %6.3f \n",x0);
  fprintf(out_file,"fwhm(ps)= %6.3f \n",pulse_len);
  fprintf(out_file,"no. of electrons= %-4.3e \n",n_e_t);
  fprintf(out_file,"pulse shape= %1d \n",pulse_shp);
  if (pulse_shp==8) fprintf(out_file,"no. of Gaussians= %1i \n", ngauss);
  else
  fprintf(out_file,"asymmetry factor= %-3.1f \n",epsl);
  fprintf(out_file,"grating period(mm)= %-5.3f \n",el);
  fprintf(out_file, "grating width(mm)= %-5.3f \n", wid);
  fprintf(out_file,"order= %1d \n",m);
  fprintf(out_file,"phi step= %3.1f \n",phistep);
  fprintf(out_file,"first column= emission angle theta \n");
  fprintf(out_file,"second column= value for phi=0 \n");
  fprintf(out_file,"third column= value for phi= %3.1f \n", phistep);
  fprintf(out_file,"ninth column= value for phi= %3.1f \n", 7*phistep);
  fprintf(out_file,"tenth column= average between phi=0 & phi= +-%3.1f \n", 3*phistep); 
  fclose(out_file);
  
            /* next differential energy in the two polarisations */

  out_file=fopen("dE_pol.txt","w");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { ii=0;
    for (ii=0; ii<10; ii++)
    {if (ii==0) fprintf(out_file, "%-5.1f  ", bin[i][ii]);
    else fprintf(out_file, " %.3e ", bin_p1[i][ii]); }
    fprintf(out_file, "\n"); fprintf(out_file, "       ");
    ii=1;
    for (ii==1; ii<10; ii++)
    { fprintf(out_file,  " %.3e ", bin_p2[i][ii]);}
    fprintf(out_file, "\n");       
  }
  fprintf(out_file,"\n");
  fprintf(out_file,"gamma= %6.1f \n",gam);
  fprintf(out_file,"x0(mm)= %6.3f \n",x0);
  fprintf(out_file,"fwhm(ps)= %6.3f \n",pulse_len);
  fprintf(out_file,"no. of electrons= %-4.3e \n",n_e_t);
  fprintf(out_file,"pulse shape= %1d \n",pulse_shp);
  if (pulse_shp==8) fprintf(out_file,"no. of Gaussians= %1i \n", ngauss);
  else
  fprintf(out_file,"asymmetry factor= %-3.1f \n",epsl);
  fprintf(out_file,"grating period(mm)= %-5.3f \n",el);
  fprintf(out_file, "grating width(mm)= %-5.3f \n", wid);
  fprintf(out_file,"order= %1d \n",m);
  fprintf(out_file,"phi step= %3.1f \n",phistep);
  fprintf(out_file,"first column= emission angle theta \n");
  fprintf(out_file,"second column= value for phi=0 \n");
  fprintf(out_file,"third column= value for phi= %3.1f \n", phistep);
  fprintf(out_file,"ninth column= value for phi= %3.1f \n", 7*phistep);
  fprintf(out_file,"tenth column= average between phi=0 & phi= +-%3.1f \n", 3*phistep);
  fclose(out_file);

                       /* diff. energy for single electron */
  out_file=fopen("sing_elec.txt","w");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { fprintf(out_file,"\n");
    ii=0;
    for (ii=0; ii<10; ii++)
    {if (ii==0) fprintf(out_file, "%-5.1f  ", bin_single[i][ii]);
    else fprintf(out_file, " %.3e ", bin_single[i][ii]); }
    fprintf(out_file, "\n"); fprintf(out_file, "       ");
    ii=1;
    for (ii==1; ii<10; ii++)
    { fprintf(out_file,  " %.3e ", bin_single_p1[i][ii]);}
    fprintf(out_file, "\n"); fprintf(out_file, "       ");
    ii=1;
    for (ii==1; ii<10; ii++)
    { fprintf(out_file,  " %.3e ", bin_single_p2[i][ii]);}       
  }
  fprintf(out_file,"\n");
  fprintf(out_file,"gamma= %6.1f \n",gam);
  fprintf(out_file,"x0(mm)= %6.3f \n",x0);
  fprintf(out_file,"grating period(mm)= %-5.3f \n",el);
  fprintf(out_file, "grating width(mm)= %-5.3f \n", wid);
  fprintf(out_file,"order= %1d \n",m);
  fprintf(out_file,"phi step= %3.1f \n",phistep);
  fprintf(out_file,"first column= emission angle theta \n");
  fprintf(out_file,"second column= value for phi=0 \n");
  fprintf(out_file,"third column= value for phi= %3.1f \n", phistep);
  fprintf(out_file,"ninth column= value for phi= %3.1f \n", 7*phistep);
  fprintf(out_file,"tenth column= average between phi=0 & phi= +-%3.1f \n", 3*phistep);
  fclose(out_file);
  
                   /* degree of polarisation derived from single electron values */
  out_file=fopen("pol.txt","w");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { fprintf(out_file,"\n");
    ii=0;
    for (ii=0; ii<9; ii++)
    {if (ii==0) fprintf(out_file, "%-5.1f  ", bin_single[i][ii]);
    else fprintf(out_file, " %.3e ", (bin_single_p1[i][ii] - bin_single_p2[i][ii])/
                                     (bin_single_p1[i][ii] + bin_single_p2[i][ii]));
    fprintf(out_file, "  "); }       
  }
  fprintf(out_file,"\n");
  fprintf(out_file,"gamma= %6.1f \n",gam);
  fprintf(out_file,"x0(mm)= %6.3f \n",x0);
  fprintf(out_file,"grating period(mm)= %-5.3f \n",el);
  fprintf(out_file, "grating width(mm)= %-5.3f \n", wid);
  fprintf(out_file,"order= %1d \n",m);
  fprintf(out_file,"phi step= %3.1f \n",phistep);
  fprintf(out_file,"first column= emission angle theta \n");
  fprintf(out_file,"second column= value for phi=0 \n");
  fprintf(out_file,"third column= value for phi= %3.1f \n", phistep);
  fprintf(out_file,"ninth column= value for phi= %3.1f \n", 7*phistep);
  fprintf(out_file,"tenth column= average between phi=0 & phi= +-%3.1f \n", 3*phistep);
  fclose(out_file);
  printf("\n");
 
                     /* Collected energy, bunch & single electron */
  out_file=fopen("ener.txt","w");   
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { fprintf(out_file,"\n"); ii=0;
    for (ii=0; ii<7; ii++)          
    {if (ii<3) fprintf(out_file, " %-5.1f \t", outbin[i][ii]);
     /*else if (ii==3) fprintf(out_file, " %-6.3e \t", outbin[i][ii]);*/
     else  fprintf(out_file, " %-6.3e \t ", outbin[i][ii]); }             
  }
  fprintf(out_file,"\n"); fprintf(out_file,"\n");
  fprintf(out_file, " single electron data follow below \n");
  fprintf(out_file,"\n");
  i=(int)(alpha1*180/pi);
  for (i=(int)(alpha1*180/pi); i<171; i++)
  { fprintf(out_file,"\n"); ii=0;
    for (ii=0; ii<7; ii++)          
    {if (ii<3) fprintf(out_file, " %-5.1f \t", outbin_single[i][ii]);
     /*else if (ii==3) fprintf(out_file, " %-6.3e \t", outbin_single[i][ii]);*/
     else  fprintf(out_file, " %-6.3e \t", outbin_single[i][ii]); }
  }   
  fprintf(out_file,"\n");
  fprintf(out_file," This case had the following parameters: \n");
  fprintf(out_file,"gamma= %6.1f \n",gam);
  fprintf(out_file,"beam height(mm)= %6.3f \n",x0);
  fprintf(out_file,"pulse length (FWHM, ps)= %6.3f \n",pulse_len);
  fprintf(out_file,"no. of electrons= %-4.3e \n",n_e_t);
  fprintf(out_file,"pulse shape code= %1d \n",pulse_shp);
  if (pulse_shp==8) fprintf(out_file,"no. of Gaussians= %1d \n",ngauss);   
  else if (pulse_shp==7)
  {fprintf(out_file,"Lorentz shape with a+b= %-5.3f ",apb);
   fprintf(out_file, "and leading fraction of bunch= %-4.2f \n",a);}
  else
  fprintf(out_file,"asymmetry factor= %-3.1f \n",epsl);
  fprintf(out_file,"grating period(mm)= %-5.3f \n",el);
  fprintf(out_file, "grating width(mm)= %-5.3f \n", wid);
  fprintf(out_file,"order= %1d \n",m);
  fprintf(out_file,"phi step= %3.1f \n",phistep);
  fprintf(out_file,"limiting phi= %3.1f ",((philmt-1)*phistep)); fprintf(out_file," deg.\n");
  fprintf(out_file,"first column= emission angle theta \n");
  fprintf(out_file,"second column= theta2 \n");
  fprintf(out_file,"third column= theta1 \n");
  fprintf(out_file,"fourth column= diff. energy per bunch J/sr/cm, total \n");
  fprintf(out_file,"fifth column= energy per bunch (J) in P1 \n");
  fprintf(out_file,"sixth column= energy per bunch (J) in P2\n");
  fprintf(out_file,"seventh column= degree of polarisation (averaged) \n");
  fprintf(out_file,"total radiated energy per bunch(J)=%-4.3e \n", sum1);
  fprintf(out_file,"total radiated energy single electron(J)=%-4.3e \n", sum2);    
  fprintf(out_file, " \n ");
  fclose(out_file);
  printf("Output of differential energy vs. theta & phi in dE.txt\n ");
  printf("Output of differential energy in p1 &p2 vs. theta & phi in dE_pol.txt\n ");
  printf("Output of single electron diff.energy (total & p1, p2) vs. theta & phi in sing_elec.txt\n ");
  printf("Degree of polarisation, based on single-electron yield, in pol.txt\n ");
  printf("Collected energy in file ener.txt\n ");
  printf("All done!!\n\n");
  
                                    /*re-entry options*/
  printf("Type:\n");
  printf("1 for new beam position \n");
  printf("2 for new sigma_x \n");
  printf("3 for new bunch shape\n");
  printf("4 for output file for CAUCHY & exit\n");
  printf("anything else to exit this programme\n");
  scanf("%d",&ans1);
  if (ans1==1) goto esc20;
  else if (ans1==2) goto esc24;
  else if (ans1==3) goto esc1;
  else if (ans1==4)
  {printf("Name of output file?? ");
   scanf("%s",&outfile_name[0]);
   out_file=fopen(outfile_name,"w");   
   i=(int)(alpha1*180/pi);
   for (i=(int)(alpha1*180/pi); i<171; i++)
   { if ( (outbin[i][0] - ((int)(outbin[i][0] / 10))*10 ) == 0 && outbin[i][0] >= 40 && outbin[i][0] <= 140 )
    { fprintf(out_file, " %-5.1f\t%-6.3e\t1\n", outbin[i][0], (outbin[i][4]+outbin[i][5]));
    }
   }
   fprintf( out_file, "%-5.3f\n", el ); 
   fclose(out_file);}
esc999:  printf(" End of programme\n");

return (0);
}

complex Ccreate(double r, double i)
{
	complex t;
	t.r = r;
	t.i = i;
	return t;
}

complex Cadd(complex a, complex b)
{
	complex t;
	t.r = a.r + b.r;
	t.i = a.i + b.i;
	return t;
}

complex Cmultiply(complex a, complex b)
{
	complex t;
	t.r = a.r * b.r - a.i * b.i;
	t.i = a.i * b.r + a.r * b.i;
	return t;
}

complex Cdivide(complex a, complex b)
{
	complex t;
	double n;
	n = b.r * b.r + b.i * b.i;
	t.r = (a.r * b.r + a.i * b.i) / n;
	t.i = (a.i * b.r - a.r * b.i) / n;
	return t;
}
	
complex Cconjugate(complex a)
{
	complex t;
	t.r = a.r;
	t.i = - a.i;
	return t;
}



vector Vcreate(double xr, double xi, double yr, double yi, double zr, double zi)
{
	vector temp;

	temp.x = Ccreate(xr, xi);
	temp.y = Ccreate(yr, yi);
	temp.z = Ccreate(zr, zi);

	return temp;
}

vector Vadd(vector v1, vector v2)
{
	vector temp;

	temp.x = Cadd(v1.x, v2.x);
	temp.y = Cadd(v1.y, v2.y);
	temp.z = Cadd(v1.z, v2.z);

	return temp;
}

vector Vmultiply(vector v, complex c)
{
	vector temp;

	temp.x = Cmultiply(v.x, c);
	temp.y = Cmultiply(v.y, c);
	temp.z = Cmultiply(v.z, c);

	return temp;
}

vector Vconjugate(vector v)
{
	vector temp;

	temp.x = Cconjugate(v.x);
	temp.y = Cconjugate(v.y);
	temp.z = Cconjugate(v.z);

	return temp;
}

complex Vdot(vector v1, vector v2)
{
	complex temp;

	temp = Cadd(Cmultiply(v1.x, v2.x), Cmultiply(v1.y, v2.y));
	temp = Cadd(temp, Cmultiply(v1.z, v2.z));

	return temp;
}


static double QU1(Integer ndim, double var[])                     /*integrand of real part of Gx, Gz*/
{  
   extern double theta,  xv, k, kx, alpha, beta, gam, D, dum, psi;
   double u, tem1, arg, y, z;
   static NagError fail;
   
   y=var[0]; z=var[1];   
   u= xv -z*tan(alpha)+ dum;
   arg= k*sqrt(u*u+y*y)/(beta*gam);
   tem1= s18adc(arg,&fail);
   return (2*k*u*cos(D*z-kx*dum+psi)*tem1*cos(ky*y)/(sqrt(u*u+y*y)*beta*gam));
}

static double QU2(Integer ndim, double var[])                     
{  
   extern double theta, xv, k, kx, alpha, beta, gam, D, dum, psi;
   double u, tem1, arg, y, z;
   static NagError fail;
      
   y=var[0]; z=var[1];   
   u= xv -z*tan(alpha)+ dum;
   arg= k*sqrt(u*u+y*y)/(beta*gam);
   tem1= s18adc(arg,&fail);                      
   return (2*k*u*sin(D*z-kx*dum+psi)*tem1*cos(ky*y)/(sqrt(u*u+y*y)*beta*gam));
}

static double QU3(Integer ndim, double var[])                    
{  
   extern double theta, xv, k, kx, alpha, beta, gam, D, dum, psi;
   double u, tem1, arg, y, z;
   static NagError fail;
   
   y=var[0]; z=var[1];   
   u= xv -z*tan(alpha)+ dum;
   arg= k*sqrt(u*u+y*y)/(beta*gam);
   tem1= s18adc(arg,&fail);
   return (2*k*y*sin(D*z-kx*dum+psi)*tem1*sin(ky*y)/(sqrt(u*u+y*y)*beta*gam));
}

static double QU4(Integer ndim, double var[])                    
{  
   extern double theta, xv, k, kx, alpha, beta, gam, D, dum, psi;
   double u, tem1, arg, y, z;
   static NagError fail;
   
   y=var[0]; z=var[1];   
   u= xv -z*tan(alpha)+ dum;
   arg= k*sqrt(u*u+y*y)/(beta*gam);
   tem1= s18adc(arg,&fail);
   return (2*k*y*cos(D*z-kx*dum+psi)*tem1*sin(ky*y)/(sqrt(u*u+y*y)*beta*gam));                               
}

static double GAU(double x)
{
   extern double x0,sigma_x;
   return (exp(-(x-x0)*(x-x0)/(2*sigma_x*sigma_x)));
}

static double FUN(double x)
{
   double kk;
   extern double x0,sigma_x,lam_e;
      kk=2/lam_e;
   return (exp(-(x-x0)*(x-x0)/(2*sigma_x*sigma_x))*exp(-kk*x));
}

static double FUN_X(double x)
{
   double kk;
   extern double x0,sigma_x,lam_e;
      kk=1/lam_e;
   return (exp(-(x-x0)*(x-x0)/(2*sigma_x*sigma_x))*exp(-kk*x));
}

static double LOR1(double t)
{ 
    extern double a;
    return (1/(a*a+t*t));
}

static double LOR2(double t)
{
    extern double b;
    return (1/(b*b+t*t));
}

    
    
    
