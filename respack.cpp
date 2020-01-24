#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm> 
#include <stdio.h>
#include <math.h>
#include <string>
#include <ctime>
#include <time.h>  
#include <sys/stat.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include "gsl/gsl_errno.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "LHAPDF/LHAPDF.h"
using namespace std;
//using namespace LHAPDF;

//--integrators
namespace QAG{
	double integrate(double (*f)(double,void *),double a,double b,double rerr,int key){
		//--key=1,6
		double result, error;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
		//gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		gsl_integration_qag (&F, a,b, 0,rerr, 50000,1,w, &result, &error);
		gsl_integration_workspace_free (w);
		return result;
		}
}
namespace QAGS{
	double integrate(double (*f)(double,void *),double a,double b,double rerr){
		double result, error;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
		//gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		gsl_integration_qags (&F, a,b, 0,rerr, 50000,w, &result, &error);
		//gsl_integration_qag (&F, a,b, 0,rerr, 50000,1,w, &result, &error);
		gsl_integration_workspace_free (w);
		return result;
		}
}
namespace CQUAD{

	double integrate(double (*f)(double,void *),double a,double b,double rerr){
		double result, abserror;
		gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (100);
		//gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		size_t neval= 100;
		gsl_integration_cquad (&F, a,b, 0,rerr,w,&result, &abserror,&neval);
		gsl_integration_cquad_workspace_free (w);
		return result;
		}
} 
namespace QAGIU{
	double integrate(double (*f)(double,void *),double a,double rerr){
		double result, error;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
		gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		gsl_integration_qagiu(&F, a, 0.0, rerr, 50000,w, &result, &error);
		gsl_integration_workspace_free (w);
		return result;
	 	}
} 
namespace QAGIL{
	double integrate(double (*f)(double,void *),double a,double rerr){
		double result, error;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (50000);
		gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		gsl_integration_qagil(&F, a, 0.0, rerr, 50000,w, &result, &error);
		gsl_integration_workspace_free (w);
		return result;
	 	}
} 
namespace GL{
	double result, error;
	double integrate(double (*f)(double,void *),double a,double b,int n){
		gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc (n);
		//gsl_set_error_handler_off ();
		gsl_function F;
		F.function = f;
		F.params = NULL;
		result=gsl_integration_glfixed (&F,a ,b ,t);
		gsl_integration_glfixed_table_free;
		return result;
	 	}

} 

//--QCD params
namespace QCD{ 
	int flag=0; 
	double CF=4.0/3.0;
	double CA=3.0;
	double NC=3.0;
	double TF=0.5;
	double pi=M_PI;
	double pi2=pi*pi;
	double Nf=0.0;
	double b0=0.0;
	double b1=0.0;

	//double get_alphaS(double muR){
	//	double alphaS=LHAPDF::alphasPDF(muR);
	//	Nf=LHAPDF::getNf();
	//	return alphaS;
	//}
}

//--fragmentation function
namespace BFG{
	//------------------------------------------------------------           
	//-- Subroutine to compute a point in x and q**2         
	//-- FOR THE FOLLOWING FRAGMENTATION FUNCTIONS:           
	//-- GLUON,  UP   ,DOWN, STRANGE , CHARM , BOTTOM , TOP .     
	//-- Range of validity :           
	//--		1.E-03 < X < .99       
	//--		LOG10 2  <LOG10 Q**2 < 10        
	//-- AUTHORS :	P CHIAPPETTA
	//--				J.PH. GUILLET                        
	//--				M.GRECO                              
	//-- CGW FOR Q02=2GEV2         METHOD P.NASON
	//-- NEXT TO LEADING EVOLUTION
	//-- TABLE FOR LHC
	//-- 	I = 1    GLUON
	//-- 	I = 2    UP ( = ANTIUP )
	//-- 	I = 3    DOWN ( = ANTIDOWN )
	//-- 	I = 4    STRANGE ( = ANTISTRANGE )
	//-- 	I = 5    CHARM ( = ANTICHARM )
	//-- 	I = 6    BOTTOM ( = ANTIBOTTOM )
	//-- 	I = 7    TOP ( = ANTITOP )
	//-- Z IS X AND QSTAR2 IS Q SQUARE         
	//-- STRFUN specifies the desired function 
	//-- 	IFUN =  1     GLUON           'GLUON' (abbrev. 'GL')                   
	//-- 	        2     UP              'UPVAL' (abbrev. 'UP')                   
	//-- 	        3     DOWN            'DOVAL' (abbrev. 'DO')                   
	//-- 	        4     STRANGE         'SBAR'  (abbrev. 'SB')                   
	//-- 	        5     CHARM           'CBAR ' (abbrev. 'CB')                   
	//-- 	        6     BOTTOM          'BBAR ' (abbrev. 'BB')                   
	//-- 	        7     TOP             'TBAR ' (abbrev. 'TB')                   
	//------------------------------------------------------------           

 	struct data{
		int NX;
		int NQ2;
		int NY;
		int IQ2MAX;
		int NF;
		double xmin,xmax;
		double Q2min,Q2max;
		vector<double>Q2;
		vector<double>X;
		vector<double>Y;
		vector<vector<vector<double> > >Z;
	};

	data D4,D5;

	const int nbx=1; //--only odd values
	const int nby=1; //--only odd values
	vector<double> bx; 
	vector<double> by; 
	vector<double> ibx; 
	vector<double> iby; 
	vector<vector<double> >bz;

	double result,err;

 	void setup(){

		const int NX=45;
		//double CC[3];
		double X[NX]={
			0.4627e-01,0.5131e-01,0.5690e-01,0.6311e-01,0.6999e-01,0.7762e-01, 
			0.8608e-01,0.9547e-01,0.1059e+00,0.1174e+00,0.1302e+00,0.1444e+00,
			0.1602e+00,0.1776e+00,0.1970e+00,0.2185e+00,0.2423e+00,0.2687e+00,
			0.2980e+00,0.3305e+00,0.3666e+00,0.4065e+00,0.4508e+00,0.5000e+00,
			0.5321e+00,0.5643e+00,0.5964e+00,0.6286e+00,0.6607e+00,0.6929e+00,
			0.7250e+00,0.7571e+00,0.7893e+00,0.8214e+00,0.8536e+00,0.8857e+00,
			0.90e0    ,0.9179e+00,0.93e0    ,0.94e0    ,0.9500e+00,0.96e0    ,
			0.9700e+00,0.98e0    , 0.99e+00 };

		//--set 4-flav data

		D4.NX=45;
		D4.NQ2=8;
		D4.NY=D4.NQ2;
		D4.IQ2MAX=22;
		D4.NF=6;

		D4.X.resize(D4.NX);
		D4.Y.resize(D4.NY);
		D4.Q2.resize(D4.NQ2);

		double Q2_4fv[]={2.0,2.8,3.9,5.4,7.5,10.5,14.5,20.24};
		for (int i=0;i<D4.NQ2;i++) D4.Q2[i]=Q2_4fv[i];
		for (int i=0;i<D4.NQ2;i++) D4.Y[i]=log10(Q2_4fv[i]);
		for (int i=0;i<D4.NX;i++) D4.X[i] =log10(X[i]);


		D4.xmin=X[0];
		D4.xmax=X[NX-1];
		D4.Q2min=D4.Q2[0];
		D4.Q2max=D4.Q2[D4.NQ2-1];

		
		//--set 5-flav data

		D5.NX=45;
		D5.NQ2=22;
		D5.NY=D5.NQ2;
		D5.IQ2MAX=22;
		D5.NF=6;

		D5.X.resize(D5.NX);
		D5.Y.resize(D5.NY);
		D5.Q2.resize(D5.NQ2);

		double Q2_5fv[]={20.26,33,55,93,156,261,435,732,1224,2048,
										3427,5733,9592,16047,26847,44915,75144,1.257e+5,
										2.1e+5,3.52e+5,5.89e+5,9.848e+5};
		for (int i=0;i<D5.NQ2;i++) D5.Q2[i]=Q2_5fv[i];
		for (int i=0;i<D5.NQ2;i++) D5.Y[i]=log10(Q2_5fv[i]);
		for (int i=0;i<D5.NX;i++) D5.X[i]=log10(X[i]);


		D5.xmin=X[0];
		D5.xmax=X[NX-1];
		D5.Q2min=D5.Q2[0];
		D5.Q2max=D5.Q2[D5.NQ2-1];


		//--load table == pert+nonpert

		D4.Z.resize(D4.NQ2);
		for (int i=0;i<D4.NQ2;i++){
			D4.Z[i].resize(6);
			}
	
		D5.Z.resize(D5.NQ2);
		for (int i=0;i<D5.NQ2;i++){
			D5.Z[i].resize(6);
			}

		ifstream table("bfg_table.dat");//this is set II
		//ifstream table("test.dat");	
		string str;
		string sD4 ("D4");
		string sD5 ("D5");
	
		char ch[2];
		int i,j,k;
		
		while (getline(table,str)){		

			if (str.find(sD4)!=std::string::npos){
				istringstream is(str);
				is>>ch>>i>>j;
				//cout<<4<<" "<<i<<" "<<j<<endl;	
				k=4;
				continue;
				}

			if (str.find(sD5)!=std::string::npos){
				istringstream is(str);
				is>>ch>>i>>j;
				//cout<<5<<" "<<i<<" "<<j<<endl;	
				k=5;
				continue;
				}

			istringstream is(str);
			double x;
			while (is>>x){
				//cout<<x<<endl;
				//cout<<x<<endl;	
				if (k==4) D4.Z[i-1][j-1].push_back(x);
				if (k==5) D5.Z[i-1][j-1].push_back(x);
				}
			}



		table.close();


		bx.resize(nbx*2+1);
		by.resize(nby*2+1);
		ibx.resize(nbx*2+1);
		iby.resize(nby*2+1);
		bz.resize(nbx*2+1);
		for (int i=0;i<nbx*2+1;i++) bz[i].resize(nby*2+1);
		} 

	int locate(vector<double> & xx,int n,double x){

		int jl=0;
		int ju=n+1;
		do{ 

			int jm=(ju+jl)/2;
			if((xx[n-1]>=xx[0])&&(x>=xx[jm])) {
				jl=jm;
				}
			else {
				ju=jm;
				}
		}while (ju-jl>1); 

		int j;
		if(x==xx[0]) j=1;
		else if(x==xx[n-1]) j=n-2;
		else j=jl;
		return j;		
	}

	void get_block(int ix,int iy,data & D,int ifun){

		//--get bx
		if (ix-nbx<=0){
			for (int i=0;i<nbx*2+1;i++){
				bx[i]=D.X[i];
				ibx[i]=i;
				}
			}
		else if (ix-D.NX>=0){
			for (int i=-1;i<nbx*2;i++){
				bx[i+1]=D.X[D.NX-nbx*2+i-1];
				ibx[i+1]=D.NX-nbx*2+i-1;
				}		
			}
		else{
			for (int i=-1;i<nbx*2;i++){
				bx[i+1]=D.X[ix-nbx+i+1];
				ibx[i+1]=ix-nbx+i+1;
				}
			}

		//--get by
		if (iy-nby<=0){
			for (int i=0;i<nby*2+1;i++){
				by[i]=D.Y[i];
				iby[i]=i;
				}
			}
		else if (iy-D.NY>=0){
			for (int i=-1;i<nby*2;i++){
				by[i+1]=D.Y[D.NY-nby*2+i-1];
				iby[i+1]=D.NY-nby*2+i-1;
				}		
			}
		else{
			for (int i=-1;i<nby*2;i++){
				by[i+1]=D.Y[iy-nby+i+1];
				iby[i+1]=iy-nby+i+1;
				}
			}

		//--get bz
		for (int i=0;i<nbx*2+1;i++){
			for (int j=0;j<nby*2+1;j++){
				bz[i][j]=D.Z[iby[j]][ifun][ibx[i]];
				}
			}
	}

	double polint(vector<double>xa,vector<double>ya,int n,double x){

		const int nmax=10;
		vector<double> c;
		vector<double> d;
		c.resize(nmax);
		d.resize(nmax);

		int ns=1;
		double dif=abs(x-xa[0]);
		

		for (int i=1;i<=n;i++){
			double dift=abs(x-xa[i-1]);
			if (dift<dif){
				ns=i;
				dif=dift;
				}
			c[i-1]=ya[i-1];
			d[i-1]=ya[i-1];
			//cout<<"--"<<c[i-1]<<" "<<d[i-1]<<endl;
			}

		double y=ya[ns-1];
		double dy=0;
		ns=ns-1;

		//cout<<"n="<<n<<endl;
		for (int m=1;m<=n-1;m++){

			for (int i=1;i<=n-m;i++){

				double ho=xa[i-1]-x;
				double hp=xa[i+m-1]-x;
				double w=c[i+1-1]-d[i-1];
				double den=ho-hp;

				den=w/den;
				d[i-1]=hp*den;
				c[i-1]=ho*den;
				//cout<<"--- "<<i-1<<"  "<<c[i-1]<<" "<<d[i-1]<<endl;
				}

			if (2*ns < n-m){
				dy=c[ns+1-1];
				//cout<<"+1+ "<<dy<<endl;
				}

			else {
				dy=d[ns-1];
				ns=ns-1;
				//cout<<"+2+ "<<dy<<endl;
				}

			y=y+dy;	
			//cout<<y<<endl;
		}
		return y;
	}

	double polint2(vector<double>&bx,vector<double>&by,vector<vector<double> >&bz, int m,int n,double x1,double x2){

		const int nmax=20;
		const int mmax=20;
		vector<double>ymtmp; ymtmp.resize(mmax);
		vector<double>yntmp; yntmp.resize(nmax);

		for (int j=0;j<m;j++){
			for (int k=0;k<n;k++){
				yntmp[k]=bz[j][k];
				//cout<<">>"<<yntmp[k]<<endl;
				}
			ymtmp[j]=polint(by,yntmp,n,x2);
			//cout<<ymtmp[j]<<endl;
			}
		return polint(bx,ymtmp,m,x1);
		}

	double interp(double z,double Q2,int ifun,data & D){                

		if (Q2 <= D.Q2min){ 
			printf("ERR: Q2 =%10.3e is too small\n ",Q2);
			printf("Q2min=%10.3e\n",D.Q2min);
			return 0;
			}

		double a=log10(z);
		double b=log10(Q2);


		int ix = locate(D.X,D.NX,a);
		int iy = locate(D.Y,D.NY,b);

		//for (int i=0;i<D.X.size();i++) printf("%3.d %10.3e %10.3e %3.d \n",i,a,D.X[i],ix);
		//for (int i=0;i<D.Y.size();i++)printf("%3.d %10.3e %10.3e %3.d \n",i,b,D.Y[i],iy);

		get_block(ix,iy,D,ifun);
	
		//cout<<"bx"<<endl;	
		//for (int i=0;i<bx.size();i++) cout<<bx[i]<<endl;
		//cout<<"by"<<endl;	
		//for (int i=0;i<by.size();i++) cout<<by[i]<<endl;
	
		//cout<<"ibx"<<endl;	
		//for (int i=0;i<ibx.size();i++) cout<<ibx[i]<<endl;
		//cout<<"iby"<<endl;	
		//for (int i=0;i<iby.size();i++) cout<<iby[i]<<endl;

		//for (int i=0;i<bx.size();i++){
		//	for (int j=0;j<by.size();j++){
		//		printf("%3d %3d %10.3e\n",i,j,bz[i][j]);
		//	
		//		}
		//	}

		//for (int i=0;i<3;i++){
		//	for (int j=0;j<3;j++){
		//		cout<<ibx[i]<<" "
		//				<<iby[j]<<" "
		//				<<D.Z[iby[j]][0][ibx[i]]
		//				<<endl;;
		//		}
		//	}

		//printf(">>>> %10.7e\n",D.Z[0][0][0]);
		//printf(">>>> %10.7e\n",D.Z[0][0][1]);
		//printf(">>>> %10.7e\n",D.Z[0][0][2]);
		//printf(">>>> %10.7e\n",D.Z[0][0][3]);

		//cout<<D.Z[0][0].size()<<endl;

		//D4.Z.resize(D4.NQ2);
		//for (int i=0;i<D4.NQ2;i++){
		//	D4.Z[i].resize(6);
		//	}

		return polint2(bx,by,bz,nbx*2+1,nby*2+1,a,b);                   
		}

	double get_ff(int ifun, double x,double Q2){
		double dum=0.0;
		//if (1e-3 < x && x < 0.99){
		if (x < 0.99){
			if (Q2<=20.26) dum=interp(x,Q2,ifun,D4);	
			else dum=interp(x,Q2,ifun,D5);
			}
		return dum/x; //--the fits are given as z*D(z)
	}

}



//--kinematics/options  
struct PARAMS{
	double rs;
	int imode;
	int iobs;
	int iord;
	int ieta;
	int ipT; //--for RES
	int irap; //--for eta or xF
	double pTmin;
	double pTmax;
	double etamin;
	double etamax;
	double zR;
	double zIF;
	double zFF;
	int isol;
	double R;
	double eps_pT;
	double eps;
 	};

//--inc_nlo
namespace INC{

	extern "C"  {
		void   setup_dir_ (double*,double*,double*,double*,double*);

		void   siglo_   	(double*,double*,double*,double*);

		double fqqb1_		  (double*);
		double fqqb1mu_		(double*);
		double fqqb2_			(double*,double*);
		double fqqb3_			(double*,double*);

		double fqg1_			(double*);
		double fqg1mu_		(double*);
		double fqg2_			(double*,double*);
		double fqg3_			(double*,double*);

		double fgq1_			(double*);
		double fgq1mu_		(double*);
		double fgq2_			(double*,double*);
		double fgq3_			(double*,double*);

		double fqqb_			(double*,double*);
		double fqg_				(double*,double*);
		double fgq_				(double*,double*);
		double fgg_				(double*,double*);
		double fqq_				(double*,double*);
		double fqb_				(double*,double*);
		double fqbs_			(double*,double*);
		double fqs_				(double*,double*,double*,double*);	

		double fqqbc_			(double*,double*);
		double fqgc_			(double*,double*);
		double fgqc_			(double*,double*);
		double fggc_			(double*,double*);
		double fqqc_			(double*,double*);
		double fqbc_			(double*,double*);
		double fqbsc_			(double*,double*);
		double fqsc_			(double*,double*,double*,double*);	

		void   setup_frag_ (double*,double*,double*,double*,double*,double*,double*,double*);
		void   fbor_  (double*,double*,double*);
		void   fmu_  (double*,double*,double*);
		void   stru_  (double*,double*,double*,double*,double*,double*,double*,double*,double*,
									double*,double*,double*,double*,double*,double*,double*,double*,double*,
									double*,double*,double*,double*,double*,double*,double*);

		double fdel1_  (double*,double*);
		double fdel2_  (double*,double*);
		double fvwpl1_  (double*,double*,double*);
		double fvwpl2_  (double*,double*,double*);
		double fvlo1_  (double*,double*,double*);
		double fvlo2_  (double*,double*,double*);
		double fresc1_  (double*,double*,double*);
		double fresc2_  (double*,double*,double*);
		double frescc1_  (double*,double*,double*);
		double frescc2_  (double*,double*,double*);

		void set_j0_  (int*);

		#define dir_aux_setup  setup_dir_
		#define siglo   siglo_   	

		#define FQQB1		fqqb1_			
		#define FQQB1MU	fqqb1mu_		
		#define FQQB2		fqqb2_			
		#define FQQB3		fqqb3_		

		#define FQG1		fqg1_			
		#define FQG1MU	fqg1mu_		
		#define FQG2		fqg2_			
		#define FQG3		fqg3_			

		#define FGQ1		fgq1_			
		#define FGQ1MU	fgq1mu_		
		#define FGQ2		fgq2_			
		#define FGQ3		fgq3_			
		
		#define FQQB	fqqb_	
		#define FQG		fqg_		
		#define FGQ		fgq_		
		#define FGG		fgg_		
		#define FQQ		fqq_		
		#define FQB		fqb_		
		#define FQBS	fqbs_	
		#define FQS		fqs_				

		#define FQQBC		fqqbc_	
		#define FQGC		fqgc_		
		#define FGQC		fgqc_		
		#define FGGC		fggc_		
		#define FQQC		fqqc_		
		#define FQBC		fqbc_		
		#define FQBSC	  fqbsc_	
		#define FQSC		fqsc_				


		#define fra_aux_setup   setup_frag_
		#define FBOR   					fbor_
		#define FMU   					fmu_
		#define STRU   					stru_

		#define FDEL1   					fdel1_
		#define FDEL2   					fdel2_
		#define FVWPL1   					fvwpl1_
		#define FVWPL2   					fvwpl2_
		#define FVLO1   					fvlo1_
		#define FVLO2   					fvlo2_
		#define FRESC1   					fresc1_
		#define FRESC2   					fresc2_
		#define FRESCC1   				frescc1_
		#define FRESCC2   				frescc2_

		#define set_j0  set_j0_

		}

	double rs;
	int iord;
	int ieta; // 0==full integration  1==partial integration 
	int irap; // 0==eta 1==xf
	int imode;
	int isol; 
	int iobs;
	double delta,epsil,R_,eps_pT,eps;
	double pTmin,pTmax;
	double etamin,etamax;
	double zR,zIF,zFF;
	LHAPDF::PDF* pdf1;
	LHAPDF::PDF* pdf2;

	void STRUCI(double X,double Q2,int j,int IH,double * U,double *UB,double *D,double *DB,double *S,double *C,double *B,double *G){

		//--NOTE: ALL DISTRIBUTIONS SHOULD BE X TIMES THE DISTRIBUTION !!!
		double Q,U0,D0,UB0;
		Q=pow(Q2,0.5);
		if (j==1){
			*G = pdf1->xfxQ(0, X, Q);
			*D = pdf1->xfxQ(1, X, Q);
			*U = pdf1->xfxQ(2, X, Q);
			*S = pdf1->xfxQ(3, X, Q);
			*C = pdf1->xfxQ(4, X, Q);
			*DB = pdf1->xfxQ(-1, X, Q);
			*UB = pdf1->xfxQ(-2, X, Q);
			*B = 0.0;

			//*G	=xfxM(j,X,Q,0);
			//*D	=xfxM(j,X,Q,1);
			//*U	=xfxM(j,X,Q,2);
			//*S	=xfxM(j,X,Q,3);
			//*C	=xfxM(j,X,Q,4);
			//*DB =xfxM(j,X,Q,-1);
			//*UB =xfxM(j,X,Q,-2);
			//*B	=0.0;//xfxM(j,X,Q,5);
		} else if (j==2) {
			*G = pdf2->xfxQ(0, X, Q);
			*D = pdf2->xfxQ(1, X, Q);
			*U = pdf2->xfxQ(2, X, Q);
			*S = pdf2->xfxQ(3, X, Q);
			*C = pdf2->xfxQ(4, X, Q);
			*DB = pdf2->xfxQ(-1, X, Q);
			*UB = pdf2->xfxQ(-2, X, Q);
			*B = 0.0;
		} else {
			std::cout << "Invalid j in method STRUCI!" << std::endl;
		}
		



		//--ppb
		if (imode == 1 && IH == 2){ 
			U0=*U;
			D0=*D;
			*U=*UB;
			*D=*DB;
			*UB=U0;
			*DB=D0;
			}
		//--pN
		else if (imode == 2 && IH == 2){
			U0 =0.5*(*U+*D);
			UB0=0.5*(*UB+*DB);
			*U=U0;
			*D=U0;
			*UB=UB0;
			*DB=UB0;
			}
		}
	void STRUCF(double z,double Q2,double *U,double *UB,double *D,double *DB,double *S,double *SB,double *C,double *CB,double *B,double *BB,double *G){
		*G =z*BFG::get_ff(0,z,Q2);
		*U =z*BFG::get_ff(1,z,Q2);
		*D =z*BFG::get_ff(2,z,Q2);
		*S =z*BFG::get_ff(3,z,Q2);
		*C =z*BFG::get_ff(4,z,Q2);
		*B =z*0.0;//BFG::get_ff(5,z,Q2);
		*UB=*U;
		*DB=*D;
		*SB=*S;
		*CB=*C;
		*BB=*B;
		}
	double integrand_dir(double * R, size_t dim, void * params){

		double alphaEM=0.0072973525376;
		double CF=4.0/3.0;
		double CA=3.0;
		double NC=3.0;
		double TF=0.5;

		double jac_pT=pTmax-pTmin;
		double pT=pTmin + R[0]*jac_pT; 
		if (isol==1){
			delta=R_;
			epsil=eps_pT/pT;
			if (epsil<eps) epsil=eps; 			
			}

		double pT2=pT*pT;
		double muR2=pow(zR*pT,2);
		double muIF2=pow(zIF*pT,2);
		double muFF2=pow(zFF*pT,2);
		double alphaS=pdf1->alphasQ2(muR2);
		double alphaSHO=alphaS*iord;
		
		double xT=2*pT/rs;
		double xT2=xT*xT;

		double etamin_,etamax_;
		if (ieta==0){
			etamax_ = log((1+sqrt(1-xT2))/xT);
			etamin_ =-etamax_;
			}
		if (ieta==1){
			etamax_ = etamax;
			etamin_ = etamin;
			}

		double jac_eta=etamax_-etamin_;
		double eta=etamin_ + R[1]*jac_eta; 
		if (ieta==1&&irap==1){
			double xF=eta;
			eta=log((xF+sqrt(xF*xF+xT2))/xT);
			}

		double vmin=xT/2.0*exp(eta);
		double vmax=1.0-xT/2.0*exp(-eta);
		double jac_v=vmax-vmin;
		double v=vmin+R[2]*jac_v;

		double wmin=xT*exp(eta)/2.0/v;
		double wmax=1.0;
		double jac_w=(wmax-wmin);
		double w=wmin+R[3]*jac_w;

		dir_aux_setup(&pT,&muIF2,&muFF2,&QCD::Nf,&delta);


		double x10=xT/2.0*exp(eta)/v;
		double x1=xT/2.0*exp(eta)/v/w;
		double x2=xT/2.0*exp(-eta)/(1.0-v);


		double GL10,U10,UB10,D10,DB10,C10,S10,B10;
		double GL1,U1,UB1,D1,DB1,C1,S1,B1;
		double GL2,U2,UB2,D2,DB2,C2,S2,B2;

		STRUCI(x10 ,muIF2,1,10,&U10 ,&UB10 ,&D10 ,&DB10 ,&S10 ,&C10 ,&B10 ,&GL10);
		//STRUCI(x1 ,muIF2,1,1,&U1 ,&UB1 ,&D1 ,&DB1 ,&S1 ,&C1 ,&B1 ,&GL1);
		STRUCI(x2 ,muIF2,2,2,&U2 ,&UB2 ,&D2 ,&DB2 ,&S2 ,&C2 ,&B2 ,&GL2);

		double SIGGQ,SIGQG,SIGQQB;
		siglo(&SIGQQB,&SIGQG,&SIGGQ,&v);

		double PQQB0,PQG0,PGQ0;
		PQQB0= 4.0/9.0*(U10*UB2+UB10*U2)	
					+1.0/9.0*(D10*DB2+DB10*D2)	
					+1.0/9.0*(S10*S2+S2*S10)		
					+4.0/9.0*(C10*C2+C2*C10);

		PQG0 =	(  4.0/9.0*(U10+UB10)			
							+1.0/9.0*(D10+DB10)			
							+1.0/9.0*(S10+S10)			
							+4.0/9.0*(C10+C10)) * GL2;

		PGQ0 =	(  4.0/9.0*(U2+UB2)			
							+1.0/9.0*(D2+DB2)
							+1.0/9.0*(S2+S2)
							+4.0/9.0*(C2+C2)) * GL10;

		double LAA=log(1.0-wmin);
		double LAA2=LAA*LAA;
		double SHAT=x10*x2*pow(rs,2);
		double XLMU=log(muR2/SHAT);


		double zero=0.0;	
		double one=1.0;	
		double FUNCHO1=
  	 PQQB0*(SIGQQB	+ alphaSHO*( FQQB1(&v) 		
															+XLMU*FQQB1MU(&v) 
															+LAA*FQQB2(&v,&one)
															+0.5*LAA2*FQQB3(&v,&one)))

		+PQG0*(SIGQG + alphaSHO*(	 FQG1(&v) 
														+XLMU*FQG1MU(&v)
														+LAA*FQG2(&v,&one)
														+0.5*LAA2*FQG3(&v,&one)))

		+PGQ0*(SIGGQ + alphaSHO*(	 FGQ1(&v) 
														+XLMU*FGQ1MU(&v)
														+LAA*FGQ2(&v,&one)
														+0.5*LAA2*FGQ3(&v,&one)));

		if (isol==1){
			double v2=v*v;
			double ELO=log(epsil/v/(1.0-wmin));
			double SUB =-alphaSHO/M_PI*pT2*pow(delta/cosh(eta),2)/SHAT
									*ELO*(SIGQQB*PQQB0*(CF-NC*v*(1.0-v))/v/(1.0-v)
									+SIGQG*PQG0*(CF*v2+NC*(1.0-v))/(1.0-v)/v
									+SIGGQ*PGQ0*(CF*pow(1.0-v,2)+NC*v)/v/(1.0-v) );
			FUNCHO1+=-SUB;
			}

		FUNCHO1*=jac_v;

		STRUCI(x1 ,muIF2,1,1,&U1 ,&UB1 ,&D1 ,&DB1 ,&S1 ,&C1 ,&B1 ,&GL1);

		double GQQB=FQQB(&v,&w);
		double GQG=FQG(&v,&w);
		double GGQ=FGQ(&v,&w);
		double GGG=FGG(&v,&w);
		double GQQ=FQQ(&v,&w);
		double GQB=FQB(&v,&w);
		double GQBS=FQBS(&v,&w);
		double GQS1=FQS(&v,&w,&one,&zero);
		double GQS2=FQS(&v,&w,&zero,&one);
		double GQS3=FQS(&v,&w,&one,&one)-GQS1-GQS2;
		double ERAT=v*(1.0-w)/(1.0-v+v*w);

		if ((isol==1)&&(ERAT>=epsil)){
			GQG=GQG-FQGC(&v,&w);
			GGQ=GGQ-FGQC(&v,&w);
			GGG=GGG-FGGC(&v,&w);
			GQQ=GQQ-FQQC(&v,&w);
			GQB=GQB-FQBC(&v,&w);
			GQBS=GQBS-FQBSC(&v,&w);
			double GQS1C=FQSC(&v,&w,&one,&zero);
			double GQS2C=FQSC(&v,&w,&zero,&one);
			double GQS3C=FQSC(&v,&w,&one,&one)-GQS1C-GQS2C;
			GQS1=GQS1-GQS1C;
			GQS2=GQS2-GQS2C;
			GQS3=GQS3-GQS3C;
			}

		double GS1=4.0/9.0 * GQS1 - 2.0/9.0 * GQS3 + 1.0/9.0 * GQS2;
		double GS2=4.0/9.0 * GQS1 + 2.0/9.0 * GQS3 + 1.0/9.0 * GQS2;
		double GS3=1.0/9.0 * GQS1 - 2.0/9.0 * GQS3 + 4.0/9.0 * GQS2;
		double GS4=1.0/9.0 * GQS1 + 2.0/9.0 * GQS3 + 4.0/9.0 * GQS2;
		double GS5=1.0/9.0 * GQS1 + 1.0/9.0 * GQS3 + 1.0/9.0 * GQS2;
		double GS6=1.0/9.0 * GQS1 - 1.0/9.0 * GQS3 + 1.0/9.0 * GQS2;
		double GS7=4.0/9.0 * GQS1 + 4.0/9.0 * GQS3 + 4.0/9.0 * GQS2;
		double GS8=4.0/9.0 * GQS1 - 4.0/9.0 * GQS3 + 4.0/9.0 * GQS2;

		double PQQB= 4.0/9.0*(U1*UB2+UB1*U2)
								+1.0/9.0*(D1*DB2+DB1*D2)
								+1.0/9.0*(S1*S2+S2*S1)
								+4.0/9.0*(C1*C2+C2*C1);
		double PQG =(	 4.0/9.0*(U1+UB1)
									+1.0/9.0*(D1+DB1)
									+1.0/9.0*(S1+S1)
									+4.0/9.0*(C1+C1)) * GL2;
		double PGQ =(	 4.0/9.0*(U2+UB2)
									+1.0/9.0*(D2+DB2)
									+1.0/9.0*(S2+S2)
									+4.0/9.0*(C2+C2) ) * GL1;
		double PGG =GL1*GL2* (4.0/9.0+1.0/9.0+1.0/9.0+4.0/9.0);
		double PQQ = 4.0/9.0*(U1*U2+UB1*UB2)
								+1.0/9.0*(D1*D2+DB1*DB2)
								+1.0/9.0*(S1*S2+S1*S2)
								+4.0/9.0*(C1*C2+C2*C1);
		double PQB = 4.0/9.0*(U1*UB2+UB1*U2)
								+1.0/9.0*(D1*DB2+DB1*D2)
								+1.0/9.0*(S1*S2+S2*S1)
								+4.0/9.0*(C1*C2+C2*C1);
		double PQBS=(U1*UB2+UB1*U2)*(1.0/9.0+1.0/9.0+4.0/9.0)
							+(D1*DB2+DB1*D2)*(1.0/9.0+4.0/9.0+4.0/9.0)
							+(S1*S2+S2*S1)*(4.0/9.0+1.0/9.0+4.0/9.0)
							+(C1*C2+C2*C1)*(4.0/9.0+1.0/9.0+1.0/9.0);
		double PQS =(U1*D2+UB1*DB2)*GS1
								+(UB1*D2+U1*DB2)*GS2
								+(D1*U2+DB1*UB2)*GS3
								+(DB1*U2+D1*UB2)*GS4
								+(U1*S2+UB1*S2+C1*D2+C1*DB2+2.0*C1*S2)*(GS1+GS2)
								+(S1*U2+S1*UB2+D1*C2+DB1*C2+2.0*S1*C2)*(GS3+GS4)
								+(S1*D2+S1*DB2+D1*S2+DB1*S2)*(GS5+GS6)
								+(U1*C2+UB1*C2+C1*U2+C1*UB2)*(GS7+GS8);

		double FUNCHO2 =  
			 (FQQB2(&v,&w)*PQQB-FQQB2(&v,&one)*PQQB0)/(1.0-w)
			+(FQQB3(&v,&w)*PQQB-FQQB3(&v,&one)*PQQB0)*log(1.0-w)/(1.0-w)
			+GQQB*PQQB
			+(FQG2(&v,&w)*PQG-FQG2(&v,&one)*PQG0)/(1.0-w)
			+(FQG3(&v,&w)*PQG-FQG3(&v,&one)*PQG0)*log(1.0-w)/(1.0-w)
			+GQG*PQG
			+(FGQ2(&v,&w)*PGQ-FGQ2(&v,&one)*PGQ0)/(1.0-w)
			+(FGQ3(&v,&w)*PGQ-FGQ3(&v,&one)*PGQ0)*log(1.0-w)/(1.0-w)
			+GGQ*PGQ
			+GGG*PGG
			+GQQ*PQQ
			+GQB*PQB
			+GQBS*PQBS
			+PQS;

		FUNCHO2*=alphaSHO*jac_v*jac_w;

		double dum=(FUNCHO1+FUNCHO2)*alphaEM*alphaS/M_PI/pow(pT,4);

		if (ieta==0) dum*=jac_eta;

		if (iobs==0) dum*=2*M_PI*pow(pT,4); //--pT^3 dsig/dpT
		if (iobs==1) dum*=1.0; 							//--E dsig/d3p
		if (iobs==2) dum*=2*M_PI*pT; 				//--dsig/dpT/deta
		if (iobs==3) dum*=2*M_PI*pT; 				//--dsig/dpT

		return dum;
 		}
	double integrand_fra(double * R, size_t dim, void * params){

		double CF=4.0/3.0;
		double CA=3.0;
		double NC=3.0;
		double TF=0.5;

		double AL=1.0;
		double CQ=0.0;
		double VC=NC*NC-1.0;
		double V1=VC*VC/NC;
		double V2=VC/NC;
		double V3=(pow(NC,4)-1.0)/2.0/(NC*NC);
		double V4=VC*VC/2.0/(NC*NC);
		double S=rs*rs;

		double jac_pT=pTmax-pTmin;
		double pT=pTmin + R[0]*jac_pT; 
		if (isol==1){
			delta=R_;
			epsil=eps_pT/pT;
			if (epsil<eps) epsil=eps; 			
			}

		double pT2=pT*pT;
		double muR2=pow(zR*pT,2);
		double muIF2=pow(zIF*pT,2);
		double muFF2=pow(zFF*pT,2);
		double alphaS=pdf1->alphasQ2(muR2);
		

		double xT=2*pT/rs;
		double xT2=xT*xT;

		if (ieta==0){
			etamax = log((1+sqrt(1-xT2))/xT);
			etamin =-etamax;
			}

		double jac_eta=etamax-etamin;
		double eta=etamin + R[1]*jac_eta; 
		if (ieta==1&&irap==1){
			double xF=eta;
			eta=log((xF+sqrt(xF*xF+xT2))/xT);
			}

		double V =1.0-pT/rs*exp(-eta);
		double W =pT2/S/V/(1.0-V);

		fra_aux_setup(&pT,&muIF2,&muFF2,&QCD::Nf,&delta,&S,&V,&W);
		
		double x3min=1.0-V+V*W;
		if (isol==1){
			double x3cut=1.0/(1.0+epsil);
			if (x3min<x3cut) x3min=x3cut;
			}
		double x3max=1.0;
		double jac_x3=x3max-x3min;
		double x3=x3min+R[2]*jac_x3;

		//x3=x3min+0.01;

		double vmin=V*W/x3;
		double vmax=1.0-(1.0-V)/x3;
		double jac_v=vmax-vmin;
		double v=vmin+R[3]*jac_v;

		//v=vmin+0.02;

		double wmin=V*W/x3/v;
		double wmax=1.0;
		double jac_w=wmax-wmin;
		double w=wmin+R[4]*jac_w;

		//w=wmin+0.03;

		double CC[16];
		for (int J0=0;J0<16;J0++){
			if (J0==15||J0==14) CC[J0]=8.0*VC*VC;
			else if(J0==13||J0==12||J0==9||J0==8||J0==7) CC[J0]=8.0*VC*NC;
			else CC[J0]=8.0*NC*NC;
			}

		double x1=V*W/v/w/x3;
		double x2=(1.0-V)/(1.0-v)/x3;
		double Bx1=V*W/v/x3;
		double Bx2=x2;
		double BXJAC=jac_x3*jac_v/Bx1/Bx2/(x3*x3);
		double XJAC=jac_x3*jac_v/x1/x2/(x3*x3);

		//if (x1>1.0 or x1<0.0) printf("ERR::x1=%10.5e\n",x1);
		//if (x2>1.0 or x2<0.0) printf("ERR::x2=%10.5e\n",x2);
		//if (Bx1>1.0 or Bx1<0.0) printf("ERR::Bx1=%10.5e\n",Bx1);
		//if (Bx2>1.0 or Bx2<0.0) printf("ERR::Bx2=%10.5e\n",Bx2);


		double GRRT[16];
		double GRRC[16];
		double GL1,UP1,UPB1,DO1,DOB1,CH1,ST1,BO1;
		double GL2,UP2,UPB2,DO2,DOB2,CH2,ST2,BO2;
		double DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,DPCB,DPB,DPBB,DPG;
		STRUCI(Bx1 ,muIF2,1,1,&UP1 ,&UPB1 ,&DO1 ,&DOB1 ,&ST1 ,&CH1 ,&BO1 ,&GL1);
		STRUCI(Bx2 ,muIF2,2,2,&UP2 ,&UPB2 ,&DO2 ,&DOB2 ,&ST2 ,&CH2 ,&BO2 ,&GL2);
		STRUCF(x3 ,muFF2,&DPU,&DPUB,&DPD,&DPDB,&DPS,&DPSB,&DPC,&DPCB,&DPB,&DPBB,&DPG);
		STRU(	&UP1,&UPB1,&DO1,&DOB1,&ST1,&CH1,&GL1,
					&UP2,&UPB2,&DO2,&DOB2,&ST2,&CH2,&GL2,
					&DPU,&DPUB,&DPD,&DPDB,&DPS,&DPSB,&DPC,&DPCB,&DPG,GRRT,GRRC);

		double s=Bx1*Bx2*S;

		double F01[16];
		double F02[16];
		double FDELMU1[16];
		double FDELMU2[16];

		double one_minus_v=1.0-v;
		FBOR(&v,&s,F01);
		FBOR(&one_minus_v,&s,F02);
		FMU(&v,&s,FDELMU1);
		FMU(&one_minus_v,&s,FDELMU2);
		double XLMU=log(s/muR2);

		double BORN[16];
		double DMU[16];
		for (int J0=0;J0<16;J0++){
			BORN[J0]=(F01[J0]*GRRT[J0] + F02[J0]*GRRC[J0])*BXJAC;
			DMU[J0]=(FDELMU1[J0]*GRRT[J0]+FDELMU2[J0]*GRRC[J0])*XLMU*BXJAC;
			}
	
		double GPPT[16];
		double GPPC[16];
		STRUCI(x1 ,muIF2,1,1,&UP1 ,&UPB1 ,&DO1 ,&DOB1 ,&ST1 ,&CH1 ,&BO1 ,&GL1);
		STRU(	&UP1,&UPB1,&DO1,&DOB1,&ST1,&CH1,&GL1,
					&UP2,&UPB2,&DO2,&DOB2,&ST2,&CH2,&GL2,
					&DPU,&DPUB,&DPD,&DPDB,&DPS,&DPSB,&DPC,&DPCB,&DPG,GPPT,GPPC);


		double UN=1.0;
		double ECOND  = 1.0-x3*(1.0+epsil)*(1.0-v+v*w);


		double FHODEL1[16];
		for (int i=0;i<16;i++) FHODEL1[i]=0.0;


		int IA1[]={0,2,4,5,10,11,12,13,14,15};
		for (int i=0;i<10;i++){
			int J0=IA1[i];
			set_j0(&J0);
			FHODEL1[J0]=	(FDEL1(&v,&x3) 
										+FVWPL1(&UN,&v,&x3)*log(1.0-wmin) 
										+FVLO1(&UN,&v,&x3)*pow(log(1.-wmin),2)/2.0)/(1.0-v) 
										-(wmax-wmin)/(1.0-v) 
										*(FVWPL1(&UN,&v,&x3)+FVLO1(&UN,&v,&x3)*log(1.0-w))/(1.0-w);
			
			//cout<<FDEL1(&v,&x3) <<endl;
			//cout<<FVWPL1(&UN,&v,&x3)<<endl;
			//cout<<log(1.0-wmin) <<endl;
			//cout<<FVLO1(&UN,&v,&x3)<<endl;
			//cout<<pow(log(1.-wmin),2)/2.0<<endl;
			//cout<<1.0/(1.0-v) <<endl;
			//cout<<(wmax-wmin)/(1.0-v) <<endl;
			//cout<<FVWPL1(&UN,&v,&x3)<<endl;
			//cout<<FVLO1(&UN,&v,&x3)<<endl;
			//cout<<log(1.0-w)/(1.0-w)<<endl;

			FHODEL1[J0]*=1.0/CC[J0]*BXJAC;

			}


		double FHODEL2[16];

		for (int i=0;i<16;i++) FHODEL2[i]=FHODEL1[i];

		int IA2[]={0,2,4,10,12,13};
		for (int i=0;i<6;i++){
			int J0=IA2[i];
			set_j0(&J0);
			FHODEL2[J0]= 	(FDEL2(&v,&x3)
										+FVWPL2(&UN,&v,&x3)*log(1.-wmin)
										+FVLO2(&UN,&v,&x3)*pow(log(1.-wmin),2)/2.0)/(1.0-v) 
										-(wmax-wmin)/(1.0-v) *
										(FVWPL2(&UN,&v,&x3)+FVLO2(&UN,&v,&x3)*log(1.-w))/(1.-w);
			FHODEL2[J0]*=1.0/CC[J0]*BXJAC;


 			}

		double FHOREST1[16];
		for (int J0=0;J0<16;J0++){
			set_j0(&J0);
			FHOREST1[J0]= (wmax-wmin)/(1.0-v) * (FVWPL1(&w,&v,&x3)/w/(1.0-w)
										+FVLO1(&w,&v,&x3)/w*log(1.0-w)/(1.0-w)+FRESC1(&w,&v,&x3)/w);

			if ((isol==1)&&(ECOND>0.0)){FHOREST1[J0]=FHOREST1[J0]-(wmax-wmin)*FRESCC1(&w,&v,&x3);	}

			//cout<<1<<"  "<<FHOREST1[J0]<<endl;
			//cout<<2<<"  "<<(wmax-wmin)*FRESCC1(&w,&v,&x3)<<endl;


			FHOREST1[J0]=FHOREST1[J0]/CC[J0]*XJAC;
			}

		double FHOREST2[16]={0};
		for (int i=0;i<16;i++) FHOREST2[i]=FHOREST1[i];

		for (int i=0;i<9;i++){
			int IA3[]={0,2,4,7,8,9,10,12,13};
			int J0=IA3[i];
			set_j0(&J0);
			FHOREST2[J0]=(wmax-wmin)/(1.0-v) * (FVWPL2(&w,&v,&x3)/w/(1.0-w)
									+FVLO2(&w,&v,&x3)/w*log(1.0-w)/(1.0-w)
									+FRESC2(&w,&v,&x3)/w);

			if ((isol==1)&&(ECOND>0.0)){
				FHOREST2[J0]=FHOREST2[J0]-(wmax-wmin)*FRESCC2(&w,&v,&x3);
				}
	
			FHOREST2[J0]=FHOREST2[J0]/CC[J0]*XJAC;
			}

		double GHD=0.0;
		double GHE=0.0;
		for (int J0=0;J0<16;J0++){
			GHD+=BORN[J0] + iord*alphaS*( FHODEL1[J0]*GRRT[J0]+FHODEL2[J0]*GRRC[J0]+DMU[J0]);
			GHE+=iord*alphaS*( FHOREST1[J0]*GPPT[J0] + FHOREST2[J0]*GPPC[J0] );
			}

		double dum=(GHD+GHE)*alphaS*alphaS/M_PI/S;

		if (ieta==0) dum*=jac_eta;
		if (iobs==0) dum*=2*M_PI*pow(pT,4); //--pT^3 dsig/dpT
		if (iobs==1) dum*=1.0; 							//--E dsig/d3p
		if (iobs==2) dum*=2*M_PI*pT; 				//--dsig/dpT/deta
		if (iobs==3) dum*=2*M_PI*pT; 				//--pT dsig/dpT

		return dum;
  	}	
	double xs_dir(){ 
		double res, err;
		double xl[4] = { 0, 0, 0 ,0};
		double xu[4] = { 1, 1, 1 ,1};

		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_monte_function G = { &integrand_dir,4, NULL};
		gsl_rng_env_setup ();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (4);
		gsl_monte_vegas_params  params;
		gsl_monte_vegas_params_get(s,&params);

		gsl_monte_vegas_integrate (&G, xl, xu, 4,1000,r, s,&res, &err);
		params.stage=2;

		//--single long evalation
		gsl_monte_vegas_integrate (&G, xl, xu, 4,10000, r, s,&res, &err);

		gsl_rng_free (r);
		return res;
 		}
	double xs_fra(){ 
		double res, err;
		double xl[5] = {0, 0, 0, 0 ,0};
		double xu[5] = {1, 1, 1, 1 ,1};

		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_monte_function G = { &integrand_fra,5, NULL};
		gsl_rng_env_setup ();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);
		gsl_monte_vegas_params  params;
		gsl_monte_vegas_params_get(s,&params);

		params.verbose=0;
		gsl_monte_vegas_integrate (&G, xl, xu, 5, 1000, r, s,&res, &err);
		params.stage=2;
		gsl_monte_vegas_integrate (&G, xl, xu, 5,10000, r, s,&res, &err);

		gsl_rng_free (r);
		return res;
  	}
	double xs(struct PARAMS * params){  
		rs=params->rs;
		imode=params->imode;
		iobs=params->iobs;
		iord=params->iord;
		ieta=params->ieta;
		irap=params->irap;
		pTmin=params->pTmin;
		pTmax=params->pTmax;
		etamin=params->etamin;
		etamax=params->etamax;

		cout<<"iobs"<<iobs<<endl;
		isol=params->isol;
		R_=params->R;
		eps_pT=params->eps_pT;
		eps=params->eps;


		zR=params->zR;
		zIF=params->zIF;
		zFF=params->zFF;	
		double dir=xs_dir();
		cout<<"dir="<<scientific<<dir<<endl;
		double fra=xs_fra();
		cout<<"fra="<<scientific<<fra<<endl;
		cout<<"tot="<<scientific<<dir+fra<<endl;

		return dir+fra;
		} 

 	}

//--main
int main( int argc, char  *argv[] ){

	if (argc!=3){
		printf("ERR: run as: $./pro card out\n");		
		return 0;
		}


	//--set parameters

	ifstream card(argv[1]);
	string pdfcol;
	int pdfID1,pdfID2,pQCD,irap;
	struct PARAMS params;
	
	string line;
	istringstream iss;

	//getline(card,line);	iss.str(line);	iss>>pdfcol;	
	getline(card,line);	iss.str(line);	iss>>pQCD;	
	getline(card,line);	iss.str(line);	iss>>pdfID1;	
	getline(card,line);	iss.str(line);	iss>>pdfID2;	
	getline(card,line);	iss.str(line);	iss>>params.rs;	
	getline(card,line);	iss.str(line);	iss>>params.imode;	
	getline(card,line);	iss.str(line);	iss>>params.iobs;	
	getline(card,line);	iss.str(line);	iss>>params.ipT;	
	getline(card,line);	iss.str(line);	iss>>params.zR;	
	getline(card,line);	iss.str(line);	iss>>params.zIF;	
	getline(card,line);	iss.str(line);	iss>>params.zFF;	
	getline(card,line);	iss.str(line);	iss>>params.etamin;	
	getline(card,line);	iss.str(line);	iss>>params.etamax;	
	getline(card,line);	iss.str(line);	iss>>params.irap;	
	getline(card,line);	iss.str(line);	iss>>params.isol;	
	getline(card,line);	iss.str(line);	iss>>params.R;	
	getline(card,line);	iss.str(line);	iss>>params.eps_pT;	
	getline(card,line);	iss.str(line);	iss>>params.eps;	

	//--load initial and final states distribution
	INC::pdf1 = LHAPDF::mkPDF(pdfID1);
	INC::pdf2 = LHAPDF::mkPDF(pdfID2);
	BFG::setup();

	ofstream file;
	file.open(argv[2]);
	
	file<<"pdfID1="<<pdfID1<<endl;
	file<<"pdfID2="<<pdfID2<<endl;	
	file<<setw(10)<<"pTmin";
	file<<setw(10)<<"pTmax";
	file<<setw(13)<<"xs";
	file<<endl<<endl;

	double xs=0.0;

	while(getline(card,line)){
		iss.str(line);	
		iss>>params.pTmin>>params.pTmax;	

		cout<<endl<<"----XS @ pT bin----"<<endl;
		cout.precision(3);		
		cout<<endl;
		cout<<setw(10)<<"pTmin="<<setw(10)<<fixed<<params.pTmin<<endl;
		cout<<setw(10)<<"pTmax="<<setw(10)<<fixed<<params.pTmax<<endl;


			params.isol=0;//-- use full z integration of isol 
			params.ieta=1;
			params.iobs=0;
			params.iord=1; 
			cout<<endl<<"computing incnlo2"<<endl;
			xs=INC::xs(&params);
			cout<<"xs="<<setw(13)<<scientific<<xs<<endl;


		file.precision(3);		
		file<<setw(10)<<fixed<<params.pTmin;
		file<<setw(10)<<fixed<<params.pTmax;
		file<<setw(13)<<scientific<<xs;
		file<<endl;

 		}

	return 0;
 	} 





















