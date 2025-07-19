//CONSISTENT TRANSPORT CLSVOF ALGORITHM USING WLIC-THINC SCHEME FOR VOLUME FRACTIONS AND ENO2 SCHEME FOR LEVEL SETS.
//BOUNDARY CONDITIONS ARE AS FOLLOWS.
//LEFT AND RIGHT BOUNDARY - FREE SLIP AND NO PENETRATION
//BOTTOM AND TOP BOUNDARY - NO SLIP AND NO PENETRATION
class CLSVOF:public MBASE
{
	int **tag;	//tags of the interfacial cells
	double **Ft;	//intermediate volume fractions
	double **Phit;	//intermediate level set functions
	double **rhot;	//intermediate density field
	double **A_xt,**A_yt;	//intermediate advection terms for X and Y momentum equations
	double beta;	//parameter to control slope and thickness of interface jump
	double mass_act;	//actual mass

	void updt_ghost(double **Phi);	//update the ghost cells of cell centered field
	double ENO2(int flag,int i,int j,double **Phi,double V);	//ENO2 scheme
	double UP1(int flag,int i,int j,double **Phi,double V);	//1st order upwind scheme
	double QUICK(int flag,int i,int j,double **Phi,double V);	//QUICK scheme
	void tag_X(double **Fa);	//tagging algorithm in X direction
	void tag_Y(double **Fa);	//tagging algorithm in Y direction
	void adv_X(int t);	//solve the advection equation(implicit discretization in X direction)
	void adv_Y(int t);	//solve the advection equation(implicit discretization in Y direction)
	double LS_corr(double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02);	//enforce mass conservation to the level set function
	void reinit(int t);	//LS reinitialization algorithm
	public:
			CLSVOF(); ~CLSVOF();
			void ini(double w,double h,double xc,double yc);
			void den_ini();	//initialize density field at n+1 time step
			void solve(int n);	//CLSVOF advection algorithm
			void prop_updt();	//update density and advection field for next time step
			void mass_err();	//calculate mass error
			void lsvf_write(int t);	//tecplot file output
			void ls_complete(int t);	//ls file output including the ghost cells
};
CLSVOF::CLSVOF()
{
	tag=new int*[J+2];
	Phit=new double*[J+2];
	Ft=new double*[J+2];
	rhot=new double*[J+1];
	A_xt=new double*[J+1];
	A_yt=new double*[J+1];
	for(int i=0;i<J+2;i++)
	{
		tag[i]=new int[I+2];
		Phit[i]=new double[I+2];
		Ft[i]=new double[I+2];
		if(i<J+1)
		{
			rhot[i]=new double[I+1];
			A_xt[i]=new double[I+1];
			A_yt[i]=new double[I+1];
		}
	}
	cout<<"CLSVOF: MEMORY ALLOCATED"<<endl;
}
CLSVOF::~CLSVOF()
{
	for(int i=0;i<J+2;i++)
	{
		delete[] tag[i];
		delete[] Phit[i];
		delete[] Ft[i];
		if(i<J+1)
		{
			delete[] rhot[i];
			delete[] A_xt[i];
			delete[] A_yt[i];
		}
	}
	delete[] tag;
	delete[] Phit;
	delete[] Ft;
	delete[] rhot;
	delete[] A_xt;
	delete[] A_yt;
	cout<<"CLSVOF: MEMORY RELEASED"<<endl;
}
void CLSVOF::ini(double w,double h,double xc,double yc)
{
	beta=2.3;	//three mesh space smoothing
	//beta=3.5;	//one mesh space smoothing
	cout<<"CLSVOF: beta = "<<beta<<endl;
	INI ms(Xm,Ym,CX,CY,F,Phi,(beta/dx),w,h,xc,yc);
	ms.LS(tag);
	ms.VF();	//initial exact volume fraction is calculated here
	updt_ghost(F);
	mass_act=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass_act+=F[j][i];
	for(int j=1;j<=J;j++)	//initialize density, viscosity, and advection field
	{
		for(int i=1;i<=I;i++)
		{
			rho_n[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			A_x[j][i]=u[j][i];
			A_y[j][i]=v[j][i];
		}
	}
	for(int j=1;j<=J;j++)	//left and right ghost cell values
	{
		rho_n[j][0]=rho_n[j][1];
		mu[j][0]=mu[j][1];
		rho_n[j][I+1]=rho_n[j][I];
		mu[j][I+1]=mu[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost cell values
	{
		rho_n[0][i]=rho_n[1][i];
		mu[0][i]=mu[1][i];
		rho_n[J+1][i]=rho_n[J][i];
		mu[J+1][i]=mu[J][i];
	}
}
void CLSVOF::den_ini()
{
	for(int j=0;j<=J+1;j++)	//calculate density field at n+1 step (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
}
void CLSVOF::prop_updt()
{
	for(int j=0;j<=J+1;j++)	//including ghost nodes
	{
		for(int i=0;i<=I+1;i++)
		{
			rho_n[j][i]=rho_np1[j][i];
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			if((i>=1)&&(i<=I)&&(j>=1)&&(j<=J))	//only inner domain
			{
				A_x[j][i]=u[j][i];
				A_y[j][i]=v[j][i];
			}
		}
	}
}
void CLSVOF::updt_ghost(double **Phia)
{
	for(int j=1;j<=J;j++)	//left and right ghost nodes (Neumann bc)
	{
		Phia[j][0]=Phia[j][1];
		Phia[j][I+1]=Phia[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost nodes (Neumann bc)
	{
		Phia[0][i]=Phia[1][i];
		Phia[J+1][i]=Phia[J][i];
	}
}
double CLSVOF::ENO2(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	double plus,minus;	//plus and minus flux
	if(flag==0)	//X direction flux
	{
		plus=Phia[j][i+1]-0.5*MINMOD((Phia[j][i+1]-Phia[j][i]),(Phia[j][i+2]-Phia[j][i+1]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j][i-1]),(Phia[j][i+1]-Phia[j][i]));
	}
	else if(flag==1)	//Y direction flux
	{
		plus=Phia[j+1][i]-0.5*MINMOD((Phia[j+1][i]-Phia[j][i]),(Phia[j+2][i]-Phia[j+1][i]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j-1][i]),(Phia[j+1][i]-Phia[j][i]));
	}
	if(V>0.0) return minus;
	else return plus;
}
double CLSVOF::UP1(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j][i+1];
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j+1][i];
	}
	else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
}
double CLSVOF::QUICK(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return (0.75*Phia[j][i]+0.375*Phia[j][i+1]-0.125*Phia[j][i-1]);
		else return (0.75*Phia[j][i+1]+0.375*Phia[j][i]-0.125*Phia[j][i+2]);
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return (0.75*Phia[j][i]+0.375*Phia[j+1][i]-0.125*Phia[j-1][i]);
		else return (0.75*Phia[j+1][i]+0.375*Phia[j][i]-0.125*Phia[j+2][i]);
	}
	else { cout<<"CLSVOF: ERROR IN QUICK SCHEME!"<<endl; return 0; }
}
void CLSVOF::tag_X(double **Fa)
{
	for(int j=0;j<=J+1;j++)	//reinitialize the cell tags (including ghost cells)
		for(int i=0;i<=I+1;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//inner domain
	{
		for(int i=1;i<=I;i++)
		{
			if((Fa[j][i]>TRUNC_l)&&(Fa[j][i]<(1.0-TRUNC_u)))	//locate interfacial cell
			{
				tag[j][i]=1;	//tag interfacial cell
				tag[j][i+1]=1; tag[j][i-1]=1;	//tag neighbouring cells in X direction
			}
		}
	}
}
void CLSVOF::tag_Y(double **Fa)
{
	for(int j=0;j<=J+1;j++)	//reinitialize the cell tags (including ghost cells)
		for(int i=0;i<=I+1;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//inner domain
	{
		for(int i=1;i<=I;i++)
		{
			if((Fa[j][i]>TRUNC_l)&&(Fa[j][i]<(1.0-TRUNC_u)))	//locate interfacial cell
			{
				tag[j][i]=1;	//tag interfacial cell
				tag[j+1][i]=1; tag[j-1][i]=1;	//tag neighbouring cells in Y direction
			}
		}
	}
}
double CLSVOF::LS_corr(double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02)
{
	double P[9],A[9];
	double beta1,om[3],eta[3];	//interface thickness parameter, Gaussian weights, and Gaussian points
	beta1=beta/dx;
	om[0]=4.0/9.0; om[1]=om[2]=5.0/18.0;	//normalized Gaussian weights
	eta[0]=0.0; eta[1]=sqrt(3.0/5.0); eta[2]=-sqrt(3.0/5.0);	//Gaussian points
	double gamma,temp,D=-1.0,C=2.0*(Fa-0.5);	//initial guess = -1.0
	double x1,y1;
	double func,func1;	//function and its derivative
	int cnt=0;	//no of iterations of NR method
	for(int l2=0;l2<3;l2++)	//calculation of gamma
	{
		for(int l1=0;l1<3;l1++)
		{
			x1=(0.5*dx)*eta[l1]; y1=(0.5*dy)*eta[l2];
			P[l2*3+l1]=(a_00+a_10*x1+a_01*y1+a_11*x1*y1+a_20*pow(x1,2.0)+a_02*pow(y1,2.0));	//calculation of interfacial polynomial
			temp=beta1*P[l2*3+l1];
			if((l1==0)&&(l2==0)) gamma=temp;	//initialize gamma in the 1st iteration
			else if(gamma>temp) gamma=temp;
		}
	}
	gamma=-gamma+EPS;	//final value of gamma
	for(int l2=0;l2<3;l2++)	//calculation of A_g
		for(int l1=0;l1<3;l1++)
			A[l2*3+l1]=tanh(beta1*P[l2*3+l1]+gamma);
	do	//Newton-Raphson method
	{
		func=func1=0.0;	//reinitialization
		temp=D;	//store value of previous iteration (variable is reused)
		for(int l2=0;l2<3;l2++)	//calculation of function and its derivatives
		{
			for(int l1=0;l1<3;l1++)
			{
				func+=om[l1]*om[l2]*(A[l2*3+l1]+D)/(1.0+A[l2*3+l1]*D);
				func1+=om[l1]*om[l2]*(1.0-pow(A[l2*3+l1],2.0))/(pow((1.0+A[l2*3+l1]*D),2.0));
			}
		}
		func-=C;
		D-=func/func1;
		cnt++;
		if(func>0.0) cout<<"CLSVOF: ERROR IN FUNC"<<endl;
		if(func1<0.0) cout<<"CLSVOF: ERROR IN FUNC1"<<endl;
	}
	while(abs(D-temp)>=EPS);
	return ((atanh(D)+gamma)/beta1);
}
void CLSVOF::reinit(int t)
{
	double h=MAX2(dx,dy);
	double temp,a,b;
	double a_00,a_01,a_10,a_11,a_20,a_02;	//interface polynomial coefficients
	//--------------------INITIALIZATION SCHEME (BASED ON ADVECTED LS FIELD)--------------------------------------------------
	for(int j=0;j<=J+1;j++)	//reinitialize the tag values and Phit values
	{
		for(int i=0;i<=I+1;i++)
		{
			tag[j][i]=0;
			Phit[j][i]=0.0;
		}
	}
	for(int j=1;j<=J;j++)	//tagging and correction based on volume fractions
	{
		for(int i=1;i<=I;i++)
		{
			if((F[j][i]>0.01)&&(F[j][i]<0.98))	//interfacial cell
			{
				tag[j][i]=1;
				a_00=Phi[j][i];
				a_10=0.5*(Phi[j][i+1]-Phi[j][i-1])/dx;
				a_01=0.5*(Phi[j+1][i]-Phi[j-1][i])/dy;
				a_11=0.25*(Phi[j+1][i+1]-Phi[j+1][i-1]-Phi[j-1][i+1]+Phi[j-1][i-1])/(dx*dy);
				a_20=0.5*(Phi[j][i+1]-2.0*Phi[j][i]+Phi[j][i-1])/(pow(dx,2.0));
				a_02=0.5*(Phi[j+1][i]-2.0*Phi[j][i]+Phi[j-1][i])/(pow(dy,2.0));
				Phit[j][i]=Phi[j][i]+LS_corr(F[j][i],a_00,a_10,a_01,a_11,a_20,a_02);
			}
		}
	}
	for(int j=1;j<=J;j++)	//reset LS in untagged cells and correct LS values in tagged cells
	{
		for(int i=1;i<=I;i++)
		{
			if(tag[j][i]==0) Phi[j][i]=SGN(Phi[j][i])*100.0;
			else if(tag[j][i]==1) Phi[j][i]=Phit[j][i];
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(free-slip boundaries)
	{
		Phi[j][0]=Phi[j][1]; tag[j][0]=tag[j][1];
		Phi[j][I+1]=Phi[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(no-slip boundaries)
	{
		Phi[0][i]=Phi[1][i]; tag[0][i]=tag[1][i];
		Phi[J+1][i]=Phi[J][i]; tag[J+1][i]=tag[J][i];
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
void CLSVOF::adv_X(int t)
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and advected velocity for X and Y momentum equation
	double Fx,Fy;	//WLIC advective flux
	double nx,ny,omx,omy;	//interface normal and weights for WLIC scheme
	int iup,jup,alp_x,alp_y,lam_x,lam_y;	//constants of THINC scheme
	double a1,a3,a32,xt,yt;	//interface parameters of THINC scheme
	double dtdx=dt/dx,dtdy=dt/dy;	//constant ratios for the advection equation
	a3=exp(beta); a32=pow(a3,2.0);
	tag_X(F);	//tag cells for advection in X direction
	for(int j=1;j<=J;j++)	//implicit discretization in X direction
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phi,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=0.0;	//no slip
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phi,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=0.0;	//no slip
			}
			else
			{
				iup=(u_EW[j][i]<0.0)?(i+1):i;
				if(abs(F[j][iup])<=TRUNC_l) F_adv[1]=0.0;	//empty cell
				else if(abs(F[j][iup]-1.0)<=TRUNC_u) F_adv[1]=u_EW[j][i];	//completely filled cell
				else	//interfacial cell (WLIC-THINC scheme)
				{
					lam_x=(u_EW[j][i]<0.0)?0:1;
					alp_x=(F[j][iup+1]>=F[j][iup-1])?1:-1;
					nx=0.5*abs(Phi[j][iup+1]-Phi[j][iup-1])/dx; ny=0.5*abs(Phi[j+1][iup]-Phi[j-1][iup])/dy;	//abs of interface normal
					omx=nx/(nx+ny); omy=1.0-omx;	//weights for WLIC
					a1=exp(beta/alp_x*(2.0*F[j][iup]-1.0));
					xt=0.5/beta*log((a32-a1*a3)/(a1*a3-1.0));
					Fx=-0.5*omx/dt*(-u_EW[j][i]*dt+alp_x*dx/beta*log((cosh(beta*(lam_x-xt-u_EW[j][i]*dt/dx)))/(cosh(beta*(lam_x-xt)))));
					Fy=omy*u_EW[j][i]*F[j][iup];
					F_adv[1]=Fx+Fy;	//WLIC-THINC scheme for inner domain
				}
				LS_adv[1]=ENO2(0,i,j,Phi,u_EW[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*u_EW[j][i];
				if((tag[j][i]==0)&&(tag[j][i+1]==0)&&(i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use QUICK scheme
				{
					u_f[1]=QUICK(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=QUICK(0,i,j,A_y,u_EW[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=UP1(0,i,j,A_y,u_EW[j][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dtdx*(F_adv[1]-F_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));	//implicit discretization in X direction
			Phit[j][i]=(Phi[j][i]-dtdx*(u_EW[j][i]*LS_adv[1]-u_EW[j][i-1]*LS_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			rhot[j][i]=(rho_n[j][i]-dtdx*(rho_adv[1]-rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dtdx*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dtdx*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
		}
	}
	updt_ghost(Ft);	//update ghost cells of Ft
	updt_ghost(Phit);	//update ghost cells of Phit
	tag_Y(Ft);	//tag cells for advection in Y direction
	for(int i=1;i<=I;i++)	//explicit discretization in Y direction
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phit,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phit,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_xt[J][i];	//outflow
				v_f[1]=A_yt[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			else
			{
				jup=(v_NS[j][i]<0.0)?(j+1):j;
				if(abs(Ft[jup][i])<=TRUNC_l) F_adv[1]=0.0;	//empty cell
				else if(abs(Ft[jup][i]-1.0)<=TRUNC_u) F_adv[1]=v_NS[j][i];	//completely filled cell
				else	//interfacial cell (WLIC-THINC scheme)
				{
					lam_y=(v_NS[j][i]<0.0)?0:1;
					alp_y=(Ft[jup+1][i]>=Ft[jup-1][i])?1:-1;
					nx=0.5*abs(Phit[jup][i+1]-Phit[jup][i-1])/dx; ny=0.5*abs(Phit[jup+1][i]-Phit[jup-1][i])/dy;	//abs of interface normal
					omx=nx/(nx+ny); omy=1.0-omx;	//weights for WLIC
					a1=exp(beta/alp_y*(2.0*Ft[jup][i]-1.0));
					yt=0.5/beta*log((a32-a1*a3)/(a1*a3-1.0));
					Fx=omx*v_NS[j][i]*Ft[jup][i];
					Fy=-0.5*omy/dt*(-v_NS[j][i]*dt+alp_y*dy/beta*log((cosh(beta*(lam_y-yt-v_NS[j][i]*dt/dy)))/(cosh(beta*(lam_y-yt)))));
					F_adv[1]=Fx+Fy;	//WLIC-THINC scheme for inner domain
				}
				LS_adv[1]=ENO2(1,i,j,Phit,v_NS[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*v_NS[j][i];
				if((tag[j][i]==0)&&(tag[j+1][i]==0)&&(i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use QUICK scheme
				{
					u_f[1]=QUICK(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=QUICK(1,i,j,A_yt,v_NS[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=UP1(1,i,j,A_yt,v_NS[j][i]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(F_adv[1]-F_adv[0]);	//explicit discretization in Y direction
			Phi[j][i]=Phit[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(rho_adv[1]-rho_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC_l) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-F[j][i])<=TRUNC_u) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::adv_Y(int t)
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and advected velocity for X and Y momentum equation
	double Fx,Fy;	//WLIC advective flux
	double nx,ny,omx,omy;	//interface normal and weights for WLIC scheme
	int iup,jup,alp_x,alp_y,lam_x,lam_y;	//constants of THINC scheme
	double a1,a3,a32,xt,yt;	//interface parameters of THINC scheme
	double dtdx=dt/dx,dtdy=dt/dy;	//constant ratios for the advection equation
	a3=exp(beta); a32=pow(a3,2.0);
	tag_Y(F);	//tag cells for advection in Y direction
	for(int i=I;i>=1;i--)	//implicit discretization in Y direction
	{
		for(int j=J;j>=1;j--)
		{
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phi,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_x[J][i];	//outflow
				v_f[1]=A_y[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phi,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			else
			{
				jup=(v_NS[j-1][i]<0.0)?j:(j-1);
				if(abs(F[jup][i])<=TRUNC_l) F_adv[0]=0.0;	//empty cell
				else if(abs(F[jup][i]-1.0)<=TRUNC_u) F_adv[0]=v_NS[j-1][i];	//completely filled cell
				else	//interfacial cell (WLIC-THINC scheme)
				{
					lam_y=(v_NS[j-1][i]<0.0)?0:1;
					alp_y=(F[jup+1][i]>=F[jup-1][i])?1:-1;
					nx=0.5*abs(Phi[jup][i+1]-Phi[jup][i-1])/dx; ny=0.5*abs(Phi[jup+1][i]-Phi[jup-1][i])/dy;	//abs of interface normal
					omx=nx/(nx+ny); omy=1.0-omx;	//weights for WLIC
					a1=exp(beta/alp_y*(2.0*F[jup][i]-1.0));
					yt=0.5/beta*log((a32-a1*a3)/(a1*a3-1.0));
					Fx=omx*v_NS[j-1][i]*F[jup][i];
					Fy=-0.5*omy/dt*(-v_NS[j-1][i]*dt+alp_y*dy/beta*log((cosh(beta*(lam_y-yt-v_NS[j-1][i]*dt/dy)))/(cosh(beta*(lam_y-yt)))));
					F_adv[0]=Fx+Fy;	//WLIC-THINC scheme for inner domain
				}
				LS_adv[0]=ENO2(1,i,j-1,Phi,v_NS[j-1][i]);	//ENO2 scheme for inner domain
				rho_adv[0]=(rho_1-rho_0)*F_adv[0]+rho_0*v_NS[j-1][i];
				if((tag[j][i]==0)&&(tag[j-1][i]==0)&&(i>2)&&(i<I-2)&&(j>2)&&(j<J-2))	//use QUICK scheme
				{
					u_f[0]=QUICK(1,i,j-1,A_x,v_NS[j-1][i]);
					v_f[0]=QUICK(1,i,j-1,A_y,v_NS[j-1][i]);
				}
				else	//use UP1 scheme
				{
					u_f[0]=UP1(1,i,j-1,A_x,v_NS[j-1][i]);
					v_f[0]=UP1(1,i,j-1,A_y,v_NS[j-1][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dtdy*(F_adv[1]-F_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));	//implicit discretization in Y direction
			Phit[j][i]=(Phi[j][i]-dtdy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			rhot[j][i]=(rho_n[j][i]-dtdy*(rho_adv[1]-rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dtdy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dtdy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			F_adv[1]=F_adv[0];	//updation for the next cell
			LS_adv[1]=LS_adv[0];
			rho_adv[1]=rho_adv[0];
			u_f[1]=u_f[0]; v_f[1]=v_f[0];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
		}
	}
	updt_ghost(Ft);	//update ghost cells of Ft
	updt_ghost(Phit);	//update ghost cells of Phit
	tag_X(Ft);	//tag cells for advection in X direction
	for(int j=J;j>=1;j--)	//explicit discretization in X direction
	{
		for(int i=I;i>=1;i--)
		{
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phit,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=0.0;	//no slip
			}
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phit,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=0.0;	//no slip
			}
			else
			{
				iup=(u_EW[j][i-1]<0.0)?i:(i-1);
				if(abs(Ft[j][iup])<=TRUNC_l) F_adv[0]=0.0;	//empty cell
				else if(abs(Ft[j][iup]-1.0)<=TRUNC_u) F_adv[0]=u_EW[j][i-1];	//completely filled cell
				else	//interfacial cell (WLIC-THINC scheme)
				{
					lam_x=(u_EW[j][i-1]<0.0)?0:1;
					alp_x=(Ft[j][iup+1]>=Ft[j][iup-1])?1:-1;
					nx=0.5*abs(Phit[j][iup+1]-Phit[j][iup-1])/dx; ny=0.5*abs(Phit[j+1][iup]-Phit[j-1][iup])/dy;	//abs of interface normal
					omx=nx/(nx+ny); omy=1.0-omx;	//weights for WLIC
					a1=exp(beta/alp_x*(2.0*Ft[j][iup]-1.0));
					xt=0.5/beta*log((a32-a1*a3)/(a1*a3-1.0));
					Fx=-0.5*omx/dt*(-u_EW[j][i-1]*dt+alp_x*dx/beta*log((cosh(beta*(lam_x-xt-u_EW[j][i-1]*dt/dx)))/(cosh(beta*(lam_x-xt)))));
					Fy=omy*u_EW[j][i-1]*Ft[j][iup];
					F_adv[0]=Fx+Fy;	//WLIC-THINC scheme for inner domain
				}
				LS_adv[0]=ENO2(0,i-1,j,Phit,u_EW[j][i-1]);	//ENO2 scheme for inner domain
				rho_adv[0]=(rho_1-rho_0)*F_adv[0]+rho_0*u_EW[j][i-1];
				if((tag[j][i]==0)&&(tag[j][i-1]==0)&&(i>2)&&(i<I-2)&&(j>2)&&(j<J-2))	//use QUICK scheme
				{
					u_f[0]=QUICK(0,i-1,j,A_xt,u_EW[j][i-1]);
					v_f[0]=QUICK(0,i-1,j,A_yt,u_EW[j][i-1]);
				}
				else	//use UP1 scheme
				{
					u_f[0]=UP1(0,i-1,j,A_xt,u_EW[j][i-1]);
					v_f[0]=UP1(0,i-1,j,A_yt,u_EW[j][i-1]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(F_adv[1]-F_adv[0]);	//explicit discretization in X direction
			Phi[j][i]=Phit[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(u_EW[j][i]*LS_adv[1]-u_EW[j][i-1]*LS_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(rho_adv[1]-rho_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[1]=F_adv[0];	//updation for the next cell
			LS_adv[1]=LS_adv[0];
			rho_adv[1]=rho_adv[0];
			u_f[1]=u_f[0]; v_f[1]=v_f[0];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC_l) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-F[j][i])<=TRUNC_u) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::solve(int n)
{
	if((n%2)==0) adv_Y(n);	//Strang splitting
	else adv_X(n);
	reinit(n);
}
void CLSVOF::lsvf_write(int t)
{
	string fname="ls_vol_frac_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::mass_err()
{
	double mass=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass+=F[j][i];
	cout<<"CLSVOF: MASS ERROR = "<<(mass-mass_act)<<endl;
}
void CLSVOF::ls_complete(int t)
{
	string fname="ls_comlpete_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS INCLUDING GHOST CELLS\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+2<<", J="<<J+2<<", DATAPACKING=BLOCK, SOLUTIONTIME="<<t*dt<<endl;
	double x=0.0-0.5*dx,y=0.0-0.5*dy;
	for(int j=0;j<=J+1;j++)	//X coordinates
	{
		for(int i=0;i<=I+1;i++)
		{
			p_out<<" "<<x;
			x+=dx;
		}
		p_out<<endl;
		x=0.0-0.5*dx;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)	//Y coordinates
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<y;
		p_out<<endl;
		y+=dy;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: COMPLETE LEVEL SET FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
