#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  PI          3.14159265358979323846
#define  N           1000000
#define  DIS_N       500000
#define  PARAMS      4
#define  DATA_BIN    20

#define  GRID_POINTS_GAMMA  5
#define  GRID_POINTS_HM     9
#define  GRID_POINTS_HM_EXT 10
#define  GRID_POINTS_TEMP   9
#define  GRID_POINTS_TAU    5

const double Gammao [GRID_POINTS_GAMMA] = {1.7, 1.5, 1.3, 1.1, 0.9};
const double Uo     [GRID_POINTS_HM]    = {20.873074,17.848053,14.544626,11.548957,8.7345047,5.8747130,3.0728977};
const double UoExt  [GRID_POINTS_HM_EXT]= {20.873074,17.848053,14.544626,11.548957,8.7345047,5.8747130,3.0728977,0};
const double LogTo  [GRID_POINTS_TEMP]  = {4.2, 4.1, 4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4};
const double Tau    [GRID_POINTS_TAU]   = {1.9171994,       1.7574327,       1.5976661,       1.4378995,       1.2781329};

//Prior variables
double glow = 0.9,  ghigh = 1.7;  //better keep 0.9 gamma
double ulow = 3.0,  uhigh = 20.873;
double tlow = 3.4, thigh =  4.2;
  
//MCMC Routines
void   MCMC(double start[PARAMS], double var[PARAMS], double gammao[N], double uo[N], double logto[N], double tau[N], double likeli[N], double proposed_model_chain_write[N][DATA_BIN]);
double BMVariate(double U, double V, double mu, double sigma);
double Likeli(double X[]);
double Posterior(double X[PARAMS]);
double Prior(double X[PARAMS]);
void   GRConv(double chains[][N], int no_of_chains);
void   Trilinear(double X[PARAMS], double proposed_model[]);
double LinearInterp(double x, double f1, double x2, double f2, double p1);
void   ExtModel();

//Files reading and writing functions
int ReadModel(char *file);
int ReadMockData(char *file);
int WriteChainModels(char *file, double proposed_model_chain_write[N][DATA_BIN]);
int WriteChainParams(char *file, double gammao[], double uo[], double logto[], double tau[], double likeli[]);
int WriteProposedData(char *file, double proposeddata[]);

double Median(double array[], int size);
double Variance (double array[], int size);
double Mean(double array[], int size);

//General handy functions
int    Search (const double array[], const int size, double search);

const int no_models = 17;
const char Models[no_models][30]={"0.3HM01","0.8HM01","1.45HM01","2.2HM01","3.1HM01","4.2HM01","5.3HM01","Dz15cool","Dz12cool","Dz9cool","Dz7cool","Dz15","Dz12","Dz9","2.2HM01_g1.4","2.2HM01_g1.2","2.2HM01_g1.0"};
double LogTo5[no_models] ={3.6764058,3.9824309,4.1587368,4.2823762,4.3797316,4.4705590,4.5280076,3.9230869,3.9281259,3.9225704,3.9297924,3.9862149,4.0122729,4.2126883,4.2781419,4.2635771,4.2541653};
double Gammao5[no_models]={1.4258432,1.4590440,1.4701708,1.4824612,1.4683188,1.4747988,1.4757064,1.4932402,1.4967602,1.4951383,1.4749911,1.5063205,1.4950317,1.5231377,1.3700523,1.0838902,0.91937702};


double data     [DATA_BIN];
double yerr     [DATA_BIN];
double cov      [DATA_BIN][DATA_BIN];
double model    [GRID_POINTS_GAMMA][GRID_POINTS_HM][GRID_POINTS_TEMP][GRID_POINTS_TAU][DATA_BIN];
double modelext [GRID_POINTS_GAMMA][GRID_POINTS_HM_EXT][GRID_POINTS_TEMP][GRID_POINTS_TAU][DATA_BIN];
double proposed_model_write[DATA_BIN];

int main()
{
  char   file[300],file_prop[300];
  int    i,j,k;
  double start[PARAMS], var[PARAMS];
  double gammao[N],uo[N],logto[N],tau[N],likeli[N];
  double proposed_model[DATA_BIN];
  double proposed_model_chain_write[N][DATA_BIN];
  double X[PARAMS],lX;   
  
  //Reading Model Data set for HM01 models
  sprintf(file,"model_SNR50_noconv_0.5pt.dat");  
  ReadModel(file);
  ExtModel();
  
  //Values from 3.1HM01 0.5 MCMC chain
  //X[0] = 1.1935356; X[1] = 15.714408;  X[2] = 4.3351302;  X[3] = 1.4662320;
  //parameters values for 3.1HM01
       
  //X[0] = 1.5231377; X[1] = 11.3;  X[2] = 4.2126883;  X[3] = 1.5275579;
  //Trilinear(X,proposed_model);
  //Reading Mock Data
  //sprintf(file,"SNR15_noconv_10qso_0.7pt/md_3.1HM01.dat");
  //ReadMockData(file);
  //printf("\n%f ",Likeli(X));
  //for (i=0;i<DATA_BIN;i++)
  //printf("%f,",proposed_model[i]);
  
  //sprintf(file_prop,"prop_mcmc_3.1HM01.dat");
  //WriteProposedData(file_prop,proposed_model);  

  for (i=3;i<4;i++)
  {
      //Reading Mock Data
      sprintf(file,"SNR50_noconv_50qso_0.5pt/bs1.3/md_%s.dat",Models[i]);
      //sprintf(file,"BeckerMD/md_Dz12cool_losn.dat");
      ReadMockData(file);
      
      for(j=0;j<DATA_BIN;j++)printf("\n%f  %f  %f",data[j],yerr[j],cov[j][j]);
  
      start[0] = 1.45;  start[1]=5.0; start[2]=4.0,  start[3]=1.52;   
      var[0]   =  0.1;  var[1]  =0.5;   var[2]=0.05,   var[3]=0.02;
      MCMC(start,var,gammao,uo,logto,tau,likeli,proposed_model_chain_write);
         
      sprintf(file,"SNR50_noconv_50qso_0.5pt/bs1.3/chain_%s.dat",Models[i]);  
      WriteChainParams(file,gammao,uo,logto,tau,likeli);
      sprintf(file,"SNR50_noconv_50qso_0.5pt/bs1.3/fits_%s.dat",Models[i]);      
      WriteChainModels(file,proposed_model_chain_write);
  }

 return 0;
}


int WriteChainParams(char *file, double gammao[], double uo[], double logto[], double tau[], double likeli[])
{  
  double temp,temp2;
  temp = N;
  temp = DIS_N;
  
  FILE *output;  
  
  if(!(output = fopen(file,"wb")))    
  {
    printf("\nFile Cant be opened");
    return 0;      
  }
  
  printf("\nWriting File %s",file);  
  //Writing header File
  //fwrite(&temp,sizeof(double),1,output);               
  //fwrite(&temp2,sizeof(double),1,output);                
  //Writing Prior values
  //fwrite(&ghigh  ,sizeof(double),1,output);                
  //fwrite(&glow   ,sizeof(double),1,output);                  
  //fwrite(&uhigh  ,sizeof(double),1,output);                  
  //fwrite(&ulow   ,sizeof(double),1,output);                  
  //fwrite(&thigh  ,sizeof(double),1,output);                  
  //fwrite(&tlow   ,sizeof(double),1,output);                  
  //fwrite(&tauhigh,sizeof(double),1,output);                
  //fwrite(&taulow ,sizeof(double),1,output); 
  
  //Writing Chains for paramters
  fwrite(gammao,sizeof(double),N,output);                
  fwrite(uo,    sizeof(double),N,output);                
  fwrite(logto, sizeof(double),N,output);                
  fwrite(tau,sizeof(double),N,output);                   
  fwrite(likeli,sizeof(double),N,output);                   

  fclose(output);  
  return 1;  
  
}


int WriteChainModels(char *file, double proposed_model_chain_write[N][DATA_BIN])
{  
  
  FILE *output; 
  int i,j;
  
  if(!(output = fopen(file,"wb")))    
  {
    printf("\nFile Cant be opened\n");
    return 0;      
  }
  
  printf("\nWriting File %s",file);              
  for (i=0;i<N;i++)
  {
    for (j=0;j<DATA_BIN;j++)
    {
    fwrite(&proposed_model_chain_write[i][j],sizeof(double),1,output); 
    }
  }
  fclose(output);  
  return 1;  
  
}


int Search (const double array[], const int size, double search)
{ 
  int i=0;  
    for (i=1;i<size;i++)             
      if(search >= array[i])    
      return i;      
}

//Linear interpolation in 1D
double LinearInterp(double x1, double x2, double f1, double f2, double p)
{
  double value;
  value = ( (p-x1)*(f2-f1)/(x2-x1) ) + f1;
  //printf("(%f %f %f) (%f %f)  (%f)\n",x1,x2,p,f1,f2, value);
  return value;
}  

void Trilinear(double X[], double proposed_model[])
{
  int i;
  int index[PARAMS]; 
  double a[8],b[4],c[2],d;
  double w;
  
  index[0] = Search(Gammao,GRID_POINTS_GAMMA ,X[0]);
  index[1] = Search(UoExt ,GRID_POINTS_HM_EXT,X[1]);  
  index[2] = Search(LogTo ,GRID_POINTS_TEMP  ,X[2]);   
  index[3] = Search(Tau   ,GRID_POINTS_TAU   ,X[3]);    

  for(i=0;i<DATA_BIN;i++)
  {     
    a[0]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]-1][index[2]-1][index[3]-1][i], modelext[index[0]][index[1]-1][index[2]-1][index[3]-1][i], X[0] ); 
    a[1]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]  ][index[2]-1][index[3]-1][i], modelext[index[0]][index[1]  ][index[2]-1][index[3]-1][i], X[0] ); 
    a[2]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]-1][index[2]  ][index[3]-1][i], modelext[index[0]][index[1]-1][index[2]  ][index[3]-1][i], X[0] ); 
    a[3]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]  ][index[2]  ][index[3]-1][i], modelext[index[0]][index[1]  ][index[2]  ][index[3]-1][i], X[0] ); 
    a[4]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]-1][index[2]-1][index[3]  ][i], modelext[index[0]][index[1]-1][index[2]-1][index[3]  ][i], X[0] ); 
    a[5]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]  ][index[2]-1][index[3]  ][i], modelext[index[0]][index[1]  ][index[2]-1][index[3]  ][i], X[0] ); 
    a[6]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]-1][index[2]  ][index[3]  ][i], modelext[index[0]][index[1]-1][index[2]  ][index[3]  ][i], X[0] ); 
    a[7]=LinearInterp(Gammao[index[0]-1],Gammao[index[0]],modelext[index[0]-1][index[1]  ][index[2]  ][index[3]  ][i], modelext[index[0]][index[1]  ][index[2]  ][index[3]  ][i], X[0] ); 
    
    b[0]=LinearInterp(UoExt[index[1]-1], UoExt[index[1]], a[0], a[1], X[1] );      
    b[1]=LinearInterp(UoExt[index[1]-1], UoExt[index[1]], a[2], a[3], X[1] );  
    b[2]=LinearInterp(UoExt[index[1]-1], UoExt[index[1]], a[4], a[5], X[1] );      
    b[3]=LinearInterp(UoExt[index[1]-1], UoExt[index[1]], a[6], a[7], X[1] );          
     
    c[0]=LinearInterp(LogTo[index[2]-1], LogTo[index[2]], b[0], b[1], X[2] );
    c[1]=LinearInterp(LogTo[index[2]-1], LogTo[index[2]], b[2], b[3], X[2] );
    
    d=LinearInterp(   Tau[index[3]-1], Tau[index[3]],   c[0], c[1], X[3] );
    proposed_model[i] = d;    
    proposed_model_write[i] = d;
  } 
}

double Likeli (double X[])
{
  double num=0, temp=0,temp2=0;
  double proposed_model[DATA_BIN],S[DATA_BIN][DATA_BIN];
  char   filename[300];
  int i,j;

  Trilinear(X,proposed_model);
  //Assuming all points are independent measurements
  for(i=0;i<DATA_BIN;i++)
  {//if yerr is small can cause possible very large negative values which gives zeros as exp(num)
   temp2  += pow(  (data[i]-pow(10,proposed_model[i])) / yerr[i],2) - log(yerr[i]) - 0.91894;
  }
  
  //Assuming they are correlated using inverse of covariance matrix
//   for(i=0;i<DATA_BIN;i++)
//   {   
//     temp = 0;
//     for(j=0;j<DATA_BIN;j++)
//     {
//       temp  += cov[i][j] * (data[j]-pow(10,proposed_model[j]));
//     }
//       temp2 += (data[i]-pow(10,proposed_model[i])) * temp; 
//   }
    
  
  temp2 *= -0.5; 
  return (temp2);
}

//Posterior function Posterior  = Likelihood x Prior
double Posterior (double X[PARAMS])
{
  if( Prior(X) != -INFINITY )
  return Likeli(X)+Prior(X);
  else 
    return -INFINITY;
}

//Prior constant probability function
double Prior (double X[PARAMS])
{ 
  double taueff_z5 = 1.5275579;
  //45% of the value at z=4.915 from becker sigma for taueff
  double sigma_taueff = 0.04*taueff_z5;
  double taueff_prior;
  
   //Gamma prior condition
   if(X[0]<= ghigh && X[0] >= glow)     
   {
     //Specific internal energy prior condition
     if(X[1]<= uhigh && X[1] >= ulow)     
     {             
       //logTo prior condition
       if(X[2]<= thigh && X[2] >= tlow)
       {
	 //Taueff prior gaussian distribution around values from becker 2010
	 //taueff_prior = pow(  (X[3]-taueff_z5) / sigma_taueff,2  ) - log(sigma_taueff) - 0.91894;
	 //return -0.5*taueff_prior;
	 //if(X[3]<= tauhigh && X[3] >= taulow)     
         //{
	 //Taueff prior gaussian distribution around values from becker 2010
	 taueff_prior = pow(  (X[3]-taueff_z5) / sigma_taueff,2  ) - log(sigma_taueff) - 0.91894;
	 return -0.5*taueff_prior;
	 //return 0;
	   	 
        //}              
        //else     
	  //return -INFINITY;	 
	   	 
       }              
       else     
	 return -INFINITY;     
     }     
     else     
       return -INFINITY;   
   }   
   else
     return -INFINITY;  
  
}



void MCMC (double start[PARAMS], double var[PARAMS], double gammao[N], double uo[N],
 double logto[N], double tau[N], double likeli[N], double proposed_model_chain_write[N][DATA_BIN])
{  
  int accept=0,trial=0, i,ii;
  double Xt[PARAMS],Y[PARAMS], accept_par[PARAMS];  
  double U,lX,lY;   
   
  srand48(time(NULL));
  accept = 0;      
  trial  = 0;
  
  //Starting paramters (gammao,uo,to)
  accept_par[0] = start[0];
  accept_par[1] = start[1];
  accept_par[2] = start[2];
  accept_par[3] = start[3];
      
  //Random variables      
  Xt[0] = BMVariate(drand48(),drand48(),accept_par[0],var[0]); 
  Xt[1] = BMVariate(drand48(),drand48(),accept_par[1],var[1]);        
  Xt[2] = BMVariate(drand48(),drand48(),accept_par[2],var[2]);        
  Xt[3] = BMVariate(drand48(),drand48(),accept_par[3],var[3]); 
  
  printf("\nRuning MCMC from (%f,%f,%f,%f):\n",start[0],start[1],start[2],start[3]);  
  lX=Posterior(Xt);  

      do
      {	
	Y[0] = BMVariate(drand48(),drand48(),accept_par[0],var[0]);
	Y[1] = BMVariate(drand48(),drand48(),accept_par[1],var[1]);
	Y[2] = BMVariate(drand48(),drand48(),accept_par[2],var[2]);
	Y[3] = BMVariate(drand48(),drand48(),accept_par[3],var[3]);
	   	     
              lY = Posterior(Y);
	      if(  log(drand48()) <= (lY-lX) )      
	      {
		//if( ( 1.4 < Y[0] && Y[0] < 1.5)  &&  ( 1.44 < Y[3] && Y[3] < 1.58) )      
	        //printf("\n%f %f %f %f ::%f ",Y[0],Y[1],Y[2],Y[3],lY); 		

		lX=lY;
		accept_par[0]   = Y[0];
		accept_par[1]   = Y[1];
		accept_par[2]   = Y[2];
		accept_par[3]   = Y[3];
		
		gammao[accept] = Y[0];
		uo    [accept] = Y[1];
		logto [accept] = Y[2];
		tau   [accept] = Y[3];
		likeli[accept] = lX;
		
		accept++;
		
		for (ii=0;ii<DATA_BIN;ii++)
		{
		  proposed_model_chain_write[accept][ii]=proposed_model_write[ii];
		}
	
	      }
	   
	      trial++;
	      
	      if(accept%100==0)	      
	      printf("\rAR %f   Completed %0.1f",(100*(double)accept/(double)trial),(double)accept/N*100.0);	      
		
      }while(accept<N);   
  
  printf("\nTrials = %d\n",N);
  printf("Acceptence Ratio = %0.3f%\n",(float)100.0*accept/trial);  
    
}


int WriteProposedData(char *file, double proposed_model[])
{  
  int i;  
  FILE *output;  
    
  if(!(output = fopen(file,"wb")))    
  {
    printf("\nFile Cant be opened");
    return 0;      
  }
  
  printf("\nWriting File %s",file);       
  fwrite(proposed_model ,sizeof(double),DATA_BIN,output);           
  return 1;
}


int ReadMockData(char *file)
{  
  int i,j;
  double temp;
  
  FILE *input;
    
  if(!(input = fopen(file,"rb")))    
  {
    printf("File (%s) Cant be opened\n",file);
    return 0;      
  }
  
  printf("\nReading File %s",file);       
  fread(data ,sizeof(double),DATA_BIN,input);         
  fread(yerr,sizeof(double),DATA_BIN,input);     
  
  for (i=0;i<DATA_BIN;i++)
  {
    for (j=0;j<DATA_BIN;j++)
    {
    fread(&temp,sizeof(double),1,input);     
    //printf("%g",temp);
    cov[i][j] = temp;
    }
  }
  
  
  fclose(input);   
  return 1;  
}


int ReadModel(char *file)
{
  
  int i,j,k,l,m;  
  FILE *input;  
  double temp;
  
  if(!(input = fopen(file,"rb")))    
  {
    printf("File (%s) Cant be opened\n",file);
    return 0;      
  }
  
  printf("Reading File %s\n",file);     
  
  for(i=0;i<GRID_POINTS_GAMMA;i++)
    for(j=0;j<GRID_POINTS_HM;j++)
      for(k=0;k<GRID_POINTS_TEMP;k++)
	for(l=0;l<GRID_POINTS_TAU;l++) 
	  for(m=0;m<DATA_BIN;m++)
	  {
	    fread(&model[i][j][k][l][m],sizeof(double),1,input);    	
	    //printf("%f\t",model[i][j][k][l][m]); 
	  }
  
  fclose(input); 	
  return 1; 
}


void ExtModel()
{  
  int i,j,k,l,m;          
  double deriv_forw, deriv_back;
  
  //copying the model array into extmodel array   
  for(i=0;i<GRID_POINTS_GAMMA;i++)
    for(j=0;j<GRID_POINTS_HM;j++)
      for(k=0;k<GRID_POINTS_TEMP;k++)
	for(l=0;l<GRID_POINTS_TAU;l++) 
	  for(m=0;m<DATA_BIN;m++)
	     modelext[i][j][k][l][m]=model[i][j][k][l][m];	          
   
    for(i=0;i<GRID_POINTS_GAMMA;i++)//loop for gammao values      
      for(k=0;k<GRID_POINTS_TEMP;k++)//loop for logto values      
	for(l=0;l<GRID_POINTS_TAU;l++)//loop for tau_eff values      
	   for(m=0;m<DATA_BIN;m++)//data grid loop
	      {	  
		deriv_forw=(model[i][GRID_POINTS_HM-2][k][l][m] - model[i][GRID_POINTS_HM-1][k][l][m])/(Uo[GRID_POINTS_HM-2]-Uo[GRID_POINTS_HM-1]);		
		//deriv_back=(model[i][1][k][l][m] - model[i][0][k][l][m])/(Uo[1]-Uo[0]);
		
		modelext[i][GRID_POINTS_HM_EXT-1][k][l][m] = model[i][GRID_POINTS_HM-1][k][l][m] + 
		                                             ( deriv_forw * (UoExt[GRID_POINTS_HM_EXT-1]-Uo[GRID_POINTS_HM-1]) );
		//modelext[i][0][k][l][m] = model[i][0][k][l][m] + 
		//                                             ( deriv_back * (UoExt[0]-Uo[0]) );
		
	      }  
      
}


//Standard normal distribution is with Mean = 0, and Variance = 1;
//We will used Box-Muller method to generate two random varaiates from standard normal distribution
//X=mean+variance*Z  where Z is standard normal distribution
double BMVariate(double U, double V, double mu, double sigma)
{return  sqrt(-2.0*log(U))*cos(2.0*PI*V)*sigma  +  mu;}



double Mean(double array[], int size) 
{
    double value=0;
    int i;
    for(i=0; i<size; i++)
    {
        value+=array[i];
    }
    return((double)value/size);
}

double Variance (double array[], int size)
{
  int i;
  double mean_data, var_data=0;  
  mean_data = Mean(array,size);  
   
  for(i=0; i<size; i++)   
  {        
    var_data+=pow( (array[i]-mean_data),2 );    
  }  
  
  return (double)var_data/(size-1);  
  
}
 
 

double Median(double array[], int size) 
{
    double temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<size; i++) 
    {
        for(j=i+1; j<size; j++) 
	{
            if(array[j] < array[i]) 
	    {
                // swap elements
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }
    }
 
    if(size%2==0) 
    {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((array[size/2] + array[size/2 - 1]) / 2.0);
    } 
    
    else 
    {
        // else return the element in the middle
        return array[size/2];
    }
}


void GRConv (double chains[][N], int no_of_chains)
{
  int i,j;    
  double temp[DIS_N];
  double avg_mean_chains;
  double chain2chain_var;
  double avg_var_chains;
  double R, over_sigma2,  V;  
  double *mean_chains = new double[no_of_chains];
  double *var_chains  = new double[no_of_chains];    
  
  for(i=0; i<no_of_chains; i++)       
  {      
    for(j=0; j<DIS_N; j++)          
    {temp[j] = chains[i][j+(N-DIS_N-1)];}       
    
    //mean of each chain            
      mean_chains[i]=Mean(temp,DIS_N);      
      var_chains[i]=Variance(temp,DIS_N);
      printf("mean = %f\tvariance = %f\n",mean_chains[i],var_chains[i]);         
  }    
  
    avg_mean_chains=Mean(mean_chains,no_of_chains);    
    avg_var_chains =Mean(var_chains,no_of_chains);    
    chain2chain_var= Variance(mean_chains,no_of_chains);    
    
    V =  (1 - 1/(double)(DIS_N)) * avg_var_chains + (1+1/(double)no_of_chains)*chain2chain_var;
    R =  V/avg_var_chains;
    
    printf("=========================================\n");    
    printf("W   = %g\n",avg_var_chains);      
    printf("V   = %g\n",V);      
    printf("R   = %g\n",R);     
}

