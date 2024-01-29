# four-params-mcmc
Provides a c++ code to do a fast MCMC search with 4 params quickly. Performs trilinear interpolation saves the chains and check for convergence. However, user has to provide the grid and grid of models as files, which is a 1D line of sight power spectrum.  


# MCMC part

```c++

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

```

# Trilinear interpolation

``` c++

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


```

# Checking convergence of MCMC changing - using mean and variance of chunks

```c++

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

```


# Results

![Screenshot 2024-01-29 at 14 55 05](https://github.com/nicenustian/four-params-mcmc/assets/111900566/b2416caf-d254-4ab8-a666-7c4e0396a52c)


# Citation

```
@ARTICLE{2016MNRAS.463.2335N,
       author = {{Nasir}, Fahad and {Bolton}, James S. and {Becker}, George D.},
        title = "{Inferring the IGM thermal history during reionization with the Lyman {\ensuremath{\alpha}} forest power spectrum at redshift z ≃ 5}",
      journal = {\mnras},
     keywords = {methods: numerical, intergalactic medium, quasars: absorption lines, dark ages, reionization, first stars, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2016,
        month = dec,
       volume = {463},
       number = {3},
        pages = {2335-2347},
          doi = {10.1093/mnras/stw2147},
archivePrefix = {arXiv},
       eprint = {1605.04155},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2335N},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```



