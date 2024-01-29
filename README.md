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



