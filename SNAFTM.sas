/*code to simulate data according to a SNAFTM*/

options nomprint nosymbolgen nomlogic nonotes;


%macro simulate(
subjects=2500,            /* number of subjects in the simulated data sets */
    tpoints= 10,            /* maximum number of measurements for each subject in the simulated data */
    
           
   psi1= 0.3,      /* true psi*/
	
  gamma0 = log(3/7), gamma1 = 2, gamma2 = log(0.5), gamma3 = log(1.5),
                              /* gamma =parameters for L*/
	
   alpha0= log(2/7) , alpha1=0.5 , alpha2 = 0.5, alpha3 = log(4),                      
                              /* alpha = parameters for A*/
	

	cut=30,

	lam=0.01,/*scale parameter for weibull*/

	shape=1,/*shape parameter for weibull*/

	dataout= Simdata_all
);


     data &dataout;
       length  A Am1 L    Am1L  censor1 censor2 Y Ym eligible  maxT T 3 ;

         

        	 sT0 = %sysevalf(&sT0) ;
             sL = %sysevalf(&sL) ;
             sA = %sysevalf(&sA) ;
             sY = %sysevalf(&sY);
                 
        
          
                              
               do id=1 to %sysevalf(&subjects); 
					call ranexp(sT0,expvar); 
                    T0=(expvar/(&lam**&shape))**(1/&shape); 

					
					if T0<&cut then IT0=1;
					else IT0=0;    
					
					maxT=&tpoints;
                    initA = .;
                    C1time = .; 
                    C2time = . ;
                    array Larray L1-L&tpoints;
                    array Aarray A1-A&tpoints;
                    array Yarray Y1-Y&tpoints;
					array Ymarray Ym1-Ym&tpoints;
                    array pAarray pA1-pA&tpoints;
                    array pLarray pL1-pL&tpoints;
					array chkarray chk1-chk&tpoints;
					array tarray t1-t&tpoints;
                    
                   
                    
                    array Censearray Cense1-Cense&tpoints;
                    array Censarray  Cens1-Cens&tpoints;

                    array pcensarray pCens1-pCens&tpoints;
                    array pcensearray pCense1-pCense&tpoints;

                    time = 1.0;
                    j = 1;

					Ym1=0;
					
                    logitL = &gamma0 + &gamma1 * IT0 ;
					
					
                    pL = 1/(1+(1/exp(logitL)));
                    pLarray(1) = pL;
                    
                    if pL=1 then do; 
                         L1 = 1;
                    end;
                    else do;                      
                         call ranbin(sL,1,pL,L1);                 
                    end;                             
                                        
                    logitA = &alpha0 + &alpha1 * L1 ;                    
                    pA1 = 1/(1+(1/exp(logitA)));
					
                    
                                    
                         call ranbin(sA,1,pA1,A1);                      
                  
                
                    
					/*generate Y's*/
					chk1 = exp(&psi1*A1);

					if T0 > chk1 then do;
						Y1=0;
						t1=.;
					end;
					else if T0 <= chk1 then do;
						Y1=1;
						t1=T0*exp(-&psi1*A1);
					end;
					
      
                    cens_1 = 0; 
                    cense = 0;
                    j = 1;                    
                   cens = 0; ;
                    pcens = 0;
                    Censearray(1) = cense;
                    Censarray(1) = cens; 
                    cens_1 = cens;
                    pcensarray(1) = pCens; 
                   
           
           
                    do j=2 to &tpoints;
                         time = j;
                         if j = C2time + 1 then leave ; 
                         if j = C1time + 1 then leave;  
						

						 Ymarray(j)=Yarray(j-1);
						 

						
                         logitL = &gamma0 + &gamma1 * IT0 + &gamma2 * Aarray(j-1) + &gamma3 * Larray(j-1);
						
						pL=exp(logitL)/(1+exp(logitL));
                         pLarray(j) = pL;
                                          
                         if pL = 1 then do; 
                             Larray(j) = 1;
                         end;
                         else do;                          
                              call ranbin(sL,1,pL,Larray(j));                            
                         end;
                     
 						logitA = &alpha0 + &alpha1 * Larray(j) + &alpha2 * Larray(j-1) + &alpha3 * Aarray(j-1);
                        
                         pAarray(j)=exp(logitA)/(1+exp(logitA));
						
                         pp= pAarray(j);                                                
                                                     
                              call ranbin(sA,1,pAarray(j),Aarray(j));                                 
                    
                         if Aarray(j)=1 then timeA=j;
                   
						 if Ymarray(j)=1 then do;
							Aarray(j)=0;
							Larray(j)=0;
						 end;
						
						/*generate Y's*/
						chkarray(j) = chkarray(j-1) + exp(&psi1*Aarray(j));

						if T0 > chkarray(j) then do;
							Yarray(j)=0;
							tarray(j)=.;
							
						end;
						else if T0 <= chkarray(j) then do;
							Yarray(j)=1;
							if Ymarray(j)=1 then do;
							tarray(j)=tarray(j-1);
							end;
							else do;
							tarray(j)=(j-1) + ((T0-chkarray(j-1))*exp(-&psi1*Aarray(j)));
							end;
							
						end;
						
                        cens = 0; ;
                        pcens = 0;
                         
                         Censearray(j) = cense; 
                         Censarray(j) = cens; 
                         pCensarray(j ) = pCens;
                         pCensearray(j) = pCense ;
                         cens_1 = cens;

					
                    end; /* end j = 2 to &tpoints */
                    if C1time = . then C1time = &tpoints ; 
                    if C2time = . then C2time = &tpoints ;
                     
                    time_used = min(C1time,C2time);
                    do tpoint=1 to  time_used; 
						  
                         tpoint2 = tpoint-1;                         
                         L=Larray(tpoint);
						 chk=chkarray(tpoint); 
                                             pL = pLarray(tpoint);                          
                                             A=Aarray(tpoint);
                         
                         _pA_t=pAarray(tpoint);                                                      
                         if tpoint=1 then Am1=0; else Am1=Aarray(tpoint-1);
                         if tpoint=1 then Lm1=0; else Lm1=Larray(tpoint-1);
                         Am1L = Am1*L;
                         
                    
                         censor1 = Censarray(tpoint);
                         pC_t = 1 - pCensarray(tpoint) ; /* want p(cens = 0 ) = 1 - p(cens = 1) */

                         censor2 = Censearray(tpoint);
                         pC2_t = 1 - pCensearray(tpoint);

                         
						 	  Ym = Ymarray(tpoint);
                              if time_used < &tpoints then do;
                                   if tpoint < time_used then Y = Yarray(tpoint);
                                   else Y = . ;
                              end;
                              else if time_used = &tpoints then Y = Yarray(tpoint) ;
							  T=tarray(maxT);
                    
                             
                        
                         output;

                    end; /* end of tpoint loop */
               end;  /* end of id loop */
          

          keep id  tpoint  tpoint2 T0 IT0 chk Y Ym  A Am1 L Lm1 _pA_t T maxT 
               pL  ;

		
               call symput('sT0',sT0);
               call symput('sA',sA);
               call symput('sL',sL);
              
               call symput('sY',sY);
               
        
 
           run;


			

		   data &dataout;
				set &dataout;
				if Ym=0;
			run;

  			



%mend simulate;



         %let sT0 = 123;
          %let sL =  1217607105;
          %let sA =  1920525422;
          %let sY = 2;
          


          %simulate();

