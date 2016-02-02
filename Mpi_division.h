




											
	int Xdv = (X_out)/nproc;				    
	int Xr = (X_out)-Xdv*nproc;			    
											
	for (int i=0; i < nproc; i++) {         
											
		if (i<Xr) {                         
											
			gstart[i] = i*(Xdv+1);			
			gend0[i] = gstart[i]+Xdv;		
											
		}									
		else {                              
											
			gstart[i] = i*Xdv+Xr;		
			gend0[i] = gstart[i]+Xdv-1;		
											
		}									
											
		gcount[i] = gend0[i]-gstart[i]+1;	

		gend[i] = gcount[i]+2;	

		//printf("%d\t%d\t%d\t%d\n",gstart[i],gcount[i],gend0[i],gend[i]);
											
	}										


	
