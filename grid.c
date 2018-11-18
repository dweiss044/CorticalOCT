#include "fdtd-macro.h"
#include "fdtd-alloc1.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// All values listed below @1550 nm wavelength

#define N_DURA 1.4022 //Dura refractive index
#define N_CSF 1.335 //Cerebral Spinal Fluid refractive index
#define N_MYELIN 1.4 //myelin refractive index
#define EPSR_WATER 1.7509 //dielectric constant of water
#define X_DURA 1000 //Dura starting point
#define X_CORTEX 2000 //Cortex starting point
#define SIZE 350 //number of nerve fibres (rows in filegrid.txt)

void gridInit(Grid *g) {
  double imp0=377.0, epsr_dura, epsr_fibres, epsr_csf, epsr_myelin, n_fibres;
	float temp;
	int mm, nn, i=0, j=0, r[SIZE], r1, xLocation, yLocation, cx[SIZE], cy[SIZE],l[SIZE], temp2;
	FILE *read, *out;
	char worg[1], activity[1], dielectric[20];

  Type    = tmZGrid;
  SizeX   = 13000; // size of domain in x with dx = 100nm, simulation domain = 1.3 mm
  SizeY   = 2000;	//size of domain in y  with dy = 100nm, simulation domain = .2 mm
  MaxTime = 60000; // duration of simulation
  Cdtds   = 1.0/sqrt(2.0); // Courant number

	//Allocate memory
  ALLOC_2D(g->hx,  SizeX,  SizeY-1, double);
  ALLOC_2D(g->chxh,SizeX,  SizeY-1, double);
  ALLOC_2D(g->chxe,SizeX,  SizeY-1, double);
  ALLOC_2D(g->hy,  SizeX-1,SizeY, double);
  ALLOC_2D(g->chyh,SizeX-1,SizeY, double);
  ALLOC_2D(g->chye,SizeX-1,SizeY, double);
  ALLOC_2D(g->ez,  SizeX,  SizeY, double);
  ALLOC_2D(g->ceze,SizeX,  SizeY, double);
  ALLOC_2D(g->cezh,SizeX,  SizeY, double);

	printf("Sample (s) or reference (r) arm? ");
  scanf(" %s", worg);
	if(strcmp(worg, "s") == 0){
		printf("You have chosen to simulate the sample arm.\n");
		sprintf(dielectric,"dielectric_%1s",worg);
		printf("Is the nerve active? (y/n)\n");
 		scanf(" %s", activity);
		if(strcmp(activity, "y") == 0)
			n_fibres = 1.33798662;
		else if(strcmp(activity, "n") == 0)
			n_fibres = 1.338;
		else{
			printf("Invalid character.\n");
			exit(-1);
		}
		//Compute dielectric constants
		epsr_dura = N_DURA*N_DURA;
		epsr_fibres = n_fibres*n_fibres;
		epsr_csf = N_CSF*N_CSF;
		epsr_myelin = N_MYELIN*N_MYELIN;
		
	
		//read file containing the positions of the nerve fibres -- center positions (x,y) and radii
		read = fopen("filegridCortex.txt","r");
		for(i = 0; i<SIZE; i++){
		  // UNCOMMENT if seg faulting at initGrid
		  char str[12];
		  sprintf(str, "%d", i);
		  fscanf(read, "%d %d %d %d",&cx[i],&cy[i],&r[i],&l[i]);
		}
		fclose(read);
	  
		out = fopen(dielectric, "wb");
 		/* set electric-field update coefficients */
 		for (mm=0; mm<SizeX; mm++){
  	  for (nn=0; nn<SizeY; nn++) {
				if( mm < X_DURA){ //water
					Ceze(mm, nn) = 1.0;
					Cezh(mm, nn) = Cdtds*imp0 /  EPSR_WATER;
				}else if(mm >= X_DURA && mm < X_CORTEX){ // Dura Mater TODO: make slightly more fine structure
  					Ceze(mm, nn) = 1.0;
		  			Cezh(mm, nn) = Cdtds*imp0 /epsr_dura;
				}else{
					for(j = 0; j<SIZE ;j++){
					  //     xLocation = mm - cx[j];
					  //	yLocation = nn - cy[j];
					  // TODO: fix pyramidal cells
						if(mm > cx[j]-l[j] && mm < cx[j]){
						  if(nn > cy[j]-r[j] && nn < cy[j]){
								Ceze(mm, nn) = 1.0;
								Cezh(mm, nn) = Cdtds*imp0 /epsr_fibres;
								break;
							}
				    
						// TODO: add myelin sheath
						                Ceze(mm, nn) = 1.0;
								Cezh(mm, nn) = Cdtds*imp0 /epsr_csf;
						}
						
					
						
						/*
						if(xLocation*xLocation + yLocation*yLocation <= r[j]*r[j]){//nerve fibres and their myelin sheath
							if(xLocation*xLocation + yLocation*yLocation <= r1*r1){
								Ceze(mm, nn) = 1.0;
								Cezh(mm, nn) = Cdtds*imp0 /epsr_fibres;
								break;
							}
							Ceze(mm, nn) = 1.0;
							Cezh(mm, nn) = Cdtds*imp0 /epsr_myelin;
							break;
							}
						*/
							else{
								Ceze(mm, nn) = 1.0;
								Cezh(mm, nn) = Cdtds*imp0 /epsr_csf;
						        }
					}
					
				}
				temp = (float)Cezh(mm, nn);
				fwrite(&temp, sizeof(float), 1, out); // write the float
		
	  }
		}
		fclose(out);
	}else if(strcmp(worg, "r") == 0){
		printf("You have chosen to simulate the reference arm.\n");
  	sprintf(dielectric,"dielectric_%1s",worg);
		out = fopen(dielectric, "wb");

  /* set electric-field update coefficients */
 		for (mm=0; mm<SizeX; mm++){
   		for (nn=0; nn<SizeY; nn++) {
				if( mm < X_DURA){ //air
					Ceze(mm, nn) = 1.0;
					Cezh(mm, nn) = Cdtds*imp0;
				}else{//mirror
					Ceze(mm,nn) = 0.0;
       		Cezh(mm,nn) = 0.0;
				}
				temp = (float)Cezh(mm, nn);
				fwrite(&temp, sizeof(float), 1, out); // write the float
	   	}
		}
	fclose(out);

	}else{
		printf("Invalid character.\n");
		exit(-1);
	}

  /* set magnetic-field update coefficients */
  for (mm=0; mm<SizeX; mm++)
    for (nn=0; nn<SizeY-1; nn++) {
      Chxh(mm,nn) = 1.0;
      Chxe(mm,nn) = Cdtds/imp0;
    }

  for (mm=0; mm<SizeX-1; mm++)
    for (nn=0; nn<SizeY; nn++) {
      Chyh(mm,nn) = 1.0;
      Chye(mm,nn) = Cdtds/imp0;
    }

  return;
}
