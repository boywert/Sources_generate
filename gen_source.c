#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "L-Galaxies.h"

float cal_H(float z, float H0, float Om) {
  return H0*sqrt(Om*(1.+z)*(1.+z)*(1.+z) + (1.-Om));
}
double delta_t(float z_max, float z_min,float Om, float H0, float h) {
  int i,n = 100000;
  double dt,dz,F1,F2,z1,z2;
  dt = 0.;
  dz = (z_max - z_min)/n;
  for(i=0;i<n;i++) {
    z1 = z_min + dz*i;
    z2 = z1+dz;
    F1 = 1./(1.+z1)/cal_H(z1,h*H0,Om);
    F2 = 1./(1.+z2)/cal_H(z2,h*H0,Om);
    dt += 0.5*(F1+F2)*dz;
  }
  return dt;
}

int main(int argc, char **argv)
{
  struct LGalaxy *lgal;
  float write_buff;
  FILE *fp,*sp;
  float *sfr_float;
  double *Sfr;
  int cubep3m_p = 1728;
  float boxsize = 47.0;
  int grid = 306;
  float gridsize = boxsize/grid;
  int cubep3m_cell = cubep3m_p*2;
  long long cell;
  int firstfile = 0;
  int lastfile = 127;
  char basename[2000];
  char outputfolder[2000];
  char zlistfile[2000];
  char outputname[2048];
  char zlist_string[100][1024];
  char filename[2048];
  char buff[1000];
  int nTrees,nGals;
  int dummy,*dummyarray;
  int i,j,k;
  int nSnaps,selected_snap;
  double gridmass;
  double total_cell,total_p;
  const float Mpc2m = 3.08567758e22;
  float m2Mpc = 1./Mpc2m;
  const float Msun2kg = 1.98855e30;
  float kg2Msun = 1./Msun2kg;
  const float year2sec = 3600.*24*365.25; 
  const float m2km = 0.001;
  const float omegam = 0.27;
  float G = 6.674e-11;   // SI
  float h = 0.7;
  float H0 = 100.0;        // km/s / (Mpc/h)
  float pi = 4.0*atan(1.0);
  double rho_crit_0;
  double gridmass_c; //to convert msun to gridmass
  double dt,z1,z2;
  if(argc == 6) {
    sscanf(argv[1],"%d",&selected_snap);
    sscanf(argv[2],"%s",basename);
    sscanf(argv[3],"%s",outputfolder);
    sscanf(argv[4],"%s",zlistfile);
    sscanf(argv[5],"%d",&grid);
  }
  else {
    printf("./gen_source <snap> <basename> <output> <zlist> <ngrid>\nExit\n");
    exit(1);
  }
  printf("argc = %d , argv[0] = %s, argv[1] = %s\n",argc,argv[0],argv[1]);
  sprintf(buff,"mkdir -p %s",outputfolder);
  system(buff);

  G = G*(m2km*m2km) * (m2Mpc) / (kg2Msun); //  (Mpc/h) (km/s)^2 / (Msun/h)
  rho_crit_0 = 3.* (H0*H0)/ (8.*pi*G); //  # (Msun/h)/(Mpc/h)^3
  printf("rho_crit_0 = %g\n",rho_crit_0);
  total_p = (double)cubep3m_p* (double)cubep3m_p* (double)cubep3m_p;
  total_cell = (double)cubep3m_cell*(double)cubep3m_cell*(double)cubep3m_cell;
  printf("total cell = %g\n",total_cell);
  gridmass = (double)omegam*(double)rho_crit_0*(double)(boxsize*boxsize*boxsize)/total_cell; // /(double)h; // Msun
  printf("G = %g, gridmass = %lf\n",G,gridmass);
  gridmass_c = 1./gridmass;
  printf("gridmass_c = %lg\n",gridmass_c);
  fp = fopen(zlistfile,"r");
  i=0;
  while (fscanf(fp, "%s", zlist_string[i]) != EOF) {
    i++;
  }
  nSnaps = i;
  fclose(fp);
  printf("Total snapshot : %d\n",nSnaps);
  for(j=0;j<nSnaps;j++) {
    if(j == selected_snap) {
      printf("Converting snapshot : %d\n",selected_snap);
      Sfr = calloc(grid*grid*grid,sizeof(double));
      if(j > 0) {
	sscanf(zlist_string[j],"%lg",&z1);
	sscanf(zlist_string[j-1],"%lg",&z2);
	dt = delta_t((float)z2,(float)z1,(float)omegam, H0, h);
      } else
	dt = 0.;
      dt *= (Mpc2m/h/1000.)/year2sec;
      printf("dt = %lg\n",dt);
	// km/s / (Mpc/h)
      // read the previous snapshot to make cumulative
      /* if(j > 0) { */
      /* 	sprintf(outputname,"%s/%s.dat",outputfolder,zlist_string[j-1]); */
      /* 	sfr_float = malloc(grid*grid*grid*sizeof(float)); */
      /* 	sp = fopen(outputname,"rb"); */
      /* 	fread(sfr_float,sizeof(float),grid*grid*grid,sp); */
      /* 	for(k=0;k<grid*grid*grid;k++) */
      /* 	  Sfr[k] += (double)sfr_float[k]; */
      /* 	fclose(sp); */
      /* 	free(sfr_float); */
      /* } */
      for (i=firstfile;i<=lastfile;i++) {
	sprintf(filename, "%s%s_%d",basename,zlist_string[j],i);
	if(i == firstfile || i == lastfile)
	  printf("Reading %s\n.....\n",filename);
	fp = fopen(filename,"rb");
	fread(&nTrees, sizeof(int), 1, fp);
	fread(&nGals, sizeof(int),1, fp);
	lgal = malloc(sizeof(struct LGalaxy)*nGals);
	fseek(fp, nTrees*sizeof(int), SEEK_CUR); // skip nGalsperTree
	fread(lgal,sizeof(struct LGalaxy),nGals,fp);
	for(k=0;k<nGals;k++) {
	  cell = (int)(lgal[k].Pos[0]/gridsize) + (int)(lgal[k].Pos[1]/gridsize)*grid + (int)(lgal[k].Pos[2]/gridsize)*grid*grid;
	  Sfr[cell] += (double)(lgal[k].Sfr*gridmass_c*dt);
	}
	free(lgal);
	fclose(fp);
      }
      sprintf(outputname,"%s/%s.dat",outputfolder,zlist_string[j]);
      printf("output to file: %s\n",outputname);
      fp = fopen(outputname,"wb+");
      fwrite(&grid,sizeof(int),1,fp);
      fwrite(&grid,sizeof(int),1,fp);
      fwrite(&grid,sizeof(int),1,fp);
      for(i=0;i<grid*grid*grid;i++) {
	write_buff = (float)Sfr[i];
	fwrite(&write_buff,sizeof(float),1,fp);
      }
      fclose(fp);
      printf("Finish converting snap :%d\n",selected_snap);
      free(Sfr);	    
    }
  }
  return 0;
}
