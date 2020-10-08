#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcd.h"

#define LL 500

int main(int argc, char *argv[])
{
	float *x,*y,*z,*x1,*y1,*z1,lx,ly,lz,timestep,q,time;
	float *avfc,dx,dy,dz,pi,*cmx,*cmy,*cmz,dt;
	double *fcx,*fcy,*fcz;
	int N,NN,i,wcell,flag,nset,tbsave,j,t,m,tmax,whichfile,jj;
	int *times,ntimes,itemp,ilast,*llist,offset,nts,ii;
	int nfiles,*olist,*res,fileconfigs,N1;
	long int pos=-1,n,gt,*norm,start;
	char file[LL],filename[500][LL],line[LL];
	FILE *input,
             *fs;
	long int read_dcd_head(FILE *input, int *N,int flag);
	int gdcdp(float *x, float *y, float *z, char file[],long int n,int flag,long int *pos, int wcell);
	long int dcd_info(char file[], int *N,int *nset, int *tbsave, float *timestep, int *wcell);
	float cm(float *x,int N);
	
        fs = fopen("fsz.dat","w");
	pi = 4*atan(1);
	lx = 28.2;
	ly = 28.2;
	lz = 28.2;
	q = 6.1;
	offset = 0;
	nfiles = 1;
	fileconfigs = 2000;
	nts = 50;
	NN = 1000;
	N1 = 800;
	start = 0;	
	if(argc > 1) nfiles = atoi(argv[1]);
	if(argc > 2) start = atoi(argv[2]);
	
	for(i=0;i<nfiles;i++){
		sprintf(filename[i],"traj%i.dcd",i+1);
	}
	fprintf(stderr,"%s\n",filename[0]);
	nset = 0;
	for(i=0;i<nfiles;i++){
		dcd_info(filename[i],&N,&j,&tbsave,&timestep,&wcell);
		if(i==0) fileconfigs = j;
		if(i>0){
			if(j != fileconfigs){
				fprintf(stderr,"Not the same number of configurations in each file.");
				exit(1);
			}
		}
		nset += j;
	}
	dt = ((float) tbsave);
	tmax = nset;
	NN = N;
	N1 = ((int) (((double) N)*0.5));
	fprintf(stderr,"%i %i %i %i\n",nset,fileconfigs,NN,N1);
	start = fileconfigs*start;
	
	itemp = 1;
	ntimes = 0;
	j = 1;
	while(itemp < tmax){
		for(i=0;i<nts;i++){
			itemp += j;
			ntimes++;
			if(itemp >= tmax) break;
		}
		j = 5*j;
	}
	
	times = (int *) malloc(ntimes*sizeof(int));
	norm = (long int *) malloc(ntimes*sizeof(long int));
	fcx = (double *) malloc(ntimes*sizeof(double));
	fcy = (double *) malloc(ntimes*sizeof(double));
	fcz = (double *) malloc(ntimes*sizeof(double));
	x = (float *) malloc(N*sizeof(float));
	y = (float *) malloc(N*sizeof(float));
	z = (float *) malloc(N*sizeof(float));
	x1 = (float *) malloc(N*sizeof(float));
	y1 = (float *) malloc(N*sizeof(float));
	z1 = (float *) malloc(N*sizeof(float));
	cmx = (float *) malloc(nset*sizeof(float));
	cmy = (float *) malloc(nset*sizeof(float));
	cmz = (float *) malloc(nset*sizeof(float));
	
	itemp = 1;
	ntimes = 0;
	j = 1;
	times[ntimes] = 0;
	ntimes++;
	while(itemp < tmax){
		for(i=0;i<nts;i++){
			times[ntimes] = times[ntimes-1] + j;
			itemp = times[ntimes];
			ntimes++;
			if(itemp >= tmax) break;
		}
		j = 5*j;
	}
	ntimes--;
	
	flag = 0;
	
	for(n=0;n<nset;n++){
		cmx[n] = cmy[n] = cmz[n] = 0;
	}
	
	for(i=0;i<ntimes;i++){
		fcx[i] = 0;
		fcy[i] = 0;
		fcz[i] = 0;
		norm[i] = 0;
	}
	
	for(n=start;n<nset;n+=10){
		whichfile = n/fileconfigs;
      		gt = n-fileconfigs*whichfile;
		gdcdp(x,y,z,filename[whichfile],gt,flag,&pos,wcell);
		if(fabs(cmx[n]) < 1e-5){
			cmx[n] = cm(x,N);
			cmy[n] = cm(y,N);
			cmz[n] = cm(z,N);
		}
		for(j=0;j<ntimes;j++){
			t = times[j] + n;
			if(t < nset){
				whichfile = t/fileconfigs;
				gt = t-fileconfigs*whichfile;
				gdcdp(x1,y1,z1,filename[whichfile],gt,flag,&pos,wcell);
				if(fabs(cmx[t]) < 1e-5){
					cmx[t] = cm(x1,N);
					cmy[t] = cm(y1,N);
					cmz[t] = cm(z1,N);
				}
				for(i=0;i<N;i++){
					//i = llist[ii];
					dx = x[i] - x1[i] - cmx[n] + cmx[t];
					dy = y[i] - y1[i] - cmy[n] + cmy[t];
					dz = z[i] - z1[i] - cmz[n] + cmz[t];
					fcx[j] += cos(q*dx);
					fcy[j] += cos(q*dy);
					fcz[j] += cos(q*dz);
				}
				norm[j]++;
			}
		}
		fprintf(stderr,"%lu\n",n);
	}
	
	printf("#q = %f\n",q);
	for(i=0;i<ntimes;i++){
		//fc[i] /= ((float) norm[i]);
		fcx[i] /= ((float) N);
		fcy[i] /= ((float) N);
		fcz[i] /= ((float) N);
		time = times[i]*dt;
		fprintf(fs,"%f   %f \n",time,fcz[i]/((float) norm[i]));//me
		printf("%f %f %f %f %lu %i\n",time,fcx[i]/((float) norm[i]),fcy[i]/((float) norm[i]),fcz[i]/((float) norm[i]),norm[i],times[i]);
	}
}

float cm(float *x, int N)
{
	int i;
	double cm;
	
	cm = 0;
	for(i=0;i<N;i++){
		cm += x[i];
	}
	
	cm /= N;
	
	return cm;
}
