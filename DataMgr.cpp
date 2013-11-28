#include "DataMgr.h"
#include <stdlib.h>

#include <stdio.h>


#if 0
d_type ***alloc3d(size_t dim1,size_t  dim2,size_t  dim3,size_t nsize)     /* allocate and return a pointer to a 3 dimensional array */
{
   size_t i,j;            /* loop counters */
   d_type	*p,             /* pointer to beginning of allocated space */
	*space,         /* pointer to beginning of allocated space */
	***begin;       /* pointer to 3 dimensional array */

   if((space = (d_type *)malloc(dim1*dim2*dim3*nsize)) == NULL){
      printf("alloc3d: can't allocate space \n");
      return(NULL);
   }
   if((begin = (d_type ***)malloc(dim1*sizeof(d_type **))) == NULL) {
      printf("alloc3d: can't allocate pointers at dim1\n");
      return(NULL);
   }
   for(i = 0; i < dim1; i++) {
      if((begin[i] = (d_type **)malloc(dim2*sizeof(d_type *))) == NULL) {
         printf("alloc3d: can't allocate pointers at dim2\n");
         return(NULL);
      }
   }

   p = space;
   for(i = 0; i < dim1; i++) {
      for(j = 0; j < dim2; j++) {
         begin[i][j] = p + (i * dim2 * dim3 * nsize) + (j*dim3*nsize);
      }
   }
   return(begin);
}


template <typename T>
int DataMgr::readBinary(const char* filename)
{
	FILE *fdata;
//	fdata = fopen("data/hipip","rb");
	fdata = fopen("data/uncBrain.vol","rb");
	if(NULL == fdata)
	{
		printf("%s","Failed to open data!");
		system("pause");
		return 0;
	}

	fscanf(fdata,"%d %d %d",&nx,&ny,&nz);					//read data head "#X,#Y,#Z"
	printf("%d %d %d\n",nx,ny,nz);
	getc(fdata);											// skip "LF"(New Line)

	f = (d_type ***)alloc3d(nx,ny,nz,sizeof(d_type));	//allocate space for a 3D array "f"
	int data_num = nx*ny*nz;
	fread(&f[0][0][0],sizeof(d_type),data_num,fdata);			//read data from file to f
	//getc(fdata);
	//if(!feof(fdata)){
	//	printf("End of File does not reach");
	//}
	for(int i=0;i<data_num;i++)
	{
		d_type tmp = *(&f[0][0][0]+i);
		if(tmp<d_min)
			d_min = tmp;
		if(tmp>d_max)
			d_max = tmp;
	}
}



template <typename T>
int DataMgr::readASCII(const char* filename)
{
	FILE *fdata;
	fdata = fopen(filename,"r");
	if(NULL == fdata)
	{
		printf("%s","Failed to open data!");
		system("pause");
		return 0;
	}

	fscanf(fdata,"%d %d %d",&nx,&ny,&nz);					//read data head "#X,#Y,#Z"
	printf("%d %d %d\n",nx,ny,nz);
	getc(fdata);											// skip "LF"(New Line)

	f = (d_type ***)alloc3d(nx,ny,nz,sizeof(d_type));	//allocate space for a 3D array "f"
	int data_num = nx*ny*nz;
	for(int i = 0;i<nx;i++)
		for(int j = 0;j<ny;j++)
			for(int k = 0; k<nz;k++)
				fscanf (fdata, "%i", &f[i][j][k]);

//	fread(&f[0][0][0],sizeof(d_type),data_num,fdata);			//read data from file to f
	//getc(fdata);
	//if(!feof(fdata)){
	//	printf("End of File does not reach");
	//}
	for(int i=0;i<data_num;i++)
	{
		d_type tmp = *(&f[0][0][0]+i);
		if(tmp<d_min)
			d_min = tmp;
		if(tmp>d_max)
			d_max = tmp;
	}
}

template <typename T>
void DataMgr::SetDataRange(d_type min, d_type max)
{
	d_min = min;
	d_max = max;
}


template <typename T>
void DataMgr::printAsciiToFile()
{
	FILE *fdataASCII;
	fdataASCII = fopen("data/myASCII","w");
	for(int i = 0;i<nx;i++){
		fprintf(fdataASCII,"**************slice x = %i****************\n",i+1);
		for(int j = 0; j < ny; j++){
			for(int k = 0; k < nz; k++){
				//printf("%i ",f[i][j][k]);
				fprintf(fdataASCII,"%i ",f[i][j][k]);			//print ascii data to a file
			}
			//printf("\n");
			fprintf(fdataASCII,"\n");
		}
		//printf("\n");
		fprintf(fdataASCII,"\n");
	}
	fclose(fdataASCII);
}

#endif