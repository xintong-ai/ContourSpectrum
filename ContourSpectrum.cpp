#include <time.h>
#include <math.h>

#include "ContourSpectrum.h"
//#include "Vector.h"

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#define eps 1.0e-06

//using namespace Math;

//int nx_t ,ny_t,nz_t;					//tetrahedra counter
//in_type area_min, area_max;					//minimun contour area and maximum contour area
//in_type *vertices;					//vertices position and area value for VBO
//GLuint VAO;
//GLuint VBO_area;
//int sample_num;
//int non_uni_num;
//const char* filename_data = "data/uncBrain.vol";
const char* filename_data = "D:/data/combustion/raw_data/jet_0032/chiContour_OHValue.vtu";
//const char* filename_data = "data/testSpectrum.vtu";
//const char* filename_data = "data/ex-blow_5_0_3.vtu";

/*
void ContourSpectrum::ComputeContourSpectrum(
	vec3<in_type> v0, vec3<in_type> v1, vec3<in_type> v2, vec3<in_type> v3, 
	in_type value[4],
	float poly0[4], float poly1[4], float poly2[4])
{
	

	if(val>=cur->f1 && val<=cur->f4)			
	{
		if (abs(val - cur->f2)<eps)
		{
			area[i] += cur->s1;
		}
		else if(abs(val - cur->f3)<eps)
		{
			area[i] += cur->s2;
		}
		else if(val<cur->f2)
			area[i] += (pow((val - cur->f1)/(cur->f2 - cur->f1),2) * cur->s1);
		else if(val>cur->f3)
			area[i] += (pow((cur->f4 - val)/(cur->f4 - cur->f3),2) * cur->s2);
		else if(val>cur->f2 && val<cur->f3)
		{
			if(abs(cur->f1 - cur->f2)<eps)
			{
				if(abs(cur->f3 - cur->f4)<eps)			//pattern 3
				{
					in_type v_1 = val - cur->f1;
					in_type v_2 = val - cur->f2;
					in_type v_3 = val - cur->f3;
					in_type v_4 = val - cur->f4;
					vec3<in_type> I = (v_1*cur->v1 - v_4*cur->v4)/(cur->f4 - cur->f1);
					vec3<in_type> J = (v_2*cur->v2 - v_4*cur->v4)/(cur->f4 - cur->f2);
					vec3<in_type> K = (v_2*cur->v2 - v_3*cur->v3)/(cur->f3 - cur->f2);
				//	in_type aa = ;
					area[i] += cross(I-J,K-J).length();
				//	printf("aa:%f \n",aa);
				}
				else							//pattern 2
				{
					in_type t1 = pow((cur->f3 - val)/ (cur->f3 - cur->f2),2);
					in_type t2 = pow((cur->f4 - val)/ (cur->f4 - cur->f3),2);
					in_type t3 = pow((cur->f4 - cur->f2)/ (cur->f4 - cur->f3),2);
				//	area[i] += ;
				//	in_type bb = ;
					area[i] += t1 * cur->s1 + (t2 - t3 * t1) * cur->s2;
				//	printf("bb:%f \n",bb);
				}
			}
			else								//pattern 1
			{
				in_type t1 = pow((val - cur->f1)/ (cur->f2 - cur->f1),2);
				in_type t2 = pow((cur->f3 - cur->f1)/ (cur->f2 - cur->f1),2);
				in_type t3 = pow((val - cur->f2)/ (cur->f3 - cur->f2),2);
			//	in_type cc = ;
				area[i] += (t1 - t2 * t3) * cur->s1 + t3 * cur->s2;
			//		printf("cc:%f \n",cc);
			}
		}
	}
	
}

*/
/*
//sorting result: *y1<*y2<*y3<*y4
void precompute(int i,int j,int k,int x_dec_1,int x_dec_2,int x_dec_3,int x_dec_4,tet *ptet)
{
	int *min1,*min2,*max1,*max2,*mid1,*mid2;
	int *y1,*y2,*y3,*y4;					//pointing to sorted vertex index
	int uf1,uf2,uf3,uf4;						//f(i) is the values of the unsorted vertices
	d_type *pmin1,*pmin2,*pmax1,*pmax2,*pmid1,*pmid2;	//intermediate variables for sorting

	//binary: 0 or 1
	int xb1[3];
	int xb2[3];
	int xb3[3];
	int xb4[3];

	xb1[0] = ((x_dec_1&4)>0?1:0);
	xb1[1] = ((x_dec_1&2)>0?1:0);
	xb1[2] = ((x_dec_1&1)>0?1:0);
	
	xb2[0] = ((x_dec_2&4)>0?1:0);
	xb2[1] = ((x_dec_2&2)>0?1:0);
	xb2[2] = ((x_dec_2&1)>0?1:0);
	
	xb3[0] = ((x_dec_3&4)>0?1:0);
	xb3[1] = ((x_dec_3&2)>0?1:0);
	xb3[2] = ((x_dec_3&1)>0?1:0);
	
	xb4[0] = ((x_dec_4&4)>0?1:0);
	xb4[1] = ((x_dec_4&2)>0?1:0);
	xb4[2] = ((x_dec_4&1)>0?1:0);
	
	//full value of the index (no order)
	int x1[3];
	int x2[3];
	int x3[3];
	int x4[3];

	x1[0] = i + xb1[0];
	x1[1] = j + xb1[1];
	x1[2] = k + xb1[2];
	
	x2[0] = i + xb2[0];
	x2[1] = j + xb2[1];
	x2[2] = k + xb2[2];
	
	x3[0] = i + xb3[0];
	x3[1] = j + xb3[1];
	x3[2] = k + xb3[2];
	
	x4[0] = i + xb4[0];
	x4[1] = j + xb4[1];
	x4[2] = k + xb4[2];

	//uf(i) is the values of the unsorted vertices
	uf1 = f[x1[0]][x1[1]][x1[2]];
	uf2 = f[x2[0]][x2[1]][x2[2]];
	uf3 = f[x3[0]][x3[1]][x3[2]];
	uf4 = f[x4[0]][x4[1]][x4[2]];

	//if the value of four vertices are equal, then no more computation needed, since no isocontour pass through.
	//AD cross g == 0
	if(uf1==uf2&&uf2==uf3&&uf3==uf4)	
	{
		ptet->uniform = true;
		non_uni_num ++;
		return;
	}
	ptet->uniform = false;

	// sort algorithm, which I made by myself and is quicker than the other sort algorithm designed for sorting large scale data
	if(uf1<uf2){
		min1 = xb1;
		max1 = xb2;
	}
	else{
		min1 = xb2;
		max1 = xb1;
	}

	if(uf3<uf4)	{
		min2 = xb3;
		max2 = xb4;
	}
	else{
		min2 = xb4;
		max2 = xb3;
	}

	if(f[i+min1[0]][j+min1[1]][k+min1[2]]<f[i+min2[0]][j+min2[1]][k+min2[2]])	{
		y1 = min1;
		mid1 = min2;
	}
	else{
		y1 = min2;
		mid1 = min1;
	}

	if(f[i+max1[0]][j+max1[1]][k+max1[2]]>f[i+max2[0]][j+max2[1]][k+max2[2]])	{
		y4 = max1;
		mid2 = max2;
	}
	else{
		y4 = max2;
		mid2 = max1;
	}

	if(f[i+mid1[0]][j+mid1[1]][k+mid1[2]]<f[i+mid2[0]][j+mid2[1]][k+mid2[2]])	{
		y2 = mid1;
		y3 = mid2;	
	}
	else{
		y2 = mid2;
		y3 = mid1;	
	}


	//printf("before:%d %d %d %d\n",(f[x1[0]][x1[1]][x1[2]]),(f[x2[0]][x2[1]][x2[2]]),(f[x3[0]][x3[1]][x3[2]]),(f[x4[0]][x4[1]][x4[2]]));
	//printf("after:%d %d %d %d\n",(f[i+y1[0]][j+y1[1]][k+y1[2]]),(f[i+y2[0]][j+y2[1]][k+y2[2]]),(f[i+y3[0]][j+y3[1]][k+y3[2]]),(f[i+y4[0]][j+y4[1]][k+y4[2]]));

	//store the sorted vertex indices into tetrahedron instance
	for(int i=0;i<3;i++)
	{
		ptet->v1[i] = y1[i];
		ptet->v2[i] = y2[i];
		ptet->v3[i] = y3[i];
		ptet->v4[i] = y4[i];
	}

	//store the sorted values into tetrahedron instance
	ptet->f1 = f[i+y1[0]][j+y1[1]][k+y1[2]];
	ptet->f2 = f[i+y2[0]][j+y2[1]][k+y2[2]];
	ptet->f3 = f[i+y3[0]][j+y3[1]][k+y3[2]];
	ptet->f4 = f[i+y4[0]][j+y4[1]][k+y4[2]];

	//A*g = B, A is a 3*3 matrix, B is a size 3 vector
	in_type B[3];
	B[0] = ptet->f2 - ptet->f1;
	B[1] = ptet->f3 - ptet->f2;
	B[2] = ptet->f4 - ptet->f3;
//	printf("begin: \n B:%f %f %f\n",B[0],B[1],B[2]);

	vec3<in_type> y12(y2[0] - y1[0],y2[1] - y1[1],y2[2] - y1[2]);
	vec3<in_type> y13(y3[0] - y1[0],y3[1] - y1[1],y3[2] - y1[2]);
	vec3<in_type> y14(y4[0] - y1[0],y4[1] - y1[1],y4[2] - y1[2]);
	vec3<in_type> y23(y3[0] - y2[0],y3[1] - y2[1],y3[2] - y2[2]);
	vec3<in_type> y24(y4[0] - y2[0],y4[1] - y2[1],y4[2] - y2[2]);
	vec3<in_type> y34(y4[0] - y3[0],y4[1] - y3[1],y4[2] - y3[2]);

	in_type A[] = {		
	y12[0],y12[1],y12[2],
	y23[0],y23[1],y23[2],
	y34[0],y34[1],y34[2]};
//	printf("A:\n");
//	for(int i = 0;i<9;i++)
//		printf("%f ",A[i]);

	//determinent of A
	in_type det_A =	A[0]*A[4]*A[8]	+A[1]*A[5]*A[6]	+A[2]*A[3]*A[7]
						-A[2]*A[4]*A[6]	-A[0]*A[5]*A[7]	-A[1]*A[3]*A[8];	
	vec3<in_type> g;		//gradient
	g[0] = ((A[4]*A[8]-A[5]*A[7])*B[0] + (A[2]*A[7]-A[1]*A[8])*B[1] + (A[1]*A[5]-A[2]*A[4])*B[2])/det_A;
	g[1] = ((A[5]*A[6]-A[3]*A[8])*B[0] + (A[0]*A[8]-A[2]*A[6])*B[1] + (A[2]*A[3]-A[0]*A[5])*B[2])/det_A;
	g[2] = ((A[3]*A[7]-A[4]*A[6])*B[0] + (A[1]*A[6]-A[0]*A[7])*B[1] + (A[0]*A[4]-A[1]*A[3])*B[2])/det_A;
	ptet->gradient = g;				//store the gradient

	// if the smaller two and larger two have the same values, s1 and s2 should be 0. like 0 0 1 1 or 2 2 5 5
	//no other possibility that s1 and s2 can both be 0
	if(ptet->f1 == ptet->f2 &&ptet->f3 == ptet->f4)
	{
		ptet->s1 = 0;
		ptet->s2 = 0;
		return;
	}

	in_type y12g = y12.dot(g);
	in_type y13g = y13.dot(g);
	in_type y14g = y14.dot(g);
	in_type y24g = y24.dot(g);
	in_type y34g = y34.dot(g);

	//AC cross g == 0
	if(abs(y13g)<eps)		//in this case point E is point A, point F is point C
	{
		ptet->s2 = (y12.cross(y23)).length()/2;
		ptet->s1 = ptet->s2;
		return;
	}

	//BD cross g == 0
	if(abs(y24g)<eps)
	{
		ptet->s1 = (y23.cross(y24)).length()/2;
		ptet->s2 = ptet->s1;
		return;
	}

	
	if(abs(y12g)<eps)		//AB cross g == 0     if y13g = 0, then y12g = 0
	{
		ptet->s1 = 0;
	}
	else
	{
		vec3<in_type> y14_y13 = y14.cross(y13);
		vec3<in_type> y14_y12 = y14.cross(y12);
		vec3<in_type> y12_y13 = y12.cross(y13);
		ptet->s1 = (pow(y12g,2)/(y14g*y13g) * y14_y13 - y12g/y14g *y14_y12 - y12g/y13g * y12_y13).length()/2;
	}
	
	
	if(abs(y34g)<eps)		//CD cross g == 0
	{
		ptet->s2 = 0;
	}
	else
	{
		vec3<in_type> y14_y24 = y14.cross(y24);
		vec3<in_type> y14_y34 = y14.cross(y34);
		vec3<in_type> y34_y24 = y34.cross(y24);
		ptet->s2 = (pow(y34g,2)/(y14g*y24g) * y14_y24 - y34g/y14g * y14_y34 -  y34g/y24g * y34_y24).length()/2;
	}
	return;
}

void spectrum(in_type *isov,in_type *area,int n,tet *cur)
{
	////convert data type to inner data type for the following calculation
	//vec3<in_type> cur_vert_1,cur_vert_2,cur_vert_3,cur_vert_4;
	//for(int i=0;i<3;i++)
	//{
	//	cur_vert_1[i] = cur->v1[i];
	//	cur_vert_2[i] = cur->v2[i];
	//	cur_vert_3[i] = cur->v3[i];
	//	cur_vert_4[i] = cur->v4[i];
	//}
	//float f1,f2,f3,f4;
	//f1 = cur->f1;
	//f2 = cur->f2;
	//f3 = cur->f3;
	//f4 = cur->f4;

	for(int i=0;i<n;i++)	//for each isovalue whose properties are needed
	{
		in_type val = isov[i];
		//if the contour of this isovalue is between the value range of the current tetrahedron, count its area
		
		//the reason why it is >= and <= here is that we can have the possible case "0 0 0 1", where s1 = the triangle area
		if(val>=cur->f1 && val<=cur->f4)			
		{
			if (abs(val - cur->f2)<eps)
			{
				area[i] += cur->s1;
			}
			else if(abs(val - cur->f3)<eps)
			{
				area[i] += cur->s2;
			}
			else if(val<cur->f2)
				area[i] += (pow((val - cur->f1)/(cur->f2 - cur->f1),2) * cur->s1);
			else if(val>cur->f3)
				area[i] += (pow((cur->f4 - val)/(cur->f4 - cur->f3),2) * cur->s2);
			else if(val>cur->f2 && val<cur->f3)
			{
				if(abs(cur->f1 - cur->f2)<eps)
				{
					if(abs(cur->f3 - cur->f4)<eps)			//pattern 3
					{
						in_type v_1 = val - cur->f1;
						in_type v_2 = val - cur->f2;
						in_type v_3 = val - cur->f3;
						in_type v_4 = val - cur->f4;
						vec3<in_type> I = (v_1*cur->v1 - v_4*cur->v4)/(cur->f4 - cur->f1);
						vec3<in_type> J = (v_2*cur->v2 - v_4*cur->v4)/(cur->f4 - cur->f2);
						vec3<in_type> K = (v_2*cur->v2 - v_3*cur->v3)/(cur->f3 - cur->f2);
					//	in_type aa = ;
						area[i] += cross(I-J,K-J).length();
					//	printf("aa:%f \n",aa);
					}
					else							//pattern 2
					{
						in_type t1 = pow((cur->f3 - val)/ (cur->f3 - cur->f2),2);
						in_type t2 = pow((cur->f4 - val)/ (cur->f4 - cur->f3),2);
						in_type t3 = pow((cur->f4 - cur->f2)/ (cur->f4 - cur->f3),2);
					//	area[i] += ;
					//	in_type bb = ;
						area[i] += t1 * cur->s1 + (t2 - t3 * t1) * cur->s2;
					//	printf("bb:%f \n",bb);
					}
				}
				else								//pattern 1
				{
					in_type t1 = pow((val - cur->f1)/ (cur->f2 - cur->f1),2);
					in_type t2 = pow((cur->f3 - cur->f1)/ (cur->f2 - cur->f1),2);
					in_type t3 = pow((val - cur->f2)/ (cur->f3 - cur->f2),2);
				//	in_type cc = ;
					area[i] += (t1 - t2 * t3) * cur->s1 + t3 * cur->s2;
				//		printf("cc:%f \n",cc);
				}
			}
		}
	//	printf("area %i:%f \n",i,area[i]);
	}
}

void procTet(int i,int j,int k,int x_dec_1,int x_dec_2,int x_dec_3,int x_dec_4,tet *ptet,in_type *isov,in_type *area,int n)
{
	precompute(i,j,k,x_dec_1,x_dec_2,x_dec_3,x_dec_4,ptet);
	if(false == ptet->uniform)
		spectrum(isov,area,n,ptet);
}

void compute()
{
	//printf("cc:%f %f %f\n",cc[0],cc[1],cc[2]);
	//counter of tetrahedra, whose size it one smaller than the size of mesh vertices
	nx_t = nx - 1;
	ny_t = ny - 1;
	nz_t = nz - 1;

	sample_num = 100;
	in_type *isovalue =(in_type *)malloc(sample_num*sizeof(in_type *));
	d_max -= eps;
	d_min += eps;
	in_type sample_step = in_type(d_max - d_min)/(sample_num-1);

	//assign the isovalues
	for(int i=0;i<sample_num;i++)
	{
		//offset(eps) makes sure the sample isovalues are not equal to the data value, so less "share face ambiguity" would happen
		isovalue[i] = (in_type)d_min + i * sample_step; 
	}

	in_type *iso_area = (in_type *)malloc(sample_num*sizeof(in_type *));
	for(int i=0;i<sample_num;i++)	//initiate the area values to be 0 for the following accumulation
	{
		iso_area[i] = 0.0;
	}

	for(int i = 0;i<nx_t;i++){
		printf("**************decompose slice x = %i****************\n",i+1);
		for(int j = 0; j < ny_t; j++){
			for(int k = 0; k < nz_t; k++){
				//tetrahedra 
			//	fprintf(fdataASCII,"%i ",f[i][j][k]);			//print ascii data to a file
				tet cur_tet;
				if (0 == (i + j + k)%2)				//two tetrahedra decomposition pattern of cube
				{
					procTet(i,j,k,0,2,3,6,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,3,5,6,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,4,5,6,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,1,3,5,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,3,5,6,&cur_tet,isovalue,iso_area,sample_num);
				}
				else
				{
					procTet(i,j,k,0,1,2,4,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,4,5,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,2,4,6,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,2,3,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,2,4,7,&cur_tet,isovalue,iso_area,sample_num);
				}
			}
	//	free(tet_slice);
		}
	}

	vertices = (in_type *)malloc(sample_num*2*sizeof(in_type *));
	for(int i = 0;i<sample_num;i++)
	{
		vertices[i*2] = isovalue[i];
		in_type tmp = iso_area[i];
		vertices[i*2 +1] = tmp;
		if(tmp<area_min)
			area_min = tmp;
		if(tmp>area_max)
			area_max = tmp;
	}
	for(int i=0;i<sample_num;i++)
		printf("area %i: %f %f\n",i,isovalue[i],iso_area[i]);
	//printf("t:%d %d %d",t[3][38][21][2].v1[0],t[3][38][21][2].v1[1],t[3][38][21][2].v1[2]);
	//printf("t:%d %d %d",t[3][38][21][2].v2[0],t[3][38][21][2].v2[1],t[3][38][21][2].v2[2]);
	//printf("t:%d %d %d",t[3][38][21][2].v3[0],t[3][38][21][2].v3[1],t[3][38][21][2].v3[2]);
	//printf("t:%d %d %d",t[3][38][21][2].v4[0],t[3][38][21][2].v4[1],t[3][38][21][2].v4[2]);
//	printf("n_tet:%d",n_tet);

	//tet *prev = p_tet_prev;
	//int i = 0;
	//while(prev != p_tet_first)
	//{
	//	printf("n:%d s1:%f  s2:%f g:%f %f %f\n",++i,prev->s1,prev->s2,prev->gradient[0],prev->gradient[1],prev->gradient[2]);
	//	prev = prev->prev;
	//	i++;
	//	//if(abs(pprev->s1) <eps && abs(pprev->s2) <eps)
	//	//{
	//	//	printf("here");
	//	//}
	//}
	//printf("number:%d",i);
}

*/



////future work: this can be a table
//void ContourSpectrum::GenLoc(int i0, int i1, int i2, int i3,
//	vec3<d_type> v0, vec3<d_type> v1, vec3<d_type> v2, vec3<d_type> v3)
//{
//
//}

//template <typename T>
//inline void ContourSpectrum::Cube2Tetra(
//	d_type v0, d_type v1, d_type v2, d_type v3,
//	d_type v4, d_type v5, d_type v6, d_type v7,
//	tetra<T> t[5])
//{
//	//t[0]
//}
/*
void ContourSpectrum::ProcessAllCubes()
{
	_dm->allPoly = new polynomial[nx_t * ny_t * nz_t * 3];
	for(int i = 0;i<nx_t;i++){
		printf("**************decompose slice x = %i****************\n",i+1);
		for(int j = 0; j < ny_t; j++){
			for(int k = 0; k < nz_t; k++){
				//tetrahedra 
			//	fprintf(fdataASCII,"%i ",f[i][j][k]);			//print ascii data to a file
				tet cur_tet;
				if (0 == (i + j + k)%2)				//two tetrahedra decomposition pattern of cube
				{

					ComputeContourSpectrum(

					procTet(i,j,k,0,2,3,6,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,3,5,6,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,4,5,6,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,1,3,5,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,0,3,5,6,&cur_tet,isovalue,iso_area,sample_num);
				}
				else
				{
					procTet(i,j,k,0,1,2,4,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,4,5,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,2,4,6,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,2,3,7,&cur_tet,isovalue,iso_area,sample_num);

					procTet(i,j,k,1,2,4,7,&cur_tet,isovalue,iso_area,sample_num);
				}
			}
	//	free(tet_slice);
		}
	}

}
*/


#if 0

void ContourSpectrum::Build()
{
	int nx_t = _dm->nx - 1;
	int ny_t = _dm->ny - 1;
	int nz_t = _dm->nz - 1;
	for(int iz = 0; iz < nz_t; iz++)
	{
		printf("**************decompose slice x = %i****************\n",iz+1);
		for(int iy = 0; iy < ny_t; iy++)
		{
			for(int ix = 0; ix < nx_t; ix++)
			{
			//	Cute2Tetra()
			}
		}
	}
}

//template <typename T>
//void ContourSpectrum::GenSpectrumOneTetra(tetra<T> *ptet)
//{
//
//}


#endif


int main(int argc, char* argv[])
{
	double tstart, tstop, ttime;
	tstart = (double)clock()/CLOCKS_PER_SEC;

	DataMgr<float> dm;
	//for calculate the minimum and maximun data (only for BYTE data)
	/*d_min = 255;
	d_max = 0;*/
	//dm.SetDataRange(0, 255);
	//dm.readBinary(filename_data);
	dm.ReadVTKUnstructuredGrid(filename_data);

	ContourSpectrum<float> cs(&dm);
	cs.GenSpectrumTriangle();
	cs.OutputSpectrum("spectrum.txt");
	//float testValue = 0.5;

//	cout<<"The length of the value "<<testValue <<" is :"<<cs.QuerySpectrumTriangle(testValue )<<endl;
	//contour spectrum value ranges:
	//area_min = 99999999;
	//area_max = -99999999;

	//non_uni_num = 0;

	printf("*********Contour Spectrum by Xin Tong*********\n");

//	readASCII();
//	readBinary();

//	printAsciiToFile();

//	compute();
	
//	free(&f[0][0][0]);

	//d_type v[] = {5,4,6,2};
	//printf("%d %d %d %d\n",v[0],v[1],v[2],v[3]);
	//d_type *p1,*p2,*p3,*p4;
	//p1 = (d_type *)malloc(sizeof(d_type *));
	//p2 = (d_type *)malloc(sizeof(d_type *));
	//p3 = (d_type *)malloc(sizeof(d_type *));
	//p4 = (d_type *)malloc(sizeof(d_type *));
	//sort(&(v[0]),&(v[1]),&(v[2]),&(v[3]),p1,p2,p3,p4);
	//printf("%d %d %d %d\n",*p1,*p2,*p3,*p4);

	tstop = (double)clock()/CLOCKS_PER_SEC;
	ttime= tstop-tstart;
	printf("Time Cost: %f seconds\n",ttime);
	//printf("non-uniform number:%d\n",non_uni_num);
	printf("\n");

//	system("pause");
	
//	initGL(argc,argv);
	return 0;
}