#ifndef DATA_MGR_H
#define DATA_MGR_H

//#include "Matrix.h"
//#include "Vector.h"
#include "my_vector_types.h"
#include "vector"
#include <vtkSmartPointer.h>
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPointData.h"

using namespace std;

typedef unsigned char d_type; 
typedef float in_type; 

template <typename T>
struct tetra
{
	float3 v[4];
	float3 gradient;
	T f1,f2,f3,f4;
	//bool uniform;
	float s1,s2;
	//tetra *prev;
};
//typedef struct tetrahedron tet;

struct polynomial
{
	in_type range[2];	// range :[min max]
	in_type coe[3];
};

template <typename T>
struct polynomial2
{
	T xmin;	// range :[min max]
	T xmax;	// range :[min max]
	float a;
	float b;  //y = a * x + b
};

template <typename T>
struct spectrum2
{
//	polynomial2<T> p[2];
	float x0, x1, x2;
	float k1, k2;
};

//template <typename T>
//struct point3
//{
//	float3 p;
//
//}

template <typename T>
struct triangle
{
	float3 p[3];	//points
	float v[3];			//values
	float3 gradient;

	triangle(	float3 p0, float3 p1, float3 p2,
				float v0, float v1, float v2)
	{
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;

		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
	}
};

//typedef struct triangle triangle;



template <typename T>
class DataMgr
{
	in_type d_min, d_max;					//minimun data value and maximum data value
	d_type ***f;							//voxel data
	d_type step_x, step_y, step_z;


public:
	int nx,ny,nz;							//voxel counter
	float valueMin, valueMax;
	polynomial* allPoly;
	vector<triangle<T>> triangles;
	vector<spectrum2<T>> allSpectrums;

	int readASCII(const char* filename);
	int readBinary(const char* filename);
	void SetDataRange(d_type min, d_type max);
	void printAsciiToFile();
	void ImportTriangles(vtkUnstructuredGrid* grid, vector<triangle<T>> &trias);

	void ReadVTKUnstructuredGrid(const char* filename);
	void OutputSpectrum();
};

template <typename T>
void DataMgr<T>::ReadVTKUnstructuredGrid(const char* filename)
{
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(filename);
    reader->Update(); // Needed because of GetScalarRange
    vtkUnstructuredGrid* grid = reader->GetOutput();
    


	ImportTriangles(grid, triangles);
}



template <typename T>
void DataMgr<T>::ImportTriangles(vtkUnstructuredGrid* grid, vector<triangle<T>> &trias)
{
	vtkPoints* vtkPts = grid->GetPoints();
    vtkCellArray* vtkCls = grid->GetCells();
    cout<<"nc="<<vtkCls->GetNumberOfCells()<<endl;
	cout<<"np="<<vtkPts->GetNumberOfPoints()<<endl;

	vtkPointData* ptsData = grid->GetPointData();
	for (int i = 0; i < ptsData->GetNumberOfArrays(); i++)
	{
		std::cout << "\tArray " << i
			<< " is named "
			<< (ptsData->GetArrayName(i) ? ptsData->GetArrayName(i) : "NULL")
			<< std::endl;
	}
	vtkDataArray* ptsValue = ptsData->GetArray(0);//3 is thickness

    vtkNew<vtkIdList> pts;
    float3 p[3];
	T v[3];
    double *coord;
    for(int c = 0; c < vtkCls->GetNumberOfCells(); c++)
    {
        vtkCls->GetCell(c * 4, pts.GetPointer());
        for(int i = 0; i < 3; i++)
        {
			int id = pts->GetId(i);
            coord = vtkPts->GetPoint(id);//();
			
			//vtkPtsData->GetField(
            p[i].x = coord[0];
            p[i].y = coord[1];
			p[i].z = coord[2];
			v[i] = ptsValue->GetVariantValue(id).ToDouble();
        }
        triangle<T> t1(p[0], p[1], p[2], v[0], v[1], v[2]);
        trias.push_back(t1);
    }

	valueMin = FLT_MAX;
	valueMax = - FLT_MAX;

	for(int ip = 0; ip < vtkPts->GetNumberOfPoints(); ip++)
	{
		float value = ptsValue->GetVariantValue(ip).ToDouble();
		if(value < valueMin)
			valueMin = value;
		if(value > valueMax)
			valueMax = value;
	}
}
#endif