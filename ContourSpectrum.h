#ifndef CONTOUR_SPECTRUM_H
#define CONTOUR_SPECTRUM_H
#include "DataMgr.h"
#include "fstream"
template <typename T>
class ContourSpectrum
{
	DataMgr<T>* _dm;
	//void ComputeContourSpectrum(
	//	vec3<in_type> v0, vec3<in_type> v1, vec3<in_type> v2, vec3<in_type> v3, 
	//	in_type value[4],
	//	float poly0[4], float poly1[4], float poly2[4]);

	void Build();
	
	//void GenLoc(int i0, int i1, int i2, int i3,
	//	float3 v0, float3 v1, float3 v2, float3 v3);

	//template <typename T>
	//inline void Cube2Tetra(
	//d_type v0, d_type v1, d_type v2, d_type v3,
	//d_type v4, d_type v5, d_type v6, d_type v7,
	//tetra<T> t[5]);

	//template <typename T>
	//void GenSpectrumOneTetra(tetra<T> *ptet);

	//template <typename T>

public:
	ContourSpectrum(DataMgr<T>* dm);
	void GenSpectrumTriangle();
	float QuerySpectrumTriangle(float v);
	void OutputSpectrum(char* filename);
};

template <typename T>
ContourSpectrum<T>::ContourSpectrum(DataMgr<T>* dm)
{
	_dm = dm;
}



template <typename T>
inline spectrum2<T> triangle2spectrum(triangle<T> tr)
{
	spectrum2<T> ret;
	int i0, i1, i2; //v[i0] < v[i1] < v[i2]
	//sort
	if(tr.v[0] < tr.v[1]) {//0<1
		if(tr.v[1] < tr.v[2]) {
			i0 = 0;
			i1 = 1;
			i2 = 2;
		}else{//2<1
			if(tr.v[0] < tr.v[2]){
				i0 = 0;
				i1 = 2;
				i2 = 1;
			}else {
				i0 = 2;
				i1 = 0;
				i2 = 1; }}
	}else{//1<0
		if(tr.v[2] < tr.v[1]){
			i0 = 2;
			i1 = 1;
			i2 = 0;
		}else{	//1<2
			if(tr.v[0] < tr.v[2]){
				i0 = 1;
				i1 = 0;
				i2 = 2;
			}else{
				i0 = 1;
				i1 = 2;
				i2 = 0;}}}
	/*cout<<tr.v[0]<<","<<tr.v[1]<<","<<tr.v[2]<<endl;
	cout<<i0<<","<<i1<<","<<i2<<endl;*/
	//
	ret.x0 = tr.v[i0];
	ret.x1 = tr.v[i1];
	ret.x2 = tr.v[i2];
	float v2_v0 = tr.v[i2] - tr.v[i0];
	if( v2_v0 <= 0)
	{
		cout<<"v2_v0:"<<v2_v0 <<endl;
	}
	else
	{
		float v1_v0 = ret.x1 - ret.x0;
		float v2_v1 = ret.x2 - ret.x1;
		//bug here!!!!
		float3 pMid = 	(v1_v0 / v2_v0) * (tr.p[i2] - tr.p[i0]) + tr.p[i0];
		float len = length(pMid - tr.p[i1]);
		if(v1_v0 > 0 )
			ret.k1 = len / v1_v0;
		else
			ret.k1 = 0;
		if(v2_v1 > 0)
			ret.k2 = len / v2_v1;
		else
			ret.k2 = 0;
	}
	return ret;
}


template <typename T>
void ContourSpectrum<T>::GenSpectrumTriangle()
{
	
	for(int i = 0; i < _dm->triangles.size(); i++)
	{
		//cout<<i <<"\t";
		_dm->allSpectrums.push_back( triangle2spectrum(_dm->triangles[i]));
	}
}

template <typename T>
void ContourSpectrum<T>::OutputSpectrum(char* filename)
{
	int nBin = 20;
	float step = (_dm->valueMax - _dm->valueMin) / 20;
	ofstream of(filename);
	for(float v = _dm->valueMin; v <= _dm->valueMax ; v+= step)
		of<<QuerySpectrumTriangle(v)<<"\t";
	of.close();
}

template<typename T>
float ContourSpectrum<T>::QuerySpectrumTriangle(float v)
{
	float lengthSum = 0;
	for(int i = 0; i < _dm->allSpectrums.size(); i++)
	{
		if(_dm->allSpectrums[i].x0 <= v)
		{
			if(v < _dm->allSpectrums[i].x1 )
			{
				lengthSum += ((v - _dm->allSpectrums[i].x0) * _dm->allSpectrums[i].k1);
			}
			else if(v < _dm->allSpectrums[i].x2 )
			{
				lengthSum += ((_dm->allSpectrums[i].x2 - v) * _dm->allSpectrums[i].k2);
			}
		}
	}
	return lengthSum;
}


#endif