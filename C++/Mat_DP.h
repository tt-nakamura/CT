#ifndef __Mat_DP_h__
#define __Mat_DP_h__

#include "Mat.h"
#include "Mat3D.h"
#include<cmath>

double max(const Mat_DP&);
double min(const Mat_DP&);

void save(const char*, Mat_DP&);
void load(Mat_DP&, const char *file_name);

void WriteBMP32(const char*, const Mat_DP&);
void ReadBMP32(Mat_DP&, const char*, int=0);
void WriteBMP32(const char*, const Mat3D_DP&);
void ReadBMP32(Mat3D_DP&, const char*);

#endif // __Mat_DP_h__
