#ifndef __Radon_h__
#define __Radon_h__

#include "nr.h"
#include "Mat_DP.h"

struct Radon {// Discrete Radon Transform
    Mat_DP d[4];
    inline int size() const { return d[0].ncols(); }
    inline Mat_DP& operator[](int i) { return d[i]; }
    inline const Mat_DP& operator[](int i) const { return d[i]; }
    void SetSize(int);
};

void scan(Radon&, const Mat_DP&);
void BackScan(Mat_DP&, const Radon&);
void reconstruct(Mat_DP&, const Radon&);
void RadonFromSinogram(Radon&, const Mat_DP&);
void SinogramFromRadon(Mat_DP&, const Radon&);
void stitch(Mat_DP&, const Radon&);
void filtering(Radon&, const Radon&);

void scan(Mat_DP&, const Mat_DP&);// slow
void BackScan(Mat_DP&, const Mat_DP&);// slow
void filtering(Mat_DP&, const Mat_DP&);
void reconstruct(Mat_DP&, const Mat_DP&);// slow

double interp2d(double, double, const Vec_DP&, const Vec_DP&, const Mat_DP&, double);
double interp(double, const Vec_DP&, const Vec_DP&, double);
void realft(Vec_IO_DP&, const int);

#endif // __Radon_h__
