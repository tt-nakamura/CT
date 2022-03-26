#include<fstream>
#include "Mat_DP.h"

double min(const Mat_DP& A)
{
    double a(A[0][0]);
    for(int i=0; i<A.nrows(); i++)
        for(int j=0; j<A.ncols(); j++)
            if(A[i][j] < a) a = A[i][j];
    return a;
}

double max(const Mat_DP& A)
{
    double a(A[0][0]);
    for(int i=0; i<A.nrows(); i++)
        for(int j=0; j<A.ncols(); j++)
            if(A[i][j] > a) a = A[i][j];
    return a;
}

void save(const char *file_name, Mat_DP& A)
{
    int i,j,m(A.nrows()),n(A.ncols());
    std::ofstream s(file_name, std::ofstream::binary);
    s.write((const char *)&m, sizeof(int));
    s.write((const char *)&n, sizeof(int));
    for(i=0; i<m; i++) for(j=0; j<n; j++)
        s.write((const char*)&A[i][j], sizeof(double));
}

void load(Mat_DP& A, const char *file_name)
{
    int i,j,m,n;
    std::ifstream s(file_name, std::ifstream::binary);
    s.read((char *)&m, sizeof(int));
    s.read((char *)&n, sizeof(int));
    A.SetDims(m,n);
    for(i=0; i<m; i++) for(j=0; j<n; j++)
        s.read((char *)&A[i][j], sizeof(double));
}

void WriteBMP32(const char *file_name, const Mat_DP& A)
// write matrix data A to bitmap file
// R,G,B are set to same value A[i,j]
//   and normalized to values between 0 and 1
// i increases from top to bottom
// j increases from left to right
{
    int i,j,k,m(A.nrows()),n(A.ncols());
    double a(min(A)),b(max(A)-a);
    Mat3D_DP B;
    B.SetDims(m,n,3);
    for(i=0; i<m; i++) for(j=0; j<n; j++)
        for(k=0; k<3; k++) B[i][j][k] = (A[i][j]-a)/b;
    WriteBMP32(file_name, B);
}

void ReadBMP32(Mat_DP& A, const char *file_name, int color)
// read matrix data A from bitmap file
// color = 0,1,2 for R,G,B
// A[i,j] are real value between 0 and 1
{
    int i,j;
    Mat3D_DP B;
    ReadBMP32(B, file_name);
    A.SetDims(B.dim1(), B.dim2());
    for(i=0; i<B.dim1(); i++)
        for(j=0; j<B.dim2(); j++)
            A[i][j] = B[i][j][color];
}
