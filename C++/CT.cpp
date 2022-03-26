#include "Radon.h"
#include<cmath>

static double PI(atan(1)*4);
static double PI2(PI/2);

void scan(Mat_DP& B, const Mat_DP& A)
// B = sinogram (Radon transform) of image A
// input:
//   A = image data (shape (M,N))
//   n = B.nrows()
//     = number of parallel X-rays (also used as
//       number of integration points along an X-ray
//   m = B.ncols() = number of directions of X-ray
//   if n==0, n is set to power of 2 (>M) and
//            m is set to 2*n
// output: B = sinogram of A (shape(n,m))
{
    if(B.nrows()==0) {
        int n(1);
        while(n <= A.nrows()) n<<=1;
        B.SetDims(n, n<<1);
    }
    int i,j,k;
    int M(A.nrows()),N(A.ncols());
    int n(B.nrows()),m(B.ncols());
    double X((M-1)/2.), Y((N-1)/2.), R(sqrt(X*X + Y*Y));
    double dr(2*R/(n-1)), dth(PI/m);
    double r,s,theta,cth,sth,x1,y1;
    Vec_DP x(M), y(N);
    for(i=0; i<M; i++) x[i] = i-X;
    for(j=0; j<N; j++) y[j] = j-Y;
    for(j=0; j<m; j++) {// 0 <= theta < pi
        theta = j*dth;
        cth = cos(theta);
        sth = sin(theta);
        for(i=0; i<n; i++) {// number of parallel X-rays
            r = i*dr - R;
            B[i][j] = 0;
            for(k=0; k<n; k++) {// integration along X-ray
                s = k*dr - R;
                x1 = r*cth - s*sth;
                y1 = r*sth + s*cth;
                B[i][j] += interp2d(x1,y1,x,y,A,0);
            }
        }
    }
}

void filtering(Mat_DP& B, const Mat_DP& A)
// B = high-pass filter applied to A along axis=0
// &A==&B is allowed
{
    int i,j,n(A.nrows()),m(A.ncols());
    if(n&(n-1)) error("n must be power of 2");
    double c(PI/n),c2(2./n);
    Vec_DP v(n);
    B.SetDims(n,m);
    for(j=0; j<m; j++) {
        for(i=0; i<n; i++) v[i] = A[i][j];
        realft(v,1);// FFT
        v[0] = 0;
        v[1] *= PI2;
        for(i=2; i<n; i++) v[i] *= (i>>1)*c;
        realft(v,-1);// inverse FFT
        for(i=0; i<n; i++) B[i][j] = v[i]*c2;
    }
}

void BackScan(Mat_DP& B, const Mat_DP& A)
// B = inverse Radon transform of sinogram A
// input:
//   A = sinogram after filtering (shape(n,m))
//   M,N = B.nrows(),B.ncols()
//       = height and width of output image
//   if M==0, M,N are both set to n/2
// ouput: B = image restored from A (shape(M,N))
{
    int n(A.nrows()), m(A.ncols());
    if(B.nrows()==0) B.SetDims(n>>1, n>>1);
    int i,j,k;
    int M(B.nrows()), N(B.ncols());
    double X((M-1)/2.), Y((N-1)/2.), R(sqrt(X*X + Y*Y));
    double dr(2*R/(n-1)), dth(PI/m);
    double x,y,r1,theta;
    Vec_DP r(n),f[m];
    for(i=0; i<n; i++) r[i] = i*dr - R;
    for(k=0; k<m; k++) {
        f[k].SetLength(n);
        for(i=0; i<n; i++) f[k][i] = A[i][k];
    }
    for(i=0; i<M; i++) {
        x = i-X;
        for(j=0; j<N; j++) {
            y = j-Y;
            B[i][j] = 0;
            for(k=0; k<m; k++) {// integration wrt theta
                theta = k*dth;
                r1 = x*cos(theta) + y*sin(theta);
                B[i][j] += interp(r1,r,f[k],0);
            }
            B[i][j] /= m;
        }
    }
}

void reconstruct(Mat_DP& B, const Mat_DP& A)
{
    Mat_DP C;
    filtering(C,A);
    BackScan(B,C);
}
