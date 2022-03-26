// Dicrete Radon Transform
// referece:
//   M. L. Brady, "A Fast Discrete Approximation Algorithm for the
//     Radon Transform" SIAM Journal on Computing 27 (1998) 107
//   W. H. Press, "Dicrete Radon Transform has an Exact, Fast Inverse..."
//     Proceedings of the National Academy of Sciences 103 (2006) 19249

#include "Radon.h"
#include<cmath>

static double PI4(atan(1)); // pi/4
static double PI2(PI4*2);   // pi/2
static double PI(PI2*2);

void Radon::SetSize(int n) {
    int i,n2(n*2);
    if(n&(n-1)) error("n must be power of 2");
    for(i=0; i<4; i++) d[i].SetDims(n2,n,0.);
}

void scan(Mat_DP& a, int h, int y)
// recursive Radon transform
// input:
//   a = image data (shape(2n,n))
//   h = width of transform region
//   y = left coordinate of the region
// output:
//   a[i,y+j] = sum_{k=0}^{h-1} a[x,y+k]
//                for 0<=i<2n and 0<=j<h where
//     (x,y+k) moves from (i,y) to (i-j,y+h-1)
{
    int i,j,k,l,m(a.ncols()+h);
    int h1(h>>1),y1(y+h1);
    double b[h];
    if(h1>1) {// divide and conquer
        scan(a,h1,y);
        scan(a,h1,y1);
    }
    for(i=m-1; i>=0; i--) {
        for(j=0; j<h; j++) {
            k = j>>1; l = (j+1)>>1;
            b[j] = a[i][y+k];
            if(i>=l) b[j] += a[i-l][y1+k];
        }
        for(j=0; j<h; j++) a[i][y+j] = b[j];
    }
}

void scan(Radon& d, const Mat_DP& A)
// d = fast Radon transform of A
// input: A = image data (shape(n,n))
// output: d(r,theta) (shape(4,2n,n))
//   d[0,i,j] = sum_{y=0}^{n-1} A[x,y] (0<=theta<=45)
//     (x,y) moves from (i,0) to (i-j,n-1)
//   d[1,i,j] = sum_{x=0}^{n-1} A[x,y] (45<=theta<=90)
//     (x,y) moves from (0,i) to (n-1,i-j)
//   d[2,i,j] = sum_{x=n-1}^0 A[x,y] (90<=theta<=135)
//     (x,y) moves from (n-1,i) to (0,i-j)
//   d[3,i,j] = sum_{y=n-1}^0 A[x,y] (135<=theta<=180)
//     (x,y) moves from (i,n-1) to (i-j,0)
{
    int n(A.nrows());
    if(A.ncols()!=n) error("image must be square");
    int i,j,k;
    d.SetSize(n);
    for(i=0; i<n; i++) for(j=0; j<n; j++) {
        d[0][i][j] = A[i][j];
        d[1][i][j] = A[j][i];
        d[2][i][j] = A[n-1-j][i];
        d[3][i][j] = A[n-1-i][j];
    }
    for(k=0; k<4; k++) scan(d[k],n,0);
}

void BackScan(Mat_DP& a, int h, int y)
// recursive inverse Radon transform
// input:
//   a = scanned data (shape(2n,n))
//   h = width of transform region
//   y = left coordinate of the region
// output:
//   a[i,y+j] = sum_{k=0}^{h-1} a[x,y+k]
//                for 0<=i<2n and 0<=j<h where
//     (x,y+k) moves from (i,y) to (i+j,y+h-1)
{
    int i,j,k,l,m(a.nrows()-h);
    int h1(h>>1),y1(y+h1);
    double b[h];
    if(h1>1) {// divide and conquer
        BackScan(a,h1,y);
        BackScan(a,h1,y1);
    }
    for(i=0; i<m; i++) {
        for(j=0; j<h; j++) {
            k = j>>1; l = (j+1)>>1;
            b[j] = a[i][y+k] + a[i+l][y1+k];
        }
        for(j=0; j<h; j++) a[i][y+j] = b[j];
    }
}

void BackScan(Mat_DP& A, const Radon& d)
// A = inverse fast Radon transform of d
// input: d = scanned and filtered data (shape(4,2n,n))
// output: A = image restored from d (shape(n,n))
//   A[i,j] = sum_{y=0}^{n-1} (
//       d[0,x,y] (x,y) moves from (i,0) to (i+j,n-1)
//     + d[1,x,y] (x,y) moves from (j,0) to (j+i,n-1)
//     + d[2,x,y] (x,y) moves from (j,0) to (j+i',n-1)
//     + d[3,x,y] (x,y) moves from (i',0) to (i'+j,n-1)
//   )/4/(n-1)    where i'=n-1-i
{
    int i,j,n(d.size());
    double c(0.25/(n-1));
    Mat_DP a[4];
    for(i=0; i<4; i++) {
        a[i] = d[i];
        // avoid double counting rays
        if(i&1) for(j=0; j<n; j++)
            a[i][j][0] = a[i][j][n-1] = 0;
        BackScan(a[i],n,0);
    }
    A.SetDims(n,n);
    for(i=0; i<n; i++) for(j=0; j<n; j++)
        A[i][j] = (a[0][i][j]
                 + a[1][j][i]
                 + a[2][j][n-1-i]
                 + a[3][n-1-i][j])*c;
}

void stitch(Mat_DP& A, const Radon& d)
//   /|\   /|\
//  / | \ / | \
// |  |  |  |  |
// |  |  |  |  |
// | / \ | / \ |
// |/   \|/   \|
// 0 45 90 135 180
{
    int i,j,k,n(d.size()),n2(n*2),n4(n*4);
    A.SetDims(n2,n4,0.);
    for(j=0; j<n; j++) {
        k = (n-j)>>1;
        for(i=0; i<n2-k; i++) {
            A[i+k][j]      = d[0][i][j];
            A[i+k][n2-1-j] = d[1][i][j];
            A[i+k][n2+j]   = d[2][i][j];
            A[i+k][n4-1-j] = d[3][i][j];
        }
    }
}

void RadonFromSinogram(Radon& d, const Mat_DP& A)
// input: A = output of scan() in CT.cpp
//        n = d.size() = image size
//        if n==0, n is set to A.nrows()/2
// output: d = input to BackScan() in FastCT.cpp
{
    if(d.size()==0) d.SetSize(A.nrows()/2);
    int i,j,n(d.size()),n2(n*2);
    int N(A.nrows()), M(A.ncols());
    double n1(n-1), R(n1/sqrt(2));
    double dr(2*R/(N-1)), dth(PI/M);
    double r1,th1,cth,jn;
    Vec_DP r(N),th(M);
    for(i=0; i<N; i++) r[i] = i*dr - R;
    for(i=0; i<M; i++) th[i] = i*dth;
    for(j=0; j<n; j++) {
        th1 = atan2(j,n1);// slope
        cth = cos(th1);
        jn = (j+n1)/2;
        for(i=0; i<n2; i++) {
            r1 = (i - jn)*cth;// transverse coordinate
            d[0][i][j] = interp2d(r1,th1,    r,th,A,0)*cth;
            d[1][i][j] = interp2d(r1,PI2-th1,r,th,A,0)*cth;
            d[2][i][j] = interp2d(r1,PI2+th1,r,th,A,0)*cth;
            d[3][i][j] = interp2d(r1,PI-th1, r,th,A,0)*cth;
        }
    }
    for(i=0; i<n; i++) d[3][i][0] = d[0][n-1-i][0];
}

void SinogramFromRadon(Mat_DP& A, const Radon& d)
// input:
//   d = output of scan() in FastCT.cpp
//   N = A.nrows() = number of parallel X-rays
//   M = A.ncols() = number of directions of X-rays
//   f N==0, (N,M) are set to (d.size()*2, d.size()*4)
// output: A = input to BackScan() in CT.cpp
{
    int n(d.size());
    if(A.nrows()==0) A.SetDims(n<<1, n<<2);
    int i,j,k,n2(n*2);
    int N(A.nrows()),M(A.ncols());
    double n1(n-1), R((n1)/sqrt(2));
    double dr(2*R/(N-1)), dth(PI/M);
    double th,x1,y1,sc,yn;
    Vec_DP x(n2),y(n);
    for(i=0; i<n2; i++) x[i] = i;
    for(j=0; j<n; j++) y[j] = j;
    for(j=0; j<M; j++) {
        th = j*dth;
        k = int(floor(th/PI4));// 0,1,2,3
        th = (th - PI2*(k+1>>1))*(k&1 ? -1:1);
        y1 = n1*tan(th);
        yn = (y1+n1)/2;
        sc = 1/cos(th);
        for(i=0; i<N; i++) {
            x1 = (i*dr - R)*sc + yn;
            A[i][j] = interp2d(x1,y1,x,y,d[k],0)*sc;
        }
    }
}

void filtering(Radon& b, const Radon& a)
// b = high-pass filter applied to a
// &b==&a is allowed
{
    int n(a.size());
    if(b.size() != n) b.SetSize(n);
    int i,j,k,n2(n*2);
    double c(PI/n2),c2(1./n);
    Vec_DP v(n2);
    for(k=0; k<4; k++) for(j=0; j<n; j++) {
        for(i=0; i<n2; i++) v[i] = a[k][i][j];
        realft(v,1);// FFT
        v[0] = 0;
        v[1] *= PI2;
        for(i=2; i<n2; i++) v[i] *= (i>>1)*c;
        realft(v,-1);// inverse FFT
        for(i=0; i<n2; i++) b[k][i][j] = v[i]*c2;
    }
}

void reconstruct(Mat_DP& A, const Radon& d)
{
    Radon a;
    filtering(a,d);
    BackScan(A,a);
}