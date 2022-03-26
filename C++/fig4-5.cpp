#include "Radon.h"

main() {
    Radon d;
    Mat_DP A;
    ReadBMP32(A, "fig1.bmp");
    scan(d,A);
    stitch(A,d);
    WriteBMP32("fig4.bmp", A);
    reconstruct(A,d);
    WriteBMP32("fig5.bmp", A);
}