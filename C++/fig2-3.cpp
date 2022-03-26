#include "Radon.h"

main() {
    Mat_DP A,B;
    ReadBMP32(A, "fig1.bmp");
    scan(B,A);
    save("fig2.dat", B);
    WriteBMP32("fig2a.bmp", B);
    filtering(B,B);
    WriteBMP32("fig2b.bmp", B);
    BackScan(A,B);
    WriteBMP32("fig3.bmp", A);
}
