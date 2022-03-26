#include "Radon.h"

main() {
    Radon d;
    Mat_DP A;
    load(A, "fig2.dat");
    RadonFromSinogram(d,A);
    reconstruct(A,d);
    WriteBMP32("fig7.bmp", A);
}