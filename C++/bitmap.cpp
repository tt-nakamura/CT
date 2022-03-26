#include<fstream>
#include<cmath>
#include "Mat3D.h"

void WriteBMP32(const char *file_name, const Mat3D_DP& M)
// write 32bit bitmap image to file
// M[i][j][k] = intensity of color k at pixel(i,j)
//   as real numbers between 0 and 1
//   i increases from top to bottom
//   j increases from left to right
//   k = 0,1,2 for R,G,B
{
    long i,j,k;
    long biHeight(M.dim1()), biWidth(M.dim2());
    unsigned char c;
    unsigned long bfOffset(54);
    unsigned long biSizeImage(biHeight * biWidth * 4);
    unsigned long bfSize(biSizeImage + bfOffset);
    unsigned long biSize(40);
    unsigned long zero(0);
    unsigned int one(1);
    unsigned int biBitCount(32);
    std::ofstream s(file_name, std::ofstream::binary);
    s.write("BM", 2);
    s.write((const char*)&bfSize, 4);
    s.write((const char*)&zero, 2);// bfReserved1
    s.write((const char*)&zero, 2);// bfReserved2
    s.write((const char*)&bfOffset, 4);
    s.write((const char*)&biSize, 4);
    s.write((const char*)&biWidth, 4);
    s.write((const char*)&biHeight, 4);
    s.write((const char*)&one, 2);// biPlanes
    s.write((const char*)&biBitCount, 2);
    s.write((const char*)&zero, 4);// biCompression
    s.write((const char*)&biSizeImage, 4);
    s.write((const char*)&zero, 4);// biXPixPerMeter
    s.write((const char*)&zero, 4);// biYPixPerMeter
    s.write((const char*)&zero, 4);// biClrUsed
    s.write((const char*)&zero, 4);// biClrImportant
    for(i=biHeight-1; i>=0; i--)
        for(j=0; j<biWidth; j++) {
            for(k=2; k>=0; k--) {
                c = (unsigned char)(round(M[i][j][k]*255));
                s.put(c);
            }
            s.put(0);
        }
}

void ReadBMP32(Mat3D_DP& M, const char *file_name)
// read 32bit bitmap image from file
{
    long i,j,k;
    long biHeight, biWidth;
    unsigned char c;
    std::ifstream s(file_name, std::ifstream::binary);
    s.seekg(18);
    s.read((char*)&biWidth, 4);
    s.read((char*)&biHeight, 4);
    M.SetDims(biHeight, biWidth, 3);
    s.seekg(28);
    s.seekg(54);
    for(i=biHeight-1; i>=0; i--)
        for(j=0; j<biWidth; j++) {
            for(k=2; k>=0; k--) {
                c = s.get();
                M[i][j][k] = c/255.;
            }
            s.get();
        }
}

void ReadBMP8(Mat3D_DP& M, const char *file_name)
// read 8bit bitmap image from file
{
    long i,j,k;
    long biHeight, biWidth;
    unsigned char c;
    unsigned char color[256][4];
    std::ifstream s(file_name, std::ifstream::binary);
    s.seekg(18);
    s.read((char*)&biWidth, 4);
    s.read((char*)&biHeight, 4);
    M.SetDims(biHeight, biWidth, 3);
    s.seekg(54);
    for(i=0; i<256; i++)
        for(j=0; j<4; j++) s.read((char*)&color[i][j], 1);
    for(i=biHeight-1; i>=0; i--)
        for(j=0; j<biWidth; j++) {
            c = s.get();
            for(k=0; k<3; k++)
                M[i][j][k] = color[c][2-k]/255.;
        }
}
