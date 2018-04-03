#include "ift.h"

iftMatrix *createRadonMatrix(iftImage *img, int theta)
{
    iftMatrix *resMatrix = NULL;

    iftVector v1 = {.x = -((float)img->xsize / 2.0), .y = -((float)img->ysize / 2.0), .z = 0.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);

    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftVector v2 = {.x = (D / 2.0), .y = (D / 2.0), .z = 0.0};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);

    resMatrix = iftMultMatricesChain(3, transMatrix2, rotMatrix, transMatrix1);

    return resMatrix;
}

iftMatrix *imagePixelToMatrix(iftImage *img, int p)
{
    iftMatrix *pixMat = iftCreateMatrix(1, 4);
    iftMatrixElem(pixMat, 0, 0) = p % img->xsize;
    iftMatrixElem(pixMat, 0, 1) = p / img->xsize;
    iftMatrixElem(pixMat, 0, 2) = 0;
    iftMatrixElem(pixMat, 0, 3) = 1;

    return pixMat;
}

iftImage *radonTransform(iftImage *img)
{
    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftImage *R = iftCreateImage(180, D, 1);

    int progress = 0;
#pragma omp parallel for shared(progress)
    for(int theta = 0; theta < 180; theta++) {
        fprintf(stdout, "Progress: %.1f %s\r", (progress/180.0*100.0), "%"); fflush(stdout);

        /* compute the translated/rotated image (pixel-wise) */
        iftMatrix *M = createRadonMatrix(img, -theta);
        iftImage *imgQ = iftCreateImage(D, D, 1);
        for(int p = 0; p < img->n; p++) {
            iftMatrix *Ip = imagePixelToMatrix(img, p);
            iftMatrix *Iq = iftMultMatrices(M, Ip);
            int x = (int)iftMatrixElem(Iq, 0, 0);
            int y = (int)iftMatrixElem(Iq, 0, 1);
            iftImgVal2D(imgQ, x, y) = img->val[p];
            iftDestroyMatrix(&Ip);
            iftDestroyMatrix(&Iq);
        }

        /* apply the Radon transform */
        for(int rho = 0; rho < D; rho++) {
            int xq = rho;
            for(int yq = 0; yq < D; yq++)
                iftImgVal2D(R, theta, rho) += iftImgVal2D(imgQ, xq, yq);
        }

        iftDestroyMatrix(&M);
        iftDestroyImage(&imgQ);
        progress++;
    }
    fprintf(stdout, "\n"); fflush(stdout);

    return R;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
        iftError("Usage: Reconstruction <input-image.png>","main");

    timer *t1 = iftTic();
    char *imgFileName = iftCopyString(argv[1]);

    /* compute the Radon transform */
    iftImage *img = iftReadImageByExt(imgFileName);
    iftImage *imgRadon = radonTransform(img);
    iftImage *imgRadonNorm = iftNormalize(imgRadon, 0.0, 255.0); // normalize only to apply the color table
    printf("Time to compute the Radon Transform: %s\n", iftFormattedTime(iftCompTime(t1, iftToc())));

    /* create a color table for the resulting image */
    iftImage *imgRadonColTab = iftCreateColorImage(imgRadon->xsize, imgRadon->ysize, 1, 8);
    iftColorTable *ctb = iftCreateHotIronColorTable(256);
    for(int p=0; p < imgRadon->n; p++) {
        iftColor YCbCr = iftRGBtoYCbCr(ctb->color[imgRadonNorm->val[p]], 255);
        imgRadonColTab->val[p] = YCbCr.val[0];
        imgRadonColTab->Cb[p] = YCbCr.val[1];
        imgRadonColTab->Cr[p] = YCbCr.val[2];
    }

    /* save the resulting images */
    char fileName[256];
    sprintf(fileName, "radon_transform_%s.png", iftFilename(imgFileName, iftFileExt(imgFileName)));
    iftWriteImageByExt(imgRadon, fileName);

    sprintf(fileName, "radon_transform_%s_colortable.png", iftFilename(imgFileName, iftFileExt(imgFileName)));
    iftWriteImageByExt(imgRadonColTab, fileName);

    iftDestroyImage(&img);
    iftDestroyImage(&imgRadon);
    iftDestroyImage(&imgRadonNorm);
    iftDestroyImage(&imgRadonColTab);
    iftDestroyColorTable(&ctb);

    return(0);
}
