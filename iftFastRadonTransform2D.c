#include "ift.h"

/* this function creas the rotation/translation matrix for the given theta */
iftMatrix *createRadonMatrix(iftImage *img, int theta)
{
    iftMatrix *resMatrix = NULL;

    iftVector v1 = {.x = (float)img->xsize / 2.0, .y = (float)img->ysize / 2.0, .z = 0.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);

    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftVector v2 = {.x = -(D / 2.0), .y = -(D / 2.0), .z = 0.0};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);

    resMatrix = iftMultMatricesChain(3, transMatrix1, rotMatrix, transMatrix2);

    return resMatrix;
}

/* this function uses de DDA algorithm to paint a line from p1 to pn in img */
void dda(iftImage *img, iftVoxel p1, iftVoxel pn)
{
    int n, k;
    float H=0.0;
    iftVoxel p;

    if (p1->x == pn->x && p1->y == pn->y)
        n=1; 
    else{
        Dx=pn->x - p1->x;
        Dy=pn->y - p1->y; 

        if( abs(Dx) >= abs(Dy) ){
            n = abs(Dx)+1;
            dx = sign(Dx);
            dy = dx * (Dy/Dx);
        }
        else{ 
            n = abs(Dy)+1;
            dy = sign(Dy);
            dx = dy * (Dx/Dy);
        }
    }

    p=p1;

    // calcular H como interpolacao

    for (k = 0; k < n; i++)
    {
        img->val[(int)(p->x)][(int)(p->y)] = H;
        p->x = p->x + dx;
        p->y = p->y + dy;
    }


}

/* this function applies the fast Radon transform (i.e. it uses the DDA algorithm) */
iftImage *fastRadonTransform(iftImage *img)
{
    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftImage *R = iftCreateImage(180, D, 1);





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
    iftImage *imgRadon = fastRadonTransform(img);
    printf("Time to compute the Radon Transform: %s\n", iftFormattedTime(iftCompTime(t1, iftToc())));

    /* save the resulting images */
    char fileName[256];
    sprintf(fileName, "fast_radon_transform_%s.png", iftFilename(imgFileName, iftFileExt(imgFileName)));
    iftWriteImageByExt(imgRadon, fileName);

    iftDestroyImage(&img);
    iftDestroyImage(&imgRadon);

    return(0);
}
