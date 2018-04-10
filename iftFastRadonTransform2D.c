#include "ift.h"

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
//#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)])

int sign( int x ){
    if(x > 0)
        return 1;
    if(x < 0) 
        return -1;
    return 0;
}


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

float* DDA(iftImage *img, iftVoxel p1, iftVoxel pn)
{
    int n, k;
    iftVoxel p;
    float * J;
    float Dx=(pn.x - p1.x),Dy=(pn.y - p1.y);
    float dx=0,dy=0;
    //iftCreateMatrix* J;


    if (p1.x == pn.x && p1.y == pn.y)
        n=1; 
    else{
        Dx=pn.x - p1.x;
        Dy=pn.y - p1.y; 

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

    J = (float*)(malloc(n*sizeof(float)));
    p=p1;

    // TODO: calcular I como interpolacao

    
    for (k = 0; k < n; k++)
    {
        //J[k] = img->val[ROUND(p.x)][ROUND(p.y)];
        J[k] =  iftImgVal2D(img, p.x, p.y);


        p.x = p.x + dx;
        p.y = p.y + dy;
    }

    return J;
}


int isValidPoint(iftImage *img, iftVoxel u)
{
    if ((u.x >= 0) && (u.x < img->xsize) &&
        (u.y >= 0) && (u.y < img->ysize)){
        return 1;
    }
    else{
        return 0;
    }
}


int findIntersection(iftMatrix *Po, iftImage *img, iftMatrix *N,int nx, int ny, iftVoxel *p1, iftVoxel *pn){
    float lamb1,lamb2,lamb3,lamb4;    
    int i;
    float max=0.0,min=999999999.0;
    iftVoxel v1, v2;
    p1->x=pn->x=p1->y=pn->y=-1;


    lamb1 = -(Po->val[1])/iftMax(N->val[1],1);
    lamb2 = (ny - Po->val[1])/iftMax(N->val[1],1);
    lamb3 = -(Po->val[1])/iftMax(N->val[0],1);
    lamb4 = (nx - Po->val[0])/iftMax(N->val[0],1);
    float lambs[] = {lamb1,lamb2,lamb3,lamb4};

    for(i=0; i < 4; i++){
        if (lambs[i]>max)
            max = lambs[i];
        if (lambs[i]<min)
            min = lambs[i];
    }

    v1.x = ROUND(Po->val[0] + min * N->val[0]);
    v1.y = ROUND(Po->val[1] + min * N->val[1]);

    if (isValidPoint(img, v1)){
        *p1=v1;
    }

    v2.x = ROUND(Po->val[0] + max * N->val[0]);
    v2.y = ROUND(Po->val[1] + max * N->val[1]);

    if (isValidPoint(img, v2)){
        *pn=v2;
        printf("%d, %d, %d, %d\n",p1->x,p1->y,pn->x,pn->y );
    }

    
    if ((p1->x != -1) && (pn->x != -1) && (p1->y != -1) && (pn->y != -1)){
        return 1;
    }
    else
        return 0;
}


int Intersection(iftMatrix *Po, iftImage *img, iftMatrix *N,int nx, int ny, iftVoxel *p1, iftVoxel *pn)
{
    // check interscetion to image 
    // 4 cases
    // y = ny-1, y=0, x=nx-1, x=0
    // p= (x0,y0) + lamb (Nx,Ny)


    // primeiro caso: lamb= -y0/Ny
    // segundo caso: lamb= (ny-1 - y0)/Ny
    // terceiro caso: lamb= -x0/Nx
    // quarto caso: lamb= (nx-1 - x0)/Nx

    float lamb;
    int found=0;
    iftVoxel v;

    // first
    // lamb= -y0/Ny

    lamb = -(Po->val[1])/iftMax(N->val[1],1);
    v.x = ROUND(Po->val[0] + lamb * N->val[0]);
    v.y = ROUND(Po->val[1] + lamb * N->val[1]);
    if (isValidPoint(img, v)){
        *p1 = v;
    }

    // sec
    // lamb= (ny-1 - y0)/Ny
    lamb = (ny - Po->val[1])/iftMax(N->val[1],1);
    v.x = ROUND(Po->val[0] + lamb * N->val[0]);
    v.y = ROUND(Po->val[1] + lamb * N->val[1]);
    if (isValidPoint(img, v)){
        found=1;
        *pn = v;
    }

    // third
    // lamb= -x0/Nx
    lamb = -(Po->val[1])/iftMax(N->val[0],1);
    v.x = ROUND(Po->val[0] + lamb * N->val[0]);
    v.y = ROUND(Po->val[1] + lamb * N->val[1]);
    if (isValidPoint(img, v)){
        found=1;
        *pn = v;

    }

    // forth
    // lamb= (nx-1 - x0)/Nx
    lamb = (nx - Po->val[0])/iftMax(N->val[0],1);

    v.x = ROUND(Po->val[0] + lamb * N->val[0]);
    v.y = ROUND(Po->val[1] + lamb * N->val[1]);
    //*pn = GetVoxelIndex(img, v);
    if (isValidPoint(img, v)){
        found=1;
        *pn = v;
    }

    

    // encontrar o mais proximo (menor lambda) e e mais longe (maior lambda)
    // maior lambda sera o pn e o menor o p1

    return found;
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

/* this function applies the fast Radon transform (i.e. it uses the DDA algorithm) */
iftImage *fastRadonTransform(iftImage *img)
{
    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    
    iftMatrix* normal,*normalVec;
    int Nx, Ny,i;
    float intensity;
    iftVoxel p1, pn;
    float* lineVals;

    Nx=Ny=D;
    iftImage *R = iftCreateImage(Nx, Ny, 1);

    for(int theta = 0; theta < 180; theta++) {

        iftMatrix *M = createRadonMatrix(img, theta);
        //iftImage *imgQ = iftCreateImage(D, D, 1);


        // normal calculation
        normalVec =  iftCreateMatrix(1, 4);
        normalVec->val[1] = 1.0;
        normalVec->val[0] = 0.0;
        normalVec->val[2] = 1.0;
        normalVec->val[3] = 0.0;

        // percorrendo a linha
        for(int p = 0; p < R->n; p++) {
            //printf("p = %d\n", p);
            
            // Po = M-1 * p
            //iftMatrix *I_ = imagePixelToMatrix(img, p);
            iftMatrix *I_ = iftCreateMatrix(1, 4);
            iftMatrixElem(I_, 0, 0) = p;
            iftMatrixElem(I_, 0, 1) = -D/2;
            iftMatrixElem(I_, 0, 2) = 0;
            iftMatrixElem(I_, 0, 3) = 1;


            // rotating line
            iftMatrix *P0_line = iftMultMatrices(M, I_);
            // iftPrintMatrix(M);

            // normal = M * normal
            //iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);
            normal = iftMultMatrices(M, normalVec);
            
            if(findIntersection(P0_line, img, normal,img->xsize,img->ysize, &p1, &pn))
            {
                //printf("%d,%d\n",pn.x, pn.y);
                // chamar o DDA
                // atribuir o valor obtido de J em Pi
                
                lineVals = DDA(img, p1, pn);

                intensity = lineVals[0];
                
                for (i = 1; i < 200; i++)
                    intensity+= lineVals[i];
                
                R->val[p] = intensity;

            }
            else{
                // sem intercessao
                
                R->val[p] = 0;
            }
        }
        
    }

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
    iftImage *normalizedImage= iftNormalize(imgRadon,0,255);

    /* save the resulting images */
    char fileName[256];
    sprintf(fileName, "fast_radon_transform_%s.png", iftFilename(imgFileName, iftFileExt(imgFileName)));
    iftWriteImageByExt(normalizedImage, fileName);

    iftDestroyImage(&img);
    iftDestroyImage(&imgRadon);

    return(0);
}
