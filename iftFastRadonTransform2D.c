#include "ift.h"

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
//#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)])

int sign( int x ){
    if(x >= 0)
        return 1;    
    return -1;
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


int LinearInterpolationValue(iftImage *img, float x, float y)
{
    iftVoxel u[4];
    float dx = 1.0;
    float dy = 1.0;
    float  P12, P34;
    int Pi;

    if ((int) (x + 1.0) == img->xsize)
        dx = 0.0;
    if ((int) (y + 1.0) == img->ysize)
        dy = 0.0;

    //closest neighbour in each direction
    u[0].x = (int)x;      u[0].y = (int)y;       
    u[1].x = (int)(x + dx); u[1].y = (int)y;       
    u[2].x = (int)x;      u[2].y = (int)(y + dy);  
    u[3].x = (int)(x + dx); u[3].y = (int)(y + dy);  


    P12 = (float)iftImgVal2D(img,u[1].x,u[1].y) * (x - u[0].x) + (float)iftImgVal2D(img,u[0].x,u[0].y) * (u[1].x - x);
    P34 = (float)iftImgVal2D(img,u[3].x,u[3].y) * (x - u[2].x) + (float)iftImgVal2D(img,u[2].x,u[2].y) * (u[3].x - x);
    Pi  = (int)P34 * (y - u[0].y) + P12 * (u[2].y - y);
    

    return Pi;
}


int DDA(iftImage *img, iftVoxel p1, iftVoxel pn)
{
    int n, k;
    //iftVoxel p;
    float px, py;
    float J=0;
    int Dx,Dy;
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
            dy = (dx * Dy)/Dx;
        }
        else{ 
            n = abs(Dy)+1;
            dy = sign(Dy);
            dx = (dy * Dx)/Dy;
        }
    }

    px = p1.x;
    py = p1.y;
    

    // TODO: calcular I como interpolacao

    for (k = 1; k < n; k++)
    {   
        J+=  (float)LinearInterpolationValue(img, px, py);
        //J+=  iftImgVal2D(img, (int)px, (int)py);
        //J+=  (float)LinearInterpolationValue(img, px, py);

        px = px + dx;
        py = py + dy;
    }

    return (int)J;
}


int isValidPoint(iftImage *img, iftVoxel u)
{
    //printf("u = %d, %d\n",u.x,u.y );
    //printf("img size = %d, %d\n",img->xsize,img->ysize );
    if ((u.x >= 0) && (u.x < img->xsize) &&
        (u.y >= 0) && (u.y < img->ysize)){
        return 1;
    }
    else{
        return 0;
    }
}


int findIntersection(iftMatrix *Po, iftImage *img, iftMatrix *N,int nx, int ny, iftVoxel *p1, iftVoxel *pn){
    float Nx, Ny;
    int x0, y0;
    float lamb;
    int found =0;

    iftVoxel v;
    p1->x=pn->x=p1->y=pn->y=-1;
    float max=-9999999.9, min=9999999.9;
    Nx = N->val[0];
    Ny = N->val[1];
    y0 = Po->val[1];
    x0 = Po->val[0];   

    if (Ny)
    {
        lamb=-y0/Ny;
        v.x = x0 + lamb*Nx;
        v.y = y0 + lamb*Ny;
        if (isValidPoint(img,v))
        {
            found+=1;
            p1->x = v.x;
            p1->y = v.y;
            max=lamb;
            min=lamb;
        }
        lamb=(ny-1-y0)/Ny;
        v.x = x0 + lamb*Nx;
        v.y = y0 + lamb*Ny;
        if (isValidPoint(img,v) && ((lamb > max) || (lamb < min)))
        {
            found+=1;
            if (p1->x != -1)
            {
                pn->x = v.x;
                pn->y = v.y;    
                max = lamb;
            }
            else{
                p1->x = v.x;
                p1->y = v.y;
                min=lamb;
            }
        }

    }

    if (Nx)
    {
        lamb=-x0/Nx;
        v.x = x0 + lamb*Nx;
        v.y = y0 + lamb*Ny;
        if (isValidPoint(img,v) && ((lamb > max) || (lamb < min)))
        {
            found+=1;
            if (p1->x != -1)
            {
                pn->x = v.x;
                pn->y = v.y;
                max = lamb;    
            }
            else{
                p1->x = v.x;
                p1->y = v.y;
                min = lamb;
            }
        }
        lamb=(nx-1-x0)/Nx;
        v.x = x0 + lamb*Nx;
        v.y = y0 + lamb*Ny;
        if (isValidPoint(img,v) && ((lamb > max) || (lamb < min)))
        {
            found+=1;
            if (p1->x != -1)
            {
                pn->x = v.x;
                pn->y = v.y;  
                max = lamb;  
            }
            else{
                p1->x = v.x;
                p1->y = v.y;
                min = lamb;
            }
        }
    }

    if (p1->x > pn->x && p1->y > pn->y){
        int auxX, auxY;
        auxX = p1->x;
        auxY = p1->y;
        p1->x = pn->x;
        p1->y = pn->y;
        pn->x = auxX;
        pn->y = auxY;
        
    }

    return (found==2);
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

    iftMatrix *normalVec;
    iftVoxel p1, pn;
    int sumIntensities;
    int progress = 0;

    iftImage *R = iftCreateImage(180, D, 1);

    for(int theta = 0; theta < 180; theta++) {


        //fprintf(stdout,"Progress: %f \r", (progress/180.0*100.0)); fflush(stdout);
        iftMatrix *M = createRadonMatrix(img, theta);

        // normal calculation
        normalVec =  iftCreateMatrix(1, 4);
        iftMatrixElem(normalVec, 0, 0) = 0;
        iftMatrixElem(normalVec, 0, 1) = 1;
        iftMatrixElem(normalVec, 0, 2) = 0;
        iftMatrixElem(normalVec, 0, 3) = 0;
        //iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);
        iftMatrix* normal = iftMultMatrices(M, normalVec);
        // percorrendo a linha
        for(int p = 0; p < D; p++) {

            // Po = M-1 * p
            //iftMatrix *I_ = imagePixelToMatrix(img, p);
            iftMatrix *I_ = iftCreateMatrix(1, 4);
            iftMatrixElem(I_, 0, 0) = p;
            iftMatrixElem(I_, 0, 1) = -D/2;
            iftMatrixElem(I_, 0, 2) = 0;
            iftMatrixElem(I_, 0, 3) = 1;

            //iftPrintMatrix(I_);
            //iftPrintMatrix(M);
            // rotating line
            //iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);
            iftMatrix *P0_line = iftMultMatrices(M,I_ );
            //iftPrintMatrix(P0_line);
            


            // normal = M * normal
            //iftMatrix *rotMatrix = iftRotationMatrix(IFT_AXIS_Z, theta);
            
            if(findIntersection(P0_line, img, normal,img->xsize,img->ysize, &p1, &pn))
            {

                // chamar o DDA
                // atribuir o valor obtido de J em Pi
                if (p1.x == pn.x && p1.y == pn.y)
                    iftImgVal2D(R, theta, p) = iftImgVal2D(img, p1.x, p1.y);
                
                

                //intensity = sumIntensities[0];
                //for (i = 1; i < line_len; i++)
                //    intensity+= sumIntensities[i];
                
                //printf("%f \n",intensity);
                else{
                    sumIntensities = DDA(img, p1, pn);
                    
                    iftImgVal2D(R, theta, p) = sumIntensities;
                }

                //R->val[p] = intensity;
            }
            else{
                // sem intercessao
                iftImgVal2D(R, theta, p) = 0;
                //R->val[p] = 0;
            }

            //printf("(p,theta) =  (%d, %d) \n",p,theta );
            iftDestroyMatrix(&P0_line);
            iftDestroyMatrix(&I_);
        }

        iftDestroyMatrix(&M);
        iftDestroyMatrix(&normalVec);
        iftDestroyMatrix(&normal);
        progress++;   
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
    iftWriteImageByExt(normalizedImage , fileName);

    iftDestroyImage(&img);
    iftDestroyImage(&imgRadon);

    return(0);
}



