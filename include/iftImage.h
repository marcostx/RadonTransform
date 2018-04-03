/**
 * @file
 * @brief Image manipulation functions.
 * 
 * @note <b>Programs:</b>
 * * @ref iftExtractROI.c = Extract ROIs of an Image
 * * @ref iftInsertROI.c = Inserts/Copies an Image (Region Of Interest (ROI)) inside a Target Image.
 * * @ref iftExtractObject.c = Extract an Object from an Image
 * * @ref iftExtractAllObjects.c = Extract All Objects (Labels) from an Image
 * * @ref iftExtractImageROI.c = Extract the Min. Bounding Box (ROI) from an Image with all objects/voxels
 *                               (non-zero) inside it.
 * *
 * * @ref iftLabelImageAreaVolume.c = Computes the area/volume from each Object from an Label Image.  
 * * @ref iftSetVoxelSize.c = Overwrites the Pixel/Voxel Sizes of an Image or a Set of Images to some given Values.
 */


#ifndef IFT_IMAGE_H
#define IFT_IMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "iftAdjacency.h"
#include "iftBMap.h"
#include "iftColor.h"
#include "iftCommon.h"
#include "iftFIFO.h"
#include "iftFile.h"
#include "iftString.h"
#include "iftSort.h"
#include "iftSet.h"



/**
 * @brief macros to access a voxel/pixel <b>p</b> in the image <b>i</b>
 * @author Peixinho
 * @date Jun, 2016
 * @{
 */
#define iftImgVoxelVal(i, v) i->val[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCb(i, v) i->Cb[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCr(i, v) i->Cr[iftGetVoxelIndex(i, v)]
/**
 * @}
 */

/**
 * @brief macros to access a pixel in coordinates <b>(x, y)</b> in the image <b>i</b>
 * @author Peixinho
 * @date Jun, 2016
 * @{
 */
#define iftImgVal2D(i, x, y) i->val[((x)+(i)->tby[(y)])]
#define iftImgCb2D(i, x, y) i->Cb[((x)+(i)->tby[(y)])]
#define iftImgCr2D(i, x, y) i->Cr[((x)+(i)->tby[(y)])]

/**
 * @}
 */


/**
 * @brief macros to access a voxel in coordinates <b>(x, y,z)</b> in the image <b>i</b>
 * @author Peixinho
 * @date Jun, 2016
 * @{
 */


#define iftImgVal(img, x, y, z) img->val[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCb(img, x, y, z) img->Cb[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCr(img, x, y, z) img->Cr[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgAlpha(img, x, y, z) img->alpha[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
/**
 * @}
 */


typedef struct ift_nifti_header iftNIfTIHeader; // NIfTI header - it comes from external/libnifit

/**
 * @brief 3D Image Plane Orientations.
 * @author Samuel Martins
 * @date Mar 1, 2016
 * @ingroup Image
 *
 * @note AXIAL    = Plane XY
 * @note CORONAL  = Plane XZ
 * @note SAGITTAL = Plane YZ
 */
typedef enum {
    IFT_AXIAL, IFT_CORONAL, IFT_SAGITTAL
} iftImagePlaneOrientation;


/**
 * @brief Image Domain (Dimensions of an Image)
 * @author Samuel Martins
 * @date May 2, 2016
 * @ingroup Image
 */
typedef struct ift_image_domain {
    int xsize;
    int ysize;
    int zsize;
} iftImageDomain;


/**
 * @brief Sizes from Image Voxels (Dimensions of a Voxel)
 * @author Samuel Martins
 * @date May 2, 2016
 * @ingroup Image
 */
typedef struct ift_voxel_size {
    float dx;
    float dy;
    float dz;
} iftVoxelSize;


/**
 * @brief iftImage definition and prototypes.
 *
 * @author Falcao
*/
typedef struct ift_image {
    /** Brightness pixels array. */
    int *val;
    /** Blue component pixels array */
    ushort *Cb;
    /** Red component pixels array */
    ushort *Cr;
    /** alpha component pixels array */
    ushort *alpha;

    /** X axis size. */
    int xsize;
    /** Y axis size. */
    int ysize;
    /** Z axis size. */
    int zsize;

    /** X axis voxel size. */
    float dx;
    /** Y axis voxel size. */
    float dy;
    /** Z axis voxel size. */
    float dz;

    /** speed-up voxel access tables */
    int *tby, *tbz;
    /** Number of pixels. */
    int n; // number of voxels

    /** NIfTI Header (only for NIfTI images) */
    iftNIfTIHeader *nii_hdr;
} iftImage;


/**
 * @brief An array of Images.
 * @author Samuel Martins
 * @date Mar 1, 2016
 * @ingroup Image
 */
typedef struct ift_image_array {
    /** Number of Images (array size) */
    size_t n;
    /** Array of Images */
    iftImage **images;
} iftImageArray;


/**
 * @brief iftImageTiles definition. It stores the coordinates for partioning an image into tiles
 * relative to a given bounding box.
 * @author Falcao
 * @ingroup Image
 */
typedef  struct ift_image_tile {
    /** Bounding box of the tiles in the original coordinate system.
        The tiles are computed only for the given bounding box, which
        encompasses the entire image when necessary.
    */
    iftBoundingBox bb;

    /** Number of tiles in the x axis. */
    int ntiles_x;
    /** Number of tiles in the y axis. */
    int ntiles_y;
    /** Number of tiles in the z axis. */
    int ntiles_z;

    /* The total number of tiles, equal to ntiles_x*ntiles_y*ntiles_z */
    int ntiles;
    /** The  bounding boxes of each tile RELATIVE TO THE START OF THE ORIGINAL BOUNDING BOX.
        That is, the first tile will have coordinates (0,0,0) relative to voxel bb.min,
        which may or may not coincide with the first voxel of the original image.
     */
    iftBoundingBox *tile_coords;
} iftImageTiles;

/**
 * @brief macros to get a specific voxel coordinate from the image/fimage/etc <b>s</b> with  voxel index <b>p</b>
 * @author Deangeli
 * @date May 12, 2016
 * @{
 */
#define iftGetXCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftGetYCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftGetZCoord(s, p) ((p) / (((s)->xsize)*((s)->ysize)))
/**
 * @}
 */


/**
 * @brief computes an index for a voxel <b>v</b> from the image/fimage/etc <b>s</b>
 * @author Deangeli
 * @date May 12, 2016
 */
#define iftGetVoxelIndex(s, v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])
#define iftDiagonalSize(s) (iftRound(sqrtf(s->xsize*s->xsize + s->ysize*s->ysize + s->zsize*s->zsize)))


/**
 * @brief Checks if the input image <b>img</b> has isotropic resolution(resolution identical in all dimensions)
 * @author Deangeli
 * @date May 19, 2016
 */
#define iftIsIsotropic(img) ( \
    ( iftIs3DImage(img) && iftAlmostZero((img)->dx-(img)->dy) && iftAlmostZero((img)->dx-(img)->dz) ) || \
    ( !iftIs3DImage(img) && iftAlmostZero((img)->dx-(img)->dy) ) \
)


/**
 * @brief Checks if the images <b>img1</b> and <b>img2</b> have the same number of voxels.
 * @author Deangeli
 * @date May 19, 2016
 */
#define iftIsVoxelSizeEqual(img1, img2) ((iftAlmostZero((img1)->dx - (img2)->dx)) && \
                                         (iftAlmostZero((img1)->dy - (img2)->dy)) && \
                                         (iftAlmostZero((img1)->dz - (img2)->dz)))

/**
 * @brief Checks if the images <b>img1</b> and <b>img2</b> are in the same domain.
 * @author Deangeli
 * @date May 19, 2016
 */
#define iftIsDomainEqual(img1, img2) (((img1)->xsize == (img2)->xsize)  && \
                                      ((img1)->ysize == (img2)->ysize) && \
                                      ((img1)->zsize == (img2)->zsize)) 

/** @brief Gets the Domain of an Image */
#define iftGetImageDomain(img) ((iftImageDomain){.xsize = (img)->xsize, .ysize = (img)->ysize, .zsize = (img)->zsize})


#define iftMaxImageRange(img_depth) ((1 << (img_depth)) - 1) // 2^img_depth -1

/**
 * @brief Image loader Prototype.
 *
 * A common interface between loading images functions. Can be used to receive any image loading function as a parameter.
 *
 */
typedef iftImage *(*iftImageLoader)(const char *filename, ...);


/**
 * @brief Verifies if two images have the same domain, aborting the program if the image domains are different.
 *
 * @param img1 First Image.
 * @param img2 Second Image.
 * @param function Name of the function that calls this function. It is used in the Error Message if the image domains are different. 
 */
void iftVerifyImageDomains(const iftImage *img1, const iftImage *img2, const char *function);


/**
 * @brief Verifies the domain and voxel sizes from two images, aborting the program in case of divergence.
 * @author Samuel Martins
 * @date Jul 14, 2014
 * 
 * @param img1              First Image.
 * @param img2              Second Image.
 * @param external_function External Function that is calling this one.
 */
void iftVerifyImages(const iftImage *img1, const iftImage *img2, const char *external_function);


/** Gets the depth (in bits) of an image */
int iftImageDepth(const iftImage* img);


/**
 * @brief Check if the image has color information.
 * @param img Image to be checked.
 * @return True for color images, False for grayscale images.
 */
bool iftIsColorImage(const iftImage *img);

/**
 * @brief Check if the image has 3D information.
 * @param img Image to be checked.
 * @return 1 for 3D images, 0 for 2D images.
 */
bool iftIs3DImage(const iftImage *img);

/**
 * @brief Image X axis size.
 * @param img Image.
 * @return X axis size.
 */
int iftXSize(iftImage *img);

/**
 * @brief Image Y axis size.
 * @param img Image.
 * @return Y axis size.
 */
int iftYSize(iftImage *img);

/**
 * @brief Image Z axis size.
 * @param img Image.
 * @return Z axis size.
 */
int iftZSize(iftImage *img);

/**
 * @brief Destroys a <b>Py</b> image.
 * @author Deangeli
 * @date May 12, 2016
 * @param img <b>Py</b> image.
 *
 */
void iftDestroyPyImage(iftImage *img);

/**
 * @brief Gets the correspondent voxel to index p.
 *
 * @param p Voxel index.
 * @return Voxel coordinates from index p.
 */
iftVoxel iftGetVoxelCoord(const iftImage *img, int p);

/**
 * @brief Creates a new iftImage, with (xsize, ysize, zsize) dimensions. This function must be a wrapper around
 * iftCreateImageFromBuffer that allocates the brightness array.
 *
 * For color images, check iftCreateColorImage()
 *
 * @param xsize X axis dimension
 * @param ysize Y axis dimension
 * @param zsize Z dimension
 * @return The new iftImage.
 */
iftImage *iftCreateImage(int xsize, int ysize, int zsize);


/**
 * @brief Creates a new iftImage with the header (domain, voxel sizes, nifti header, ...) from
 * the image src.
 *
 * @note If the source is image is colored, the output image will be too (the Cb and Cr are also allocated).
 * 
 * @param  src Source Image used to create a new image.
 * @return     The new Image.
 */
iftImage *iftCreateImageFromImage(const iftImage *src);


/**
 * @brief Creates an image from a pre-allocated buffer whose size must be xsize*ysize*zsize
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param xsize X axis dimension
 * @param ysize Y axis dimension
 * @param zsize Z dimension
 * @param val   A buffer with size xsize*ysize*zsize
 * @return The new iftImage.
 */
iftImage  *iftCreateImageFromBuffer(int xsize,int ysize,int zsize, int *val);


/**
 * @brief Creates a new colored iftImage, with (xsize, ysize, zsize) dimensions.
 *
 * @param xsize X axis dimension
 * @param ysize Y axis dimension
 * @param zsize Z dimension
 * @param depth image depth
 * @return The new iftImage.
 */
  iftImage *iftCreateColorImage(int xsize, int ysize, int zsize, int depth);


/**
 * @brief Creates an array of Images.
 * @author Samuel Martins
 * @date Mar 1, 2016
 * @ingroup Image
 * 
 * @param n The number of Images (array size).
 * @return The array of Images
 */
iftImageArray *iftCreateImageArray(size_t n);


/**
 * @brief Destroys the image.
 * @warning Receives a pointer to an image, not the image itself.
 *
 * @param img A pointer to the image to be destroyed.
 */
void iftDestroyImage(iftImage **img);


/**
 * @brief Destroys an Image Array deallocating their images too.
 * @author Samuel Martins 
 * @date Mar 1, 2016
 * @ingroup Image
 * @warning If some image from the array is actually NULL, it is ignored by the Destroyer.
 */
void iftDestroyImageArray(iftImageArray **img_arr);


/**
 * @brief Copies the voxel size from <b>src</b> to <b>dst</b>.
 *
 * @param src Source image.
 * @param dst Destination image.
 */
#define iftCopyVoxelSize(src, dst) (dst)->dx = (src)->dx; (dst)->dy = (src)->dy; (dst)->dz = (src)->dz;

/**
 * @brief Copies the blue and red (YCbCr) components from source image <src> to destinantion image <dst>.
 *
 * It only copies if src is a colored image.
 * 
 * @param src Source image.
 * @param dst Destination image.
 */
void iftCopyCbCr(const iftImage *src, iftImage *dst);

/**
 * @brief Set the image blue and red (YCbCr) components to a certain value.
 *
 * @param img Target image.
 * @param value blue and red component value.
 *
 * @warning Set Cb and Cr with 128 to ensure that the image will have black background (i.e., 128 is equivalent to 0 in RGB)
 */
void iftSetCbCr(iftImage *img, ushort value);

/**
 * @brief Set the image blue,red and alpha (YCbCrA) components to a certain value.
 *
 * @param img Target image.
 * @param value blue, red and alpha component value.
 *
 * @warning Set Cb and Cr with 128 to ensure that the image will have black background (i.e., 128 is equivalent to 0 in RGB)
 */
void  iftSetCbCrAlpha(iftImage *img, ushort value);
void  iftSetAlpha(iftImage *img, ushort value);

/**
 * @brief Set the image blue (YCbCr) components to a certain value.
 *
 * @param img Target image.
 * @param value blue component value.
 */
void iftSetCb(iftImage *img, ushort value);

/**
 * @brief Set the image red (YCbCr) components to a certain value.
 *
 * @param img Target image.
 * @param value red component value.
 */
void iftSetCr(iftImage *img, ushort value);

/**
 * @brief Check if the voxel belongs to the image domain.
 *
 * @param img Target image.
 * @param v Voxel coordinates.
 */
#define iftValidVoxel(img, v)  ((v.x >= 0) && (v.x < img->xsize) && (v.y >= 0) && (v.y < img->ysize) && (v.z >= 0) && (v.z < img->zsize))


/**
 * @brief Gets the maximum intensity in a region of the image.
 * @author Samuel Martins
 * @date Jun 16, 2016
 * 
 * @param img Target image.
 * @param bb  Target region (Bounding Box).
 * @return    Maximum intensity in the region.
 */
int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb);

/**
 * @brief Gets the maximum intensity in a region of the image. If the
 * region is not provided (mask is NULL), then it uses the entire
 * image domain.
 * @author Alexandre Falcao
 * @date Nov 16th, 2016
 * 
 * @param  img   Target image.
 * @param  mask  Target region.
 * @return Maximum intensity in the region.
 */
int iftMaximumValueInMask(const iftImage *img, iftImage *mask);

/**
 * @brief Gets the minimum intensity in a region of the image. If the
 * region is not provided (mask is NULL), then it uses the entire
 * image domain.
 * @author Alexandre Falcao
 * @date Nov 16th, 2016
 * 
 * @param  img   Target image.
 * @param  mask  Target region.
 * @return Minimum intensity in the region.
 */
int iftMinimumValueInMask(const iftImage *img, iftImage *mask);

/**
 * Gets the maximum intensity in the image.
 * @param img Target image.
 * @return Maximum intensity.
 */
int iftMaximumValue(const iftImage *img);

/**
 * @brief Computes the mean intensity in a given region of the image, as
 * provided by a mask. It uses the entire image when the mask is not
 * provided.
 * @author Alexandre Falcao   
 * @date   Nov 10th, 2016
 * @param  img  Input image.  
 * @param  mask Input region (or NULL).  
 * @return Mean intensity.
 */
  float iftMeanValue(const iftImage *img, iftImage *mask);

/**
 * @brief Computes the most frequent intensity (mode) in a given
 * region of the image, as provided by a mask. It uses the entire
 * image when the mask is not provided.  
 * @author Alexandre Falcao
 * @date   Nov 10th, 2016
 * @param  img Input image.  
 * @param  mask Input region (or NULL).
 * @return The most frequent intensity.
 */
  int iftModeValue(const iftImage *img, iftImage *mask);

/**
 * Gets the maximum blue (YCbCr) value in the image.
 * @param img Target image.
 * @return Maximum blue value.
 */
int iftMaximumCb(const iftImage *img);

/**
 * Gets the maximum red (YCbCr) value in the image.
 * @param img Target image.
 * @return Maximum red value.
 */
int iftMaximumCr(const iftImage *img);

/**
 * Gets the minimum intensity in the image.
 * @param img Target image.
 * @return Minimum intensity.
 */
int iftMinimumValue(const iftImage *img);

/**
 * Gets the minimum blue (YCbCr) value in the image.
 * @param img Target image.
 * @return Minimum blue value.
 */
int iftMinimumCb(const iftImage *img);

/**
 * Gets the minimum red (YCbCr) value in the image.
 * @param img Target image.
 * @return Minimum red value.
 */
int iftMinimumCr(const iftImage *img);

/**
 * Computes the maximum image value around a voxel <p> within a given adjacency, which includes <p>'s value.
 *
 * @author Thiago Vallin Spina
 * @date Jan 25, 2016
 *
 * @param img The input image.
 * @param p The center voxel.
 * @param A The adjacency relation.
 *
 * @return The maximum image value.
 */
int iftMaximumValueInAdjacency(iftImage *img, int p, iftAdjRel *A);

/**
 * Computes the minimum image value around a voxel <p> within a given adjacency, which includes <p>'s value.
 *
 * @author Thiago Vallin Spina
 * @date Jan 25, 2016
 *
 * @param img The input image.
 * @param p The center voxel.
 * @param A The adjacency relation.
 *
 * @return The minimum image value.
 */
int iftMinimumValueInAdjacency(iftImage *img, int p, iftAdjRel *A);

/**
 * Computes the median image value around a voxel <p> within a given adjacency, which includes <p>'s value.
 *
 * @author Thiago Vallin Spina
 * @date Jan 26, 2016
 *
 * @param img The input image.
 * @param p The center voxel.
 * @param A The adjacency relation.
 *
 * @return The median image value.
 */
int iftMedianValueInAdjacency(iftImage *img, int p, iftAdjRel *A);

/**
 * @brief Gets an Extension for the input image: ".scn" for 3D, ".ppm" for 2D color, and ".pgm" for 2D gray
 * @param  img Input Image
 * @return     Image Extension.
 *
 * @author samuka
 * @date Dec 28, 2016
 */
char *iftImageExt(const iftImage *img);


/**
 * @brief Reads a SCN image from disk. Implements iftImageLoader()
 *
 * @param filename The image path.
 * @return The loaded image.
 */
iftImage *iftReadImage(const char *filename, ...);

/**
 * @brief Reads a NIfTI image from disk.
 *
 * @param filename The image path.
 * @return The loaded image.
 *
 * @author Samuka Martins
 * @date Jun 20, 2017
 */
iftImage *iftReadImageNIfTI(const char *filename, ...);


/**
 * @brief Reads an Analyze (.hdr or .img) image from disk.
 *
 * Analyze adopts the coordinate system convention of Computer Graphics, that is,
 * the coordinate origin is in the lower left corner.
 * Thus, with the subject lying supine, the coordinate origin is on the right side of the body (x),
 * at the back (y), and at the feet (z).
 * This is the radiological convention LAS+: Right to Left, Posterior to Anterior, Inferior to Superior.
 * Moreover, the header file (.hdr) stores a field which shows the orientation (hdr->hist.orient)
 * in which the data are stored in the data file (*.img).
 *
 * libIFT adopts the coordinate system convention of Image Processing, that is,
 * the coordinate origin is in the upper left corner.
 * It also considers the subject lying on a stretcher, but the coordinate origin has a difference:
 * it is at the face (y).
 * Its convention is LPS+: Right to Left, Anterior to Posterior, Inferior to Superior.
 *
 * @note This functions was based on description: http://eeg.sourceforge.net/AnalyzeDirect_ANALYZEformat.pdf
 *
 * @param filename The image path.
 * @return The loaded image.
 *
 * @author Samuka Martins
 * @date Jul 18, 2017
 */
iftImage *iftReadImageAnalyze(const char *filename, ...);


/**
 * @brief Read an image from disk according to the extension type (Accepts PGM and PPM). Implements iftImageLoader()
 * @param filename Image path.
 * @return The loaded image.
 */
iftImage *iftReadImageByExt(const char *filename, ...);

/**
 * @brief Writes an image into disk according to the extension type (Accepts PGM and PPM).
 * @param filename Image path.
 * @param The image to be written.
 *
 * @ingroup Image
 */
void iftWriteImageByExt(const iftImage *img, const char *filename, ...);

/**
 * @brief Writes an image into disk.
 * @param img The image to be written.
 * @param filename to image path.
 *
 * @ingroup Image
 */
void iftWriteImage(const iftImage *img, const char *filename, ...);

/**
 * @brief Writes a NIfTI image into disk.
 * @param img The image to be written.
 * @param filename Image Filename (it must have .nii)
 * @ingroup Image
 *
 * @author Samuka Martins
 * @date Jun 21, 2017
 *
 * @exception If the field nii_hdr from the img is NULL, nothing is done
 */
void iftWriteImageNIfTI(const iftImage *img, const char *filename, ...);

/**
 * @brief Writes an Analyze image (.hdr and .img) into disk.
 * 
 * @param img The image to be written.
 * @param filename Image Filename (.hdr or .img)
 * @ingroup Image
 *
 * @author Samuka Martins
 * @date Jul 18, 2017
 */
void iftWriteImageAnalyze(const iftImage *img, const char *filename, ...);

/**
 * @brief Reads a compressed SCN image in GZip format from disk. The expected extensions are .gz and .zscn
 *
 * @author Thiago Vallin Spina
 * @date Mar 31, 2016
 *
 * @param filename The image path.
 * @return The loaded image
 *
 * @ingroup Image
 */
iftImage *iftReadImageGZip(const char *filename, ...);

/**
 * @brief Writes a 3D image into disk using GZip compression. The standard extensions should be .gz and .zscn.
 *
 * @author Thiago Vallin Spina
 * @date Mar 31, 2016
 *
 * @param img The image to be written
 * @param filename Image path
 * @note The image header is saved by placing all info in a single line to facilitate backwards compatibility with iftReadImage,
 * for the case when the image is decompressed with an external tool.
 *
 * @ingroup Image
 */
void iftWriteImageGZip(const iftImage *img, const char *filename, ...);

/**
 * @brief Read a PGM (P5 format) image from disk. Implements iftImageLoader()
 * @param filename Image path.
 * @return The loaded image.
 *
 * @ingroup Image
 */
iftImage *iftReadImageP5(const char *filename, ...);

/**
 * @brief Writes a PGM image into disk (P5 format).
 * @param Filename to image path.
 * @param The image to be written.
 */
void iftWriteImageP5(const iftImage *img, const char *filename, ...);

/**
 * @brief Converts the image to a PGM (P5 format) file and reads the image from disk. Implements iftImageLoader()
 * @warning This function is not thread safe, be carefull.
 * @param filename Image path.
 * @return The loaded image.
 */
iftImage *iftReadImageAsP5(const char *filename, ...);

/**
 * @brief Writes an image into disk according to the extension format
 * @author Deangeli
 * @date May 12, 2016
 * @param filename Image path.
 */
void iftWriteImageExtFormat(iftImage *img, const char *filename, ...);

/**
 * @brief Converts the image to a PGM (P6 format) file and reads the image from disk. Implements iftImageLoader()
 * @warning This function is not thread safe, be carefull.
 * @param filename Image path.
 * @return The loaded image.
 */
iftImage *iftReadImageAsP6(const char *filename, ...);

/**
 * @brief Read a PGM (P6 format) image from disk. Implements iftImageLoader()
 * @param filename Image path.
 * @return The loaded image.
 */
iftImage *iftReadImageP6(const char *filename, ...);

/**
 * @brief Writes a PGM image into disk (P6 format).
 * @param Filename to image path.
 * @param The image to be written.
 */
void iftWriteImageP6(const iftImage *img, const char *filename, ...);

/**
 * @brief Read a PGM (P2 format) image from disk. Implements iftImageLoader()
 * @param filename Image path.
 * @return The loaded image.
 */
iftImage *iftReadImageP2(const char *filename, ...);

/**
 * @brief Writes a PGM image into disk (P2 format).
 * @param Filename to image path.
 * @param The image to be written.
 */
void iftWriteImageP2(const iftImage *img, const char *filename, ...);

/**
 * @brief Read PNG image from the disk.
 * @param format Path to PNG image.
 * @author Peixinho
 */
iftImage* iftReadImagePNG(const char* format, ...);


/**
 * @brief Read JPEG image from the disk.
 * @param format Path to JPEG image.
 * @author Deangeli
 */
iftImage* iftReadImageJPEG(const char* format, ...);

/**
 * @brief Read TIFF image from the disk.
 * @param format Path to TIFF image.
 * @author Deangeli
 */
//iftImage* iftReadImageTIFF(const char* format, ...);

/**
 * @brief Write a PNG image on the disk
 * @param format Path to store file
 * @author Peixinho
 */
void iftWriteImagePNG(const iftImage* img, const char* format, ...);

/**
 * @brief Write a JPEG image on the disk
 * @param format Path to store file
 * @author Deangeli
 */
void iftWriteImageJPEG(const iftImage* img, const char* format, ...);

/**
 * @brief Write a TIFF image on the disk
 * @param format Path to store file
 * @author Deangeli
 */
//void iftWriteImageTIFF(const iftImage* image, const char* format, ...);
/**
 * @brief Write a PNG image on the disk with alpha channel
 * @param format Path to store file
 * @date Dec 1, 2016
 * @author Thiago Vallin Spina
 */
void iftWriteImagePNGWithAlpha(const iftImage* img, const iftImage *alpha, const char* format, ...);

/**
 * @brief Extracts a Slice of a given Plane (AXIAL, CORONAL, SAGITTAL) from a Volumetric Image.
 * @author Samuel Martins
 * @date Mar 1, 2016
 * @ingroup Image
 *
 * @note See a demo in demo/Miscellaneous/iftExtractSliceFrom3DImage.c
 * 
 * AXIAL    = Plane XY \n
 * CORONAL  = Plane XZ \n
 * SAGITTAL = Plane YZ \n
 * 
 * @param  vol_img           Volumetric Image
 * @param  plane_orientation The Image Plane Orientation.
 * @param  slice             The (number of the) required slice.
 * @return                   The 2D Image with the extracted Slice.
 * 
 * @warning The function works with Colored Images.
 * @exception Image is not volumetric (3D).
 * @exception Slice is < 0 or > that the last slice of the plane orientation.
 */
iftImage *iftExtractSlice(const iftImage *vol_img, iftImagePlaneOrientation plane_orientation, long slice);


/**
 * @brief Creates a copy of an image.
 *
 * @param Source image.
 * @return The copy of the image.
 */
iftImage *iftCopyImage(const iftImage *img);


/**
 * @brief Copy the NIFfTIHeader from a src image to a dst image.
 * @param src Source Image.
 * @param dst Destination Image.
 *
 * @author Samuka Martins
 * @date Jun 21, 2017.
 */
void iftCopyNIfTIHeader(const iftImage *src, iftImage *dst);


/**
 * @brief Creates a copy of an image to another pre-allocated image.
 *
 * @param src Source image.
 * @param dest Destination image
 */
void iftCopyImageInplace(const iftImage *src, iftImage *dst);


/**
 * @brief Checks if two images are equals (in voxel sizes, domain, and voxel values)
 * @author Samuka
 * @date Jan 3, 2017
 */
bool iftIsImageEqual(const iftImage *img1, const iftImage *img2);


/**
 * @brief Creates an image with a solid cube inside.
 *
 * Creates a cuboid with value <b>val</b> inside an image with <perc> times the images
 * coordinates.
 *
 * @author Samuel Martins
 * @date Jun 15, 2016
 * 
 * @param xsize X axis size.
 * @param ysize Y axis size.
 * @param zsize Z axis size.
 * @param perc  Percentage of the Image that will correspond the cuboid's volume.
 * @param val   Value of the Cuboid.
 * @return The cuboid container image.
 */
iftImage *iftCreateCuboid(int xsize, int ysize, int zsize, float perc, int val);


/**
 * @brief Converts a CSV file to image file.
 * @author Deangeli
 * @date May 12, 2016
 *
 * @param filename CSV file path.
 * @return Image described by CSV file.
 *
 */
iftImage *iftCSVtoImage(const char *filename, ...);

/**
 * @brief Check if two voxels are adjacent according to a given adjacency.
 *
 * @param img The domain image.
 * @param A The considered adjacency.
 * @param u First image voxel.
 * @param v Second image voxel.
 * @return 1 if the voxels are adjacent, 0 otherwise.
 */
char iftAdjacentVoxels(iftImage *img, iftAdjRel *A, iftVoxel u, iftVoxel v);

/**
 * @brief Gets the image gradient magnitude (ignoring the gradient direction).
 *
 * @param img The target image.
 * @param A The adjacency relationship.
 * @return The gradient magnitude image.
 */
iftImage *iftImageGradientMagnitude(const iftImage *img, iftAdjRel *Ain);


/**
 * @brief Adds a frame to image
 *
 *
 * The <b>sz</b> indicates the thickness of frame in all axes (x,y,z). For example, if we have
 * a 2D image of 5x10 pixels, and  <b>sz</b> is equal to 2, it will results
 * in an output image of 9x14 pixels. The sketch below illustrates the output
 * image generated. The "x" are the pixels values of the image without frame and "V"
 * are the pixels values on frame.
 *
 *
 *
 *  Input                           Output
 *                                 VVVVVVVVVVVVVV
 *                                 VVVVVVVVVVVVVV
 *  xxxxxxxxxx                     VVxxxxxxxxxxVV
 *  xxxxxxxxxx                     VVxxxxxxxxxxVV
 *  xxxxxxxxxx     ----->          VVxxxxxxxxxxVV
 *  xxxxxxxxxx                     VVxxxxxxxxxxVV
 *  xxxxxxxxxx                     VVxxxxxxxxxxVV
 *                                 VVVVVVVVVVVVVV
 *                                 VVVVVVVVVVVVVV
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img Image/fimage/etc where the frame will be added
 * @param sz thickness of frame
 * @param value color component value for all pixel on frame
 *
 * @return Original image on the frame of sz x sz x sz   .
 *
 */
iftImage *iftAddFrame(iftImage *img, int sz, int value);

/**
 * @brief Removes a frame from  an image
 *
 * The <b>sz</b> indicates the thickness of frame (x,y,z). For example, if we have
 * an image of 9x14 pixels, and  <b>sz</b> is equal to 2, it will results
 * in an output image of 5x10 pixels. The sketch below illustrates the output
 * image generated. The "V" are the pixels values on the frame and, "x" are the non-removed  pixels
 * values.
 *
 *  Input                           Output
 * VVVVVVVVVVVVVV
 * VVVVVVVVVVVVVV
 * VVxxxxxxxxxxVV                   xxxxxxxxxx
 * VVxxxxxxxxxxVV                   xxxxxxxxxx
 * VVxxxxxxxxxxVV   ----->          xxxxxxxxxx
 * VVxxxxxxxxxxVV                   xxxxxxxxxx
 * VVxxxxxxxxxxVV                   xxxxxxxxxx
 * VVVVVVVVVVVVVV
 * VVVVVVVVVVVVVV
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param fimg float Image/fimage/etc where the frame will be added
 * @param sz thickness of frame
 * @param value color component value for all pixel on frame
 *
 * @return Image generated by removing a rectangle frame of  sz x sz x sz
 *
 */
iftImage *iftRemFrame(const iftImage *fimg, int sz);

/**
 * @brief Adds a 3D frame to an image
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img Image/fimage/etc where the frame will be added
 * @param sx thickness of frame on x-axis
 * @param sy thickness of frame on y-axis
 * @param sz thickness of frame on z-axis
 * @param value color component value for all pixel on frame
 *
 * @return Original image on the frame of sx x sy x sz
 */
iftImage *iftAddRectangularBoxFrame(iftImage *img, int sx, int sy, int sz, int value);

/**
 * @brief Removes a 3D frame to an image
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param fimg float Image/fimage/etc where the frame will be added
 * @param sx thickness of frame on x-axis
 * @param sy thickness of frame on y-axis
 * @param sz thickness of frame on z-axis
 *
 * @return Image generated by removing a rectangle frame of  sx x sy x sz
 */
iftImage *iftRemRectangularBoxFrame(const iftImage *fimg, int sx, int sy, int sz);

/**
 * @brief Set the image brightness according to value.
 *
 * @param img Target image.
 * @param value intensity.
 */
void iftSetImage(iftImage *img, int value);

/**
 * @brief Gets the slice XY from  a 3D image, given a Z coordinate.
 *
 * @param img Target image
 * @param zcoord Z axis coordinate.
 * @return The 2D image slice.
 */
iftImage *iftGetXYSlice(const iftImage *img, int zcoord);

/**
 * @brief Gets the slice ZX from  a 3D image, given a Y coordinate.
 *
 * @param img Target image
 * @param ycoord Y axis coordinate.
 * @return The 2D image slice.
 */
iftImage *iftGetZXSlice(iftImage *img, int ycoord);

/**
 * @brief  Gets the slice YZ from  a 3D image, given a X coordinate.
 *
 * @param img Target image
 * @param xcoord X axis coordinate.
 * @return The 2D image slice.
 */
iftImage *iftGetYZSlice(iftImage *img, int xcoord);

/**
 * @brief Inserts a 2D image as a XY slice in a 3D image in a specified Z coordinate.
 *
 * @param img Target image
 * @param slice Inserted 2D image slice
 * @param zcoord Z axis coord.
 */
void iftPutXYSlice(iftImage *img, iftImage *slice, int zcoord);

/**
 * @brief Inserts a 2D image as a ZX slice in a 3D image in a specified Y coordinate.
 *
 * @param img Target image
 * @param slice Inserted 2D image slice
 * @param ycoord Y axis coord.
 */
void iftPutZXSlice(iftImage *img, iftImage *slice, int ycoord);

/**
 * @brief Inserts a 2D image as a YZ slice in a 3D image in a specified X coordinate.
 *
 * @param img Target image
 * @param slice Inserted 2D image slice
 * @param xcoord X axis coord.
 */
void iftPutYZSlice(iftImage *img, iftImage *slice, int xcoord);

/**
 * @brief Switch data from <b>X</b> axis to <b>Z</b> axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the swapped data
 *
 */
iftImage *iftSwitchXByZ(iftImage *img);

/**
 * @brief Switch data from <b>Y</b> axis to <b>Z</b> axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the swapped data
 *
 */
iftImage *iftSwitchYByZ(iftImage *img);

/**
 * @brief Switch data from <b>X</b> axis to <b>Y</b> axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the swapped data
 *
 */
iftImage *iftSwitchXByY(iftImage *img);


/**
 * @brief Reflects the image data about x-axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the reflected data
 *
 */

iftImage *iftInvertX(iftImage *img);

/**
 * @brief Reflects the image data about y-axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the reflected data
 *
 */
iftImage *iftInvertY(iftImage *img);

/**
 * @brief Reflects the image data about z-axis
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img 3D reference image
 * @return image that contains the reflected data
 *
 */
iftImage *iftInvertZ(iftImage *img);

/**
 * @brief Get the closest voxel to the image geometric center. See also iftGeometricCenter()
 *
 * @param obj Target image.
 * @return The center voxel.
 */
iftVoxel iftGeometricCenterVoxel(iftImage *obj);

/**
 * @brief Get the image geometric center. See also iftGeometricCenterVoxel()
 *
 * @param obj Target image.
 * @return The center point
 */
iftPoint iftGeometricCenter(iftImage *obj);

/**
 * @brief Computes the diagonal size of an object (non zero voxels) inside the image.
 *
 * @param obj Object image
 */
int iftObjectDiagonal(iftImage *obj);

/**
 *
 * @author Deangeli
 * @date May 13, 2016
 *
 *
 * @warning provavelmente existe algum erro nesta função. Razão: hist[img->val[p]+img_min_val] pode acessar uma
 * posição de mémoria não valida.
 */
void iftGetDisplayRange(iftImage *img, int *lower, int *higher);

/**
 * @brief Gets an image composed by the Cb band.
 *
 * @param img Target image
 * @return The chroma component image.
 */
iftImage *iftImageCb(const iftImage *img);

/**
 * @brief Gets an image composed by the Cr band.
 *
 * @param img Target image
 * @return The red component image.
 */
iftImage *iftImageCr(const iftImage *img);

/**
 * @brief Gets an image composed by the red (RGB) component of another.
 *
 * First of all, the function converts the pixels values from YCbCr to RGB, then
 * the red component is collected.
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img target image
 * @return The red (RGB) component image.
 */
iftImage *iftImageRed(iftImage *img);

/**
 * @brief Gets an image composed by the green (RGB) component of another.
 *
 * First of all, the function converts the pixels values from YCbCr to RGB, then
 * the green component is collected.
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img target image
 * @return The green (RGB) component image.
 */
iftImage *iftImageGreen(iftImage *img);

/**
 * @brief Gets an image composed by the blue (RGB) component of another.
 *
 * First of all, the function converts the pixels values from YCbCr to RGB, then
 * the blue component is collected.
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img target image
 * @return The blue (RGB) component image.
 */
iftImage *iftImageBlue(iftImage *img);

/**
 * @brief Generates a greyscale image from colored image.
 *
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img target image
 * @return A greyscale image
 */
iftImage *iftImageGray(iftImage *img);

/*
 *@brief gets the hue channel component (HSV color space) from a image
 *
 * First of all, the function converts the pixels values from YCbCr to RGB color space, then
 * RGB to HSV color space, and finally the Hue component is collected.
 *
 * @author Deangeli
 * @date May 12, 2016
 * @param img target image
 * @return a vector that contains Hue channel components.
 * */
iftImage *iftImageHue(iftImage *img);


/*
 *@brief gets the saturation channel component (HSV color space) from a image
 *
 * First of all, the function converts the pixels values from YCbCr to RGB color space, then
 * RGB to HSV color space, and finally the Saturation component is collected.
 *
 * @author Deangeli
 * @date May 19, 2016
 * @param img target image
 * @return a vector that contains Saturation channel components.
 *
 * */
iftImage *iftImageSaturation(iftImage *img);

/*
 *@brief Creates a image which its pixels values are described by a gaussian distribution
 *
 *
 * The pixels values are described as
 * \f[
 * maxval * e^{ \frac{D(x,y,z)^2}{2 \sigma^2} }
 * \f]
 *
 * where \f$ D(x,y,z) \f$ is the distance between the image voxel and the mean
 *
 * \f[
 * D(x,y) =  \parallel  (x,y,z) - (\mu_{x}, \mu_{y}, \mu_{z} \parallel
 * \f]
 *
 *
 * @author Deangeli
 * @date May 19, 2016
 *
 * @param xsize number of pixels on x-axis
 * @param ysize number of pixels on y-axis
 * @param zsize number of pixels on z-axis
 * @param mean vector mean
 * @param stdev standard deviation
 * @param maxval scaling factor that multiplies the gaussian distribution
 *
 * @return Gaussian image
 *
 * */
iftImage *iftCreateGaussian(int xsize, int ysize, int zsize, iftVoxel mean, float stdev, int maxval);

iftImage *iftRegionBorders(iftImage *label, int value);


/*
 *@brief Creates an image that its histogram can be described by a Gaussian distribution
 *
 *
 * The histogram distribution it is described as
 * \f[
 * e^{-\frac{(I_i - \mu)^2}{2 \sigma^2} }
 * \f]
 *
 * where \f$ I_i \f$ is the range where the pixels values are contained
 *
 * @author Deangeli
 * @date May 19, 2016
 *
 * @param xsize number of pixels on x-axis
 * @param ysize number of pixels on y-axis
 * @param zsize number of pixels on z-axis
 * @param mean value where the Gaussian distributuin is centred.
 * @param stdev standard deviation
 * @param maxval maximum pixel value in the distribution
 *
 * @return Image which has a Gaussian histogram
 *
 * */
iftImage *iftCreateImageWithGaussianHistogram(int xsize, int ysize, int zsize, float mean, float stdev, int maxval);

/*
 *@brief Creates an image that its histogram can be described by a combination of two Gaussian distributions
 *
 * @author Deangeli
 * @date May 19, 2016
 *
 * The histogram distribution it is described as
 * \f[
 * e^{-\frac{(I_i - \mu_{1})^2}{2 \sigma_{1}^2} } + e^{-\frac{(I_i - \mu_{2})^2}{2 \sigma_{2}^2} }
 * \f]
 *
 * @param xsize number of pixels on x-axis
 * @param ysize number of pixels on y-axis
 * @param zsize number of pixels on z-axis
 * @param mean1 first distribution  mean
 * @param mean2 second distribution  mean
 * @param stdev1 first distribution standard deviation
 * @param stdev2 second distribution standard deviation
 * @param stdev standard deviation
 * @param maxval maximum pixel value in the distribution
 *
 * @return Image which has a histogram described a combination of two Gaussian distribution
 *
 * */

iftImage *iftCreateImageWithTwoGaussiansHistogram(int xsize, int ysize, int zsize, float mean1, float stdev1,
                                                  float mean2, float stdev2, int maxval);

/**
 * @brief Reads all slices between <b>first</b> and <b>last</b> given a 3D image basename.
 * @author Deangeli
 * @date jun 17, 2016
 *
 * Gets the slices XY of 3D image based on the Z coordinate, from <b>first</b> and <b>last</b>.
 *
 * @param basename basename
 * @param first number of the first desired slice in Z-axis
 * @param last number of the last desired slice in Z-axis.
 * @param xsize number of voxels in X-axis
 * @param ysize number of voxels in Y-axis
 * @param bits_per_voxe
 * @return a 3D image that contains the desired slices.
 */
iftImage *iftReadRawSlices(char *basename, int first, int last, int xsize, int ysize, int bits_per_voxel);



void iftWriteRawSlices(iftImage *img, char *basename);

/**
 * @brief Writes a 3D image as a raw scene.
 *
 * @author Thiago Vallin Spina@
 * @param img The input image
 * @param format Formatting string to compose the filename.
 */
void iftWriteRawScene(iftImage *img, const char *format, ...);

/**
 * @brief Reads a raw file as a 3D image.
 *
 * @author Thiago Vallin Spina
 * @param xsize number of voxels in X-axis
 * @param ysize number of voxels in Y-axis
 * @param zsize number of voxels in Z-axis
 * @param format Formatting string to compose the filename
 */
iftImage *iftReadRawScene(int xsize, int ysize, int zsize, int bits_per_voxel, const char *format, ...);

/**
 * @brief Writes a 3D image as a raw scene and appeands to the filename the image's header information.
 *
 * @author Thiago Vallin Spina@
 * @param img The input image
 * @param format Formatting string to compose the filename.
 */
void iftWriteRawSceneWithInfoOnFilename(iftImage *img, const char *format, ...);

/**
 * @brief Reads a raw file as a 3D image by parsing the image's header from its filename.
 *
 * @author Thiago Vallin Spina
 * @param format Formatting string to compose the filename
 */
iftImage *iftReadRawSceneWithInfoOnFilename(const char *format, ...);

/**
 * @brief Reads a set of 2D standard images from disk as a 3D volume, by placing the images as XY slices.
 *
 * @author Thiago Vallin Spina
 * @date May 30, 2016
 *
 * @param files A set of image files.
 *
 * @return A 3D volume with the slices.
 */
iftImage *iftReadSlices(const iftFileSet *files);

iftImage *iftExtractGreyObject(iftImage *image);

iftImage *iftExtractGreyObjectPos(iftImage *image, iftVoxel *pos);

iftImage *iftCrop2DImageByBorder(iftImage *img, int border);

iftImage *iftCropImage(iftImage *img, int dx, int dy, int dz);

/**
 * @brief Gets the interpolated voxel value in a given position.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 */
int iftImageValueAtPoint(const iftImage *img, iftPoint P);

/**
 * @brief Gets the 2D image interpolated pixel value in a given position.
 *
 * @warning This function is for 2D images.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 */
int iftImageValueAtPoint2D(const iftImage *img, iftPoint P);

/**
 * @brief Gets the 2D image interpolated pixel value in a given position according to its nearest neighbor.
 *
 * @warning This function is for 2D images.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 */
int iftImageValueAtPoint2DNearestNeighbor(iftImage *img, iftPoint P);

/**
 * @brief Gets the 3D image interpolated pixel value in a given position according to its nearest neighbor.
 *
 * @warning This function is for 3D images.
 *
 * @param img Target image.
 * @param p position to be interpolated.
 * @return Interpolated intensity.
 */
int iftImageValueAtPointNearestNeighbor(iftImage *img, iftPoint P);

/**
 * @brief Gets the object variance in each axis.
 *
 * Gets the object variance for all the axes (non zero voxels).
 *
 * @param img Target image
 * @return The object variance in each axis.
 */
iftVector iftObjectAxesVariance(iftImage *img);

/* Simple voxel sampling methods. See iftMImage.h(c) for more
   complex ones. They return a binary mask with the selected
   voxels. */

int iftNumberOfElements(iftImage *mask);


iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize);

iftImage *iftSelectRegionOfInterest(int xsize, int ysize, int zsize, iftVoxel uo, iftVoxel uf);

/**
 * Set an image voxel to a specified color (YCbCr).
 *
 * @param img Target image
 * @param p Voxel index
 * @param YCbCr color
 */
static inline void iftSetYCbCr(iftImage *img, int p, iftColor YCbCr) {
    img->val[p] = YCbCr.val[0];
    img->Cb[p] = YCbCr.val[1];
    img->Cr[p] = YCbCr.val[2];
}

// Fast functions for accessing/putting RGB values into images.
// For efficiency, they do not check if @param img is colored.
static inline void iftSetRGB(iftImage *img, int p, int R, int G, int B, int normalization_value) {
    iftColor RGB;
    RGB.val[0] = R;
    RGB.val[1] = G;
    RGB.val[2] = B;

    iftSetYCbCr(img, p, iftRGBtoYCbCr(RGB, normalization_value));
}

// Fast functions for accessing/putting RGB values into images.
// For efficiency, they do not check if @param img is colored.
static inline void iftSetRGB2(iftImage *img, int p, iftColor rgb, int normalization_value) {
    iftSetRGB(img, p, rgb.val[0], rgb.val[1], rgb.val[2], normalization_value);
}

static inline iftColor iftGetRGB(iftImage *img, int p, int normalization_value) {
    iftColor YCbCr;

    YCbCr.val[0] = img->val[p];
    YCbCr.val[1] = img->Cb[p];
    YCbCr.val[2] = img->Cr[p];

    return iftYCbCrtoRGB(YCbCr, normalization_value);
}

//Converts between color spaces. This function quantizes the LAB cspace between 0 and 255.
iftImage *iftConvertColorSpace(iftImage *image, char origin_cspace, char dest_cspace);


iftImage *iftLooseImage(iftImage *image, int xsize, int ysize, int zsize);

void iftCenterImages(iftImage *image1, iftImage *image2, iftImage **centeredImage1, iftImage **centeredImage2);

iftImage *iftColorTableToImage(iftColorTable *ct, int xsize, int ysize);

void iftTickColorTableImage(iftImage *img, float minval, float maxval, int nticks, const char *filename);

// By TVS
char *iftGetImageType(const char *filename, ...);


/**
 * @brief Shifts the values of all spels according to the equation: @p img->val[i] + shift_val.
 *
 * @note shift_val: positive --> Right Shift. \n
 * shift_val: negative --> Left Shift.
 *
 * @author Samuel Martins
 * @date November 17, 2015
 * @ingroup ImageConverters
 *
 * @param[in/out] image Image to be shiffted.
 * @param[in] shift_val Shift Value.
 */
void iftShiftImageInPlace(iftImage *img, int shift_val);


/**
 * @brief Fills the bounding box in the image <img> with a given value.
 * @author Samuka
 * @date Aug 15, 2017
 */
void iftFillBoundingBoxInImage(iftImage *img, iftBoundingBox bb, int value);


/**
 * @brief Returns the Minimum (Image) Bounding Box necessary to cover all Image Objects (non-zero pixels/voxels),
 * and returns its Geometric Center (if <b>gc_out != NULL</b>).
 * @author Samuel Martins
 * @date Mar 2, 2016
 * @ingroup Image
 * @note An example can be found in @ref iftExtractImageROI.c
 *
 * E.x: The internal contour is the resulting Minimum Bounding Box, and the character <b>X</b> is its
 * Geometric Center: \n
 * <pre>
 *  ________________ 
 * |                | 
 * |  ___________   | 
 * | |        111|  | 
 * | |        222|  | 
 * | |           |  | 
 * | |      X 333|  | 
 * | |222_____333|  | 
 * |                | 
 * |____________ ___|
 * </pre>
 *
 * @param img Target Image used to extract the Bounding Box.
 * @param gc_out Return the Geometric Center from the Bounding Box to this Reference. If it is NULL, nothing is returned to it.
 * @return The Minimum (Image) Bounding Box.
 *
 * @warning If <b>gc_out</b> is <b>NULL</b>, the Geometric Center is not Returned.
 * 
 * @sa iftMinObjectBoundingBox(), iftMinLabelsBoundingBox()
 */
iftBoundingBox iftMinBoundingBox(const iftImage *img, iftPoint *gc_out);



/**
 * @brief Returns the Minimum Bounding Box necessary to cover a given Object with label <b>obj_label</b>,
 * and returns its Geometric Center (if <b>gc_out != NULL</b>).
 * @author Samuel Martins
 * @date Mar 2, 2016
 * @ingroup Image
 * @note An example can be found in @ref iftExtractObject.c
 *
 * E.x: Gets the Minimum Bounding Box from the Object 2 \n
 * <pre>
 *  ________________ 
 * |  111           | 
 * |  111 _____     | 
 * |     | 222 |    | 
 * |     |22222|    |
 * |     |_222_|    |
 * |             555| 
 * |  33333         | 
 * |   333   44444  | 
 * |         44444  | 
 * |_________44444__|
 * </pre>
 *
 * @param img Target Image used to extract the Bounding Box.
 * @param obj_label The label (code) from the required object.
 * @param gc_out Return the Geometric Center from the Bounding Box to this Reference. If it is NULL, nothing is returned to it.
 * @return The Minimum Bounding Box for the required Object.
 *
 * @warning The Min. Bounding Box from the required Object can contain other objects inside.
 * @warning If <b>gc_out</b> is <b>NULL</b>, the Geometric Center is not Returned.
 * 
 * @sa iftMinBoundingBox(), iftMinLabelsBoundingBox()
 */
iftBoundingBox iftMinObjectBoundingBox(const iftImage *img, int obj_label, iftPoint *gc_out);


/**
 * @brief Returns an array with the Minimum Bounding Boxes of all objects of a Label Array,
 * and returns all their Geometric Centers (if gcs_out != NULL).\n
 * If an object <b>DOESN'T EXIST</b>, its Bounding Box and Geometric Center will have coordinates (-1,-1,-1).
 *
 * @param label_img Label Image.
 * @param labels Input Array with the considered labels.
 * @param gcs_out Return the Geometric Center from the Bounding Boxes. If NULL, nothing is returned to it.
 * @return An array with the Minimum Bounding Box for each considered label.
 * 
 * @author Samuel Martins
 * @date Dec 16, 2016
 * @ingroup Image
 * @note An example can be found in @ref iftExtractAllObjects.c
 */
iftBoundingBox *iftMinLabelsBoundingBox(const iftImage *img, const iftIntArray *labels, iftVoxel **gcs_out);



/**
 * @brief Gets all Minimum (Image) Bounding Boxes necessary to cover all Image Objects (non-zero pixels/voxels)
 * from a set of Images, and returns their Geometric Centers (if <b>gcs_out != NULL</b>).
 * @author Samuel Martins
 * @date Apr 5, 2016
 * @ingroup Image
 *
 * @param img_paths Pathnames from the Target Images used to extract the Bounding Boxes.
 * @param gcs_out Returns the Geometric Center from the Bounding Boxes. If it is NULL, nothing is returned to it.
 * @return The Minimum (Image) Bounding Boxes.
 *
 * @warning If <b>gcs_out</b> is <b>NULL</b>, the Geometric Centers are not Returned.
 * @sa iftMinBoundingBox()
 */
iftBoundingBox *iftAllMinBoundingBox(const iftFileSet *img_paths, iftPoint **gcs_out);


/**
 * @brief Gets all Minimum Bounding Boxes necessary to cover a given Object with label <b>obj_label</b>
 * from a set of Images, and returns their Geometric Centers (if <b>gcs_out != NULL</b>).
 * @author Samuel Martins
 * @date Apr 5, 2016
 * @ingroup Image
 *
 * @param img_paths Pathnames from the Target Images used to extract the Bounding Boxes.
 * @param obj_label The label (code) from the required object.
 * @param gcs_out Returns the Geometric Center from the Bounding Boxes. If it is NULL, nothing is returned to it.
 * @return Minimum (Image) Bounding Boxes.
 *
 * @warning If <b>gcs_out</b> is <b>NULL</b>, the Geometric Centers are not Returned.
 * @sa iftMinObjectBoundingBox()
 */
iftBoundingBox *iftAllMinObjectBoundingBox(const iftFileSet *img_paths, int obj_label, iftPoint **gcs_out);


/**
 * @brief Gets all Minimum Bounding Boxes for all objects from label array <b>obj_label</b>
 * 
 * Gets the Min. Bounding Box and Geo Center for each label from each Image.
 * It returns a matrix n_labels x n_imgs, such that mbbs[o][i] and gcs[o][i] is the bounding box and
 * geometric center, respectively, of the image i for the object o.
 * This matrix configuration will make it easier its use in other functions.
 * 
 * @param  label_paths Pathnames from the Target Images used to extract the Bounding Boxes.
 * @param  labels      Array with the labels from each object.
 * @param  gcs_out     Returns a matrix n_labels x n_imgs, such that gcs[o][i] is the geometric center of
 *                     the image i for the object o.
 * @return             Matrix n_labels x n_imgs, such that mbbs[o][i] is the bounding box of the image i for the object o.
 */
iftBoundingBox **iftAllMinLabelsBoundingBox(const iftFileSet *label_paths, const iftIntArray *labels,
                                            iftVoxel ***gcs_out);

/**
 * @brief Gets the a Bounding Box that fits all Bounding Boxes from the Input Array, that is,
 * the Maximum Bounding Box.
 * @author Samuel Martins
 * @date Apr 4, 2016
 * @ingroup Image
 * 
 * @param  boxes Array of Bounding Boxes.
 * @param  n_boxes Number of Bounding Boxes.
 * @return The bounding box that fits all input bounding boxes. 
 */
iftBoundingBox iftMaxBoundingBox(iftBoundingBox *bbs, size_t n_boxes);


/**
 * @brief Gets the central voxel from a bounding box. 
 * @author Samuka
 * @date Jan 5, 2017
 */
iftVoxel iftBoundingBoxCenterVoxel(iftBoundingBox bb);


/**
 * @brief Extract a Region Of Interest (ROI) of an Image from a Bounding Box.
 * @date Mar 2, 2016
 * @ingroup Image
 * @note An example is @ref iftExtractROI.c
 * 
 * @param  img Image whose ROI will be extracted.
 * @param  bb  Bounding Box of the ROI.
 * @return     An ROI Image.
 *
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function extracts what it is possible.
 * @warning The Size of the ROI Image is the the same of the Input Bounding
 */
iftImage *iftExtractROI(const iftImage *img, iftBoundingBox bb);

/**
 * @brief Extract the voxels with a given label inside a Region Of Interest (ROI) of an Image.
 * @date Mar 2, 2016
 * @ingroup Image
 * @note An example is @ref iftExtractROI.c
 *
 * @param  src_img   Source Image whose Object will be extract from a ROI.
 * @param  bb        Bounding Box of the ROI.
 * @param  obj_label Label of the Target Object.
 * @return           An ROI Image with only the Object with Label <b>obj_label</b>.
 *
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function extracts what it is possible.
 * @warning The Size of the ROI Image is the the same of the Input Bounding
 */
iftImage *iftExtractObjectInsideROI(const iftImage *src_img, iftBoundingBox bb, int obj_label);


/**
 * @brief Extract the voxels of the objects with labels <labels int array> inside a Region Of Interest (ROI) of an Image.
 * Voxel with other labels will have label 0.
 * @date Jul 31, 2017
 * @ingroup Image
 *
 * @param  src_img   Source Image whose Labels will be extract from a ROI.
 * @param  bb        Bounding Box of the ROI.
 * @param  labels    Array with the labels to be extracted.
 * @return           An ROI Image with only the Object with Label <b>obj_label</b>.
 */
iftImage *iftExtractLabelsInsideROI(const iftImage *src_img, iftBoundingBox bb, const iftIntArray *labels);


/**
 * @brief Inserts an Image (Region Of Interest ROI) inside a Target Image from an initial position.
 * @date Mar 3, 2016
 * @ingroup Image
 * @note An example is @ref iftInsertROI.c
 * 
 * @param  roi Image (ROI) to be inserted.
 * @param  target Target Image where the ROI will be inserted.
 * @param  begin  Initial (coordinate) position where the ROI will be inserted..
 *
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function inserts what it is possible.
 */
void iftInsertROI(const iftImage *roi, iftImage *target, iftVoxel begin);


/**
 * @brief Inserts an Image (Region Of Interest ROI) inside a Target Image by aligning their Center Points.
 * @date Mar 14, 2016
 * @ingroup Image
 * 
 * @param  roi Image (ROI) to be inserted.
 * @param  target Target Image where the ROI will be inserted.
 *
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function inserts what it is possible.
 */
void iftInsertROIByCenter(const iftImage *roi, iftImage *target);

/**
 * @brief Inserts an Image (Region Of Interest ROI) inside a Target Image by aligning the given positions.
 * @date Nov 16, 2016
 * @author Thiago Vallin Spina
 * @ingroup Image
 *
 * @param  roi Image (ROI) to be inserted.
 * @param  roi_pos Position in the ROI image used as reference for pasting on to target.
 * @param  target Target Image where the ROI will be inserted.
 * @param  target_pos Position in the Target image used as reference for pasting ROI.
 *
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function inserts what it is possible.
 */
void iftInsertROIByPosition(const iftImage *roi, const iftVoxel roi_pos, iftImage *target, const iftVoxel target_pos);

/**
 * @brief Copies only the values from the Source Image ROI returning an Image with same Dimensions.
 * @author Samuel Martins 
 * @date Mar 5, 2016
 * @ingroup Image
 * 
 * @param  src Source Image where the ROI will be copied.
 * @param  bb  Bounding Box of the ROI to be copied.
 * @return     The Image (with same dimensions of the Source) with only values of ROI copied.
 * 
 * @warning If the ROI Image doesn't fit entirely inside the Target Image, the function copies what it is possible.
 */
iftImage *iftCopyImageROI(const iftImage *src, iftBoundingBox bb);


/**
 * @brief Extracts the object of label <b>obj_label</b> of an Image.
 * @date May 17, 2016
 * @ingroup Image
 *
 * @param  src Source Image.
 * @param  obj_label Label of the Target Object.
 * @return           A Label Image with only the required object.
 *
 * @warning If the Source Image does not have the required image, an empty image (only with the background)
 * will be returned.
 */
iftImage *iftExtractObject(const iftImage *src, int obj_label);


/**
 * @brief Extracts the labels from int array <labels> from an Image.
 * @date Jul 31, 2017
 * @ingroup Image
 *
 * @param  src_img Source Image.
 * @param  labels  Label of the Target Object.
 * @return         A Label Image with the required objects (labels).
 */
iftImage *iftExtractLabels(const iftImage *src_img, const iftIntArray *labels);


/**
 * @brief Counts the Number of Spels (Pixels/Voxels) of a given Object with label <b>obj_label</b>.
 * @author Samuel Martins
 * @date Mar 5, 2016
 * @ingroup Image
 * @note See a demo in @ref iftLabelImageAreaVolume.c
 *  
 * @param  label     Input Label Image.
 * @param  obj_label Label of the required Object.
 * @return           The number of spels (pixels/voxels) of the require Object found in the Input Label Image.
 * 
 * @exception <b>obj_label</b> is < 0 or > its Number of Labels (maximum value of the Label Image)
 */
int iftCountObjectSpels(const iftImage *label, int obj_label);


/**
 * @brief Counts the Number of Spels (Pixels/Voxels) of a given Object with label <b>obj_label</b> 
 * inside a Bounding Box.
 * @author Samuel Martins
 * @date Apr 17, 2016
 * @ingroup Image
 *  
 * @param  label_img Input Label Image.
 * @param  obj_label Label of the required Object.
 * @param bb Bounding Box.
 * @return           The number of spels (pixels/voxels) of the require Object found in the Input Label Image.
 */
int iftCountObjectSpelsFromBoundingBox(const iftImage *label_img, int obj_label, iftBoundingBox bb);


/**
 * @brief Counts the Number of Spels (Pixels/Voxels) of each Object (Label) in the Label Image with 
 * labels in the range [0..n_objects].
 * @author Samuel Martins
 * @date Mar 5, 2016
 * @ingroup Image
 * @note See a demo in @ref iftLabelImageAreaVolume.c
 * 
 * @param  label Input Label Image.
 * @return An int array with positions [0..n_labels] of the number of spels (pixels/voxels) from each Object (Label) in the Input Label Image.
 *
 * @warning The size of the array is n_labels+1, where the position 0 is the background.
 */
iftIntArray *iftCountLabelSpels(const iftImage *label);


/**
 * @brief Returns the Area/Volume (in mm²/mm³) of a given Object with label <b>obj_label</b> from a Label Image.
 * @author Samuel Martins
 * @date Mar 5, 2016
 * @ingroup Image
 * @note See a demo in @ref iftLabelImageAreaVolume.c
 * 
 * @param  label     Input Label Image.
 * @param  obj_label Label of the required Object.
 * @return           The area/volume (in mm²/mm³) of the require Object found in the Label Image.
 *
 * @note To compute the pixel/voxel in mm²/mm³, the function uses the spel displacement values from the Label Image.
 * @exception <b>obj_label</b> is < 0 or > its Number of Labels (maximum value of the Label Image)
 */
float iftAreaVolumeOfObject(const iftImage *label, int obj_label);


/**
 * @brief Returns the Area/Volume (in mm²/mm³) of a given Object with label <b>obj_label</b> from a
 * Bounding Box.
 * @author Samuel Martins
 * @date Apr 17, 2016
 * @ingroup Image
 * 
 * @param  label_img     Input Label Image.
 * @param  obj_label Label of the required Object.
 * @param  bb Bounding Box where the volume will be computed.
 * @return           The area/volume (in mm²/mm³) of the require Object found inside the Bounding Box.
 *
 * @note To compute the pixel/voxel in mm²/mm³, the function uses the spel displacement values from the Label Image.
 * @exception <b>obj_label</b> is < 0 or > its Number of Labels (maximum value of the Label Image)
 */
double iftAreaVolumeOfObjectFromBoundingBox(const iftImage *label_img, int obj_label, iftBoundingBox bb);


/**
 * @brief Returns the Area/Volume (in mm²/mm³) of each Object (Label) from a Label Image.
 * @author Samuel Martins
 * @date Mar 5, 2016
 * @ingroup Image
 * @note See a demo in @ref iftLabelImageAreaVolume.c
 * 
 * @param  label Input Label Image.
 * @return       A float array with positions [0..n_labels] of the volumes in mm (millimeters) of each Object (Label) found in the Label Image.
 *
 * @note To compute the pixel/voxel in mm²/mm³, the function uses the spel displacement values from the Label Image.
 * @warning The size of the array is n_labels+1, where the position 0 is the background.
 */
iftFloatArray *iftAreaVolumeOfLabels(const iftImage *label);


void iftInsertObject(iftImage *bin, iftImage *label, int obj_code, iftVoxel pos);


/**
 * @brief Shifts an Image according to given displacements.
 * @author Samuel Martins
 * @date Jul 5, 2016
 * 
 * @param  img Image to be Shifted.
 * @param  dx  Displacement in x-axis.
 * @param  dy  Displacement in y-axis.
 * @param  dz  Displacement in z-axis.
 * @return     The resulting Shifted Image.
 */
iftImage *iftShiftImage(const iftImage *img, int dx, int dy, int dz);


/**
 * @brief Gets the X tile coordinate from the tile at given index.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles An iftImageTiles record.
 * @param i     The index of the tile.
 * @return The X tile coordinate.
 */
#define iftGetXTileCoord(tiles, i) (((i) % (((tiles)->ntiles_x)*((tiles)->ntiles_y))) % (tiles)->ntiles_x)
/**
 * @brief Gets the Y tile coordinate from the tile at given index.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles An iftImageTiles record.
 * @param i     The index of the tile.
 * @return The Y tile coordinate.
 */
#define iftGetYTileCoord(tiles, i) (((i) % (((tiles)->ntiles_x)*((tiles)->ntiles_y))) / (tiles)->ntiles_x)

/**
 * @brief Gets the Z tile coordinate from the tile at given index.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles An iftImageTiles record.
 * @param i     The index of the tile.
 * @return The Z tile coordinate.
 */
#define iftGetZTileCoord(tiles, i) ((i) / (((tiles)->ntiles_x)*((tiles)->ntiles_y)))

/**
 * @brief Gets the index of the tile given a set of x,y,z coordinates passed as an iftVoxel.
 *
 * @param tiles An iftImageTiles record.
 * @param x     The X coordinate of the tile.
 * @param y     The Y coordinate of the tile.
 * @param z     The Z coordinate of the tile.
 * @return The index of the tile.
 */
#define iftGetTileIndex(tiles, x, y, z) ((x) + (y)*(tiles)->ntiles_x + (z)*(tiles)->ntiles_x*(tiles)->ntiles_y)

/**
 * @brief Allocates an iftImageTiles record setting all of its members accordingly.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param ntiles_x The number of tiles in the X axis.
 * @param ntiles_y The number of tiles in the Y axis.
 * @param ntiles_z The number of tiles in the Z axis.
 * @param coords   The coordinates of the first voxel for each tile.
 * @return The allocated record.
 */
iftImageTiles *iftCreateImageTiles(iftBoundingBox *tile_coords, int ntiles_x, int ntiles_y, int ntiles_z,
                                   iftBoundingBox bb);

/**
 * @brief Deallocates an iftImageTiles record and nullifies the pointer.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles The address of a pointer to an iftImageTiles record.
 */
void iftDestroyImageTiles(iftImageTiles **tiles);

/**
 * @brief Reads an iftImageTiles record from a file where it was saved in binary form.
 *
 * @author Thiago V. Spina
 * @date June 04, 2016
 *
 * @param filename The path to the file to be read.
 * @return The read record.
 */
iftImageTiles* iftReadImageTiles(const char *filename);

/**
 * @brief Writes an iftImageTiles record to a file in binary form.
 *
 * @author Thiago V. Spina
 * @date June 04, 2016
 *
 * @param tiles    A pointer to an iftImageTiles record.
 * @param filename The path to the file to be written.
 */
void iftWriteImageTiles(const iftImageTiles *tiles, const char *filename);


/**
 * @brief Copies an Image Tiles struct.
 * @author Samuel Martins
 * @date Jun, 2, 2016
 */
iftImageTiles *iftCopyImageTiles(const iftImageTiles *tiles);


/**
 * @brief Computes the number of tiles that may be used to divide an image into roughly equally sized patches
 * with a maximum number of pixels per tile.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param xsize         The number of pixels in the original image in the X axis.
 * @param ysize         The number of pixels in the original image in the Y axis.
 * @param zsize         The number of pixels in the original image in the Y axis.
 * @param max_tile_size Maximum number of pixels that each tile may have.
 * @param ntiles_x      Pointer that returns the number of tiles in the X axis.
 * @param ntiles_y      Pointer that returns the number of tiles in the Y axis.
 * @param ntiles_z      Pointer that returns the number of tiles in the Z axis.
 */
void iftNumberOfEquallyDimensionedTilesWithGivenMaximumSize(int xsize, int ysize, int zsize,
                                                            unsigned long max_tile_size, int *ntiles_x,
                                                            int *ntiles_y,
                                                            int *ntiles_z);
/**
 * @brief Computes an iftImageTiles record given an image and a number of tiles for each image dimension.
 *
 * The dimensions of the image tiles are equally sized to approximately divide the image into regular patches.
 * In practice, this function uses iftComputeBoundingBoxImageTiles with the entire image as a bounding box.
 *
 * @author Thiago V. Spina
 * @date Feb 26, 2015
 *
 * @param img    An image that must be tiled..
 * @param ntiles_x The number of tiles in the X axis.
 * @param ntiles_y The number of tiles in the Y axis.
 * @param ntiles_z The number of tiles in the Z axis.
 * @return The allocated record.
 *
 * @warning The number of tiles may actually be less than ntiles_x*ntiles_y*ntiles_z since we prioritize the fact that
 * the patches must have approximately the same size. REFER TO tiles->ntiles_x/y/z when accessing the actual number of
 * computed tiles.
 *
 * @sa iftBoundingBoxImageTilesByEquallyDividingAxes
 */
iftImageTiles *iftImageTilesByEquallyDividingAxes(const iftImage *img, int ntiles_x, int ntiles_y, int ntiles_z);

/**
 * @brief Computes an iftImageTiles record given a bounding box and a number of tiles for each image dimension.
 *
 * The dimensions of the image tiles are equally sized to approximately divide the bounding box into regular patches.
 * The number of tiles in each dimension is adjusted to ensure that the images have the same size approximately.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param bb    A bounding box for a given image, which need not be specified.
 * @param ntiles_x The number of tiles in the X axis.
 * @param ntiles_y The number of tiles in the Y axis.
 * @param ntiles_z The number of tiles in the Z axis.
 *
 * @warning The number of tiles may actually be less than ntiles_x*ntiles_y*ntiles_z since we prioritize the fact that
 * the patches must have approximately the same size. REFER TO tiles->ntiles_x/y/z when accessing the actual number of
 * computed tiles.
 *
 * @return The allocated record.
 */
iftImageTiles *iftBoundingBoxImageTilesByEquallyDividingAxes(iftBoundingBox bb, int ntiles_x, int ntiles_y,
                                                             int ntiles_z);

/**
 * @brief This function computes a series of image tiles with a given size and stride steps, which may overlap.
 *
 * Each tile will be have a size exactly D = (tile_xsize, tile_ysize, tile_zsize) or less for those around the image's
 * border. The tiles overlap by D - S, where S is the stride step given by S = (xstride, ystride, zstride). If S > D then
 * the tiles will not overlap but will leave a space of D - S voxels among them. S must be at least (1,1,1) or the function
 * issues an error. This function uses iftStridedBoundingBoxImageTiles.
 *
 * @param tile_xsize The tile's size in the X axis.
 * @param tile_ysize The tile's size in the Y axis.
 * @param tile_zsize The tile's size in the Z axis.
 * @param xstride    The stride in the X axis.
 * @param ystride    The stride in the Y axis.
 * @param zstride    The stride in the Z axis.
 * @return The allocated record.
 * @sa iftBoundingBoxImageTilesByStriding
 */
iftImageTiles *iftImageTilesByStriding(iftImage *img, int tile_xsize, int tile_ysize, int tile_zsize,
                                       int xstride, int ystride, int zstride);

/**
 * @brief This function computes a series of tiles with a given size and stride steps, which may overlap, inside the
 * bounding box.
 *
 * Each tile will be have a size exactly D = (tile_xsize, tile_ysize, tile_zsize) or less for those around the bounding
 * box's border. The tiles overlap by D - S voxels, where S is the stride step given by S = (xstride, ystride, zstride).
 * If S > D then the tiles will not overlap but will leave a space of D - S voxels among them. S must be at least
 * (1,1,1) or the function issues an error. This function uses iftBoundingBoxImageTilesByStriding.
 *
 * @param tile_xsize The tile's size in the X axis.
 * @param tile_ysize The tile's size in the Y axis.
 * @param tile_zsize The tile's size in the Z axis.
 * @param xstride    The stride in the X axis.
 * @param ystride    The stride in the Y axis.
 * @param zstride    The stride in the Z axis.
 * @return The allocated record.
 */
iftImageTiles *iftBoundingBoxImageTilesByStriding(iftBoundingBox bb, int tile_xsize, int tile_ysize,
                                                  int tile_zsize, int xstride, int ystride, int zstride);


/**
 * @brief Returns the tile index of a voxel from the original image domain.
 *
 * The tile corresponds to the first one that contains the voxel of interest, since image tiles may overlap themselves
 * causing regions with voxels intersecting multiple tiles at once.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles A pointer to an iftImageTiles record.
 * @param v     A voxel in the coordinate space of the original image.
 * @return The tile index of voxel v or NIL if the voxel is outside the original bounding box.
 */
int iftGetIndexFromFirstTileIntersectingVoxel(iftImageTiles *tiles, iftVoxel v);


/**
 * @brief Returns the indices of all tiles that contain a voxel from the original image domain.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles A pointer to an iftImageTiles record.
 * @param v     A voxel in the coordinate space of the original image.
 * @return The tile indices intersecting voxel v or NULL if the voxel is outside the original bounding box.
 */
iftSet*iftGetIndicesFromAllTilesIntersectingVoxel(iftImageTiles *tiles, iftVoxel v);

/**
 * @brief Computes the tile extremities in the original image coordinates.
 *
 * Recall that the tile coordinates are relative to the specified bounding box. These relative coordinates are only
 * absolute when the selected bounding box corresponds to the entire original image.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param tiles A pointer to an iftImageTiles record.
 * @param tile  The tile index.
 * @return The bounding box of the tile.
 */
iftBoundingBox iftGetTileInOriginalImageCoordinates(iftImageTiles *tiles, int tile);

/**
 * @brief Extracts an specific tile from an image given pre-computed tile coordinates.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img   The image from which the tile will be extracted.
 * @param tiles A pointer to an iftImageTiles record.
 * @param tile  The tile index.
 * @return The extracted tile.
 */
iftImage *iftExtractTile(iftImage *img, iftImageTiles *tiles, int tile);

/**
 * @brief Splits the entire image into a number of tiles specified for each dimension.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img      The image from which the tile will be extracted.
 * @param ntiles_x The number of tiles in the X axis.
 * @param ntiles_y The number of tiles in the Y axis.
 * @param ntiles_z The number of tiles in the Z axis.
 * @param tiles    Returns a pointer to an iftImageTiles record. It may not be NULL since it returns the data structure
 * with the actual number of computed tiles.
 * @return The extracted tiles.
 *
 * @sa iftImageTilesByEquallyDividingAxes
 */
iftImage **iftSplitImageIntoTilesByEquallyDividingAxes(iftImage *img, int ntiles_x, int ntiles_y, int ntiles_z,
                                                       iftImageTiles **tiles);

/**
 * @brief Splits the entire image into a series of tiles specified for each dimension.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img      The image from which the tile will be extracted.
 * @param ntiles_x The number of tiles in the X axis.
 * @param ntiles_y The number of tiles in the Y axis.
 * @param ntiles_z The number of tiles in the Z axis.
 * @param tiles    Returns a pointer to an iftImageTiles record. It may not be NULL since it returns the data structure
 * with the actual number of computed tiles.
 *
 * @return The extracted tiles.
 */
iftImage **iftSplitImageIntoTilesByStriding(iftImage *img, int tile_xsize, int tile_ysize, int tile_zsize, int xstride,
                                            int ystride, int zstride, iftImageTiles **tiles);
/**
 * @brief Recomposes a series of tiles into the original image.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img_tiles The image tiles.   
 * @param tiles     A pointer to an iftImageTiles record
 * @return The image of the recomposed tiles.
 *
 * @note If the image tiles overlap, this is NOT checked by this function and the result will simply overlap as
 * well.
 */
iftImage *iftRecomposeImageFromTiles(iftImage **img_tiles, iftImageTiles *tiles);


/**
 * @brief  Similarly to iftAddFrame, iftAddPadding pads the image with the given value but only after the
 * end of the image.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img   The original image
 * @param sx    The padding for the X axis.
 * @param sy    The padding for the Y axis.
 * @param sz    The padding for the Z axis.
 * @param value The value for padding.
 * @return The padded image.
 */
iftImage *iftAddPadding(iftImage *img, int sx, int sy, int sz, int value);


/**
 * @brief  Similarly to iftRemFrame, iftRemPadding removes padding from the end of the image only.
 *
 * @author Thiago V. Spina
 * @date October 23, 2015
 *
 * @param img   The padded image
 * @param sx    The padding for the X axis.
 * @param sy    The padding for the Y axis.
 * @param sz    The padding for the Z axis.
 * @return The unpadded image.
 */
iftImage *iftRemPadding(iftImage *fimg, int sx, int sy, int sz);

/**
 * @brief Checks if the image is a label image contained ALL labels from 0 to n_objects.
 * @author Samuel Martins
 * @date Feb 10, 2016
 * @ingroup Image
 * 
 * @param  label_img  Image to be checked.
 * @param  n_objects Number of objects from the image except the background (label 0).
 * @return True, if it is a label image and it contains all labels. False, otherwise.
 * 
 * @note The background (label 0) is actually ignored.
 * @warning The image must contain all labels
 */
bool iftIsLabelImage(const iftImage *label_img, int n_objects);


/**
 * @brief Gets the number of objects from a label image.
 * @author Samuel Martins
 * @date Feb 20, 2016
 * @ingroup Image
 * 
 * Since it is expected that the label image has voxel values in the range
 * [0..n_objects], where the value 0 corresponds to the background and the
 * ith value correponds to the object i, the number of objects of the image is its maximum value. 
 * 
 * @param  label_img  Label Image.
 * @return The number of objects from the input image.
 *
 * @warning It does not check if really the label image only has voxel values in the range
 * 0..n_objects. To ckeck that, use the function iftIsLabelImage()
 * @warning It does not check if the image is in grayscale or not.
 * 
 * @sa iftGetNumberOfObjectsFromImagePathname()
 */
int iftGetNumberOfObjectsFromImage(const iftImage *label_img);


/**
 * @brief Gets the number of objects from a label image pathname.
 * @author Samuel Martins
 * @date Feb 20, 2016
 * @ingroup Image
 *
 * @note See the documentation of iftGetNumberOfObjectsFromImage().
 * @sa iftGetNumberOfObjectsFromImage()
 */
int iftGetNumberOfObjectsFromImagePathname(const char *label_img_pathname);


/**
 * @brief Converts a Bit Map into a Binary Image.
 * @author Samuel Martins
 * @date Apr 7, 2016
 * @ingroup Image
 *
 * @param bmap Bit Map to be converted.
 * @param xsize Size of the X axis from the Image Domain.
 * @param ysize Size of the Y axis from the Image Domain.
 * @param zsize Size of the Z axis from the Image Domain.
 * 
 * @warning The Input Image Domain (sizes of the axes) must have the same size of the Bit Map.
 */
iftImage *iftBMapToBinImage(const iftBMap *bmap, int xsize, int ysize, int zsize);


/**
 * @brief Converts a Bin Image into a Bit Map.
 * @author Samuel Martins
 * @date Apr 12, 2016
 * @ingroup Image
 */
iftBMap *iftBinImageToBMap(const iftImage *bin_img);


/**
 * @brief Get the Labels from a Label Image.
 * @author Samuel Martins
 * @date Apr 17, 2016
 * 2ingroup Image
 */
iftIntArray *iftGetLabels(const iftImage *label_img);


/**
 * @brief Gets the most frequent label from a Label Image by Majority Voting.
 * @author Samuel Martins
 * @date May 31, 2016
 *
 * @param  label_img         Label Image.
 * @param  bool exclude_zero Tells if the Zero value (label) must be ignored in the majority voting.
 * @return                   The most frequent label from the Label Image.
 *
 * @note For ties, the first label checked is chosen.
 * @warning Label Image must have values >= 0 (No check for this in the function).
 */
int iftGetLabelByMajorityVoting(const iftImage *label_img, bool ignore_zero);


/**
 * @brief Gets the True Labels for given Superpixels/Supervoxel using Majority Voting.
 * @author Samuel Martins
 * @date May 31, 2016
 *
 * For each supervoxel from <b>super_img</b>, a majority voting over their pixels on Label Image
 * <b>label_img</b> determines its True Label. \n
 * * The Supervoxel Image MUST HAVE ONLY labels from 1 to n. Zero labels will be ignored. \n
 * Thus, for example, the supervoxel <b>i</b> will be the sample <b>i-1</b> in the resulting True Label Array. \n\n
 *
 * @param  label_img Label Image.
 * @param  super_img Superpixel Image.
 * @return           Array of size (n_supervoxels+1) with the true labels from each supervoxel from [0, n_supervoxels].
 *
 * @note For ties, the first label checked is chosen.
 * @warning Label Image and Superpixel Image must have values >= 0 (No check for this in the function).
 */
iftIntArray *iftGetSupervoxelTrueLabelByMajorityVoting(const iftImage *label_img, const iftImage *super_img);


/**
 * @brief Converts an Adjacency Relation into a Binary Image, where each reached voxel by the adjacency
 * has value 1
 * @author Samuel Martins
 * @date Apr 29, 2016
 * @ingroup Image
 *
 * The Adjacency relation corresponds to the displacements from the reference voxel coordinate <b>ref_voxel</b>. \n
 * Such reference voxel should belong to the domain of the image to be constructed, which corresponds
 * to the sizes passed as parameters.
 * 
 * @param adj       Adjacency Relation to be converted into a Binary Image.
 * @param xsize     Size of the X axis from the resulting Binary Image.
 * @param ysize     Size of the Y axis from the resulting Binary Image.
 * @param zsize     Size of the Z axis from the resulting Binary Image.
 * @param ref_voxel Reference voxel coordinate (from the resulting Binary Image) where the displacements
 *                  are computed.
 * @return Converted Binary Image 
 * 
 * @warning The reference voxel should belong to the domain of the Resulting Binary Image
 */
iftImage *iftAdjRelToImage(const iftAdjRel *adj, int xsize, int ysize, int zsize, iftVoxel ref_voxel);

/**
 * @brief Assigns colors to the luminance values of an input image
 * according to the blue to red (rainbow) color map 
 * 
 * @author Alexandre Falcao
 * @date Aug, 5th, 2016
 * @ingroup Image
 *
 */
iftImage *iftColorCoding(iftImage* img);


/**
 * @brief Gets the Orientation from a string: ["AXIAL", "CORONAL", "SAGITTAL"] 
 * 
 * @author Samuel Martins
 * @date Nov 1, 2016
 * @ingroup Image
 *
 */
iftImagePlaneOrientation iftGetImagePlaneOrientation(const char *orientation);

/**
 * @brief Read integer image in vtk file format
 * 
 * @author Alexandre Falcao
 * @date Nov 6th, 2016
 * @ingroup Image
 * @param  filename  Filename of the input integer image in vtk file format 
 * @return integer image
 *
 */
iftImage *iftReadVTKImage(const char *filename);

/**
 * @brief Write integer image in vtk file format
 * 
 * @author Alexandre Falcao
 * @date Nov 6th, 2016
 * @ingroup Image
 * @param img       Input integer image
 * @param filename  Filename of the output integer image in vtk file format 
 *
 */
void iftWriteVTKImage(iftImage *img, char *filename);


/**
 * @brief Translate the content of an image according to a displacement vector.
 * 
 * If a given voxel is moved out of the image domain, it is ignored. \n
 * Voxels without new values have zero-value.
 * 
 * @param  img  Image whose content will be translated.
 * @param  disp Displacement vector used to translated the image's content.
 * @return      Translated Image.
 */
iftImage *iftTranslateImageContent(const iftImage *img, iftVector disp_vec);

void iftSetImageBoundingBoxValue(iftImage *img, iftBoundingBox bb, int value );

/**
 * @brief Relabels a grayscale image, the first value in the image will be one, and the values will be consecutives. This function serves to correct label images
 * @param img An image
 * @author Adan Echemendia Montero
 * @return A relabeled image
 */
iftImage *iftRelabelGrayScaleImage(const iftImage *img);
/** @} */

#ifdef __cplusplus
}
#endif


#endif


