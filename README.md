# TIFF Image Filters

The following programs implement various transform an spatial filters on TIFF images.

## Spatial filters implemented

* Min Filter
* Max Filter
* Mean Filte
* Median Filter  
* Average Filter
* Laplacican Fil
* Sobel Horz & Vert Filter
* Robert's Cross Gradient Filter

## Transform filters implemented

* Fourier Transform
* Inverse Fourier Transform
* Notch Filter
* Low Pass Filter
* High Pass Filter
* Laplacican Filter
* Homomorphic Filter

## Compiling

```bash
cd image-filters
make
```

## Usage

After compiling, to use spatial filters run:-

```bash
./FiltersSpatial
```

To use transform filters run:-

```bash
./FiltersTransform
```

## Sample I/O

Input image

![input](./lena50.tif)

Output image after fourier transformation

![output](./lena50.tiff)
