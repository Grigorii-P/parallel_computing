#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

//TODO сделать другой путь /Users/grigoriipogorelov/Desktop/pc hw/lena512.bmp

#define NUM_THREADS 64

typedef struct
{
    int threadId;
    unsigned char *src1, *src2, *src3;
    float *dst1, *dst2, *dst3;
    int width;
    int height;
    float sigma;
    int ksize;
} parameters;

//todo можно тоже попробовать распараллелить
int read_BMP(char* filename, unsigned char *info, unsigned char **dataR, unsigned char **dataG, unsigned char **dataB, int *size, int *width, int *height, int *offset, int *row_padded)
{
    int i = 0, j, k, read_bytes, h, w, o, p;
    unsigned char *data;
    
    FILE* f = fopen(filename, "rb");
    
    if(f == NULL)
    {
        printf ("Invalid filename: %s\n", filename);
        return -1;
    }
    
    
    read_bytes = fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
    if (read_bytes != 54)
    {
        printf ("Error at read: %d instead of 54 bytes", read_bytes);
        return -1;
    }
    
    
    // extract image data from header
    *width = *(int*)&info[18];
    *height = *(int*)&info[22];
    *size = *(int*)&info[2];
    *offset = *(int*)&info[10];
    *row_padded = (*width*3 + 3) & (~3);
    
    
    //printf ("Filename: %s, Width: %d, Row_padded: %d, Height: %d, Size:  %d, Offset: %d\n", filename, *width, *row_padded, *height, *size, *offset);
    w = *width;
    p = *row_padded;
    h = *height;
    o = *offset;
    
    data = (unsigned char*) malloc (p * h);
    *dataR = (unsigned char*) malloc (w * h);
    *dataG = (unsigned char*) malloc (w * h);
    *dataB = (unsigned char*) malloc (w * h);
    
    fseek(f, sizeof(unsigned char) * o, SEEK_SET);
    read_bytes = fread(data, sizeof(unsigned char), p * h, f);
    if (read_bytes != p * h)
    {
        printf ("Error at read: %d\n", read_bytes);
        free (data);
        return -1;
    }
    for (k = 0; k < h; k++)
    {
        i = k * p;
        for (j = 0; j < w; j++)
        {
            (*dataB)[k*w + j] = data[i];
            (*dataG)[k*w + j] = data[i + 1];
            (*dataR)[k*w + j] = data[i + 2];
            
            //printf ("BGR %d %d i= %d: %d %d %d\n", k, j, i, data[i], data[i+1], data[i+2]);
            i+= 3;
        }
    }
    
    free (data);
    fclose(f);
    return 0;
}

int write_BMP(char* filename, float *dataB, float *dataG, float *dataR, unsigned char *header, int offset, int width,  int row_padded, int height)
{
    int write_bytes = 0, i, pad_size;
    FILE* f = fopen(filename, "wb");
    unsigned char null_byte = 0, valR, valB, valG;
    
    write_bytes = fwrite (header, sizeof(unsigned char), offset, f);
    if (write_bytes < offset)
    {
        printf( "Error at writing the header\n");
        return -1;
    }
    
    
    for (i = 0; i< width*height; i++)
    {
        if ( dataB[i] > 256.0f || dataR[i] > 256.0f || dataG[i] > 256.0f ){
            printf( "Error: invalid value %f %f %f", dataB[i], dataG[i], dataR[i]);
            return -1;
        }
        
        valB = dataB[i];
        valG = dataG[i];
        valR = dataR[i];
        write_bytes = fwrite(&valB, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valB);
            return -1;
        }
        write_bytes = fwrite(&valG, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valG);
            return -1;
        }
        write_bytes = fwrite(&valR, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valR);
            return -1;
        }
        
        if ((i + 1) % width == 0 ) {
            pad_size = row_padded - width *3 ;
            while( pad_size-- > 0 ) {
                fwrite(&null_byte, sizeof( unsigned char), 1, f );
            }
        }
    }
    
    fclose (f);
    return 0;
}



float convolve(const float *kernel, const float *buffer, const int ksize) {
    float sum = 0.0f;
    int i;
    for(i=0; i<ksize; i++)
    {
        sum += kernel[i]*buffer[i];
    }
    return sum;
}


void gaussian_blur_row(int threadId, unsigned char *src, float *dst, int width, int height, float sigma, int ksize)
{
    int x, y, i, x1, y1;
    
    int halfksize = ksize / 2;
    float sum = 0.f, t;
    float *kernel, *buffer;
    
    
    
    // create Gaussian kernel
    kernel = (float*)malloc(ksize * sizeof(float));
    buffer = (float*)malloc(ksize * sizeof(float));
    
    if (!kernel || !buffer)
    {
        printf ("Error in memory allocation!\n");
        return;
    }
    
    // if sigma too small, just copy src to dst
    if (ksize <= 1)
    {
        for (y = 0; y < height; y++)
            for (x = 0; x < width; x++)
                dst[y*width + x] = src[y*width + x];
        return;
    }
    
    
    //compute the Gaussian kernel values
    for (i = 0; i < ksize; i++)
    {
        x = i - halfksize;
        t = expf(- x * x/ (2.0f * sigma * sigma)) / (sqrt(2.0f * M_PI) * sigma);
        kernel[i] = t;
        sum += t;
    }
    for (i = 0; i < ksize; i++)
    {
        kernel[i] /= sum;
        //printf ("Kernel [%d] = %f\n", i, kernel[i]);
    }
    
    
    // blur each row
    int p = threadId;
    while (p < height) {
        y = p;
        for (x1 = 0; x1 < halfksize; x1++) {
            buffer[x1] = (float) src[y * width];
        }
        
        for (x1 = halfksize; x1 < ksize - 1; x1++) {
            buffer[x1] = (float) src[y * width + x1 - halfksize];
        }
        
        for (x1 = 0; x1 < width; x1++) {
            i = (x1 + ksize - 1) % ksize;
            
            if (x1 < width - halfksize) {
                buffer[i] = (float) src[y * width + x1 + halfksize];
            } else {
                buffer[i] = (float) src[y * width + width - 1];
            }
            
            dst[y * width + x1] = convolve(kernel, buffer, ksize);
        }
        p += NUM_THREADS;
    }
    
    // clean up
    free(kernel);
    free(buffer);
}

void gaussian_blur_column(int threadId, unsigned char *src, float *dst, int width, int height, float sigma, int ksize)
{
    int x, y, i, x1, y1;
    
    int halfksize = ksize / 2;
    float sum = 0.f, t;
    float *kernel, *buffer;
    
    
    
    // create Gaussian kernel
    kernel = (float*)malloc(ksize * sizeof(float));
    buffer = (float*)malloc(ksize * sizeof(float));
    
    if (!kernel || !buffer)
    {
        printf ("Error in memory allocation!\n");
        return;
    }
    
    // if sigma too small, just copy src to dst
    if (ksize <= 1)
    {
        for (y = 0; y < height; y++)
            for (x = 0; x < width; x++)
                dst[y*width + x] = src[y*width + x];
        return;
    }
    
    
    //compute the Gaussian kernel values
    for (i = 0; i < ksize; i++)
    {
        x = i - halfksize;
        t = expf(- x * x/ (2.0f * sigma * sigma)) / (sqrt(2.0f * M_PI) * sigma);
        kernel[i] = t;
        sum += t;
    }
    for (i = 0; i < ksize; i++)
    {
        kernel[i] /= sum;
        //printf ("Kernel [%d] = %f\n", i, kernel[i]);
    }
    
    
    // blur each column
    int p = threadId;
    while (p < height) {
        x = p;
        for (y1 = 0; y1 < halfksize; y1++) {
            buffer[y1] = dst[0 * width + x];
        }
        for (y1 = halfksize; y1 < ksize - 1; y1++) {
            buffer[y1] = dst[(y1 - halfksize) * width + x];
        }
        
        for (y1 = 0; y1 < height; y1++) {
            i = (y1 + ksize - 1) % ksize;
            if (y1 < height - halfksize) {
                buffer[i] = dst[(y1 + halfksize) * width + x];
            } else {
                buffer[i] = dst[(height - 1) * width + x];
            }
            
            dst[y1 * width + x] = convolve(kernel, buffer, ksize);
        }
        p += NUM_THREADS;
    }
    
    
    // clean up
    free(kernel);
    free(buffer);
}

void * executor_row(void *par) {
    parameters *p = (parameters *)par;
    //    printf("Thread ID : %d\n", p->threadId);
    gaussian_blur_row(p->threadId, p->src1, p->dst1, p->width, p->height, p->sigma, p->ksize);
    gaussian_blur_row(p->threadId, p->src2, p->dst2, p->width, p->height, p->sigma, p->ksize);
    gaussian_blur_row(p->threadId, p->src3, p->dst3, p->width, p->height, p->sigma, p->ksize);
    return NULL;
}

void * executor_column(void *par) {
    parameters *p = (parameters *)par;
    //    printf("Thread ID : %d\n", p->threadId);
    gaussian_blur_column(p->threadId, p->src1, p->dst1, p->width, p->height, p->sigma, p->ksize);
    gaussian_blur_column(p->threadId, p->src2, p->dst2, p->width, p->height, p->sigma, p->ksize);
    gaussian_blur_column(p->threadId, p->src3, p->dst3, p->width, p->height, p->sigma, p->ksize);
    return NULL;
}


int main(int argc, char ** argv)
{
    //    pthread_t *threads;
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    parameters *arg;
    
    unsigned char info[54], *dataR = NULL, *dataG = NULL, *dataB = NULL;
    int blur_size, ret_code = 0, size, width, height, offset, row_padded;
    char *in_filename, *out_filename;
    float* dstB, *dstR, *dstG, sigma;
    
    
    if (argc != 5)
    {
        printf ("Usage: %s <filename.bmp> <sigma> <blur_size> <output_filename.bmp>", argv[0]);
        return -1;
    }
    in_filename = argv[1];
    out_filename = argv[4];
    blur_size = atoi (argv[3]);
    sigma = atof (argv[2]);
    ret_code = read_BMP(in_filename, info, &dataR, &dataG, &dataB, &size, &width, &height, &offset, &row_padded);
    if (ret_code < 0)
    {
        free (dataB);
        free (dataR);
        free (dataG);
        return -1;
    }
    
    dstB = (float*)malloc (width*height* sizeof(float));
    dstR = (float*)malloc (width*height* sizeof(float));
    dstG = (float*)malloc (width*height* sizeof(float));
    
    //    threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    
    //    arg = (parameters *)malloc(sizeof(parameters));
    arg = (parameters *)malloc(sizeof(parameters) * NUM_THREADS);
    
    
    for (int i = 0; i < NUM_THREADS; i++) {
        arg[i].threadId = i;
        arg[i].src1 = dataB;
        arg[i].src2 = dataR;
        arg[i].src3 = dataG;
        arg[i].dst1 = dstB;
        arg[i].dst2 = dstR;
        arg[i].dst3 = dstG;
        arg[i].width = width;
        arg[i].height = height;
        arg[i].sigma = sigma;
        arg[i].ksize = blur_size;
        //        printf("arg->threadId : %d\n", arg[i].threadId);
        pthread_create(&threads[i], NULL, executor_row, (void *)&arg[i]);
    }
    
    for (int i= 0;i < NUM_THREADS;i++){
        pthread_join(threads[i],NULL);
    }
    
    for (int i = 0; i < NUM_THREADS; i++) {
        arg[i].threadId = i;
        arg[i].src1 = dataB;
        arg[i].src2 = dataR;
        arg[i].src3 = dataG;
        arg[i].dst1 = dstB;
        arg[i].dst2 = dstR;
        arg[i].dst3 = dstG;
        arg[i].width = width;
        arg[i].height = height;
        arg[i].sigma = sigma;
        arg[i].ksize = blur_size;
        //        printf("arg->threadId : %d\n", arg[i].threadId);
        pthread_create(&threads[i], NULL, executor_column, (void *)&arg[i]);
    }
    
    
    for (int i= 0;i < NUM_THREADS;i++){
        pthread_join(threads[i],NULL);
    }
    
    ret_code = write_BMP (out_filename, dstB, dstG, dstR, info, offset, width, row_padded, height);
    
    free (dstB);
    free (dstR);
    free (dstG);
    free (dataB);
    free (dataR);
    free (dataG);
    free(arg);
    
    //    pthread_exit(NULL); //todo pthread_exit(NULL);
    return ret_code;
}

