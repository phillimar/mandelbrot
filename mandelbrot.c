#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define DIM_X (3840 * 4)
#define DIM_Y (2160 * 4)
#define MAX_ITER 1000
#define MAX_BOUND 1000000

void map_colour(unsigned short index, unsigned short *rgb)
{
    float normalised_index = 1 - ((float)index / (float)MAX_ITER);

    unsigned short value = 65535 * (normalised_index * normalised_index * normalised_index * normalised_index * normalised_index * normalised_index * normalised_index);

    /* flip endian-ness */
    unsigned short high = (value & 0xff00) >> 8;
    unsigned short low = (value & 0x00ff) << 8;
    value = high | low;

    rgb[0] = value; /* red */
    rgb[1] = value; /* green */
    rgb[2] = value; /* blue */
}

void write_ppm(char *filename, unsigned int dimx, unsigned int dimy, unsigned short *img_data)
{
    int i, j;
    static unsigned short colour[3];

    FILE *fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%d %d\n65535\n", dimx, dimy); /* Magic PPM file header - max colour value of 65535 */
    for (j = 0; j < dimy; ++j)
    {
        for (i = 0; i < dimx; ++i)
        {
            map_colour(img_data[i + j * dimx], colour);
            fwrite(colour, 2, 3, fp);
        }
    }
    fclose(fp);
}

/* single threaded version */
void mandelbrot(double minx, double miny, double xsize, unsigned int dimx, unsigned int dimy, unsigned short max_iter, unsigned short *img_data)
{
    unsigned int i, j;
    double increment = xsize / (double)dimx;
    double c_real, c_img;
    double z_real, z_img, z_real_temp;
    unsigned short iter = 0;

    c_img = miny;
    for (j = 0; j < dimy; ++j)
    {
        c_real = minx;
        for (i = 0; i < dimx; ++i)
        {
            z_real = c_real;
            z_img = c_img;
            iter = 0;
            while (((z_real * z_real + z_img * z_img) < MAX_BOUND) && (iter < max_iter))
            {
                /* iterate Z -> z^2 + c */
                z_real_temp = z_real * z_real - z_img * z_img + c_real;
                z_img = 2 * z_real * z_img + c_img;
                z_real = z_real_temp;
                ++iter;
            }
            img_data[i + j * dimx] = iter;
            c_real += increment;
        }
        c_img += increment;
    }
}

typedef struct mandelbrot_state
{
    double minx;
    double miny;
    double increment;
    unsigned int dimx;
    unsigned int dimy;
    unsigned short max_iter;
    unsigned short *img_data;

    /* line num reads/ writes protected by line_num_mutex */
    unsigned int line_num;
    pthread_mutex_t line_num_mutex;

} mandelbrot_state;

void *mandelbrot_core(void *ptr)
{
    mandelbrot_state *state = (mandelbrot_state *)ptr;
    unsigned short max_iter = state->max_iter;
    double increment = state->increment;
    unsigned int dimx = state->dimx;
    unsigned int dimy = state->dimy;
    double minx = state->minx;
    double miny = state->miny;

    double c_real, c_img, z_real, z_img, z_real_temp;
    unsigned int iter;
    unsigned int i, j;

    while (1)
    {
        pthread_mutex_lock(&state->line_num_mutex);
        j = state->line_num++;
        pthread_mutex_unlock(&state->line_num_mutex);
        if (j >= dimy)
            break;

        c_img = miny + increment * j;
        c_real = minx;

        for (i = 0; i < dimx; ++i)
        {
            z_real = c_real;
            z_img = c_img;
            iter = 0;
            while (((z_real * z_real + z_img * z_img) < MAX_BOUND) && (iter < max_iter))
            {
                /* iterate Z -> z^2 + c */
                z_real_temp = z_real * z_real - z_img * z_img + c_real;
                z_img = 2 * z_real * z_img + c_img;
                z_real = z_real_temp;
                ++iter;
            }
            state->img_data[i + j * dimx] = iter;
            c_real += increment;
        }
    }

    return NULL;
}

/* multi threaded mandelbrot... */
void mandelbrot_multi(double minx, double miny, double xsize, unsigned int dimx, unsigned int dimy, unsigned short max_iter, unsigned short *img_data, unsigned short num_threads)
{
    unsigned int i;
    mandelbrot_state state;

    state.minx = minx;
    state.miny = miny;
    state.dimx = dimx;
    state.dimy = dimy;
    state.max_iter = max_iter;
    state.img_data = img_data;

    pthread_mutex_init(&state.line_num_mutex, NULL);
    state.line_num = 0;

    state.increment = xsize / (double)dimx;

    pthread_t *thread;
    thread = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

    for (i = 0; i < num_threads; ++i)
        pthread_create(&thread[i], NULL, mandelbrot_core, &state);

    for (i = 0; i < num_threads; ++i)
        pthread_join(thread[i], NULL);
}

int main(void)
{
    unsigned short *img_data = (unsigned short *)malloc(DIM_X * DIM_Y * sizeof(short));

    mandelbrot_multi(-2.5, -1, 3.5, DIM_X, DIM_Y, MAX_ITER, img_data, 64);

    write_ppm("mb.ppm", DIM_X, DIM_Y, img_data);

    return EXIT_SUCCESS;
}