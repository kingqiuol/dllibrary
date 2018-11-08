#include<stdio.h>
#include<stdlib.h>
#include<math.h>   
#include<assert.h>
#include<string.h>

#include"spatial_filter.h"


#define FOR_EACH_PIXEL(W, H, FUNC)  for(int j = 0; j < H; ++j){for(int i = 0; i < W; ++i)FUNC}
#define SQUARE(x) ((x)*(x))
#define ABS(x) (fabs(x))
#define TWO_PI 6.2831853
#define FLAG 2

/******************************************
         Generating filter  
*******************************************/
Image make_box_filter(int w) {
    Image filter = make_image(w, w, 1);
    FOR_EACH_PIXEL(filter.w, filter.h,
                   set_pixel(filter, i, j, 0, 1.0 / w / w);)
    return filter;
}

Image make_highpass_filter() {
    Image hi_filter = make_image(3, 3, 1);
    float data[9] = {0, -1.0, 0,
                     -1.0, 4.0, -1.0,
                     0, -1.0, 0};
    memcpy(hi_filter.data, data, sizeof(data));

    return hi_filter;
}

Image make_sharpen_filter() {
    Image sharpen_filter = make_image(3, 3, 1);

    float data[9] = {0, -1, 0,
                     -1, 5, -1,
                     0, -1, 0};
    memcpy(sharpen_filter.data, data, sizeof(data));

    return sharpen_filter;
}

Image make_emboss_filter() {
    Image img = make_image(3, 3, 1);
    float data[9] = {-2, -1, 0,
                     -1, 1, 1,
                     0, 1, 2};
    memcpy(img.data, data, sizeof(data));

    return img;
}

Image make_gaussian_filter(float sigma) {
    int w = (int) ceil(6 * sigma);
    w = w & 1 ? w : w + 1;

    int center = w / 2;
    Image filter = make_image(w, w, 1);

    int i, j;
    FOR_EACH_PIXEL(w, w, set_pixel(filter, i, j, 0, 1 / (TWO_PI * SQUARE(sigma)) *
                                                    exp(-(SQUARE(i - center) + SQUARE(j - center)) /
                                                        (2 * SQUARE(sigma))));
    );

    return filter;
}

Image make_gx_filter() {
    Image img = make_image(3, 3, 1);
    float data[9] = {-1, 0, 1,
                     -2, 0, 2,
                     -1, 0, 1};
    memcpy(img.data, data, sizeof(data));

    return img;
}

Image make_gy_filter() {
    Image img = make_image(3, 3, 1);
    float data[9] = {-1, -2, -1,
                     0, 0, 0,
                     1, 2, 1};
    memcpy(img.data, data, sizeof(data));

    return img;
}


/******************************************
              spatial filter
*******************************************/
Image convolve_image(Image im, Image filter, int preserve, int flag) {
    assert(im.c == filter.c || filter.c == 1);

    Image new_img = make_image(im.w, im.h, preserve ? im.c : 1);
    assert(filter.w % 2);
    int center_x = filter.w / 2, center_y = filter.h / 2;
    for (int c = 0; c < im.c; ++c) {
        FOR_EACH_PIXEL(new_img.w, new_img.h, {
            int cur_i = i;
            int cur_j = j;

            //边界位置跳过
            if (cur_i < center_x || cur_j < center_y || cur_i >= new_img.w - center_x || \
                    cur_j >= new_img.h - center_y) {
                continue;
            }

            float q = preserve ? 0 : get_pixel(new_img, cur_i, cur_j, 0);
            FOR_EACH_PIXEL(filter.w, filter.h,
                           q += get_pixel(im, cur_i - (center_x - i), cur_j - (center_y - j), c) *
                                get_pixel(filter, i, j, filter.c == 1 ? 0 : c);
            );
            set_pixel(new_img, cur_i, cur_j, preserve ? c : 0, q);
        });
    }

    constrain_image(new_img);
    if (flag != 0) {
        new_img = crop_image(new_img, center_x, center_y, new_img.w - 2 * center_x, new_img.h - 2 * center_y);
    }

    return new_img;
}

//边界扩充
//共有三种模式：0：不扩充；1：扩充零；2：扩充为最邻近像素值
Image boundary_expansion(Image im, int filter_size, int flag) {
    assert(flag == 0 || flag == 1 || flag == 2);

    if (flag == 0) {
        return im;
    }
    int fill = filter_size / 2;
    Image fill_image = make_image(im.w + 2 * fill, im.h + 2 * fill, im.c);
    int i, j, k;
    for (k = 0; k < fill_image.c; k++) {
        for (j = 0; j < fill_image.h; j++) {
            for (i = 0; i < fill_image.w; i++) {
                float val = 0.0;
                if (i < fill || j < fill || i >= (fill_image.w - fill) || j >= (fill_image.h - fill)) {
                    if (flag == 1) {
                        val = 0.0;
                    } else if (flag == 2) {
                        if (i < fill && j < fill) {
                            val = get_pixel(im, 0, 0, k);
                        } else if (i < fill && j >= fill && j < (fill_image.h - fill)) {
                            val = get_pixel(im, 0, j - fill, k);
                        } else if (i < fill && j >= (fill_image.h - fill)) {
                            val = get_pixel(im, 0, im.h - 1, k);
                        } else if (i >= fill && i < fill_image.w - fill && j >= (fill_image.h - fill)) {
                            val = get_pixel(im, i - fill, im.h - 1, k);
                        } else if (i >= (fill_image.w - fill) && j >= (fill_image.h - fill)) {
                            val = get_pixel(im, im.w - 1, im.h - 1, k);
                        } else if (i >= (fill_image.w - fill) && j >= fill && j < (fill_image.h - fill)) {
                            val = get_pixel(im, im.w - 1, j - fill, k);
                        } else if (i >= (fill_image.w - fill) && j < fill) {
                            val = get_pixel(im, im.w - 1, 0, k);
                        } else {
                            get_pixel(im, i - fill, 0, k);
                        }
                    }
                } else {
                    val = get_pixel(im, i - fill, j - fill, k);
                }
                set_pixel(fill_image, i, j, k, val);
            }
        }
    }

    return fill_image;
}

Image blur(Image im, int w, int flag) {
    Image filter = make_box_filter(w);
    Image b_im = boundary_expansion(im, filter.w, flag);
    Image blur_image = convolve_image(b_im, filter, 1, flag);

    free_image(filter);
    return blur_image;
}

Image GaussianBlur(Image im, float sigma, int flag) {
    Image gauss_filter = make_gaussian_filter(sigma);
    Image b_im = boundary_expansion(im, gauss_filter.w, flag);
    Image gauss_image = convolve_image(b_im, gauss_filter, 1, flag);

    free_image(gauss_filter);
    return gauss_image;
}

Image HighpassBlur(Image im, int flag) {
    Image filter = make_highpass_filter();
    Image b_im = boundary_expansion(im, filter.w, flag);
    Image highpass_image = convolve_image(b_im, filter, 0, flag);

    free_image(filter);
    return highpass_image;
}

Image SharpenBlur(Image im, float flag) {
    Image filter = make_sharpen_filter();
    Image b_im = boundary_expansion(im, filter.w, flag);
    Image sharpen_image = convolve_image(b_im, filter, 1, flag);

    free_image(filter);
    return sharpen_image;
}

Image EmbossBlur(Image im, int flag) {
    Image filter = make_emboss_filter();
    Image b_im = boundary_expansion(im, filter.w, flag);
    Image emboss_image = convolve_image(b_im, filter, 1, flag);

    free_image(filter);
    return emboss_image;
}

void l1_normalize(Image im) {
    for (int c = 0; c < im.c; ++c) {
        float sum = 0;
        FOR_EACH_PIXEL(im.w, im.h,
                       sum += get_pixel(im, i, j, c);
        )
        if (sum != 0) {
            FOR_EACH_PIXEL(im.w, im.h,
                           set_pixel(im, i, j, c, get_pixel(im, i, j, c) / sum);
            )
        }
    }
}

void feature_normalize(Image im) {
    float min, max;
    min = max = get_pixel(im, 0, 0, 0);
    for (int c = 0; c < im.c; ++c) {
        FOR_EACH_PIXEL(im.w, im.h, {
            float p = get_pixel(im, i, j, c);
            if (p < min) min = p; else if (p > max) max = p;
        });
    }
    float range = max - min;
    for (int c = 0; c < im.c; ++c) {
        FOR_EACH_PIXEL(im.w, im.h, {
            float p = get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, range == 0 ? 0 : (p - min) / range);
        });
    }
}

Image *Sobel(Image im, int flag) {
    Image gx_filter = make_gx_filter();
    Image gy_filter = make_gy_filter();

    Image b_im = boundary_expansion(im, gx_filter.w, flag);
    Image gx = convolve_image(b_im, gx_filter, 0, flag);
    Image gy = convolve_image(b_im, gy_filter, 0, flag);
    free_image(gx_filter);
    free_image(gy_filter);

    Image *result = calloc(2, sizeof(Image));
    result[0] = make_image(gx.w, gx.h, 1);
    result[1] = make_image(gy.w, gy.h, 1);

    FOR_EACH_PIXEL(gx.w, gx.h, {
        set_pixel(result[0], i, j, 0, sqrt(SQUARE(get_pixel(gx, i, j, 0)) + SQUARE(get_pixel(gy, i, j, 0))));
        float angle = atan2(get_pixel(gy, i, j, 0), get_pixel(gx, i, j, 0));
        set_pixel(result[1], i, j, 0, angle);
    });
    free_image(gx);
    free_image(gy);

    constrain_image(result[0]);
    constrain_image(result[1]);
    return result;
}
