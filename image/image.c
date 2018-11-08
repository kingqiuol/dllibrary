#include<stdio.h>
#include<stdlib.h>
#include<math.h>   
#include<time.h>
#include<string.h>

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"

#include"image.h"

#define FOR_EACH_PIXEL(W, H, FUNC) for(int j=0;j<H;j++){for(int i=0;i<W;i++){FUNC}}

/******************************************
         Basic operations 
*******************************************/
float get_pixel(Image im, int x, int y, int c) {
    assert(x < im.w && x >= 0 && y < im.h && y >= 0 &&
           c < im.c && c >= 0);
    return im.data[x + y * im.w + c * im.h * im.w];
}

void set_pixel(Image im, int x, int y, int c, float v) {
    assert(x < im.w && x >= 0 && y < im.h && y >= 0 &&
           c < im.c && c >= 0);
    im.data[x + y * im.w + c * im.h * im.w] = v;
}

void add_pixel(Image im, int x, int y, int c, float v) {
    assert(x < im.w && x >= 0 && y < im.h && y >= 0 &&
           c < im.c && c >= 0);
    im.data[x + y * im.w + c * im.h * im.w] += v;
}

void fill_image(Image im, float s) {
    int i;
    for (i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] = s;
    }
}

void constrain_image(Image im) {
    int i = 0;
    for (i = 0; i < im.w * im.h * im.c; i++) {
        if (im.data[i] < 0) im.data[i] = 0;
        if (im.data[i] > 1) im.data[i] = 1;
    }
}

Image add_image(Image a, Image b) {
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    Image img = make_image(a.w, a.h, a.c);
    for (int c = 0; c < a.c; c++) {
        FOR_EACH_PIXEL(a.w, a.h,
                       set_pixel(img, i, j, c, (float) get_pixel(a, i, j, c) + get_pixel(b, i, j, c));
        );
    }
    constrain_image(img);

    return img;
}

Image sub_image(Image a, Image b) {
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    Image img = make_image(a.w, a.h, a.c);
    for (int c = 0; c < a.c; c++) {
        FOR_EACH_PIXEL(a.w, a.h,
                       set_pixel(img, i, j, c, (float) get_pixel(a, i, j, c) - get_pixel(b, i, j, c));
        );
    }
    constrain_image(img);

    return img;

}

Image copy_image(Image p) {
    Image copy;
    copy.h = p.h;
    copy.w = p.w;
    copy.c = p.c;
    copy.data = calloc(p.h * p.w * p.c, sizeof(float));
    memcpy(copy.data, p.data, p.h * p.w * p.c * sizeof(float));

    return copy;
}

void print_image(Image im) {
    int i, j, k;
    for (k = 0; k < im.c; k++) {
        for (i = 0; i < im.h; i++) {
            for (j = 0; j < im.w; j++) {
                printf("%.21f, ", im.data[k * im.w * im.h + i * im.w + j]);
                if (j > 30) break;
            }
            printf("\n");
            if (i > 30) break;
        }
        printf("\n");
    }
}

/******************************************
        Image foundation function 
*******************************************/
Image make_empty_image(int w, int h, int c) {
    Image out;
    out.h = h;
    out.w = w;
    out.c = c;
    out.data = 0;

    return out;
}

Image make_image(int w, int h, int c) {
    Image out = make_empty_image(w, h, c);
    out.data = calloc(w * h * c, sizeof(float));

    return out;
}

Image load_image_stb(char *filename, int channels) {
    int w, h, c;
    unsigned char *data = stbi_load(filename, &w, &h, &c, channels);

    if (!data) {
        fprintf(stderr, "Cannot load image \"%s\"\nSTB Reason: %s\n",
                filename, stbi_failure_reason());
        exit(0);
    }

    if (channels) c = channels;
    Image im = make_image(w, h, c);

    int i, j, k;
    for (k = 0; k < c; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                int dest_index = i + w * j + w * h * k;
                int src_index = k + c * i + c * w * j;
                im.data[dest_index] = (float) data[src_index] / 255.;
            }
        }
    }
    if (im.c == 4) im.c = 3;

    free(data);
    return im;
}

void save_image_stb(Image im, const char *name) {
    char buff[256];
    sprintf(buff, "./%s.jpg", name);
    unsigned char *data = calloc(im.w * im.h * im.c, sizeof(char));
    int i, k;
    for (k = 0; k < im.c; ++k) {
        for (i = 0; i < im.w * im.h; ++i) {
            data[i * im.c + k] = (unsigned char) roundf((255 * im.data[i + k * im.w * im.h]));
        }
    }
    //int success = stbi_write_png(buff, im.w, im.h, im.c, data, im.w*im.c);
    int success = stbi_write_jpg(buff, im.w, im.h, im.c, data, 100);
    free(data);
    if (!success) fprintf(stderr, "Failed to write image %s\n", buff);
}

void save_image(Image im, const char *name) {
    save_image_stb(im, name);
}

Image load_image(char *filename) {
    Image out = load_image_stb(filename, 0);
    return out;
}

void free_image(Image im) {
    if (im.data) {
        free(im.data);
    }
}

void embed_image(Image source, Image dst, int dx, int dy) {
    int x, y, k;
    for (k = 0; k < source.c; k++) {
        for (y = 0; y < source.h; k++) {
            for (x = 0; x < source.w; x++) {
                float Val = get_pixel(source, x, y, k);
                set_pixel(dst, dx + x, dy + y, k, Val);
            }
        }
    }
}

Image get_image_layer(Image m, int l) {
    Image out = make_image(m.w, m.h, 1);
    int i;
    for (i = 0; i < m.h * m.w; i++) {
        out.data[i] = m.data[i + l * m.h * m.w];
    }
    return out;
}

Image collape_images_vert(Image *ims, int n) {
    int color = 1;
    int border = 1;
    int h, w, c;
    int size = ims[0].h;
    h = size;
    w = (ims[0].w + border) * n - border;
    c = ims[0].c;
    if (c != 3 || !color) {
        h = (h + border) * c - border;
        c = 1;
    }

    Image filters = make_image(w, h, c);
    int i, j;
    for (i = 0; i < n; i++) {
        int w_offset = i * (size + border);
        Image copy = copy_image(ims[i]);
        if (c == 3 && color) {
            embed_image(copy, filters, w_offset, 0);
        } else {
            for (j = 0; j < copy.c; j++) {
                int h_offset = j * (size + border);
                Image layer = get_image_layer(copy, j);
                embed_image(layer, filters, w_offset, h_offset);
                free_image(layer);
            }
        }
        free_image(copy);
    }
    return filters;
}

void normalize_image(Image img) {
    int i;
    float min = 9999999;
    float max = -9999999;

    for (i = 0; i < img.h * img.w * img.c; i++) {
        float v = img.data[i];
        if (v < min) min = v;
        if (v > max) max = v;
    }
    if (max - min < .000000001) {
        min = 0;
        max = 1;
    }
    for (i = 0; i < img.h * img.w * img.c; i++) {
        img.data[i] = (img.data[i] - min) / (max - min);
    }
}

void rgbgr_image(Image im) {
    int i;
    for (i = 0; i < im.w * im.h; ++i) {
        float swap = im.data[i];
        im.data[i] = im.data[i + im.w * im.h * 2];
        im.data[i + im.w * im.h * 2] = swap;
    }
}

#ifdef OPENCV
void show_image_cv(image p, const char *name, IplImage *disp)
{
    int x,y,k;
    if(p.c == 3) rgbgr_image(p);
    //normalize_image(copy);
    
    char buff[256];
    //sprintf(buff, "%s (%d)", name, windows);
    sprintf(buff, "%s", name);
    
    int step = disp->widthStep;
    cvNamedWindow(buff, CV_WINDOW_NORMAL); 
    //cVMoVeWindow(buff, 100*(windows%10) + 200*(windows/10), 100*(windows%10));
    ++windows;
    for(y = 0; y < p.h; ++y){
        for(x = 0; x < p.w; ++x){
            for(k= 0; k < p.c; ++k){
            disp->imageData[y*step + x*p.c + k] = (unsigned char)(get_pixel(p,x,y,k)*255);
           }
        }
    }
    if(0){
        int w = 448;
        int h = w*p.h/p.w;
        if(h > 1000){  
            h = 1000;
            w=h*p.w/p.h
        }
        IplImage *buffer = disp;
        disp = cvCreateImage(cvSize(w, h), buffer->depth, buffer->nChannels);
        cvResize(buffer, disp, CV_INTER_LINEAR);
        cvReleaseImage(&buffer);
    }
    cvShowImage(buff, disp);
}
#endif

int show_image(Image img, const char *name, int ms) {
#ifdef OPENCV
    IplImage* disp=cvCreateImage(cvSize(img.w,img.h),IPL_DEPTH_8U,img.c);
    Image copy=copy_image(img);
    constrain_image(copy);
    show_image_cv(copy,name,disp);
    free_image(copy);
    cvReleaseImage(&disp);
    int c=cvWaitKey(ms);
    if(c!=-1) c=c%256;
    return c;
#else
    fprintf(stderr, "Not compiled with OpenCV, saVing to %s.png instead\n", name);
    save_image(img, name);
    return 0;
#endif
}

void show_images(Image *ims, int n, char *window) {
    Image img = collape_images_vert(ims, n);
    normalize_image(img);
    save_image(img, window);
    show_image(img, window, 1);
    free_image(img);
}

/******************************************          
        color space transformation 
*******************************************/
Image rgb_to_gray(Image im) {
    assert(im.c == 3);
    Image gray = make_image(im.w, im.h, 1);

    int i, j;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            gray.data[i + j * gray.w] = 0.299 * im.data[i + j * im.w] + 0.587 * im.data[i + j * im.w + im.w * im.h] +
                                        0.114 * im.data[i + j * im.w + im.w * im.h * 2];
        }
    }

    return gray;
}

float three_way_max(float a, float b, float c) {
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

float three_way_min(float a, float b, float c) {
    return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

Image rgb_to_hsv(Image im) {
    assert(im.c == 3);
    Image hsv = make_image(im.w, im.h, im.c);

    int i, j;
    float R, G, B;
    float H, S, V;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            R = get_pixel(im, i, j, 0);
            G = get_pixel(im, i, j, 1);
            B = get_pixel(im, i, j, 2);

            float c_max = three_way_max(R, G, B);
            float c_min = three_way_min(R, G, B);
            float delta = c_max - c_min;

            V = c_max;
            if (delta == 0) {
                H = 0;
                S = 0;
            } else {
                S = delta / c_max;
                if (c_max == R) H = (G - B) / delta;
                else if (c_max == G) H = (B - R) / delta + 2;
                else H = (R - G) / delta + 4;

                if (H < 0) H += 6;
                H = H / 6;
            }
            set_pixel(hsv, i, j, 0, H);
            set_pixel(hsv, i, j, 1, S);
            set_pixel(hsv, i, j, 2, V);
        }
    }

    return hsv;
}

Image rgb_to_yuv(Image im) {
    assert(im.c == 3);
    Image yuv = make_image(im.w, im.h, im.c);

    int i, j;
    float B, G, R;
    float Y, U, V;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            R = get_pixel(im, i, j, 0);
            G = get_pixel(im, i, j, 1);
            B = get_pixel(im, i, j, 2);

            Y = 0.299 * R + 0.587 * G + 0.114 * B;
            U = -0.147 * R - 0.289 * G - 0.436 * B;
            V = 0.615 * R - 0.515 * G - 0.100 * B;

            set_pixel(yuv, i, j, 0, Y);
            set_pixel(yuv, i, j, 1, U);
            set_pixel(yuv, i, j, 2, V);
        }
    }

    return yuv;
}

Image hsv_to_rgb(Image im) {
    assert(im.c == 3);
    Image rgb = make_image(im.w, im.h, im.c);

    int i, j;
    float H, S, V;
    float R, G, B;
    float f, p, q, t;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            H = 6 * get_pixel(im, i, j, 0);
            S = get_pixel(im, i, j, 1);
            V = get_pixel(im, i, j, 2);

            if (S == 0) {
                R = G = B = V;
            } else {
                int index = floor(H);
                f = H - index;
                p = V * (1 - S);
                q = V * (1 - S * f);
                t = V * (1 - S * (1 - f));
                if (index == 0) {
                    R = V;
                    G = t;
                    B = p;
                } else if (index == 1) {
                    R = q;
                    G = V;
                    B = p;
                } else if (index == 2) {
                    R = p;
                    G = V;
                    B = t;
                } else if (index == 3) {
                    R = p;
                    G = q;
                    B = V;
                } else if (index == 4) {
                    R = t;
                    G = p;
                    B = V;
                } else {
                    R = V;
                    G = p;
                    B = q;
                }
            }

            set_pixel(rgb, i, j, 0, R);
            set_pixel(rgb, i, j, 1, G);
            set_pixel(rgb, i, j, 2, B);
        }
    }

    return rgb;
}

Image yuv_to_rgb(Image im) {
    assert(im.c == 3);
    Image rgb = make_image(im.w, im.h, im.c);

    int i, j;
    float Y, U, V;
    float R, G, B;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            Y = get_pixel(im, i, j, 0);
            U = get_pixel(im, i, j, 1);
            V = get_pixel(im, i, j, 2);

            R = Y + 1.13983 * V;
            G = Y + -.39465 * U + -.58060 * V;
            B = Y + 2.03211 * U;

            set_pixel(rgb, i, j, 0, R);
            set_pixel(rgb, i, j, 1, G);
            set_pixel(rgb, i, j, 2, B);
        }
    }

    return rgb;
}

/******************************************          
            image enhancement 
*******************************************/
void scale_image(Image im, float s) {
    int i;
    for (i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] *= s;
    }
    constrain_image(im);
}

void translate_image(Image im, float s) {
    int i;
    for (i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] += s;
    }
}

/******************************************          
            image pretreatment
*******************************************/
Image random_augment_image(Image im, float angle, float aspect, int low, int high, int w, int h);//随机增强图像

/**旋转和剪切图像**/
Image crop_image(Image im, int dx, int dy, int w, int h) {
    assert(dx >= 0 && dy >= 0 && w > 0 && h > 0);
    assert((dx + w) <= im.w && (dy + h) <= im.h);

    Image crop = make_image(w, h, im.c);
    int i, j, k;
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                float value = get_pixel(im, dx + i, dy + j, k);
                set_pixel(crop, i, j, k, value);
            }
        }
    }

    return crop;
}

Image center_crop_image(Image im, int w, int h) {
    assert(w <= im.w && h <= im.h);
    //生成偏移量
    int dx = (im.w - w) / 2;
    int dy = (im.h - h) / 2;

    Image center_crop = make_image(w, h, im.c);
    center_crop = crop_image(im, dx, dy, w, h);

    return center_crop;
}

Image random_crop_image(Image im, int w, int h) {
    int offset_x = im.w - w;
    int offset_y = im.h - h;
    srand((int) time(0));
    int random_x = rand() % (offset_x - 1);
    int random_y = rand() % (offset_y - 1);

    int i, j, k;
    int index = 0;
    Image crop_image = make_image(w, h, im.c);

    for (k = 0; k < im.c; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                float val = get_pixel(im, random_x + i, random_y + j, k);
                set_pixel(crop_image, i, j, k, val);
            }
        }
    }

    return crop_image;
}

float get_pixel_extend(Image im, int x, int y, int c) {
    if (x >= im.w)x = im.w - 1;
    if (x < 0)x = 0;
    if (y >= im.h)y = im.h - 1;
    if (y < 0) y = 0;
    if (c < 0 || c >= im.c) return 0;

    return get_pixel(im, x, y, c);
}

float bilinear_interpolate(Image im, float rx, float ry, int c) {
    int ix = (int) floorf(rx);
    int iy = (int) floorf(ry);

    float dx = rx - ix;
    float dy = ry - iy;

    float val = (1 - dy) * (1 - dx) * get_pixel_extend(im, ix, iy, c) + \
        dy * (1 - dx) * get_pixel_extend(im, ix, iy + 1, c) + \
        (1 - dy) * dx * get_pixel_extend(im, ix + 1, iy, c) + \
        dy * dx * get_pixel_extend(im, ix + 1, iy + 1, c);

    return val;
}

Image rotate_image(Image im, float rad) {
    float cx = im.w / 2;
    float cy = im.h / 2;
    Image rot_image = make_image(im.w, im.h, im.c);

    int x, y, c;
    for (c = 0; c < im.c; c++) {
        for (y = 0; y < im.h; y++) {
            for (x = 0; x < im.w; x++) {
                float rx = cos(rad) * (x - cx) - sin(rad) * (y - cy) + cx;
                float ry = sin(rad) * (x - cx) + cos(rad) * (y - cy) + cy;
                float val = bilinear_interpolate(im, rx, ry, c);
                set_pixel(rot_image, x, y, c, val);
            }
        }
    }

    return rot_image;
}

void rotate_image_cw(Image im, int times) {
    assert(im.w == im.h);
    times = (times + 400) % 4;
    int i, x, y, c;
    int n = im.w;
    for (i = 0; i < times; ++i) {
        for (c = 0; c < im.c; ++c) {
            for (x = 0; x < n / 2; ++x) {
                for (y = 0; y < (n - 1) / 2 + 1; ++y) {
                    float temp = im.data[y + im.w * (x + im.h * c)];
                    im.data[y + im.w * (x + im.h * c)] = im.data[n - 1 - x + im.w * (y + im.h * c)];
                    im.data[n - 1 - x + im.w * (y + im.h * c)] = im.data[n - 1 - y + im.w * (n - 1 - x + im.h * c)];
                    im.data[n - 1 - y + im.w * (n - 1 - x + im.h * c)] = im.data[x + im.w * (n - 1 - y + im.h * c)];
                    im.data[x + im.w * (n - 1 - y + im.h * c)] = temp;
                }
            }
        }
    }
}

Image flip_image(Image im) {
    Image flip = make_image(im.w, im.h, im.c);
    int i, j, k;
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < im.h; j++) {
            for (i = 0; i < im.w; i++) {
                float val = get_pixel(im, im.w - i - 1, j, k);
                set_pixel(flip, i, j, k, val);
            }
        }
    }

    return flip;
}

Image rotate_crop_image(Image im, float rad, float s, int w, int h, float dx, float dy, float aspect) {
    int x, y, c;
    float cx = im.w / 2.;
    float cy = im.h / 2.;
    Image rot = make_image(w, h, im.c);
    for (c = 0; c < im.c; ++c) {
        for (y = 0; y < h; ++y) {
            for (x = 0; x < w; ++x) {
                float rx = cos(rad) * ((x - w / 2.) / s * aspect + dx / s * aspect) -
                           sin(rad) * ((y - h / 2.) / s + dy / s) + cx;
                float ry = sin(rad) * ((x - w / 2.) / s * aspect + dx / s * aspect) +
                           cos(rad) * ((y - h / 2.) / s + dy / s) + cy;
                float val = bilinear_interpolate(im, rx, ry, c);
                set_pixel(rot, x, y, c, val);
            }
        }
    }
    return rot;
}

/**改变图像曝光，亮度和饱和度等操作**/
float rand_uniform(float min, float max) {
    srand((int) time(0));
    return min + (float) rand() * (max - min);
}

float rand_scale(float rate) {
    srand((int) time(0));
    return rand() * rate;
}

void scale_image_channel(Image im, int c, float v) {
    int i, j;
    for (j = 0; j < im.h; j++) {
        for (i = 0; i < im.w; i++) {
            float pixel = get_pixel(im, i, j, c);
            pixel = pixel * v;
            set_pixel(im, i, j, c, pixel);
        }
    }
}

Image distort_image(Image im, float hue, float sat, float val) {
    Image hsv = rgb_to_hsv(im);
    scale_image_channel(hsv, 1, sat);
    scale_image_channel(hsv, 2, val);

    int i;
    for (i = 0; i < hsv.w * hsv.h; i++) {
        hsv.data[i] = hsv.data[i] + hue;
        if (hsv.data[i] > 1) hsv.data[i] -= 1;
        if (hsv.data[i] < 0) hsv.data[i] += 1;
    }

    im = hsv_to_rgb(hsv);
    constrain_image(im);

    if (hsv.data) {
        free_image(hsv);
    }

    return im;
}

Image random_distort_image(Image im, float hue, float saturation, float exposure) {
    float dhue = rand_uniform(-hue, hue);
    float dsat = rand_scale(saturation);
    float dexp = rand_scale(exposure);

    Image dist_img = distort_image(im, dhue, dsat, dexp);

    return dist_img;

}

Image saturate_image(Image im, float sat) {
    Image hsv = rgb_to_hsv(im);
    scale_image_channel(hsv, 1, sat);
    im = hsv_to_rgb(hsv);
    constrain_image(im);

    if (hsv.data) {
        free_image(hsv);
    }

    return im;
}

Image exposure_image(Image im, float val) {
    Image hsv = rgb_to_hsv(im);
    scale_image_channel(hsv, 2, val);
    im = hsv_to_rgb(hsv);
    constrain_image(im);

    if (hsv.data) {
        free_image(hsv);
    }

    return im;
}

Image saturate_exposure_image(Image im, float sat, float exposure) {
    Image hsv = rgb_to_hsv(im);
    scale_image_channel(hsv, 1, sat);
    scale_image_channel(hsv, 2, exposure);
    im = hsv_to_rgb(hsv);
    constrain_image(im);

    if (hsv.data) {
        free_image(hsv);
    }

    return im;
}

/**调整图片尺寸**/
Image resize_image(Image im, int w, int h) {
    Image resized = make_image(w, h, im.c);//目标图像
    Image part = make_image(w, im.h, im.c);//中间推昂
    float w_scale = (float) (im.w - 1) / (w - 1);
    float h_scale = (float) (im.h - 1) / (h - 1);

    int i, j, k;
    //缩放宽度，高度保持不变
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < im.h; j++) {
            for (i = 0; i < w; i++) {
                float val = 0;
                if (i == w - 1 || im.w == 1) {
                    val = get_pixel(im, im.w - 1, j, k);
                } else {
                    float x = i * w_scale;
                    int ix = (int) x;
                    float dx = x - ix;
                    val = (1 - dx) * get_pixel(im, ix, j, k) + dx * get_pixel(im, ix + 1, j, k);
                }
                set_pixel(part, i, j, k, val);
            }
        }
    }

    //缩放高度
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < h; j++) {
            float y = j * h_scale;
            int jy = (int) y;
            float dy = y - jy;

            for (i = 0; i < w; i++) {
                float val = 0;
                if (j == h - 1 || im.h == 1) {
                    val = get_pixel(part, i, part.h - 1, k);
                } else {
                    val = (1 - dy) * get_pixel(part, i, jy, k) + dy * get_pixel(part, i, jy + 1, k);
                }
                set_pixel(resized, i, j, k, val);
            }
        }
    }

    free_image(part);
    return resized;
}

Image resize_min(Image im, int min) {
    int w = im.w, h = im.h;
    if (im.w < im.h) {
        w = min;
        h = (int) (im.h * min) / im.w;
    } else {
        h = min;
        w = (int) (w * min) / im.h;
    }
    if (w == im.w && h == im.h) {
        return im;
    }
    Image resized_min = resize_image(im, w, h);

    return resized_min;
}

Image resize_max(Image im, int max) {
    int w = im.w;
    int h = im.h;
    if (w > h) {
        h = (h * max) / w;
        w = max;
    } else {
        w = (w * max) / h;
        h = max;
    }
    if (w == im.w && h == im.h) return im;
    Image resized = resize_image(im, w, h);

    return resized;
}

/**平移增强**/
void place_image(Image im, int w, int h, int dx, int dy, Image canvas) {
    int i, j, k;
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                int rx = ((float) i / w) * im.w;
                int ry = ((float) j / h) * im.h;

                float val = bilinear_interpolate(im, rx, ry, k);
                set_pixel(canvas, i + dx, j + dy, k, val);
            }
        }
    }
}
