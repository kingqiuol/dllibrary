#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

#include"image.h"
#include"spatial_filter.h"
#include"panorama_image.h"
#include "matrix.h"

#define FOR_EACH_PIXEL(W, H, FUNC) for(int j=0;j<H;j++){for(int i=0;i<W;i++)FUNC}

#define SQUARE(x) ((x)*(x))

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define ABS(X) ((X) < 0 ? -(X) : (X))


/******************************************
         Harris corner detection 
*******************************************/
Image structure_matrix(Image im, float sigma) {
    Image S = make_image(im.w, im.h, 3);

    Image gx = make_gx_filter();
    Image gy = make_gy_filter();
    Image ix = convolve_image(im, gx, 0, 0);
    Image iy = convolve_image(im, gy, 0, 0);
    free_image(gx);
    free_image(gy);

    FOR_EACH_PIXEL(im.w, im.h, {
        float x = get_pixel(ix, i, j, 0);
        float y = get_pixel(iy, i, j, 0);
        set_pixel(S, i, j, 0, x * x);
        set_pixel(S, i, j, 1, y * y);
        set_pixel(S, i, j, 2, x * y);
    });
    free_image(ix);
    free_image(iy);

    Image result = GaussianBlur(S, sigma, 0);//高斯滤波
    free_image(S);

    return result;
}

Image cornerness_response(Image S) {
    Image R = make_image(S.w, S.h, 1);

    float alpha = 0.06;
    FOR_EACH_PIXEL(S.w, S.h, {
        float xx = get_pixel(S, i, j, 0);
        float yy = get_pixel(S, i, j, 1);
        float xy = get_pixel(S, i, j, 2);

        float det = xx * yy - xy * xy;
        float trace = xx + yy;
        set_pixel(R, i, j, 0, det - alpha * trace * trace);
    });

    return R;
}

Image nms_image(Image im, int w) {
    Image r = copy_image(im);
    int center_x = w / 2;
    int center_y = w / 2;

    FOR_EACH_PIXEL(r.w, r.h, {
        int cur_i = i;
        int cur_j = j;

        if (cur_i < center_x || cur_j < center_y || cur_i >= (r.w - center_x) || cur_j >= (r.h - center_y)) {
            continue;
        }
        float pixel = get_pixel(im, cur_i, cur_j, 0);
        FOR_EACH_PIXEL(w, w, {
                if (get_pixel(im, cur_i-center_x+i, cur_j-center_y+j, 0)>pixel){
                    set_pixel(r, cur_i, cur_j, 0, -999999);
                    break;
                }
        });
    });

    return r;
}

//生成描述子
Descriptor describe_index(Image im, int i) {
    int w = 5;
    Descriptor d;
    d.p.x = i % im.w;
    d.p.y = i / im.w;
    d.data = calloc(w * w * im.c, sizeof(float));
    d.n = w * w * im.c;
    int c, dx, dy;
    int count = 0;
    for (c = 0; c < im.c; ++c) {
        float cval = im.data[c * im.w * im.h + i];
        if (d.p.x < 2 || d.p.y < 2 || d.p.x >= im.w - 2 || d.p.y >= im.h - 2) {
            continue;
        }
        for (dx = -w / 2; dx < (w + 1) / 2; ++dx) {
            for (dy = -w / 2; dy < (w + 1) / 2; ++dy) {
                float val = get_pixel(im, i % im.w + dx, i / im.w + dy, c);
                d.data[count++] = cval - val;
            }
        }
    }

    return d;
}

//harris角点检测和特征提取
//Image im:输入图像
//float sigma:高斯平滑时的参数
//float thresh:角点的临界值
//int nms:进行非极大值抑制的窗口大小
//int *n:
//返回特征描述子的数组
Descriptor *harris_corner_detector(Image im, float sigma, float thresh, int nms, int *n) {
    //计算结构矩阵
    Image S = structure_matrix(im, sigma);

    //计算角点响应
    Image R = cornerness_response(S);

    //对角点进行非极大值抑制
    Image Rnums = nms_image(R, nms);

    int count = 0;
    int capacity = 64;

    int *arr = (int *) malloc(capacity * sizeof(int));
    FOR_EACH_PIXEL(Rnums.w, Rnums.h, {
        if (get_pixel(Rnums, i, j, 0) > thresh) {
            ++count;
            if (count > capacity) {
                capacity <<= 1;
                int *new_arr = (int *) malloc(capacity * sizeof(int));
                memcpy(new_arr, arr, (capacity >> 1) * sizeof(int));
                free(arr);
                arr = new_arr;
            }
            arr[count - 1] = i + j * Rnums.w;
        }
    });

    *n = count;
    Descriptor *discriptor = calloc(count, sizeof(Descriptor));
    Descriptor *ptr = discriptor;
    for (int i = 0; i < count; i++) {

        *ptr++ = describe_index(im, arr[i]);
    }

    free(arr);
    free_image(S);
    free_image(R);
    free_image(Rnums);

    return discriptor;
}

void free_descriptors(Descriptor *descriptor, int n) {
    int i;
    for (i = 0; i < n; i++) {
        free(descriptor[i].data);
    }
}

void draw_cross_point(Image im, Point2d p) {
    int x = p.x;
    int y = p.y;
    int c, i;
    for (c = 0; c < im.c; c++) {
        float pixel = 1.0;
        if (c == 1) { pixel = 0.0; }
        for (i = -9; i < 10; i++) {
            //横线
            int x_0 = x + i;
            int y_0 = y;
            if (x_0 <= 0) {
                set_pixel(im, 0, y_0, c, pixel);
            } else if (x_0 >= im.w - 1) {
                set_pixel(im, im.w - 1, y_0, c, pixel);
            } else { set_pixel(im, x_0, y_0, c, pixel); }

            //竖线
            int x_1 = x;
            int y_1 = y + i;
            if (y_1 <= 0) {
                set_pixel(im, x_1, 0, c, pixel);
            } else if (y_1 >= im.h - 1) {
                set_pixel(im, x_1, im.h - 1, c, pixel);
            } else { set_pixel(im, x_1, y_1, c, pixel); }
        }
    }
}

void draw_descriptors(Image im, Descriptor *descriptor, int n) {
    int i;
    for (i = 0; i < n; i++) {
        Point2d p = descriptor[i].p;
        draw_cross_point(im, p);
    }
}

/******************************************
             Feature matching
*******************************************/
//计算两个描述子之间的距离
float l1_distance(float *a, float *b, int n) {
    float dist = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        dist += fabs(a[i] - b[i]);
    }

    return dist;
}

int match_compare(const void *a, const void *b) {
    Match *ra = (Match *) a;
    Match *rb = (Match *) b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return 1;
    else return 0;
}

//对两幅图像中的特征点进行匹配
//Descriptor *a,*b:两张图片中的特征点描述子
//int an,bn:两张图片中各自描述子的数量
//int *mn:匹配的特征点的数量
//return:找到匹配的特征点
Match *match_descriptors(Descriptor *a, int an, Descriptor *b, int bn, int *mn) {
    int i, j;

    *mn = an; //匹配的点最多只有an个
    Match *match = calloc(an, sizeof(Match));
    for (j = 0; j < an; j++) {
        int best_index = 0;
        float best_distance = l1_distance(a[j].data, b[0].data, a[j].n);
        for (i = 1; i < bn; i++) {
            float d = 0.0;
            if (d = l1_distance(a[j].data, b[i].data, a[j].n) < best_distance) {
                best_distance = d;
                best_index = i;
            }
        }
        match[j].p = a[j].p;
        match[j].q = b[best_index].p;
        match[j].ai = j;
        match[j].bi = best_index;
        match[j].distance = best_distance;
    }

    qsort(match, an, sizeof(*match), &match_compare);
    int count = 0;
    int *seen = calloc(bn, sizeof(int));
    int *remove = calloc(an, sizeof(int));
    int m, n;
    for (n = 0; n < an; n++) {
        if (!seen[match[n].bi]) {
            count++;
            seen[match[n].bi] = 1;
        } else remove[n] = 1;
    }
    int k = -1;
    for (m = 0; m < count; m++) {
        while (remove[++k]);
        match[m] = match[k];
    }

    free(remove);
    *mn = count;
    free(seen);

    return match;
}

//将两张图片放到一张图片上便于绘制匹配点
Image both_images(Image a, Image b) {
    Image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);

    int i, j, k;
    for (k = 0; k < a.c; k++) {
        FOR_EACH_PIXEL(a.w, a.h, {
            set_pixel(both, i, j, k, get_pixel(a, i, j, k));
        });
    }

    for (k = 0; k < b.c; k++) {
        FOR_EACH_PIXEL(b.w, b.h, {
            set_pixel(both, i + a.w, j, k, get_pixel(b, i, j, k));
        });
    }

    return both;
}

Image draw_matches(Image a, Image b, Match *matches, int n, int inliers) {
    int i, j, k;
    for (k = 0; k < n; k++) {
        draw_cross_point(a, matches[k].p);
        draw_cross_point(b, matches[k].q);
    }

    Image both = both_images(a, b);

    for (i = 0; i < n; i++) {
        int ax = matches[i].p.x;
        int ay = matches[i].p.y;
        int bx = matches[i].q.x;
        int by = matches[i].q.y;

        for (j = ax; j < bx + a.w; j++) {
            int r = (float) (j - ax) * (by - ay) / (bx + a.w - ax) + ay;
            set_pixel(both, j, r, 0, i < inliers ? 0 : 1);
            set_pixel(both, j, r, 1, i < inliers ? 1 : 0);
            set_pixel(both, j, r, 2, 0);
        }
    }

    return both;
}

/******************************************
              Image stitching
*******************************************/
static Point2d project_point_fast(Matrix H, Point2d p, Matrix m) {
    assert(m.cols == 1 && m.rows == 3);
    assert(H.cols == 3 && H.rows == 3);

    m.data[0][0] = p.x;
    m.data[1][0] = p.y;
    m.data[2][0] = 1;
    Matrix n = matrix_mult_matrix(H, m);
    Point2d r = {
            .x = n.data[0][0] / n.data[2][0],
            .y = n.data[1][0] / n.data[2][0]
    };
    free_matrix(n);

    return r;
}

Point2d project_point(Matrix H, Point2d p) {
    double arr[3];
    double *data[3] = {&arr[0], &arr[1], &arr[2]};
    Matrix mat = {
            .rows = 3,
            .cols = 1,
            .data = data
    };
    Point2d r = project_point_fast(H, p, mat);
    return r;
}

float point_distance(Point2d p, Point2d q) {
    return SQUARE(p.x - q.x) + SQUARE(p.y - q.y);
}

//计算匹配中内点的个数
int model_inliers(Matrix H, Match *m, int n, float thresh) {
    int count = 0;
    Matrix mat = make_matrix(3, 1);

    int j = n;
    while (count < j) {
        if (point_distance(project_point_fast(H, m[count].p, mat), m[count].q) < thresh) {
            count++;
        } else {
            Match temp = m[--j];
            m[j] = m[count];
            m[count] = temp;
        }
    }

    return count;
}

//计算两张图片的单应性矩阵
//Match *matches：两张图片的匹配点对
//int n:前n个内点的个数
Matrix compute_homography(Match *matches, int n) {
    Matrix M = make_matrix(n * 2, 8);
    Matrix b = make_matrix(n * 2, 1);

    int i;
    for (i = 0; i < n; i++) {
        double x = matches[i].p.x;
        double xp = matches[i].q.x;
        double y = matches[i].p.y;
        double yp = matches[i].q.y;

        double arr1[8] = {x, y, 1, 0, 0, 0, -x * xp, -y * xp};
        double arr2[8] = {0, 0, 0, x, y, 1, -xp * yp, -y * yp};
        memcpy(M.data[i * 2], arr1, sizeof(arr1));
        memcpy(M.data[i * 2 + 1], arr2, sizeof(arr2));

        b.data[i * 2][0] = xp;
        b.data[i * 2 + 1][0] = yp;
    }
    Matrix a = solve_system(M, b);
    free_matrix(M);
    free_matrix(b);

    Matrix none = {0};
    if (!a.data) return none;

    Matrix H = make_matrix(3, 3);
    assert(a.cols == 1 && a.rows == 8);
    for (i = 0; i < 8; ++i) {
        H.data[i / 3][i % 3] = a.data[i][0];
    }
    H.data[2][2] = 1;
    free_matrix(a);

    return H;
}

void randomize_matches(Match *m, int n) {
    for (int i = n - 1; i >= 1; --i) {
        int j = (int) (rand() / (RAND_MAX / (i + 1.0)));
        j = j > i ? i : j;
        Match tmp = m[i];
        m[i] = m[j];
        m[j] = tmp;
    }
}

// 执行随机抽样一致来计算有噪声的匹配点的单应矩阵
// match *m：匹配点
// int n：匹配点的个数
// float thresh：内点的距离阈值
// int k：迭代次数
// int cutoff：内点计算提前中止条件
Matrix RANSAC(Match *m, int n, float thresh, int k, int cutoff) {
    int best = 0;
    Matrix Hb = make_translation_homography(256, 0);

    int fit_num = 4;
    best = model_inliers(Hb, m, n, thresh);
    if (best > cutoff) {
        free_matrix(Hb);
        return compute_homography(m, best);
    }
    for (int i = 0; i < k; i++) {
        randomize_matches(m, n);
        Matrix H = compute_homography(m, fit_num);
        int inliners = model_inliers(H, m, n, thresh);
        if (inliners > best) {
            free_matrix(Hb);
            Hb = compute_homography(m, inliners);
            best = inliners;
            if (best > cutoff) {
                return Hb;
            }
        } else {
            free_matrix(H);
        }
    }

    return Hb;
}

Point2d make_point(float x, float y) {
    Point2d p;
    p.x = x;
    p.y = y;
    return p;
}

//用投影变换将两幅图像拼接在一起
//Image a,b:图像拼接
//Matrix H:从图像a坐标到图像b坐标的单应性矩阵。
//return:拼接在一起的组合图像。
Image combine_images(Image a, Image b, Matrix H) {

    Matrix Hinv = matrix_invert(H);

    //将图像b的角投影到图像a坐标中
    Point2d c1 = project_point(Hinv, make_point(0, 0));
    Point2d c2 = project_point(Hinv, make_point(b.w - 1, 0));
    Point2d c3 = project_point(Hinv, make_point(0, b.h - 1));
    Point2d c4 = project_point(Hinv, make_point(b.w - 1, b.h - 1));

    //计算拼接后的图像大小
    Point2d topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    int i, j, k;
    Image c = make_image(w, h, a.c);

    for (k = 0; k < a.c; ++k) {
        for (j = 0; j < a.h; ++j) {
            for (i = 0; i < a.w; ++i) {
                set_pixel(c, i - dx, j - dy, k, get_pixel(a, i, j, k));
            }
        }
    }

    double arr[3];
    double *data[3] = {&arr[0], &arr[1], &arr[2]};
    Matrix mat = {
            .rows = 3,
            .cols = 1,
            .data = data
    };

    for (k = 0; k < a.c; ++k) {
        for (j = topleft.y; j < botright.y; ++j) {
            for (i = topleft.x; i < botright.x; ++i) {
                if (i - dx >= c.w || i - dy < 0 || j - dy >= c.h || j - dy < 0) {
                    continue;
                } else {
                    Point2d p = project_point_fast(Hinv, make_point(i, j), mat);
                    if (p.x >= 0 && p.x < b.w && p.y >= 0 && p.y < b.h) {
                        set_pixel(c, i - dx, j - dy, k, bilinear_interpolate(b, p.x, p.y, k));
                    }
                }
            }
        }
    }
    return c;
}

// 对两幅图像进行拼接
// mage a, b: 输入图像
// float sigma: 高斯平滑参数。一般取2
// float thresh: 是否为角点的阈值。一般取1-5
// int nms: 非极大值抑制遍历窗口大小. 一般取3
// float inlier_thresh: RANSAC内点的阈值。一般取2-5
// int iters: RANSAC迭代次数。一般取1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
Image panorama_image(Image a, Image b, float sigma, float thresh,
                     int nms, float inlier_thresh, int iters, int cutoff) {
    srand(10);
    int an = 0, bn = 0, mn = 0;

    //计算角点和描述子
    Descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    Descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    //进行角点匹配
    Match *m = match_descriptors(ad, an, bd, bn, &mn);

    //执行RANSAC寻找单应行矩阵
    Matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);

    //图像拼接
    Image comb = combine_images(a, b, H);

    return comb;
}

static Point2d* cylinder_point(float xc, float yc, float f, Point2d* p) {
    double theta = ((double)p->x - xc) / f;
    double xp = sin(theta);
    double yp = ((double)p->y - yc) / f;
    double zp = cos(theta);
    p->x = f * xp / zp + xc;
    p->y = f * yp / zp + yc;
    return p;
}

// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
Image cylindrical_project(Image im, float f)
{
    //TODO: project image onto a cylinder
    int xc = im.w / 2, yc = im.h / 2;
    int lx = -1, rx = im.w, px;
    do {
        double theta = ((double)++lx - xc) / f;
        px = f * sin(theta) / cos(theta) + xc;
    } while(px < 0);
    do {
        double theta = ((double)--rx - xc) / f;
        px = f * sin(theta) / cos(theta) + xc;
    } while(px >= im.w);
    int ty = -1, by = im.h, py;
    do {
        double yp = ((double)++ty - yc) / f;
        py = f * yp + yc;
    } while(py < 0);
    do {
        double yp = ((double)--by - yc) / f;
        py = f * yp + yc;
    } while(py >= im.h);
    Image c = make_image(rx - lx + 1, by - ty + 1, im.c);
    Point2d tmp;
    for (int k = 0; k < c.c; ++k) {
        FOR_EACH_PIXEL(c.w,c.h, {
                tmp.x = i + lx;
                tmp.y = j + ty;
                cylinder_point(xc, yc, f, &tmp);
                set_pixel(c, i, j, k, bilinear_interpolate(im, tmp.x, tmp.y, k));
        });
    }
    return c;
}

