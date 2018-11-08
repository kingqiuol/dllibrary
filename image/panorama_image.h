#ifndef PANORAMA_IMAGE_H
#define PANORAMA_IMAGE_H

#include "matrix.h"

//2d points
typedef struct {
    float x, y;
} Point2d;

//特征描述子
typedef struct {
    Point2d p;
    int n;
    float *data;
} Descriptor;

//匹配点
//Point2d p q:两张图片匹配点的坐标
//int ai,bi:描述子中的索引
//float distance:匹配点描述子之间的距离
typedef struct {
    Point2d p, q;
    int ai, bi;
    float distance;
} Match;

//特征点的检测
Image structure_matrix(Image im, float sigma);//计算图像的结构矩阵
Image cornerness_response(Image S);//计算角点
Image nms_image(Image im, int w);//对角点进行非极大值抑制
Descriptor *harris_corner_detector(Image im, float sigma, float thresh, int nms, int *n);//Harris角点检测
void free_descriptors(Descriptor *descriptor, int n);//释放角点空间
void draw_descriptors(Image im, Descriptor *descriptor, int n);//绘制角点


//特征点的匹配
Match *match_descriptors(Descriptor *a, int an, Descriptor *b, int bn, int *mn);//特征匹配
Image draw_matches(Image a, Image b, Match *matches, int n, int inliers);//绘制匹配点连线

//图像拼接
Matrix RANSAC(Match *m, int n, float thresh, int k, int cutoff);
Image panorama_image(Image a, Image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff);
Image cylindrical_project(Image im, float f);
#endif

