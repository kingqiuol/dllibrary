#ifndef PANORAMA_IMAGE_H
#define PANORAMA_IMAGE_H

#include "matrix.h"

//2d points
typedef struct {
    float x, y;
} Point2d;

//����������
typedef struct {
    Point2d p;
    int n;
    float *data;
} Descriptor;

//ƥ���
//Point2d p q:����ͼƬƥ��������
//int ai,bi:�������е�����
//float distance:ƥ���������֮��ľ���
typedef struct {
    Point2d p, q;
    int ai, bi;
    float distance;
} Match;

//������ļ��
Image structure_matrix(Image im, float sigma);//����ͼ��Ľṹ����
Image cornerness_response(Image S);//����ǵ�
Image nms_image(Image im, int w);//�Խǵ���зǼ���ֵ����
Descriptor *harris_corner_detector(Image im, float sigma, float thresh, int nms, int *n);//Harris�ǵ���
void free_descriptors(Descriptor *descriptor, int n);//�ͷŽǵ�ռ�
void draw_descriptors(Image im, Descriptor *descriptor, int n);//���ƽǵ�


//�������ƥ��
Match *match_descriptors(Descriptor *a, int an, Descriptor *b, int bn, int *mn);//����ƥ��
Image draw_matches(Image a, Image b, Match *matches, int n, int inliers);//����ƥ�������

//ͼ��ƴ��
Matrix RANSAC(Match *m, int n, float thresh, int k, int cutoff);
Image panorama_image(Image a, Image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff);
Image cylindrical_project(Image im, float f);
#endif

