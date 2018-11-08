#ifndef _SPATIAL_FILTER_H_
#define _SPATIAL_FILTER_H_

#include"image.h"

//�����˲���
Image make_box_filter(int w);//��ֵ�˲���
Image make_highpass_filter();//��ͨ�˲���
Image make_sharpen_filter();//���˲��� 
Image make_emboss_filter();

Image make_gaussian_filter(float sigma);//��˹�˲���
//sobel�˲���
Image make_gx_filter();

Image make_gy_filter();

//ʵ��ͼ����ǿ֮�ռ��˲�
Image boundary_expansion(Image im, int filter_size, int flag);//�߽����
Image convolve_image(Image im, Image filter, int preserve, int flag);//�������
Image blur(Image im, int w, int flag);//��ֵ�˲�
Image GaussianBlur(Image im, float sigma, int flag);//��˹�˲�
Image HighpassBlur(Image im, int flag);//��ͨ�˲�
Image SharpenBlur(Image im, float flag);//���˲�
Image EmbossBlur(Image im, int flag);//emboss�˲�

Image *Sobel(Image im, int flag);//sobel���

#endif
