#ifndef _SPATIAL_FILTER_H_
#define _SPATIAL_FILTER_H_

#include"image.h"

//生成滤波器
Image make_box_filter(int w);//均值滤波器
Image make_highpass_filter();//高通滤波器
Image make_sharpen_filter();//锐化滤波器 
Image make_emboss_filter();

Image make_gaussian_filter(float sigma);//高斯滤波器
//sobel滤波器
Image make_gx_filter();

Image make_gy_filter();

//实现图像增强之空间滤波
Image boundary_expansion(Image im, int filter_size, int flag);//边界填充
Image convolve_image(Image im, Image filter, int preserve, int flag);//卷积操作
Image blur(Image im, int w, int flag);//均值滤波
Image GaussianBlur(Image im, float sigma, int flag);//高斯滤波
Image HighpassBlur(Image im, int flag);//高通滤波
Image SharpenBlur(Image im, float flag);//锐化滤波
Image EmbossBlur(Image im, int flag);//emboss滤波

Image *Sobel(Image im, int flag);//sobel检测

#endif
