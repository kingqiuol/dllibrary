#ifndef IMAGE_H
#define IMAGE_H

//����ͼ��ṹ��
typedef struct {
    int h, w, c;
    float *data;
} Image;

//Basic operations
float get_pixel(Image im, int x, int y, int c);//��ȡ����ֵ
void set_pixel(Image im, int x, int y, int c, float v);//��������ֵ
void add_pixel(Image im, int x, int y, int c, float v);//Ԫ��ֵ���
void fill_image(Image im, float s);//���Ԫ��
Image add_image(Image a, Image b);//ͼ�����
Image sub_image(Image a, Image b);//ͼ�����
Image copy_image(Image p);//����ͼ��
void constrain_image(Image im);//����ֵԽ��Լ��������1������Ϊ1��С��0����Ϊ0
void print_image(Image im);//��ӡͼƬ

//loading image and saving and showing
Image make_image(int w, int h, int c);//����w*h*c��ͼ��
Image load_image(char *filename);//����ͼ��
void show_images(Image *im, int n, char *window);//��ʾͼ��
void save_image(Image im, const char *name);//����ͼ��
void free_image(Image im);//�ͷ�ͼ��

//color space transformation 
Image rgb_to_gray(Image im);//rgbת�Ҷ�ͼ��
Image rgb_to_hsv(Image im);//rgbתhsvͼ��
Image rgb_to_yuv(Image im);//rgbתyuvͼ��
Image hsv_to_rgb(Image im);//hsvתrgbͼ��
Image yuv_to_rgb(Image im);//yuvתrgbͼ��

//image enhancement 
void scale_image(Image im, float s);//���ԶԱȶ���ǿ
void translate_image(Image m, float s);//��������ƫ��s(����ֵ��s)

//image pretreatment
Image random_augment_image(Image im, float angle, float aspect, int low, int high, int w, int h);//���ͼ����ǿ

Image crop_image(Image im, int dx, int dy, int w, int h);//����ͼ��
Image center_crop_image(Image im, int w, int h);//���ļ���
Image random_crop_image(Image im, int w, int h);//�������ͼ��
Image rotate_image(Image im, float rad);//��תͼ��
void rotate_image_cw(Image im, int times);//��תͼ��
Image flip_image(Image im);//���ҷ�תͼ��
Image rotate_crop_image(Image im, float rad, float s, int w, int h, float dx, float dy, float aspect);//��ת����

Image saturate_image(Image im, float sat);//����ͼ��
Image exposure_image(Image im, float val);//�ع�ͼ��
Image saturate_exposure_image(Image im, float sat, float exposure);//�����ع�ͼ��
Image random_distort_image(Image im, float hue, float saturation, float exposure);//���Ť��ͼ��

float bilinear_interpolate(Image im, float rx, float ry, int c);
Image resize_image(Image im, int w, int h);//����ͼ���С
Image resize_min(Image im, int min);//������С�ߵĴ�С������һ���ȱ�������
Image resize_max(Image im, int max);//������ߴ�С������һ���ߵȱ�������

void place_image(Image im, int w, int h, int dx, int dy, Image canvas);//ͼ��ƽ����ǿ

#endif
