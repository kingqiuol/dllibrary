#ifndef IMAGE_H
#define IMAGE_H

//定义图像结构体
typedef struct {
    int h, w, c;
    float *data;
} Image;

//Basic operations
float get_pixel(Image im, int x, int y, int c);//获取像素值
void set_pixel(Image im, int x, int y, int c, float v);//设置像素值
void add_pixel(Image im, int x, int y, int c, float v);//元素值相加
void fill_image(Image im, float s);//填充元素
Image add_image(Image a, Image b);//图像相加
Image sub_image(Image a, Image b);//图像相减
Image copy_image(Image p);//拷贝图像
void constrain_image(Image im);//像素值越界约束；大于1的设置为1，小于0设置为0
void print_image(Image im);//打印图片

//loading image and saving and showing
Image make_image(int w, int h, int c);//生成w*h*c的图像
Image load_image(char *filename);//加载图像
void show_images(Image *im, int n, char *window);//显示图像
void save_image(Image im, const char *name);//保存图像
void free_image(Image im);//释放图像

//color space transformation 
Image rgb_to_gray(Image im);//rgb转灰度图像
Image rgb_to_hsv(Image im);//rgb转hsv图像
Image rgb_to_yuv(Image im);//rgb转yuv图像
Image hsv_to_rgb(Image im);//hsv转rgb图像
Image yuv_to_rgb(Image im);//yuv转rgb图像

//image enhancement 
void scale_image(Image im, float s);//线性对比度增强
void translate_image(Image m, float s);//像素整体偏移s(像素值加s)

//image pretreatment
Image random_augment_image(Image im, float angle, float aspect, int low, int high, int w, int h);//随机图像增强

Image crop_image(Image im, int dx, int dy, int w, int h);//剪切图像
Image center_crop_image(Image im, int w, int h);//中心剪切
Image random_crop_image(Image im, int w, int h);//随机剪切图像
Image rotate_image(Image im, float rad);//旋转图像
void rotate_image_cw(Image im, int times);//旋转图像
Image flip_image(Image im);//左右翻转图像
Image rotate_crop_image(Image im, float rad, float s, int w, int h, float dx, float dy, float aspect);//旋转剪切

Image saturate_image(Image im, float sat);//饱和图像
Image exposure_image(Image im, float val);//曝光图像
Image saturate_exposure_image(Image im, float sat, float exposure);//饱和曝光图像
Image random_distort_image(Image im, float hue, float saturation, float exposure);//随机扭曲图像

float bilinear_interpolate(Image im, float rx, float ry, int c);
Image resize_image(Image im, int w, int h);//调整图像大小
Image resize_min(Image im, int min);//调整最小边的大小，另外一个等比例缩放
Image resize_max(Image im, int max);//调整最长边大小，另外一个边等比例缩放

void place_image(Image im, int w, int h, int dx, int dy, Image canvas);//图像平移增强

#endif
