#include<stdio.h>
#include<stdlib.h>

#include"image.h"
#include"spatial_filter.h"
#include"panorama_image.h"

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("请输入图片名字");
        return -1;
    }
    Image img = load_image(argv[1]);
    printf("h= %d,w=%d,c=%d\n", img.h, img.w, img.c);
    //print_image(img);
    //save_image(img,"temp");
    //Image gray=rgb_to_gray(img);
    //save_image(gray,"gray");
    //Image rgb_hsv=rgb_to_hsv(img);
    //save_image(rgb_hsv,"rgb_to_hsv");
    //Image rgb_yuv=rgb_to_yuv(img);
    //save_image(rgb_yuv,"rgb_to_yuv");
    //Image yuv_rgb=yuv_to_rgb(rgb_yuv);
    //save_image(yuv_rgb,"yuv_to_rgb");
    //Image hsv_rgb=hsv_to_rgb(rgb_hsv);
    //save_image(hsv_rgb,"hsv_to_rgb");

    //scale_image(gray,0.5);
    //save_image(gray,"scale_image");

    //Image random_crop_img=random_crop_image(img,150,250);
    //save_image(random_crop_img,"random_crop_image");

    //Image dist_img=random_distort_image(img, 1,2,2);
    //save_image(dist_img,"random_disttort_image");

    //Image sat_img=saturate_image(img,0.5);
    //save_image(sat_img,"saturate_image");
    //Image val_img=exposure_image(img,0.5);
    //save_image(val_img,"exposure_image");
    //Image s_e_img=saturate_exposure_image(img,2, 2);
    //save_image(s_e_img,"saturate_exposure_image");

    //Image crop_img=crop_image(img, 0, 0, 50, 50);
    //save_image(crop_img,"crop_image");
    //Image center_crop_img=center_crop_image(img,250, 250);
    //save_image(center_crop_img,"center_crop_image");

    //#define PI 3.14159265
    //Image rot_img=rotate_image(img, (90*PI)/180);
    //save_image(rot_img,"rotate_image");

    //Image img_1=load_image("./center_crop_image.jpg");
    //rotate_image_cw(img_1, 10);
    //save_image(img_1,"rotate_image_cw");

    //Image flip_img=flip_image(img);
    //save_image(flip_img,"flip_image");

    //Image rc_img=rotate_crop_image(img, (60*PI)/180, 5.0, 250, 250, 10,10, 0.5);
    //save_image(rc_img,"rc_image");

    //Image resized_img=resize_image(img,500,350);
    //save_image(resized_img,"resized_image");

    //Image canvas=make_image(100,100,3);
    //place_image(img,50,50,10,10,canvas);
    //save_image(canvas,"place_image");

    //Image blur_img=blur(img,11,2);//均值滤波
    //save_image(blur_img,"blur");
    //Image gauss_image=GaussianBlur(img,3,2);//高斯滤波
    //save_image(gauss_image,"GaussianBlur");
    //Image hi_image=HighpassBlur(img,2);//高通滤波
    // save_image(hi_image,"HighpassBlur");
    //Image sharpen_img=SharpenBlur(img,2);
    //save_image(sharpen_img,"SharpenBlur");
    //Image emboss_img=EmbossBlur(img,2);
    //save_image(emboss_img,"EmbossBlur");
    //Image* sobel_img=Sobel(img,2);
    //save_image(sobel_img[0],"Sobel");
    //save_image(sobel_img[1],"Sobel_angle");

//    Image smat_img=structure_matrix(img, 1.5);
//    save_image(smat_img,"struct_matrix_img");
//    Image r_img=cornerness_response(smat_img);
//    save_image(r_img,"cor_response");
//    Image nms_img=nms_image(r_img, 3);
//    print_image(nms_img);
//    save_image(nms_img,"nms_image");

//    int n = 0;
//    Descriptor *p = harris_corner_detector(img, 2, 0.01, 3, &n);
//
//    Image img2 = load_image(argv[2]);
//    printf("h= %d,w=%d,c=%d\n", img2.h, img2.w, img2.c);
//    int m = 0;
//    Descriptor *q = harris_corner_detector(img2, 2, 0.0001, 7, &m);
//
//    int mn;
//    Match *matches = match_descriptors(p, n, q, m, &mn);
//
//    Image match_img = draw_matches(img, img2, matches, mn, 0);
//    save_image(match_img, "match_image");
//    Image img3=cylindrical_project(img,1200);
//    Image img4=cylindrical_project(img2,1200);
//    Image pan_img=panorama_image(img3, img4, 2, 0.0001, 7, 2, 500, 30);
//    save_image(pan_img,"panorama_image");




    return 0;
}
