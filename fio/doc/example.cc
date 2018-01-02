// `包含` fio `头文件`
#include <fio.h>
#include <iostream>
#include <cmath>

// `命名空间`
using namespace std;
using namespace blitz;

// `主程序开始`
int main(int argc, char *argv[])
// argc `和` argv `用于从命令行获取参数`
// argc `是命令行参数的数目，包含程序名称`
// argv[0], argv[1], ..., argv[argc-1] `是各个参数字符串`
{
    if (3 != argc) {        // `判断输入参数个数是否符合要求`
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "     %s in.fits out.fits\n", argv[0]);
        exit(-1);           // `强制程序结束并返回错误代码` －1
    }

    cfitsfile ff1;          // `声明一个用于操作` fits `文件的对象`
    ff1.open(argv[1]);      // `打开一个已经存在的` fits `文件`
    Array<double,2> img1;   // `声明一个用于存放图像的矩阵`
    ff1 >> img1;            // `将数据从` fits `文件导入矩阵`

    // `然后可以对这个矩阵中的数据进行各种处理`
    // `比如说，这里计算输入图像的泊松误差`
    Array<double,2> img2(img1.shape());
    int i, j;
    for (i=0; i<img2.extent(0); ++i) {
        for (j=0; j<img2.extent(1); ++j) {
            img2(i,j) = sqrt(img1(i,j));
        }
    }

    cfitsfile ff2;          // `声明另一个用于操作` fits `文件的对象`
    ff2.create(argv[2]);    // `新建一个` fits `文件`
                            // `如果文件存在，则先删除再新建`
    ff2 << img2;            // `将数据从矩阵导出到` fits `文件`

    return 0;               // `程序正常结束`
}

