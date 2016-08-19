/** @file
********************************************************************************
<PRE>
模块名: 基二时域快速傅立叶变换
文件名: FFT2.h
相关文件: FFT2.c
文件实现功能: 声明基二时域快速傅立叶变换相关函数
作者: Dake
版本: V2010.09.01
编程方式: ANSI C语言编程
授权方式: Copyright(C) Dake
联系方式: chen423@yeah.net
生成日期: 2010-06-22
--------------------------------------------------------------------------------
多线程安全性: <是/否>[，说明]
异常时安全性: <是/否>[，说明]
--------------------------------------------------------------------------------
备注: <其它说明>
--------------------------------------------------------------------------------
修改记录:
日 期        版本     修改人              修改内容
YYYY/MM/DD   X.Y      <作者或修改者名>    <修改内容>
</PRE>
*******************************************************************************/
#ifndef _FFT2_
#define _FFT2_

#include <cmath>
#include <ctime>
#include <iostream>
#include <complex>

using namespace std;
template<typename _Tp = float>
class FFT2 {
    public:
        FFT2(int N);
        ~FFT2();
        void doFFT2(complex<_Tp>  xin[]);
        void doFFT2(_Tp *xin,complex<_Tp>  *xout);
    private:
        int lgn;
        int dataNumber;
        float *SinTable;
        float *CosTable;
};

//static float *SinTable;
//static float *CosTable;
#ifdef PI
    #undef PI
#endif

#define PI 3.1415926

template<typename _Tp> 
FFT2<_Tp>::FFT2(int N)
{
    SinTable = new _Tp[N];
    CosTable = new _Tp[N];
    for(int i=0;i<N;i++) {
        SinTable[i] = sin(2.0*i*PI/N);
        CosTable[i] = cos(2.0*i*PI/N);
    }
    switch(N) {
        case 256:
        dataNumber = 256;
        lgn = 8;
        break;
        case 512:
        dataNumber = 512;
        lgn = 9;
        break;
    }
}
template<typename _Tp> 
FFT2<_Tp>::~FFT2() {
    delete SinTable;
    delete CosTable;
}


template<typename _Tp> 
void FFT2<_Tp>::doFFT2(complex<_Tp>  xin[])
{
    int LH, nm, k, j;
    complex<_Tp> w, t;
    
    LH = dataNumber / 2;
    nm = dataNumber-2;
    j = LH;
    /* 变址运算 */
    for (int i = 1; i <= nm; ++i)
    {
        if (i < j)
        {
            t = xin[j];
            xin[j] = xin[i];
            xin[i] = t;
        }
        k = LH;
        while (j >= k)
        {
            j = j - k;
            k = k / 2;
        }
        j = j + k;
    }

    {
        int B, ip;
        int le = 1;
        int pe = dataNumber;
        for (int L = 1; L<=lgn; L++)
        {
            B = le;
            le <<= 1;
            pe >>= 1;
            int p;
            for (j = 0; j < B; j++)
            {
                p = pe * j;
                w.real(CosTable[p]);
                w.imag(SinTable[p]);
                for (int i = j; i <= dataNumber-1; i = i + le)
                {
                    ip = i + B;
                    t = xin[ip] * w;
                    xin[ip] = xin[i] - t;
                    xin[i]  = xin[i] + t;
                }
            }
        }
    }
    return;
}

template<typename _Tp> 
void FFT2<_Tp>::doFFT2(_Tp *xin,complex<_Tp>  *xout)
{
    complex<_Tp>  *t = xout;
    for(int i=0;i<dataNumber;i++) {
        t->real(xin[i]);
        t->imag(0);
        t++;
    }
    doFFT2(xout);
}


//void FFT2(COMPX * xin, int N);

#endif    //_FFT2_
