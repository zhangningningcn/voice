#ifndef  MFCC_H  
#define  MFCC_H  
#include <stdio.h>  
#include "FFT2.hpp"
  
typedef struct AudioInfo_  
{  
    int sample_rate;  
    int frame_len;  
    int frame_shift;  
    float* data_in;  
    int data_len;     
    //AudioInfo():sample_rate(0),frame_len(0),frame_shift(0),data_in(NULL),data_len(0){};
    //构造函数注释掉，变成纯C版本  
}AudioInfo;  
typedef enum MFCC_TYPE_ 
{  
    MFCC_STD,  
    MFCC_DIFF_1,  
    MFCC_DIFF_2,  
}MFCC_TYPE;  
class MFCC {
    public:
            // (AudioInfo audioctx, int low, int high, int nfilters, int ndcts) 
        MFCC(AudioInfo &audioctx, int low, int high, int nfilters, int ndcts); 
        int Mfcc_Frame_std(float *in,float *out);//输出mfcc，任意帧输出  
        //int Mfcc_Frame_diff1(int iframe, float *out, int len);  
        //int Mfcc_Frame_diff2(int iframe, float *out, int len);  
        ~MFCC();  
    private:
        typedef struct MelBankInfo_  
        {  
            float **filter;  
            int     nfilters;  
            int     nfft;  
            int     low;  
            int     high;  
            //MelBankInfo():filter(NULL),nfilters(0),nfft(0),low(0),high(0){};  
        }MelBankInfo;  
        typedef struct DctInfo_  
        {  
            float  **coeff;  
            int     dctlen;  
            //DctInfo():coeff(NULL),dctlen(0){};  
        }DctInfo;  
        typedef struct MfccInfo_  
        {       
            MelBankInfo   melbank;  
            DctInfo       dct;  
          
            int     nframes;  
            int     out_nframes;//可输出的特征数  
            float *frame_data;  
            float *data_in;  
            float  *window;  
            float  *window_data;
            float  *lift_window;      
            int     frame_len;  
            int     frame_shift;  
            float preF;
            #if 0
            float  *pre1;  
            float  *pre2;  
            float  *cur;  
            float  *next1;  
            float  *next2;  
          
            float  *diff_pre1;  
            float  *diff_pre2;  
            float  *diff_cur;  
            float  *diff_next1;  
            float  *diff_next2;  
            #endif
            //MFCC_TYPE m_type;  
            //MfccInfo():nframes(0),out_nframes(0),window(NULL),lift_window(NULL),frame_len(0),
            //fft(NULL),window_data(NULL),frame_data(NULL),frame_shift(0){};
            //pre1(NULL),pre2(NULL),cur(NULL),next1(NULL),next2(NULL),m_type(MFCC_STD),;  
        }MfccInfo;  
    private:
        MfccInfo info;
        FFT2<float> *fft;
        complex<float> *fft_outBuffer; // = new [info.frame_len];
        float *mal_Buffer;  //[info.frame_len];
        void lift_window(float* p, int m);
        void DctCoeff( int m, int n, float **coeff );
        void PreEmphasise(const float *data, int len, float *out, float preF);//预加重
        float HzToMel(float  f);
        float MelToHz(float data);
        int  HzToN(float f, int fs, int nfft);
        void MelBank( int fs, int nfft, int low, int high, int nfilters, float **coeff );
        float Product(float *data1, float* data2, int len);
        float **MallocMatrix(int m, int n);
        void FreeMatrix(float **in);
        void hamming( float *win, int N);
        void hanning( float *win, int N);
        void apply_window(float* out,float* in, float* window, int window_len);
};
  
  
  
  
  
  
//MfccInfo*  MfccInit(AudioInfo audioctx, int nfft, int low, int high, int nfilters, int ndcts, MFCC_TYPE type);  
  
#endif  