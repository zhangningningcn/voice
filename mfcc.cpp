/*
 * 在别人代码基础上修改，原作者未知。
 * 一阶差分和二阶差分未实现
 */

#include  "mfcc.hpp"  
//#include  "Spl.h"  
//#include  "fftw3.h"  
#include  <malloc.h>  
#include  <math.h>  
#include  <string.h>  
#include  <assert.h>  
  
#ifndef TWOPI
#define  PI   3.1415926  
#endif
#define  EPS  0.0000001  
#ifndef TWOPI
#define  TWOPI   6.283185307179586  
#endif
#include  "FFT2.hpp"  
#include  <vector>  
#include <complex>
//#pragma comment(lib, "libfftw3f-3.lib")  

  
using namespace std;  
void MFCC::PreEmphasise(const float *data, int len, float *out, float preF)//预加重  
{
    static float pre_last = 0.f;
    float temp,temp2;
    temp= pre_last;
    pre_last = data[len-1];
    for(int i = 0; i < len; i++)  
    {
        temp2 = data[i];
        out[i] = temp2 - preF * temp;//preF * data[i-1];  
        temp = temp2;
    }  
}
  
float MFCC::HzToMel(float  f)  
{  
    return 1127*log(1.0 + f/700);  
}
float MFCC::MelToHz(float data)  
{  
    return  700 * (exp(data/1127) - 1);  
}  
int  MFCC::HzToN(float f, int fs, int nfft)  
{  
    return  f/fs *nfft+1;  
}  
void MFCC::MelBank( int fs, int nfft, int low, int high, int nfilters, float **coeff )//三角滤波器组。  
{  
    float  fre_bin = (float)fs / nfft;  
    float  low_mel = HzToMel(low);  
    float  high_mel = HzToMel(high);  
    float  mel_bw  = (high_mel - low_mel)/(nfilters + 1);  
    int  valid_nfft = nfft/2 + 1;  
      
  
    for(int j = 0; j < nfilters; j++)  
    {     
        float  mel_cent  = (j+1) * mel_bw + low_mel;  
        float  mel_left  = mel_cent - mel_bw;  
        float  mel_right = mel_cent + mel_bw;  
        float  freq_cent =  MelToHz(mel_cent);  
        float  freq_left =  MelToHz(mel_left);  
        float  freq_bw_left = freq_cent - freq_left;  
        float  freq_right = MelToHz(mel_right);  
        float  freq_bw_right = freq_right - freq_cent;  
        for(int i = 0; i < valid_nfft; i++)  
        {             
            float freq = (i) * fre_bin ;  
            if( freq > freq_left && freq < freq_right )  
            {  
                if( freq <= freq_cent)  
                {  
                    //为什么不是 
                    coeff[j][i] = 2*(freq - freq_left) / ((freq_right - freq_left)*freq_bw_left);// ？
                    //coeff[j][i] = (freq - freq_left) / freq_bw_left;  
                }  
                else  
                {  
                    coeff[j][i] = 2*(freq_right - freq) / ((freq_right - freq_left)*freq_bw_right);// ？
                    //coeff[j][i] = (freq_right - freq) / freq_bw_right;  
                }  
  
            }
            else {
                coeff[j][i] = 0;
            }
              
  
        }     
          
    }     
      
}  
  
void MFCC::DctCoeff( int m, int n,float **coeff )//标准DCT变换。  
{  
    printf("DctCoeff:\n");
    for( int i = 0; i < m; i++)  
    {  
        for(int j = 0; j < n; j++)  
        {  
            //coeff[i][j] = cos( (2*j + 1) * i *PI / (2 * n));  
            coeff[i][j] = cos( (2*j + 1) * (i+1) *PI / (2 * n));  
            printf("%f,",coeff[i][j]);
        }  
        printf("\n");
    }  
}  
  
void MFCC::lift_window(float* p, int m)//倒谱提升窗归一化。  
{  
    printf("lift_window\n");
    float  max_value = 0.0f;  
    for(int i = 0; i < m; i++)  
    {  
        p[i] = 1+ 0.5 * m * sin( PI * (i+1)/m );  
        if( p[i] > max_value)  
        {  
            max_value = p[i];  
        }  
    }  
    for(int i = 0; i < m; i++)  
    {  
        p[i] /= max_value;  
        printf("%f,",p[i]);
    }  
    printf("\n");
}  
  
float MFCC::Product(float *data1, float* data2, int len)  
{  
    float result = 0.0;  
    for(int i = 0; i < len; i++)  
    {  
        result += data1[i] * data2[i];  
    }  
    return result;  
}  
  
float **MFCC::MallocMatrix(int m, int n)
{  
    float **in = new float *[m];//(float**)malloc(m * sizeof(float*));  
    float* data = new float[m*n];//(float*)malloc( m*n*sizeof(float));  
    memset( data, 0, sizeof(float)*m*n );  
    for(int i = 0; i < m; i++)  
    {  
        in[i] = &data[i*n];  
    }  
    return in;  
}  

void MFCC::FreeMatrix(float **in)
{    
    float *data = *in;  
    if(data != NULL)  
    {  
        delete data;  
    }  
    if(in != NULL)  
    {  
        delete in;  
    }  
  
  
}  
//int Mfcc_Frame_diff1_temp(MfccInfo *p, int iframe, float *out, int len);  
  
//初始化，预加重，获取滤波器组系数，DCT系数，倒谱提升窗系数等。  
MFCC::MFCC(AudioInfo &audioctx, int low, int high, int nfilters, int ndcts)  
{           
    //MfccInfo*  info = (MfccInfo*)malloc(sizeof(MfccInfo));  
    if(audioctx.frame_len != 256 && audioctx.frame_len != 512) {
        return;
    }
    info.melbank.nfft = audioctx.frame_len;  
    info.melbank.low  = low;  
    info.melbank.high = high;  
    info.melbank.nfilters = nfilters;  
    info.dct.dctlen = ndcts;  
    // info.pre1 = NULL;  
    // info.pre2 = NULL;  
    // info.cur  = NULL;  
    // info.next1 = NULL;  
    // info.next2 = NULL;  
    // info.m_type = type;  
    //info.data_in = audioctx.data_in;//整段语音的数据流  
    info.frame_shift = audioctx.frame_shift;  
    int valid_nfft = audioctx.frame_len/2 + 1;  
    info.melbank.filter = MallocMatrix( nfilters, valid_nfft);  
    MelBank( audioctx.sample_rate, audioctx.frame_len, low, high, nfilters, info.melbank.filter);//Mel滤波器系数       
    info.dct.coeff = MallocMatrix( ndcts, nfilters);  
    DctCoeff( ndcts, nfilters, info.dct.coeff );//DCT系数  
  
    info.preF = 0.9375;  
    //整段语音高通滤波，预加重   
    //PreEmphasise( audioctx.data_in, audioctx.data_len, audioctx.data_in, preF);  
    //int nframes = (audioctx.data_len - audioctx.frame_len)/audioctx.frame_shift + 1;      
    //info.nframes = nframes;  
    //info.out_nframes = nframes;  
    info.frame_len = audioctx.frame_len;  
    info.window  = new float[audioctx.frame_len];//(float*) malloc( audioctx.frame_len * sizeof(float));  
    hamming( info.window, audioctx.frame_len);//加窗  
    info.lift_window = new float[ndcts];   //(float*)malloc( ndcts * sizeof(float));  
    lift_window(info.lift_window, ndcts);  //倒谱提升窗  
    //int  buffer_len = audioctx.frame_len > nfft ? audioctx.frame_len:nfft;  
    info.frame_data = new float[audioctx.frame_len]; //(float*) malloc( buffer_len * sizeof(float));    
    info.window_data = new float[audioctx.frame_len]; //  
    //return p;  
    fft = new FFT2<float>(info.frame_len);
    fft_outBuffer = new complex<float>[info.frame_len];
    mal_Buffer = new float[info.frame_len];
    //return;
}  
int MFCC::Mfcc_Frame_std(float *in,float *out)//输出MFCC，任意帧输出  
{
    float eng;
    float *p_newdata = info.frame_data+(info.frame_len - info.frame_shift);
    memcpy(info.frame_data,info.frame_data+info.frame_shift,(info.frame_len - info.frame_shift)*sizeof(float));
    memcpy(p_newdata,in,info.frame_shift*sizeof(float));

    PreEmphasise(p_newdata, info.frame_shift, p_newdata, info.preF);  
    apply_window(info.window_data,info.frame_data, info.window, info.frame_len);  
  

    int  nfft = info.melbank.nfft;  
    int  valid_nfft = nfft/2 + 1;  
  
    // 用new 防止栈溢出
    //complex<float> *temp;
    //temp = new complex<float>[info.frame_len];
    
    //fftwf_plan r2cP;  
    //fftwf_complex* temp = (fftwf_complex*)fftwf_malloc(sizeof( fftwf_complex ) * valid_nfft);         
    //r2cP = fftwf_plan_dft_r2c_1d( info.frame_len, info.frame_data, temp, FFTW_ESTIMATE ); //完成FFT运算  
    //ftwf_execute( r2cP ); 
    fft->doFFT2(info.window_data,fft_outBuffer);
    
    eng = 0;
    float *pw = info.window_data;
    complex<float> *t2 = fft_outBuffer;
    for (int j = 0; j < valid_nfft; ++j)  
    {  
        //info.frame_data[j] = pow( temp[j][0], 2 ) + pow( temp[j][1], 2 );//平方能量值，也可以用谱幅度值  
        *pw = t2->real()*t2->real() + t2->imag()*t2->imag();
        eng += *pw;
        pw ++;
        t2 ++;
    }

    *out = eng;
    out++;
    //fftwf_destroy_plan( r2cP );  
      
    for(int j = 0; j < info.melbank.nfilters; j++)  
    {  
        //解卷积  
        mal_Buffer[j] = log(Product(info.window_data, info.melbank.filter[j], valid_nfft)+ EPS)/log(10);             
    }
    for(int i = 0; i < info.dct.dctlen; i++)  
    {  
        float temp = 0.0;  
        for(int j = 0; j < info.melbank.nfilters; j++)  
        {  
            //DCT变换
            temp += info.dct.coeff[i][j] * mal_Buffer[j];             
        }
        *out = temp;// * info.lift_window[i];//倒谱提升  
        out++;
    }  
      
    //delete temp;
  
    return 0;  
}  
#if 0
int MFCC::Mfcc_Frame_diff1(MfccInfo *p, int iframe, float *out, int len)//标准一阶差分，输出 MFCC + 一阶差分。 逐帧输出  
{  
   assert(p->nframes >= 5 && iframe <= p->nframes -4 && p->m_type == MFCC_DIFF_1);  
   int ret = Mfcc_Frame_std(p, iframe + 4, p->next2, len);     
   int dctlen = p->dct.dctlen;  
   memcpy( out, p->cur, sizeof(float)* dctlen);//MFCC  
   float  factor = sqrt(10.0);  
   for(int i = 0; i < dctlen; i++)  
   {  
       out[i + dctlen] = (2 * p->next2[i] + p->next1[i] - 2*p->pre1[i] - p->pre2[i])/factor ;//一阶差分  
   }  
  
   float *temp = p->pre1;  
   p->pre1 = p->pre2;     
   p->pre2 = p->cur;  
   p->cur  = p->next1;  
   p->next1 = p->next2;  
   p->next2 = temp;  
   return ret;  
}  
  
int MFCC::Mfcc_Frame_diff1_temp(MfccInfo *p, int iframe, float *out, int len)//输出一阶差分  
{  
    int ret = Mfcc_Frame_std(p, iframe + 4, p->next2, len);     
    int dctlen = p->dct.dctlen;  
    float  factor = sqrt(10.0);  
    for(int i = 0; i < dctlen; i++)  
    {  
        out[i] = (2 * p->next2[i] + p->next1[i] - 2*p->pre1[i] - p->pre2[i])/factor ;//一阶差分  
    }  
  
    float *temp = p->pre1;  
    p->pre1 = p->pre2;     
    p->pre2 = p->cur;  
    p->cur  = p->next1;  
    p->next1 = p->next2;  
    p->next2 = temp;  
    return ret;  
}  
  
int MFCC::Mfcc_Frame_diff2(MfccInfo *p, int iframe, float *out, int len)//输出MFCC+1+2  
{  
    assert(p->nframes >= 9 && iframe <= p->nframes -8 && p->m_type == MFCC_DIFF_2);  
  
    int ret = Mfcc_Frame_diff1_temp(p, iframe + 8, p->diff_next2, len);    
  
    int dctlen = p->dct.dctlen;  
    memcpy( out, p->next2, sizeof(float)* dctlen);//MFCC  
    memcpy( out + dctlen, p->diff_cur, sizeof(float)* dctlen);//一阶差分  
    float  factor = sqrt(10.0);  
    for(int i = 0; i < dctlen; i++)  
    {  
        out[i + 2*dctlen] = (2 * p->diff_next2[i] + p->diff_next1[i] - 2*p->diff_pre1[i] - p->diff_pre2[i])/factor ;//二阶差分  
    }  
  
    float *temp = p->diff_pre1;  
    p->diff_pre1 = p->diff_pre2;     
    p->diff_pre2 = p->diff_cur;  
    p->diff_cur  = p->diff_next1;  
    p->diff_next1 = p->diff_next2;  
    p->diff_next2 = temp;  
    return ret;  
}  
#endif 

void MFCC::hanning( float *win, int N)  
{     
    int half = 0;  
    if ( N % 2 == 0 )  
    {  
        half = N / 2;  
  
        for (int i = 1; i <= half; ++i)  
        {  
            win[i - 1] = 0.5 - 0.5*cos(TWOPI*i / (N + 1.0));  
        }  
  
        int index = half + 1;  
        for (int i = half; i >= 1; i--)  
        {  
            win[index - 1] = win[i - 1];  
            index++;  
        }  
  
    }  
    else  
    {  
        half = (N + 1) / 2;  
  
        for (int i = 1; i <= half; ++i)  
        {  
            win[i - 1] = 0.5 - 0.5*cos(TWOPI*i / (N + 1.0));  
        }  
  
        int index = half + 1;  
        for (int i = half-1; i >= 1; i--)  
        {  
            win[index - 1] = win[i - 1];  
            index++;  
        }  
  
    }  
}  
void MFCC::hamming( float *win, int N)  
{  
    int half = 0;  
    if ( N % 2 == 0 )  
    {  
        half = N / 2;  
  
        for (int i = 1; i <= half; ++i)  
        {  
            win[i - 1] = 0.54 - 0.46*cos(TWOPI*i / (N + 1.0));  
        }  
  
        int index = half + 1;  
        for (int i = half; i >= 1; i--)  
        {  
            win[index - 1] = win[i - 1];  
            index++;  
        }  
  
    }  
    else  
    {  
        half = (N + 1) / 2;  
  
        for (int i = 1; i <= half; ++i)  
        {  
            win[i - 1] = 0.54 - 0.46*cos(TWOPI*i / (N + 1.0));  
        }  
  
        int index = half + 1;  
        for (int i = half-1; i >= 1; i--)  
        {  
            win[index - 1] = win[i - 1];  
            index++;  
        }  
  
    }  
  
  
  
}  
void MFCC::apply_window(float* out,float* in, float* window, int window_len)  
{  
    for(int i = 0; i< window_len; i++) {  
        out[i] = in[i] * window[i];  
    }
}
  
MFCC::~MFCC()  
{  
    FreeMatrix(info.melbank.filter);  
    FreeMatrix(info.dct.coeff);  
    if(info.window)  
    {  
        delete info.window;  
        info.window = NULL;  
    }  
    if(info.lift_window)  
    {  
        delete info.lift_window;  
        info.lift_window = NULL;  
    }  
    if(info.frame_data)  
    {  
        delete info.frame_data;  
        info.frame_data = NULL;  
    }  
    if(fft)
    {
        delete fft;
        fft = NULL;  
    }
    if(fft_outBuffer) {
        delete fft_outBuffer;
        fft_outBuffer = NULL;
    }
    if(mal_Buffer) {
        delete mal_Buffer;
        mal_Buffer = NULL;
    }
    #if 0
    if(info.pre1)  
    {  
        free(info.pre1);  
        info.pre1 = NULL;  
    }  
    if(info.pre2)  
    {  
        free(info.pre2);  
        info.pre2 = NULL;  
    }  
    if(info.cur)  
    {  
        free(info.cur);  
        info.cur = NULL;  
    }  
    if(info.next1)  
    {  
        free(info.next1);  
        info.next1 = NULL;  
    }  
    if(info.next2)  
    {  
        free(info.next2);  
        info.next2 = NULL;  
    }  
  
    if(info.diff_pre1)  
    {  
        free(info.pre1);  
        info.pre1 = NULL;  
    }  
    if(info.diff_pre2)  
    {  
        free(info.pre2);  
        info.pre2 = NULL;  
    }  
  
    if(info.diff_cur)  
    {  
        free(info.cur);  
        info.cur = NULL;  
    }  
    if(info.diff_next1)  
    {  
        free(info.next1);  
        info.next1 = NULL;  
    }  
    if(info.diff_next2)  
    {  
        free(info.next2);  
        info.next2 = NULL;  
    }  
    #endif
}  