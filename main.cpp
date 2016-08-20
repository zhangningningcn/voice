#include "wav.hpp"
#include  "mfcc.hpp" 

int main(void) {
    WAV wav;
    
    int framlen = wav.Open("12.wav");
    
    if(framlen < 0) {
        return framlen;
    }
    printf("fram len %d",framlen);
    float *fram;
    fram = new float[framlen];
    
    
    AudioInfo audioctx;  
    audioctx.data_in = fram;  
    //audioctx.data_len = samples;  
    audioctx.frame_len = framlen;  
    audioctx.frame_shift = framlen/2;  
    audioctx.sample_rate = wav.getSampleRate();  
    
    int frame_shift = audioctx.frame_shift;
    MFCC mfcc(audioctx, 0, audioctx.sample_rate/2, 24, 12);  
    

    int len = 13;  
    float *out = new float[len];
    FILE *fp;
    fp = fopen("mfcc.txt","w");
    if(fp == NULL) {
        return -2;
    }
    //int out_nframes = (audioctx.data_len - audioctx.frame_len)/audioctx.frame_shift + 1; 
    clock_t start,finish;
    start=clock();
    while(true)  
    {  
        if(wav.ReadBlock(fram,0) <= 0) {
            break;
        }
        //frame_shift = 0;
        for(int i=0;i<framlen;i+=frame_shift){
            mfcc.Mfcc_Frame_std(fram+i, out); 
            for(int i=0;i<13;i++) {
                fprintf(fp,"%f,",out[i]);
            }
            fprintf(fp,"\n");
        }
    }  
    
    finish=clock();
    cout << (double)(finish-start)/CLOCKS_PER_SEC << endl;
    fclose(fp);
  
}