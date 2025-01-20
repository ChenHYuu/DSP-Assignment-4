#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#define FS 48000.0f
#define FL 1000.0f
#define FH 4000.0f
#define M 128
#define PI 3.141592653589793f

#define L 480
#define P 2*M+1 
#define N 1024 // N >= L+P-1 (480+257-1=736)

#define FRAME_BASED
#define FRAME_BASED_WITH_TIME_DOMAIN
//#define USE_FFT


typedef struct _wav {
	int fs;
	char header[44];
	size_t length;
	short *LChannel;
	short *RChannel;
} wav;

int wav_read_fn(char *fn, wav *p_wav)
{
	//char header[44];
	short temp = 0;
	size_t i = 0;

	FILE *fp = fopen(fn, "rb");
	if(fp==NULL) {
		fprintf(stderr, "cannot read %s\n", fn);
		return 0;
	}
	fread(p_wav->header, sizeof(char), 44, fp);
	while( !feof(fp) ) {
		fread(&temp, sizeof(short), 1, fp);
		i++;
	}
	p_wav->length = i / 2;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		fclose(fp);
		return 0;
	}
	fseek(fp, 44, SEEK_SET);
	for(i=0;i<p_wav->length;i++) {
		fread(p_wav->LChannel+i, sizeof(short), 1, fp);
		fread(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_save_fn(char *fn, wav *p_wav)
{
	FILE *fp = fopen(fn, "wb");
	size_t i;
	if(fp==NULL) {
		fprintf(stderr, "cannot save %s\n", fn);
		return 0;
	}
	fwrite(p_wav->header, sizeof(char), 44, fp);
	for(i=0;i<p_wav->length;i++) {
		fwrite(p_wav->LChannel+i, sizeof(short), 1, fp);
		fwrite(p_wav->RChannel+i, sizeof(short), 1, fp);
	}
	fclose(fp);
	return 1;
}

int wav_init(size_t length, wav *p_wav)
{
	p_wav->length = length;
	p_wav->LChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->LChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for LChannel in wav_read_fn\n");
		return 0;
	}
	p_wav->RChannel = (short *) calloc(p_wav->length, sizeof(short));
	if( p_wav->RChannel==NULL ) {
		fprintf(stderr, "cannot allocate memory for RChannel in wav_read_fn\n");
		return 0;
	}
	return 1;
}

void wav_free(wav *p_wav)
{
	free(p_wav->LChannel);
	free(p_wav->RChannel);
}

/* hamming: for n=0,1,2,...N, length of N+1 */
float hamming(int tmpN, int n)
{
	return 0.54 - 0.46 * cosf(2*PI*((float)(n))/((float)tmpN));
}

/* low-pass filter coef: n=0,1,2...,2M */
float low_pass(int m, int n)
{
	float wc = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return wc/PI;
	}
	else {
		return sinf(wc*((float)(n-m)))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}

float band_pass(int m, int n)
{
	float wh = 2*PI*FH/FS;
    float wl = 2*PI*FL/FS;
	if(n==m) {// L'Hopital's Rule
		return 2.0*(wh/PI - wl/PI);
	}
	else {
		return 2.0*(sinf(wh*((float)(n-m)))-sinf(wl*((float)(n-m))))/PI/((float)(n-m)) * hamming(2*m+1, n);
	}
}

void conv(float x[N], float h[N], float y[N])
{
	int n, k;
	for(n=0;n<(L+P-1);n++) {
		y[n] = 0;
		for(k=0;k<P;k++) {
			if((n-k)>=0) {
				y[n] = y[n] + h[k] * x[n-k];
			}
		}
	}
}

float cosbv[N][N], sinbv[N][N];

void DFT_basis_forming(void)
{
	int n, k;
	float w;
	for(k=0;k<N;k++) {
		for(n=0;n<N;n++) {
			w = -2.0 * PI * ((float)k) * ((float)n) / ((float)N);
			printf("k=%d, n=%d, w=%f\n", k, n, w);
			cosbv[k][n] = cosf(w);
			sinbv[k][n] = sinf(w);
		}
	}
}

void DFT(float x[N], float Xre[N], float Xim[N])
{
	// X[k] = \sum_{n=0}^{N-1} x[n]\exp(-\frac{2 \pi kn}{N})
	// Xre[k] = \sum_{n=0}^{N-1} x[n]\cos(-\frac{2 \pi kn}{N})
	// Xim[k] = \sum_{n=0}^{N-1} x[n]\sin(-\frac{2 \pi kn}{N})
	int n, k;
	float w;
	for(k=0;k<N;k++) {
		Xre[k] = 0.0f; 
		Xim[k] = 0.0f; 
		for(n=0;n<N;n++) {
			Xre[k] += x[n]*cosbv[k][n];
			Xim[k] += x[n]*sinbv[k][n];
		}
	}
}

void IDFT(float x[N], float Xre[N], float Xim[N])
{	
	// x[n] = 1/N \sum_{k=0}^{N-1} X[k] \exp(\frac{2 \pi kn}{N})
	// X[k]                     = Xre[k]                   + j Xim[k];
	// \exp(\frac{2 \pi kn}{N}) = \cos(\frac{2 \pi kn}{N}) + j \sin(\frac{2 \pi kn}{N})
	
	int n, k;
	float w;
	for(n=0;n<N;n++) {
		x[n] = 0; 
		for(k=0;k<N;k++) {
			x[n] += ( Xre[k]*cosbv[k][n] + Xim[k]*sinbv[k][n] );			
		}
		x[n] = x[n]/((float)N);
	}

}


void conv_by_DFT(float x[N], float Hre[N], float Him[N], float y[N])
{
	float Xre[N], Xim[N];
	float Yre[N], Yim[N];
	int k;
	// x --> X: X=DFT(x)
	DFT(x, Xre, Xim);
	// h --> H: H=DFT(h)
	
	// Y = X H: Y = HX = (Hre + j Him) (Xre + j Xim) = Hre Xre - Him Xim + j (Hre Xim + Him Xre)
	for(k=0;k<N;k++) {
		Yre[k] = Hre[k]*Xre[k] - Him[k]*Xim[k];
		Yim[k] = Hre[k]*Xim[k] + Him[k]*Xre[k];
	}

	// Y --> y: y = IDFT(Y)
	IDFT(y, Yre, Yim);
}

int main(int argc, char **argv)
{
	wav wavin;
	wav wavout;
	char fn_in[1024] = {"blue_giant_fragment.wav"};// please find this file at https://github.com/cychiang-ntpu/ntpu-ce-mmsp-2023/blob/master/Chapter-4/blue_giant_fragment.wav
#ifndef FRAME_BASED
	char fn_out[1024] = {"out_sample_based_time_domain_filtering.wav"};
#else
    #ifdef FRAME_BASED_WITH_TIME_DOMAIN
        char fn_out[1024] = {"out_frame_based_time_domain_filtering.wav"};
    #else
        #ifndef USE_FFT
            char fn_out[1024] = {"out_frame_based_freq_domain_filtering_with_DFT.wav"};
        #else
            char fn_out[1024] = {"out_frame_based_freq_domain_filtering_with_FFT.wav"};
        #endif
    #endif
#endif
	float h_L[2*M+1] = {0};
    float h_R[2*M+1] = {0};
	int n = 0;
	float y = 0;
	int k;
	float x_r[N]; // input of each frame
	float y_r[N]; // output of each frame
	float h[N]; // inpulse response
	int r = 0;// frame index
	int R = 0;// number of frames
	float past_y[N];
	float Hre[N], Him[N];

	DFT_basis_forming();


	// read wav
	if( wav_read_fn(fn_in, &wavin) == 0 ) {
		fprintf(stderr, "cannot read wav file %s\n", fn_in);
		exit(1);
	}


	// construct low-pass filter
	for(n=0;n<(2*M+1);n++) {
		h_L[n] = band_pass(M, n);
        h_R[n] = band_pass(M, n);
	}

    /*
	for(n=0;n<(2*M+1);n++) {
		fprintf(stdout, "%.15f\n", h[n]);
	}
    */

	
	if( wav_init(wavin.length, &wavout)==0 ) {
		exit(1);
	}


// filtering (convolution)
#ifdef FRAME_BASED
	memset(past_y, 0, sizeof(float)*N);
	memset(h, 0, sizeof(float)*N);
	for(n=0;n<P;n++) {
		h[n] = h_L[n];
	}
	DFT(h, Hre, Him);

    clock_t begin = clock();
    clock_t begin_excluded;
	R = (int)ceilf(((float)(wavin.length)) / ((float)L)); // number of frames
	for(r=0;r<(R-1);r++) { // for frame 0 to frame R-2
        begin_excluded = clock();
		printf("processing frame %d/%d\n", (r+1), R);
		// left channel		
		// framing
		memset(x_r, 0, sizeof(float)*N);
		for(n=0;n<L;n++) {
			x_r[n] = (float)(wavin.LChannel[r*L+n]);
		}

#ifdef FRAME_BASED_WITH_TIME_DOMAIN
		// conv: y_r = x_r * h; in time domain
		conv(x_r, h, y_r);
#else
		// conv in freq.
    #ifdef USE_FFT
        // Please add code here
    #else
		conv_by_DFT(x_r, Hre, Him, y_r);
    #endif
#endif

		// overlap and add		
		for(n=0;n<(L+P-1);n++) {
			wavout.LChannel[r*L+n] += y_r[n];// precision problem
		}


		// right channel		
		// framing
		memset(x_r, 0, sizeof(float)*N);
		for(n=0;n<L;n++) {
			x_r[n] = (float)(wavin.RChannel[r*L+n]);
		}	

#ifdef FRAME_BASED_WITH_TIME_DOMAIN
		// conv: y_r = x_r * h; in time domain
		conv(x_r, h, y_r);
#else
        // conv in freq.
    #ifdef USE_FFT
        // Please add code here
    #else
		conv_by_DFT(x_r, Hre, Him, y_r);
    #endif
#endif

		// overlap and add		
		for(n=0;n<(L+P-1);n++) {
			wavout.RChannel[r*L+n] += y_r[n];// precision problem
		}		
	}
	// for frame R-1
	// left channel
	memset(x_r, 0, sizeof(float)*N);
	for(n=0;n<L;n++) {
		if( (r*L+n)< wavin.length ) {
			x_r[n] = (float)(wavin.LChannel[r*L+n]);
		}
	}

#ifdef FRAME_BASED_WITH_TIME_DOMAIN
    // conv: y_r = x_r * h; in time domain
    conv(x_r, h, y_r);
#else
    // conv in freq.
    #ifdef USE_FFT
        // Please add code here
    #else
		conv_by_DFT(x_r, Hre, Him, y_r);
    #endif
#endif

	// overlap and add		
	for(n=0;(r*L+n)<wavin.length;n++) {
		wavout.LChannel[r*L+n] += y_r[n];// precision problem
	}
	// right channel
	memset(x_r, 0, sizeof(float)*N);
	for(n=0;n<L;n++) {
		if( (r*L+n)< wavin.length ) {
			x_r[n] = (float)(wavin.RChannel[r*L+n]);
		}
	}

#ifdef FRAME_BASED_WITH_TIME_DOMAIN
    // conv: y_r = x_r * h; in time domain
    conv(x_r, h, y_r);
#else
    // conv in freq.
    #ifdef USE_FFT
        // Please add code here
    #else
		conv_by_DFT(x_r, Hre, Him, y_r);
    #endif
#endif

	// overlap and add		
	for(n=0;(r*L+n)<wavin.length;n++) {
		wavout.RChannel[r*L+n] += y_r[n];// precision problem
	}

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("It took %.15e msec to generate each sample with frame-based filtering\n",
        time_spent/((double)(wavin.length))/2.0);

#else /* sample-based processing */	
    clock_t begin = clock();
	for(n=0;n<wavin.length;n++) {
		y = 0;
		for(k=0;k<(2*M+1);k++) {
			if( (n-k)>=0 )
				y = y + h_L[k] * ((float)(wavin.LChannel[n-k]));
		}
		wavout.LChannel[n] = (short)(roundf(y));

		y = 0;
		for(k=0;k<(2*M+1);k++) {
			if( (n-k)>=0 )
				y = y + h_R[k] * ((float)(wavin.RChannel[n-k]));
		}
		wavout.RChannel[n] = (short)(roundf(y));
	}
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("It took %.15e msec to generate each sample with time domain linear filtering\n",
        time_spent/((double)(wavin.length))/2.0);
#endif

	memcpy(wavout.header, wavin.header, 44);
	// save wav
	if( wav_save_fn(fn_out, &wavout)==0) {
		fprintf(stderr, "cannot save %s\n", fn_out);
		exit(1);

	}
	wav_free(&wavin);
	wav_free(&wavout);
}