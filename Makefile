CC = gcc
LIBS = -lfftw3f -lm

SRC = dft_filter.c
OUT = dft_filter

# Default build target
.PHONY: all clean

# Build all configurations
all: frame_based_with_dft_time_domain frame_based_with_fft_time_domain dft_based fft_based no_define

# 1. FRAME_BASED & FRAME_BASED_WITH_TIME_DOMAIN
frame_based_with_dft_time_domain:
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)_frame_based_with_dft_time_domain.exe $(LIBS) -DFRAME_BASED -DFRAME_BASED_WITH_TIME_DOMAIN
	./dft_filter_frame_based_with_dft_time_domain.exe

# 2. FRAME_BASED & USE_FFT & FRAME_BASED_WITH_TIME_DOMAIN
frame_based_with_fft_time_domain:
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)_frame_based_with_fft_time_domain.exe $(LIBS) -DFRAME_BASED -DUSE_FFT -DFRAME_BASED_WITH_TIME_DOMAIN
	./dft_filter_frame_based_with_fft_time_domain.exe

# 3. FRAME_BASED
dft_based:
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)_dft_based.exe $(LIBS) -DFRAME_BASED
	./dft_filter_dft_based.exe

# 4. FRAME_BASED & USE_FFT
fft_based:
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)_fft_based.exe $(LIBS) -DFRAME_BASED -DUSE_FFT
	./dft_filter_fft_based.exe

# 5. No defines
no_define:
	$(CC) $(CFLAGS) $(SRC) -o $(OUT)_no_define.exe $(LIBS)
	./dft_filter_no_define.exe

# Clean all generated files
clean:
	rm -f $(OUT)_frame_based_with_fft_time_domain.exe \
	      $(OUT)_frame_based_with_dft_time_domain.exe \
	      $(OUT)_dft_based.exe \
	      $(OUT)_fft_based.exe \
	      $(OUT)_no_define.exe
	@echo "Cleaned all build files"