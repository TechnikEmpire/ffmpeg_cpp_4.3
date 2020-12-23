#pragma once

#define COMPILER_VC_WINDOWS				1

#define CONFIG_AVUTIL					1
#define CONFIG_AVFORMAT					1
#define CONFIG_AVCODEC					1
#define CONFIG_AVFILTER					1
#define CONFIG_SWSCALE					1
#define CONFIG_SWRESAMPLE				1
#define CONFIG_AVDEVICE					1
#define CONFIG_POSTPROC					1
#define CONFIG_AVRESAMPLE				0
#ifndef CONFIG_NETWORK
#define CONFIG_NETWORK					0
#endif
#define HAVE_ALTIVEC_H					0

#define FFMPEG_CONFIGURATION			""

#define CONFIG_NONFREE					0
#define CONFIG_GPL						0
#define CONFIG_GPLV3					0
#define CONFIG_LGPLV3					1

//#define AV_PIX_FMT_ABI_GIT_MASTER		0
#define HAVE_GETRUSAGE					0
#define HAVE_SYS_RESOURCE_H				0
#define HAVE_GETPROCESSTIMES			1
#define HAVE_TERMIOS_H					0
#define HAVE_KBHIT						1
#define HAVE_PEEKNAMEDPIPE				0
#define HAVE_COMMANDLINETOARGVW			1

#define HAVE_IO							1
#define HAVE_THREADS					1

#define CONFIG_RTSP_DEMUXER				1
#define CONFIG_MMSH_PROTOCOL			0
#define CONFIG_OPENCL					0

#if (COMPILER_VC_WINDOWS)

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stddef.h>
#include <float.h>
#include <windows.h>
#include <math.h>
#include <inttypes.h>

#ifndef usleep
#define 	usleep			Sleep
#endif

#define		snprintf		_snprintf
#define		strcasecmp		_stricmp
#define		strncasecmp		strncmp 
#define		pid_t			DWORD
#define		CC_TYPE			"vc"
#define		CC_VERSION		"vc10"
#define     CC_IDENT		"ms-cl"

#define HAVE_GETPROCESSMEMORYINFO		0		//pro info
#if HAVE_GETPROCESSMEMORYINFO
#pragma comment(lib,"Psapi.lib")
#endif

//math function
#define		HAVE_ATANF		1
#define		HAVE_POWF		1
#define		HAVE_ATAN2F		1
#define		HAVE_COSF		1
#define		HAVE_EXPF		1
#define		HAVE_LDEXPF		1
#define		HAVE_SINF		1

	int isinf_i386(double x);

#include <time.h>
#if !WINVER
#include <sys/time.h>
#endif

	int gettimeofday(struct timeval* tp, void* tzp);

#ifdef __cplusplus
}
#endif

#elif (COMPILER_GCC_LINUX)
#include <inttypes.h>
#include <stdint.h>
#define CC_TYPE "gcc"
#define CC_VERSION __VERSION__
#define restrict restrict
#elif (COMPILER_XCODE_IOS)
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stddef.h>
#elif (COMPILER_GCC_ANDROID)
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stddef.h>
#define CC_TYPE "gcc"
#define CC_VERSION __VERSION__
#define restrict __restrict
#endif

#define ff_dlog(a, ...) while(0)

#pragma comment(lib, "avcodec.lib")
#pragma comment(lib, "avdevice.lib")
#pragma comment(lib, "avfilter.lib")
#pragma comment(lib, "avformat.lib")
#pragma comment(lib, "avutil.lib")
#pragma comment(lib, "swresample.lib")
#pragma comment(lib, "swscale.lib")
#pragma comment(lib, "Shell32.lib")
