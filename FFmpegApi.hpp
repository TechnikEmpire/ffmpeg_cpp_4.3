#pragma once

#pragma once

#include "Native/Config.hpp"
#include <cstdint>

#ifdef FFMPEG_EXPORT
#ifdef _MSC_VER
#define FFMPEG_API __declspec(dllexport)
#else
#define FFMPEG_API __attribute__((visibility("default")))
#endif // #ifdef _MSC_VER
#else
#ifdef _MSC_VER
#define FFMPEG_API __declspec(dllimport)
#else
#define FFMPEG_API
#endif
#endif // #ifdef FFMPEG_EXPORT

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

	extern FFMPEG_API int run_ffmpeg(int argc, char** argv, int (*read_packet)(void* opaque, uint8_t* buf, int buf_size), int (*write_packet)(void* opaque, uint8_t* buf, int buf_size));

#ifdef __cplusplus
};
#endif // __cplusplus