#include "FFmpegApi.hpp"

#include "Native/ffmpeg.hpp"

#include <memory>
#include <iostream>

int run_ffmpeg(int argc, char** argv, int (*read_packet)(void* opaque, uint8_t* buf, int buf_size), int (*write_packet)(void* opaque, uint8_t* buf, int buf_size))
{
	int retVal = -9001;
	try
	{
		auto shared = std::make_shared<FFmpegNative::FFmpeg>();
		retVal = shared->FFmpegMain(argc, argv, read_packet, write_packet);
	}
	catch(std::exception& e){
		std::cout << e.what() << std::endl;
	}

	return retVal;
}
