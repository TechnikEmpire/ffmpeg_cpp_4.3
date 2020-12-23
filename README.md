# ffmpeg_cpp_4.3
Really ugly conversion of ffmpeg.exe into a C++ DLL.

The sole purpose of this project was to get FFMPEG.exe and all of its components into a structure that could be allocated and freed so it could be invoked on a per-instance basis. The code is ugly because ffmpeg's code is ugly and we just made it worse by slapping it all into structs.

Released to comply with requirements of original LGPL as it was build from a LGPL 4.3.x version.
