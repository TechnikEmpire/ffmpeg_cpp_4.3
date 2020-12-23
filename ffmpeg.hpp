/*
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#pragma once

#include "cmdutils.hpp"
#include "Config.hpp"
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdio.h>
#include <signal.h>

#include <libavformat/avformat.h>
#include <libavformat/avio.h>

#include <libavcodec/avcodec.h>

#include <libavfilter/avfilter.h>

#include <libavutil/avutil.h>
#include <libavutil/dict.h>
#include <libavutil/eval.h>
#include <libavutil/fifo.h>
#include <libavutil/hwcontext.h>
#include <libavutil/pixfmt.h>
#include <libavutil/rational.h>
#include <libavutil/threadmessage.h>

#include <libswresample/swresample.h>

#ifdef __cplusplus
}
#endif

#include <vector>

#define VSYNC_AUTO       -1
#define VSYNC_PASSTHROUGH 0
#define VSYNC_CFR         1
#define VSYNC_VFR         2
#define VSYNC_VSCFR       0xfe
#define VSYNC_DROP        0xff

#define MAX_STREAMS 1024    /* arbitrary sanity check value */

namespace FFmpegNative
{
    enum HWAccelID {
        HWACCEL_NONE = 0,
        HWACCEL_AUTO,
        HWACCEL_GENERIC,
        HWACCEL_VIDEOTOOLBOX,
        HWACCEL_QSV,
    };

    typedef struct HWAccel {
        const char* name;
        int (*init)(AVCodecContext* s);
        enum HWAccelID id;
        enum AVPixelFormat pix_fmt;
    } HWAccel;

    typedef struct HWDevice {
        const char* name;
        enum AVHWDeviceType type;
        AVBufferRef* device_ref;
    } HWDevice;

    /* select an input stream for an output stream */
    typedef struct StreamMap {
        int disabled;           /* 1 is this mapping is disabled by a negative map */
        int file_index;
        int stream_index;
        int sync_file_index;
        int sync_stream_index;
        char* linklabel;       /* name of an output link, for mapping lavfi outputs */
    } StreamMap;

    typedef struct {
        int  file_idx, stream_idx, channel_idx; // input
        int ofile_idx, ostream_idx;               // output
    } AudioChannelMap;

    typedef struct OptionsContext {
        OptionGroup* g;

        /* input/output options */
        int64_t start_time;
        int64_t start_time_eof;
        int seek_timestamp;
        const char* format;

        SpecifierOpt* codec_names;
        int        nb_codec_names;
        SpecifierOpt* audio_channels;
        int        nb_audio_channels;
        SpecifierOpt* audio_sample_rate;
        int        nb_audio_sample_rate;
        SpecifierOpt* frame_rates;
        int        nb_frame_rates;
        SpecifierOpt* frame_sizes;
        int        nb_frame_sizes;
        SpecifierOpt* frame_pix_fmts;
        int        nb_frame_pix_fmts;

        /* input options */
        int64_t input_ts_offset;
        int loop;
        int rate_emu;
        int accurate_seek;
        int thread_queue_size;

        SpecifierOpt* ts_scale;
        int        nb_ts_scale;
        SpecifierOpt* dump_attachment;
        int        nb_dump_attachment;
        SpecifierOpt* hwaccels;
        int        nb_hwaccels;
        SpecifierOpt* hwaccel_devices;
        int        nb_hwaccel_devices;
        SpecifierOpt* hwaccel_output_formats;
        int        nb_hwaccel_output_formats;
        SpecifierOpt* autorotate;
        int        nb_autorotate;

        /* output options */
        StreamMap* stream_maps;
        int     nb_stream_maps;
        AudioChannelMap* audio_channel_maps; /* one info entry per -map_channel */
        int           nb_audio_channel_maps; /* number of (valid) -map_channel settings */
        int metadata_global_manual;
        int metadata_streams_manual;
        int metadata_chapters_manual;
        const char** attachments;
        int       nb_attachments;

        int chapters_input_file;

        int64_t recording_time;
        int64_t stop_time;
        uint64_t limit_filesize;
        float mux_preload;
        float mux_max_delay;
        int shortest;
        int bitexact;

        int video_disable;
        int audio_disable;
        int subtitle_disable;
        int data_disable;

        /* indexed by output file stream index */
        int* streamid_map;
        int nb_streamid_map;

        SpecifierOpt* metadata;
        int        nb_metadata;
        SpecifierOpt* max_frames;
        int        nb_max_frames;
        SpecifierOpt* bitstream_filters;
        int        nb_bitstream_filters;
        SpecifierOpt* codec_tags;
        int        nb_codec_tags;
        SpecifierOpt* sample_fmts;
        int        nb_sample_fmts;
        SpecifierOpt* qscale;
        int        nb_qscale;
        SpecifierOpt* forced_key_frames;
        int        nb_forced_key_frames;
        SpecifierOpt* force_fps;
        int        nb_force_fps;
        SpecifierOpt* frame_aspect_ratios;
        int        nb_frame_aspect_ratios;
        SpecifierOpt* rc_overrides;
        int        nb_rc_overrides;
        SpecifierOpt* intra_matrices;
        int        nb_intra_matrices;
        SpecifierOpt* inter_matrices;
        int        nb_inter_matrices;
        SpecifierOpt* chroma_intra_matrices;
        int        nb_chroma_intra_matrices;
        SpecifierOpt* top_field_first;
        int        nb_top_field_first;
        SpecifierOpt* metadata_map;
        int        nb_metadata_map;
        SpecifierOpt* presets;
        int        nb_presets;
        SpecifierOpt* copy_initial_nonkeyframes;
        int        nb_copy_initial_nonkeyframes;
        SpecifierOpt* copy_prior_start;
        int        nb_copy_prior_start;
        SpecifierOpt* filters;
        int        nb_filters;
        SpecifierOpt* filter_scripts;
        int        nb_filter_scripts;
        SpecifierOpt* reinit_filters;
        int        nb_reinit_filters;
        SpecifierOpt* fix_sub_duration;
        int        nb_fix_sub_duration;
        SpecifierOpt* canvas_sizes;
        int        nb_canvas_sizes;
        SpecifierOpt* pass;
        int        nb_pass;
        SpecifierOpt* passlogfiles;
        int        nb_passlogfiles;
        SpecifierOpt* max_muxing_queue_size;
        int        nb_max_muxing_queue_size;
        SpecifierOpt* muxing_queue_data_threshold;
        int        nb_muxing_queue_data_threshold;
        SpecifierOpt* guess_layout_max;
        int        nb_guess_layout_max;
        SpecifierOpt* apad;
        int        nb_apad;
        SpecifierOpt* discard;
        int        nb_discard;
        SpecifierOpt* disposition;
        int        nb_disposition;
        SpecifierOpt* program;
        int        nb_program;
        SpecifierOpt* time_bases;
        int        nb_time_bases;
        SpecifierOpt* enc_time_bases;
        int        nb_enc_time_bases;
        SpecifierOpt* autoscale;
        int        nb_autoscale;
    } OptionsContext;

    typedef struct InputFilter {
        AVFilterContext* filter;
        struct InputStream* ist;
        struct FilterGraph* graph;
        uint8_t* name;
        enum AVMediaType    type;   // AVMEDIA_TYPE_SUBTITLE for sub2video

        AVFifoBuffer* frame_queue;

        // parameters configured for this input
        int format;

        int width, height;
        AVRational sample_aspect_ratio;

        int sample_rate;
        int channels;
        uint64_t channel_layout;

        AVBufferRef* hw_frames_ctx;

        int eof;
    } InputFilter;

    typedef struct OutputFilter {
        AVFilterContext* filter;
        struct OutputStream* ost;
        struct FilterGraph* graph;
        uint8_t* name;

        /* temporary storage until stream maps are processed */
        AVFilterInOut* out_tmp;
        enum AVMediaType     type;

        /* desired output stream properties */
        int width, height;
        AVRational frame_rate;
        int format;
        int sample_rate;
        uint64_t channel_layout;

        // those are only set if no format is specified and the encoder gives us multiple options
        int* formats;
        uint64_t* channel_layouts;
        int* sample_rates;
    } OutputFilter;

    typedef struct FilterGraph {
        int            index;
        const char* graph_desc;

        AVFilterGraph* graph;
        int reconfiguration;

        InputFilter** inputs;
        int          nb_inputs;
        OutputFilter** outputs;
        int         nb_outputs;
    } FilterGraph;

    typedef struct InputStream {
        int file_index;
        AVStream* st;
        int discard;             /* true if stream data should be discarded */
        int user_set_discard;
        int decoding_needed;     /* non zero if the packets must be decoded in 'raw_fifo', see DECODING_FOR_* */
#define DECODING_FOR_OST    1
#define DECODING_FOR_FILTER 2

        AVCodecContext* dec_ctx;
        AVCodec* dec;
        AVFrame* decoded_frame;
        AVFrame* filter_frame; /* a ref of decoded_frame, to be sent to filters */

        int64_t       start;     /* time when read started */
        /* predicted dts of the next packet read for this stream or (when there are
         * several frames in a packet) of the next frame in current packet (in AV_TIME_BASE units) */
        int64_t       next_dts;
        int64_t       dts;       ///< dts of the last packet read for this stream (in AV_TIME_BASE units)

        int64_t       next_pts;  ///< synthetic pts for the next decode frame (in AV_TIME_BASE units)
        int64_t       pts;       ///< current pts of the decoded frame  (in AV_TIME_BASE units)
        int           wrap_correction_done;

        int64_t filter_in_rescale_delta_last;

        int64_t min_pts; /* pts with the smallest value in a current stream */
        int64_t max_pts; /* pts with the higher value in a current stream */

        // when forcing constant input framerate through -r,
        // this contains the pts that will be given to the next decoded frame
        int64_t cfr_next_pts;

        int64_t nb_samples; /* number of samples in the last decoded audio frame before looping */

        double ts_scale;
        int saw_first_ts;
        AVDictionary* decoder_opts;
        AVRational framerate;               /* framerate forced with -r */
        int top_field_first;
        int guess_layout_max;

        int autorotate;

        int fix_sub_duration;
        struct { /* previous decoded subtitle and related variables */
            int got_output;
            int ret;
            AVSubtitle subtitle;
        } prev_sub;

        struct sub2video {
            int64_t last_pts;
            int64_t end_pts;
            AVFifoBuffer* sub_queue;    ///< queue of AVSubtitle* before filter init
            AVFrame* frame;
            int w, h;
            unsigned int initialize; ///< marks if sub2video_update should force an initialization
        } sub2video;

        int dr1;

        /* decoded data from this stream goes into all those filters
         * currently video and audio only */
        InputFilter** filters;
        int        nb_filters;

        int reinit_filters;

        /* hwaccel options */
        enum HWAccelID hwaccel_id;
        enum AVHWDeviceType hwaccel_device_type;
        char* hwaccel_device;
        enum AVPixelFormat hwaccel_output_format;

        /* hwaccel context */
        void* hwaccel_ctx;
        void (*hwaccel_uninit)(AVCodecContext* s);
        int  (*hwaccel_get_buffer)(AVCodecContext* s, AVFrame* frame, int flags);
        int  (*hwaccel_retrieve_data)(AVCodecContext* s, AVFrame* frame);
        enum AVPixelFormat hwaccel_pix_fmt;
        enum AVPixelFormat hwaccel_retrieved_pix_fmt;
        AVBufferRef* hw_frames_ctx;

        /* stats */
        // combined size of all the packets read
        uint64_t data_size;
        /* number of packets successfully read for this stream */
        uint64_t nb_packets;
        // number of frames/samples retrieved from the decoder
        uint64_t frames_decoded;
        uint64_t samples_decoded;

        int64_t* dts_buffer;
        int nb_dts_buffer;

        int got_output;
    } InputStream;

    typedef struct InputFile {
        AVFormatContext* ctx;
        int eof_reached;      /* true if eof reached */
        int eagain;           /* true if last read attempt returned EAGAIN */
        int ist_index;        /* index of first stream in input_streams */
        int loop;             /* set number of times input stream should be looped */
        int64_t duration;     /* actual duration of the longest stream in a file
                                 at the moment when looping happens */
        AVRational time_base; /* time base of the duration */
        int64_t input_ts_offset;

        int64_t ts_offset;
        int64_t last_ts;
        int64_t start_time;   /* user-specified start time in AV_TIME_BASE or AV_NOPTS_VALUE */
        int seek_timestamp;
        int64_t recording_time;
        int nb_streams;       /* number of stream that ffmpeg is aware of; may be different
                                 from ctx.nb_streams if new streams appear during av_read_frame() */
        int nb_streams_warn;  /* number of streams that the user was warned of */
        int rate_emu;
        int accurate_seek;

#if HAVE_THREADS
        AVThreadMessageQueue* in_thread_queue;
        HANDLE thread;           /* thread reading from this file */
        int non_blocking;           /* reading packets from the thread should not block */
        int joined;                 /* the thread has been joined */
        int thread_queue_size;      /* maximum number of queued packets */
#endif
    } InputFile;

    enum forced_keyframes_const {
        FKF_N,
        FKF_N_FORCED,
        FKF_PREV_FORCED_N,
        FKF_PREV_FORCED_T,
        FKF_T,
        FKF_NB
    };

#define ABORT_ON_FLAG_EMPTY_OUTPUT        (1 <<  0)
#define ABORT_ON_FLAG_EMPTY_OUTPUT_STREAM (1 <<  1)

    enum OSTFinished {
        ENCODER_FINISHED = 1,
        MUXER_FINISHED = 2,
    } ;

    typedef struct OutputStream {
        int file_index;          /* file index */
        int index;               /* stream index in the output file */
        int source_index;        /* InputStream index */
        AVStream* st;            /* stream in the output file */
        int encoding_needed;     /* true if encoding needed for this stream */
        int frame_number;
        /* input pts and corresponding output pts
           for A/V sync */
        struct InputStream* sync_ist; /* input stream to sync against */
        int64_t sync_opts;       /* output frame counter, could be changed to some true timestamp */ // FIXME look at frame_number
        /* pts of the first frame encoded for this stream, used for limiting
         * recording time */
        int64_t first_pts;
        /* dts of the last packet sent to the muxer */
        int64_t last_mux_dts;
        // the timebase of the packets sent to the muxer
        AVRational mux_timebase;
        AVRational enc_timebase;

        AVBSFContext* bsf_ctx;

        AVCodecContext* enc_ctx;
        AVCodecParameters* ref_par; /* associated input codec parameters with encoders options applied */
        AVCodec* enc;
        int64_t max_frames;
        AVFrame* filtered_frame;
        AVFrame* last_frame;
        int last_dropped;
        int last_nb0_frames[3];

        void* hwaccel_ctx;

        /* video only */
        AVRational frame_rate;
        int is_cfr;
        int force_fps;
        int top_field_first;
        int rotate_overridden;
        int autoscale;
        double rotate_override_value;

        AVRational frame_aspect_ratio;

        /* forced key frames */
        int64_t forced_kf_ref_pts;
        int64_t* forced_kf_pts;
        int forced_kf_count;
        int forced_kf_index;
        char* forced_keyframes;
        AVExpr* forced_keyframes_pexpr;
        double forced_keyframes_expr_const_values[FKF_NB];

        /* audio only */
        int* audio_channels_map;             /* list of the channels id to pick from the source stream */
        int audio_channels_mapped;           /* number of channels in audio_channels_map */

        char* logfile_prefix;
        FILE* logfile;

        OutputFilter* filter;
        char* avfilter;
        char* filters;         ///< filtergraph associated to the -filter option
        char* filters_script;  ///< filtergraph script associated to the -filter_script option

        AVDictionary* encoder_opts;
        AVDictionary* sws_dict;
        AVDictionary* swr_opts;
        AVDictionary* resample_opts;
        char* apad;
        OSTFinished finished;        /* no more packets should be written for this stream */
        int unavailable;                     /* true if the steram is unavailable (possibly temporarily) */
        int stream_copy;

        // init_output_stream() has been called for this stream
        // The encoder and the bitstream filters have been initialized and the stream
        // parameters are set in the AVStream.
        int initialized;

        int inputs_done;

        const char* attachment_filename;
        int copy_initial_nonkeyframes;
        int copy_prior_start;
        char* disposition;

        int keep_pix_fmt;

        /* stats */
        // combined size of all the packets written
        uint64_t data_size;
        // number of packets send to the muxer
        uint64_t packets_written;
        // number of frames/samples sent to the encoder
        uint64_t frames_encoded;
        uint64_t samples_encoded;

        /* packet quality factor */
        int quality;

        int max_muxing_queue_size;

        /* the packets are buffered here until the muxer is ready to be initialized */
        AVFifoBuffer* muxing_queue;

        /*
         * The size of the AVPackets' buffers in queue.
         * Updated when a packet is either pushed or pulled from the queue.
         */
        size_t muxing_queue_data_size;

        /* Threshold after which max_muxing_queue_size will be in effect */
        size_t muxing_queue_data_threshold;

        /* packet picture type */
        int pict_type;

        /* frame encode sum of squared error values */
        int64_t error[4];
    } OutputStream;

    typedef struct OutputFile {
        AVFormatContext* ctx;
        AVDictionary* opts;
        int ost_index;       /* index of the first stream in output_streams */
        int64_t recording_time;  ///< desired length of the resulting file in microseconds == AV_TIME_BASE units
        int64_t start_time;      ///< start time in microseconds == AV_TIME_BASE units
        uint64_t limit_filesize; /* filesize limit expressed in bytes */

        int shortest;

        int header_written;
    } OutputFile;

    typedef struct BenchmarkTimeStamps {
        int64_t real_usec;
        int64_t user_usec;
        int64_t sys_usec;
    } BenchmarkTimeStamps;

    enum OptGroup {
        GROUP_OUTFILE,
        GROUP_INFILE,
    };

    struct FFmpeg
    {   
        int run_as_daemon = 0;
        int nb_frames_dup = 0;
        unsigned dup_warning = 1000;
        int nb_frames_drop = 0;
        int64_t decode_error_stat[2];

        int want_sdp = 1;

        BenchmarkTimeStamps current_time;
        AVIOContext* progress_avio = NULL;

        uint8_t* subtitle_out;

        InputStream** input_streams = NULL;
        int        nb_input_streams = 0;
        InputFile** input_files = NULL;
        int        nb_input_files = 0;

        OutputStream** output_streams = NULL;
        int         nb_output_streams = 0;
        OutputFile** output_files = NULL;
        int         nb_output_files = 0;

        FilterGraph** filtergraphs;
        int        nb_filtergraphs;

#if HAVE_TERMIOS_H
        /* init terminal so that we can grab keys */
        struct termios oldtty;
        int restore_tty;
#endif

        const char* const forced_keyframes_const_names[6] = {
            "n",
            "n_forced",
            "prev_forced_n",
            "prev_forced_t",
            "t",
            NULL
        };

        FILE* vstats_file;

        char* vstats_filename;
        char* sdp_filename;

        char* videotoolbox_pixfmt;

        static const std::vector<HWAccel> hwaccels;
#if CONFIG_QSV
        char* qsv_device;
#endif
        HWDevice* filter_hw_device;

        volatile int received_sigterm = 0;
        volatile int received_nb_signals = 0;
        volatile LONG transcode_init_done = 0;
        volatile int ffmpeg_exited = 0;
        int main_return_code = 0;
        int64_t copy_ts_first_pts = AV_NOPTS_VALUE;

        const AVIOInterruptCB int_cb = { decode_interrupt_cb, this };

        int64_t last_time = -1;
        int qp_histogram[52];
        
        float audio_drift_threshold = 0.1;
        float dts_delta_threshold = 10;
        float dts_error_threshold = 3600 * 30;

        int audio_volume = 256;
        int audio_sync_method = 0;
        int video_sync_method = VSYNC_AUTO;
        float frame_drop_threshold = 0;
        int do_deinterlace = 0;
        int do_benchmark = 0;
        int do_benchmark_all = 0;
        int do_hex_dump = 0;
        int do_pkt_dump = 0;
        int copy_ts = 0;
        int start_at_zero = 0;
        int copy_tb = -1;
        int debug_ts = 0;
        int exit_on_error = 0;
        int abort_on_flags = 0;
        int print_stats = -1;
        int qp_hist = 0;
        int stdin_interaction = 1;
        int frame_bits_per_raw_sample = 0;
        float max_error_rate = 2.0 / 3;
        int filter_nbthreads = 0;
        int filter_complex_nbthreads = 0;
        int vstats_version = 2;
        int auto_conversion_filters = 1;


        int intra_only = 0;
        int file_overwrite = 0;
        int no_file_overwrite = 0;
        int do_psnr = 0;
        int input_sync;
        int input_stream_potentially_available = 0;
        int ignore_unknown_streams = 0;
        int copy_unknown_streams = 0;
        int find_stream_info = 1;

        int nb_hw_devices;
        HWDevice** hw_devices;

        static const OptionGroupDef groups[];

        //int (*read_packet)(void* opaque, uint8_t* buf, int buf_size), int (*write_packet)(void* opaque, uint8_t* buf, int buf_size)

        //std::function<int(void* opaque, uint8_t* buf, int buf_size)> _readPacket;
        //std::function<int(void* opaque, uint8_t* buf, int buf_size)> _writePacket;

        typedef int(*ffmpeg_custom_io_func)(void* opaque, uint8_t* buf, int buf_size);

        ffmpeg_custom_io_func _readPacket;
        ffmpeg_custom_io_func _writePacket;

        static constexpr size_t ReadBufferSize = 4096;
        static constexpr size_t WriteBufferSize = 4096;

        uint8_t* _readBuffer;

        uint8_t* _writeBuffer;

        CmdUtils _cmdUtils;


        void term_init(void);
        void term_exit(void);

        static void show_usage(void);

        void remove_avoptions(AVDictionary** a, AVDictionary* b);
        void assert_avoptions(AVDictionary* m);

        int guess_input_channel_layout(InputStream* ist);

        enum AVPixelFormat choose_pixel_fmt(AVStream* st, AVCodecContext* avctx,
            const AVCodec* codec, enum AVPixelFormat target);
        void choose_sample_fmt(AVStream* st, const AVCodec* codec);

        int configure_filtergraph(FilterGraph* fg);
        int configure_output_filter(FilterGraph* fg, OutputFilter* ofilter, AVFilterInOut* out);
        void check_filter_outputs(void);
        int ist_in_filtergraph(FilterGraph* fg, InputStream* ist);
        int filtergraph_is_simple(FilterGraph* fg);
        int init_simple_filtergraph(InputStream* ist, OutputStream* ost);
        int init_complex_filtergraph(FilterGraph* fg);

        void sub2video_update(InputStream* ist, int64_t heartbeat_pts, AVSubtitle* sub);

        int ifilter_parameters_from_frame(InputFilter* ifilter, const AVFrame* frame);

        int ffmpeg_parse_options(int argc, char** argv);

        HWDevice* hw_device_get_by_name(const char* name);
        int hw_device_init_from_string(const char* arg, HWDevice** dev);
        void hw_device_free_all(void);

        int hw_device_setup_for_decode(InputStream* ist);
        int hw_device_setup_for_encode(OutputStream* ost);
        int hw_device_setup_for_filter(FilterGraph* fg);

        static int hwaccel_decode_init(AVCodecContext* avctx);

        int sub2video_get_blank_frame(InputStream* ist);

        void sub2video_copy_rect(uint8_t* dst, int dst_linesize, int w, int h, AVSubtitleRect* r);

        void sub2video_push_ref(InputStream* ist, int64_t pts);

        void sub2video_heartbeat(InputStream* ist, int64_t pts);

        void sub2video_flush(InputStream* ist);

        void term_exit_sigsafe(void);

        void sigterm_handler(int sig);

        int read_key(void);

        static int decode_interrupt_cb(void* ctx);

        void ffmpeg_cleanup(int ret);

        void abort_codec_experimental(AVCodec* c, int encoder);

        void update_benchmark(const char* fmt, ...);

        void close_all_output_streams(OutputStream* ost, OSTFinished this_stream, OSTFinished others);

        void write_packet(OutputFile* of, AVPacket* pkt, OutputStream* ost, int unqueue);

        void close_output_stream(OutputStream* ost);

        void output_packet(OutputFile* of, AVPacket* pkt, OutputStream* ost, int eof);

        int check_recording_time(OutputStream* ost);

        double adjust_frame_pts_to_encoder_tb(OutputFile* of, OutputStream* ost, AVFrame* frame);

        int init_output_stream(OutputStream* ost, AVFrame* frame, char* error, int error_len);

        int init_output_stream_wrapper(OutputStream* ost, AVFrame* frame, unsigned int fatal);

        void do_audio_out(OutputFile* of, OutputStream* ost, AVFrame* frame);

        void do_subtitle_out(OutputFile* of, OutputStream* ost, AVSubtitle* sub);

        void do_video_out(OutputFile* of, OutputStream* ost, AVFrame* next_picture);

        static double psnr(double d);

        void do_video_stats(OutputStream* ost, int frame_size);

        void finish_output_stream(OutputStream* ost);

        int reap_filters(int flush);

        void print_final_stats(int64_t total_size);

        void print_report(int is_last_report, int64_t timer_start, int64_t cur_time);

        static void ifilter_parameters_from_codecpar(InputFilter* ifilter, AVCodecParameters* par);

        void flush_encoders(void);

        int check_output_constraints(InputStream* ist, OutputStream* ost);

        void do_streamcopy(InputStream* ist, OutputStream* ost, const AVPacket* pkt);

        void check_decode_result(InputStream* ist, int* got_output, int ret);

        static int ifilter_has_all_input_formats(FilterGraph* fg);

        int ifilter_send_frame(InputFilter* ifilter, AVFrame* frame);

        int ifilter_send_eof(InputFilter* ifilter, int64_t pts);

        int decode(AVCodecContext* avctx, AVFrame* frame, int* got_frame, AVPacket* pkt);

        int send_frame_to_filters(InputStream* ist, AVFrame* decoded_frame);

        int decode_audio(InputStream* ist, AVPacket* pkt, int* got_output, int* decode_failed);

        int decode_video(InputStream* ist, AVPacket* pkt, int* got_output, int64_t* duration_pts, int eof, int* decode_failed);

        int transcode_subtitles(InputStream* ist, AVPacket* pkt, int* got_output, int* decode_failed);

        int send_filter_eof(InputStream* ist);

        int process_input_packet(InputStream* ist, const AVPacket* pkt, int no_eof);

        void print_sdp(void);

        static AVPixelFormat get_format(AVCodecContext* s, const enum AVPixelFormat* pix_fmts);

        static int get_buffer(AVCodecContext* s, AVFrame* frame, int flags);

        int init_input_stream(int ist_index, char* error, int error_len);

        InputStream* get_input_stream(OutputStream* ost);

        static int compare_int64(const void* a, const void* b);

        int check_init_output_file(OutputFile* of, int file_index);

        int init_output_bsfs(OutputStream* ost);

        int init_output_stream_streamcopy(OutputStream* ost);

        void set_encoder_id(OutputFile* of, OutputStream* ost);

        void parse_forced_key_frames(char* kf, OutputStream* ost, AVCodecContext* avctx);

        void init_encoder_time_base(OutputStream* ost, AVRational default_time_base);

        int init_output_stream_encode(OutputStream* ost, AVFrame* frame);

        void report_new_stream(int input_index, AVPacket* pkt);

        int transcode_init(void);

        int need_output(void);

        OutputStream* choose_output(void);

        void set_tty_echo(int on);

        int check_keyboard_interaction(int64_t cur_time);

#if HAVE_THREADS
        static DWORD input_thread(__in LPVOID lParam);

        void free_input_thread(int i);

        void free_input_threads(void);

        int init_input_thread(int i);

        int init_input_threads(void);

        int get_input_packet_mt(InputFile* f, AVPacket* pkt);
#endif

        int get_input_packet(InputFile* f, AVPacket* pkt);

        int got_eagain(void);

        void reset_eagain(void);

        static AVRational duration_max(int64_t tmp, int64_t* duration, AVRational tmp_time_base, AVRational time_base);

        int seek_to_start(InputFile* ifile, AVFormatContext* is);

        int process_input(int file_index);

        int transcode_from_filter(FilterGraph* graph, InputStream** best_ist);

        int transcode_step(void);

        int transcode(void);

        BenchmarkTimeStamps get_benchmark_time_stamps(void);

        static int64_t getmaxrss(void);

        static void log_callback_null(void* ptr, int level, const char* fmt, va_list vl);

        int FFmpegMain(int argc, char** argv, int (*read_packet)(void* opaque, uint8_t* buf, int buf_size), int (*write_packet)(void* opaque, uint8_t* buf, int buf_size));

        static const enum AVPixelFormat* get_compliance_unofficial_pix_fmts(enum AVCodecID codec_id, const enum AVPixelFormat default_formats[]);

        char* choose_pix_fmts(OutputFilter* ofilter);

        char* describe_filter_link(FilterGraph* fg, AVFilterInOut* inout, int in);

        void init_input_filter(FilterGraph* fg, AVFilterInOut* in);

        int insert_trim(int64_t start_time, int64_t duration,
            AVFilterContext** last_filter, int* pad_idx,
            const char* filter_name);

        int insert_filter(AVFilterContext** last_filter, int* pad_idx,
            const char* filter_name, const char* args);

        int configure_output_video_filter(FilterGraph* fg, OutputFilter* ofilter, AVFilterInOut* out);

        int configure_output_audio_filter(FilterGraph* fg, OutputFilter* ofilter, AVFilterInOut* out);

        int sub2video_prepare(InputStream* ist, InputFilter* ifilter);

        int configure_input_video_filter(FilterGraph* fg, InputFilter* ifilter,
            AVFilterInOut* in);

        int configure_input_audio_filter(FilterGraph* fg, InputFilter* ifilter,
            AVFilterInOut* in);

        int configure_input_filter(FilterGraph* fg, InputFilter* ifilter,
            AVFilterInOut* in);

        static void cleanup_filtergraph(FilterGraph* fg);

        /* Define a function for building a string containing a list of
         * allowed formats. */    
        char* choose_sample_fmts(OutputFilter* ofilter);

        char* choose_sample_rates(OutputFilter* ofilter);

        char* choose_channel_layouts(OutputFilter* ofilter);

        void uninit_options(OptionsContext* o);

        void init_options(OptionsContext* o);

        int show_hwaccels(void* optctx, const char* opt, const char* arg);

        AVDictionary* strip_specifiers(AVDictionary* dict);

        int opt_abort_on(void* optctx, const char* opt, const char* arg);

        int opt_sameq(void* optctx, const char* opt, const char* arg);

        int opt_video_channel(void* optctx, const char* opt, const char* arg);

        int opt_video_standard(void* optctx, const char* opt, const char* arg);

        int opt_audio_codec(void* optctx, const char* opt, const char* arg);

        int opt_video_codec(void* optctx, const char* opt, const char* arg);

        int opt_subtitle_codec(void* optctx, const char* opt, const char* arg);

        int opt_data_codec(void* optctx, const char* opt, const char* arg);

        int opt_map(void* optctx, const char* opt, const char* arg);

        int opt_attach(void* optctx, const char* opt, const char* arg);

        int opt_map_channel(void* optctx, const char* opt, const char* arg);

        int opt_sdp_file(void* optctx, const char* opt, const char* arg);

        int opt_init_hw_device(void* optctx, const char* opt, const char* arg);

        int opt_filter_hw_device(void* optctx, const char* opt, const char* arg);

        void parse_meta_type(char* arg, char* type, int* index, const char** stream_spec);

        int copy_metadata(char* outspec, char* inspec, AVFormatContext* oc, AVFormatContext* ic, OptionsContext* o);

        int opt_recording_timestamp(void* optctx, const char* opt, const char* arg);

        AVCodec* find_codec_or_die(const char* name, enum AVMediaType type, int encoder);

        AVCodec* choose_decoder(OptionsContext* o, AVFormatContext* s, AVStream* st);

        void add_input_streams(OptionsContext* o, AVFormatContext* ic);

        void assert_file_overwrite(const char* filename);

        void dump_attachment(AVStream* st, const char* filename);

        int open_input_file(OptionsContext* o, const char* filename);

        uint8_t* get_line(AVIOContext* s);

        int get_preset_file_2(const char* preset_name, const char* codec_name, AVIOContext** s);

        int choose_encoder(OptionsContext* o, AVFormatContext* s, OutputStream* ost);

        OutputStream* new_output_stream(OptionsContext* o, AVFormatContext* oc, enum AVMediaType type, int source_index);

        void parse_matrix_coeffs(uint16_t* dest, const char* str);

        uint8_t* read_file(const char* filename);

        char* get_ost_filters(OptionsContext* o, AVFormatContext* oc,
            OutputStream* ost);

        void check_streamcopy_filters(OptionsContext* o, AVFormatContext* oc,
            const OutputStream* ost, enum AVMediaType type);

        OutputStream* new_video_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        OutputStream* new_audio_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        OutputStream* new_data_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        OutputStream* new_unknown_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        OutputStream* new_attachment_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        OutputStream* new_subtitle_stream(OptionsContext* o, AVFormatContext* oc, int source_index);

        int opt_streamid(void* optctx, const char* opt, const char* arg);

        int copy_chapters(InputFile* ifile, OutputFile* ofile, int copy_metadata);

        void init_output_filter(OutputFilter* ofilter, OptionsContext* o,
            AVFormatContext* oc);

        int init_complex_filters(void);

        int open_output_file(OptionsContext* o, const char* filename);

        int opt_target(void* optctx, const char* opt, const char* arg);

        int opt_vstats_file(void* optctx, const char* opt, const char* arg);

        int opt_vstats(void* optctx, const char* opt, const char* arg);

        int opt_video_frames(void* optctx, const char* opt, const char* arg);

        int opt_audio_frames(void* optctx, const char* opt, const char* arg);

        int opt_data_frames(void* optctx, const char* opt, const char* arg);

        int opt_default_new(OptionsContext* o, const char* opt, const char* arg);

        int opt_preset(void* optctx, const char* opt, const char* arg);

        int opt_old2new(void* optctx, const char* opt, const char* arg);

        static int opt_bitrate(void* optctx, const char* opt, const char* arg);

        int opt_qscale(void* optctx, const char* opt, const char* arg);

        static int opt_profile(void* optctx, const char* opt, const char* arg);

        int opt_video_filters(void* optctx, const char* opt, const char* arg);

        int opt_audio_filters(void* optctx, const char* opt, const char* arg);

        int opt_vsync(void* optctx, const char* opt, const char* arg);

        int opt_timecode(void* optctx, const char* opt, const char* arg);

        int opt_channel_layout(void* optctx, const char* opt, const char* arg);

        int opt_audio_qscale(void* optctx, const char* opt, const char* arg);

        int opt_filter_complex(void* optctx, const char* opt, const char* arg);

        int opt_filter_complex_script(void* optctx, const char* opt, const char* arg);

        int open_files(OptionGroupList* l, const char* inout, std::function<int (OptionsContext* o, const char* filename)> openFn);

        int opt_progress(void* optctx, const char* opt, const char* arg);

        HWDevice* hw_device_get_by_type(enum AVHWDeviceType type);

        HWDevice* hw_device_add(void);

        char* hw_device_default_name(enum AVHWDeviceType type);

        int hw_device_init_from_type(enum AVHWDeviceType type,
            const char* device,
            HWDevice** dev_out);

        HWDevice* hw_device_match_by_codec(const AVCodec* codec);

        static int hwaccel_retrieve_data(AVCodecContext* avctx, AVFrame* input);

        /**
         * Per-fftool specific help handler. Implemented in each
         * fftool, called by show_help().
         */
        void show_help_default(const char* opt, const char* arg);

#define OFFSET(x) offsetof(OptionsContext, x)
        const OptionDef options[336] = {
            /* main options */
            CMDUTILS_COMMON_OPTIONS
            { "f",              HAS_ARG | OPT_STRING | OPT_OFFSET |
                                OPT_INPUT | OPT_OUTPUT,                      {.off = OFFSET(format) },
                "force format", "fmt" },
            { "y",              OPT_BOOL,                                    {              &file_overwrite },
                "overwrite output files" },
            { "n",              OPT_BOOL,                                    {              &no_file_overwrite },
                "never overwrite output files" },
            { "ignore_unknown", OPT_BOOL,                                    {              &ignore_unknown_streams },
                "Ignore unknown stream types" },
            { "copy_unknown",   OPT_BOOL | OPT_EXPERT,                       {              &copy_unknown_streams },
                "Copy unknown stream types" },
            { "c",              HAS_ARG | OPT_STRING | OPT_SPEC |
                                OPT_INPUT | OPT_OUTPUT,                      {.off = OFFSET(codec_names) },
                "codec name", "codec" },
            { "codec",          HAS_ARG | OPT_STRING | OPT_SPEC |
                                OPT_INPUT | OPT_OUTPUT,                      {.off = OFFSET(codec_names) },
                "codec name", "codec" },
            { "pre",            HAS_ARG | OPT_STRING | OPT_SPEC |
                                OPT_OUTPUT,                                  {.off = OFFSET(presets) },
                "preset name", "preset" },
            { "map",            HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                                OPT_OUTPUT,                                  {.func_arg = std::bind(&FFmpeg::opt_map, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set input stream mapping",
                "[-]input_file_id[:stream_specifier][,sync_file_id[:stream_specifier]]" },
            { "map_channel",    HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_map_channel, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "map an audio channel from one stream to another", "file.stream.channel[:syncfile.syncstream]" },
            { "map_metadata",   HAS_ARG | OPT_STRING | OPT_SPEC |
                                OPT_OUTPUT,                                  {.off = OFFSET(metadata_map) },
                "set metadata information of outfile from infile",
                "outfile[,metadata]:infile[,metadata]" },
            { "map_chapters",   HAS_ARG | OPT_INT | OPT_EXPERT | OPT_OFFSET |
                                OPT_OUTPUT,                                  {.off = OFFSET(chapters_input_file) },
                "set chapters mapping", "input_file_index" },
            { "t",              HAS_ARG | OPT_TIME | OPT_OFFSET |
                                OPT_INPUT | OPT_OUTPUT,                      {.off = OFFSET(recording_time) },
                "record or transcode \"duration\" seconds of audio/video",
                "duration" },
            { "to",             HAS_ARG | OPT_TIME | OPT_OFFSET | OPT_INPUT | OPT_OUTPUT,  {.off = OFFSET(stop_time) },
                "record or transcode stop time", "time_stop" },
            { "fs",             HAS_ARG | OPT_INT64 | OPT_OFFSET | OPT_OUTPUT, {.off = OFFSET(limit_filesize) },
                "set the limit file size in bytes", "limit_size" },
            { "ss",             HAS_ARG | OPT_TIME | OPT_OFFSET |
                                OPT_INPUT | OPT_OUTPUT,                      {.off = OFFSET(start_time) },
                "set the start time offset", "time_off" },
            { "sseof",          HAS_ARG | OPT_TIME | OPT_OFFSET |
                                OPT_INPUT,                                   {.off = OFFSET(start_time_eof) },
                "set the start time offset relative to EOF", "time_off" },
            { "seek_timestamp", HAS_ARG | OPT_INT | OPT_OFFSET |
                                OPT_INPUT,                                   {.off = OFFSET(seek_timestamp) },
                "enable/disable seeking by timestamp with -ss" },
            { "accurate_seek",  OPT_BOOL | OPT_OFFSET | OPT_EXPERT |
                                OPT_INPUT,                                   {.off = OFFSET(accurate_seek) },
                "enable/disable accurate seeking with -ss" },
            { "itsoffset",      HAS_ARG | OPT_TIME | OPT_OFFSET |
                                OPT_EXPERT | OPT_INPUT,                      {.off = OFFSET(input_ts_offset) },
                "set the input ts offset", "time_off" },
            { "itsscale",       HAS_ARG | OPT_DOUBLE | OPT_SPEC |
                                OPT_EXPERT | OPT_INPUT,                      {.off = OFFSET(ts_scale) },
                "set the input ts scale", "scale" },
            { "timestamp",      HAS_ARG | OPT_PERFILE | OPT_OUTPUT,          {.func_arg = std::bind(&FFmpeg::opt_recording_timestamp, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the recording timestamp ('now' to set the current time)", "time" },
            { "metadata",       HAS_ARG | OPT_STRING | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(metadata) },
                "add metadata", "string=string" },
            { "program",        HAS_ARG | OPT_STRING | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(program) },
                "add program with specified streams", "title=string:st=number..." },
            { "dframes",        HAS_ARG | OPT_PERFILE | OPT_EXPERT |
                                OPT_OUTPUT,                                  {.func_arg = std::bind(&FFmpeg::opt_data_frames, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the number of data frames to output", "number" },
            { "benchmark",      OPT_BOOL | OPT_EXPERT,                       { &do_benchmark },
                "add timings for benchmarking" },
            { "benchmark_all",  OPT_BOOL | OPT_EXPERT,                       { &do_benchmark_all },
              "add timings for each task" },
            { "progress",       HAS_ARG | OPT_EXPERT,                        {.func_arg = std::bind(&FFmpeg::opt_progress, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
              "write program-readable progress information", "url" },
            { "stdin",          OPT_BOOL | OPT_EXPERT,                       { &stdin_interaction },
              "enable or disable interaction on standard input" },            
            { "dump",           OPT_BOOL | OPT_EXPERT,                       { &do_pkt_dump },
                "dump each input packet" },
            { "hex",            OPT_BOOL | OPT_EXPERT,                       { &do_hex_dump },
                "when dumping packets, also dump the payload" },
            { "re",             OPT_BOOL | OPT_EXPERT | OPT_OFFSET |
                                OPT_INPUT,                                   {.off = OFFSET(rate_emu) },
                "read input at native frame rate", "" },
            { "target",         HAS_ARG | OPT_PERFILE | OPT_OUTPUT,          {.func_arg = std::bind(&FFmpeg::opt_target, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "specify target file type (\"vcd\", \"svcd\", \"dvd\", \"dv\" or \"dv50\" "
                "with optional prefixes \"pal-\", \"ntsc-\" or \"film-\")", "type" },
            { "vsync",          HAS_ARG | OPT_EXPERT,                        {.func_arg = std::bind(&FFmpeg::opt_vsync, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "video sync method", "" },
            { "frame_drop_threshold", HAS_ARG | OPT_FLOAT | OPT_EXPERT,      { &frame_drop_threshold },
                "frame drop threshold", "" },
            { "async",          HAS_ARG | OPT_INT | OPT_EXPERT,              { &audio_sync_method },
                "audio sync method", "" },
            { "adrift_threshold", HAS_ARG | OPT_FLOAT | OPT_EXPERT,          { &audio_drift_threshold },
                "audio drift threshold", "threshold" },
            { "copyts",         OPT_BOOL | OPT_EXPERT,                       { &copy_ts },
                "copy timestamps" },
            { "start_at_zero",  OPT_BOOL | OPT_EXPERT,                       { &start_at_zero },
                "shift input timestamps to start at 0 when using copyts" },
            { "copytb",         HAS_ARG | OPT_INT | OPT_EXPERT,              { &copy_tb },
                "copy input stream time base when stream copying", "mode" },
            { "shortest",       OPT_BOOL | OPT_EXPERT | OPT_OFFSET |
                                OPT_OUTPUT,                                  {.off = OFFSET(shortest) },
                "finish encoding within shortest input" },
            { "bitexact",       OPT_BOOL | OPT_EXPERT | OPT_OFFSET |
                                OPT_OUTPUT | OPT_INPUT,                      {.off = OFFSET(bitexact) },
                "bitexact mode" },
            { "apad",           OPT_STRING | HAS_ARG | OPT_SPEC |
                                OPT_OUTPUT,                                  {.off = OFFSET(apad) },
                "audio pad", "" },
            { "dts_delta_threshold", HAS_ARG | OPT_FLOAT | OPT_EXPERT,       { &dts_delta_threshold },
                "timestamp discontinuity delta threshold", "threshold" },
            { "dts_error_threshold", HAS_ARG | OPT_FLOAT | OPT_EXPERT,       { &dts_error_threshold },
                "timestamp error delta threshold", "threshold" },
            { "xerror",         OPT_BOOL | OPT_EXPERT,                       { &exit_on_error },
                "exit on error", "error" },
            { "abort_on",       HAS_ARG | OPT_EXPERT,                        {.func_arg = std::bind(&FFmpeg::opt_abort_on, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "abort on the specified condition flags", "flags" },
            { "copyinkf",       OPT_BOOL | OPT_EXPERT | OPT_SPEC |
                                OPT_OUTPUT,                                  {.off = OFFSET(copy_initial_nonkeyframes) },
                "copy initial non-keyframes" },
            { "copypriorss",    OPT_INT | HAS_ARG | OPT_EXPERT | OPT_SPEC | OPT_OUTPUT,   {.off = OFFSET(copy_prior_start) },
                "copy or discard frames before start time" },
            { "frames",         OPT_INT64 | HAS_ARG | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(max_frames) },
                "set the number of frames to output", "number" },
            { "tag",            OPT_STRING | HAS_ARG | OPT_SPEC |
                                OPT_EXPERT | OPT_OUTPUT | OPT_INPUT,         {.off = OFFSET(codec_tags) },
                "force codec tag/fourcc", "fourcc/tag" },
            { "q",              HAS_ARG | OPT_EXPERT | OPT_DOUBLE |
                                OPT_SPEC | OPT_OUTPUT,                       {.off = OFFSET(qscale) },
                "use fixed quality scale (VBR)", "q" },
            { "qscale",         HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                                OPT_OUTPUT,                                  {.func_arg = std::bind(&FFmpeg::opt_qscale, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "use fixed quality scale (VBR)", "q" },
            { "profile",        HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = opt_profile },
                "set profile", "profile" },
            { "filter",         HAS_ARG | OPT_STRING | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(filters) },
                "set stream filtergraph", "filter_graph" },
            { "filter_threads",  HAS_ARG | OPT_INT,                          { &filter_nbthreads },
                "number of non-complex filter threads" },
            { "filter_script",  HAS_ARG | OPT_STRING | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(filter_scripts) },
                "read stream filtergraph description from a file", "filename" },
            { "reinit_filter",  HAS_ARG | OPT_INT | OPT_SPEC | OPT_INPUT,    {.off = OFFSET(reinit_filters) },
                "reinit filtergraph on input parameter changes", "" },
            { "filter_complex", HAS_ARG | OPT_EXPERT,                        {.func_arg = std::bind(&FFmpeg::opt_filter_complex, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "create a complex filtergraph", "graph_description" },
            { "filter_complex_threads", HAS_ARG | OPT_INT,                   { &filter_complex_nbthreads },
                "number of threads for -filter_complex" },
            { "lavfi",          HAS_ARG | OPT_EXPERT,                        {.func_arg = std::bind(&FFmpeg::opt_filter_complex, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "create a complex filtergraph", "graph_description" },
            { "filter_complex_script", HAS_ARG | OPT_EXPERT,                 {.func_arg = std::bind(&FFmpeg::opt_filter_complex_script, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "read complex filtergraph description from a file", "filename" },
            { "auto_conversion_filters", OPT_BOOL | OPT_EXPERT,              { &auto_conversion_filters },
                "enable automatic conversion filters globally" },
            { "stats",          OPT_BOOL,                                    { &print_stats },
                "print progress report during encoding", },
            { "attach",         HAS_ARG | OPT_PERFILE | OPT_EXPERT |
                                OPT_OUTPUT,                                  {.func_arg = std::bind(&FFmpeg::opt_attach, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "add an attachment to the output file", "filename" },
            { "dump_attachment", HAS_ARG | OPT_STRING | OPT_SPEC |
                                 OPT_EXPERT | OPT_INPUT,                     {.off = OFFSET(dump_attachment) },
                "extract an attachment into a file", "filename" },
            { "stream_loop", OPT_INT | HAS_ARG | OPT_EXPERT | OPT_INPUT |
                                OPT_OFFSET,                                  {.off = OFFSET(loop) }, "set number of times input stream shall be looped", "loop count" },
            { "debug_ts",       OPT_BOOL | OPT_EXPERT,                       { &debug_ts },
                "print timestamp debugging info" },
            { "max_error_rate",  HAS_ARG | OPT_FLOAT,                        { &max_error_rate },
                "ratio of errors (0.0: no errors, 1.0: 100% errors) above which ffmpeg returns an error instead of success.", "maximum error rate" },
            { "discard",        OPT_STRING | HAS_ARG | OPT_SPEC |
                                OPT_INPUT,                                   {.off = OFFSET(discard) },
                "discard", "" },
            { "disposition",    OPT_STRING | HAS_ARG | OPT_SPEC |
                                OPT_OUTPUT,                                  {.off = OFFSET(disposition) },
                "disposition", "" },
            { "thread_queue_size", HAS_ARG | OPT_INT | OPT_OFFSET | OPT_EXPERT | OPT_INPUT,
                                                                             {.off = OFFSET(thread_queue_size) },
                "set the maximum number of queued packets from the demuxer" },
            { "find_stream_info", OPT_BOOL | OPT_PERFILE | OPT_INPUT | OPT_EXPERT, { &find_stream_info },
                "read and decode the streams to fill missing information with heuristics" },

            /* video options */
            { "vframes",      OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,           {.func_arg = std::bind(&FFmpeg::opt_video_frames, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the number of video frames to output", "number" },
            { "r",            OPT_VIDEO | HAS_ARG | OPT_STRING | OPT_SPEC |
                              OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(frame_rates) },
                "set frame rate (Hz value, fraction or abbreviation)", "rate" },
            { "s",            OPT_VIDEO | HAS_ARG | OPT_SUBTITLE | OPT_STRING | OPT_SPEC |
                              OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(frame_sizes) },
                "set frame size (WxH or abbreviation)", "size" },
            { "aspect",       OPT_VIDEO | HAS_ARG | OPT_STRING | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(frame_aspect_ratios) },
                "set aspect ratio (4:3, 16:9 or 1.3333, 1.7777)", "aspect" },
            { "pix_fmt",      OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_STRING | OPT_SPEC |
                              OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(frame_pix_fmts) },
                "set pixel format", "format" },
            { "bits_per_raw_sample", OPT_VIDEO | OPT_INT | HAS_ARG,                      { &frame_bits_per_raw_sample },
                "set the number of bits per raw sample", "number" },
            { "intra",        OPT_VIDEO | OPT_BOOL | OPT_EXPERT,                         { &intra_only },
                "deprecated use -g 1" },
            { "vn",           OPT_VIDEO | OPT_BOOL | OPT_OFFSET | OPT_INPUT | OPT_OUTPUT,{.off = OFFSET(video_disable) },
                "disable video" },
            { "rc_override",  OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_STRING | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(rc_overrides) },
                "rate control override for specific intervals", "override" },
            { "vcodec",       OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_INPUT |
                              OPT_OUTPUT,                                                {.func_arg = std::bind(&FFmpeg::opt_video_codec, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force video codec ('copy' to copy stream)", "codec" },
            { "sameq",        OPT_VIDEO | OPT_EXPERT ,                                   {.func_arg = std::bind(&FFmpeg::opt_sameq, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "Removed" },
            { "same_quant",   OPT_VIDEO | OPT_EXPERT ,                                   {.func_arg = std::bind(&FFmpeg::opt_sameq, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "Removed" },
            { "timecode",     OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,            {.func_arg = std::bind(&FFmpeg::opt_timecode, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set initial TimeCode value.", "hh:mm:ss[:;.]ff" },
            { "pass",         OPT_VIDEO | HAS_ARG | OPT_SPEC | OPT_INT | OPT_OUTPUT,     {.off = OFFSET(pass) },
                "select the pass number (1 to 3)", "n" },
            { "passlogfile",  OPT_VIDEO | HAS_ARG | OPT_STRING | OPT_EXPERT | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(passlogfiles) },
                "select two pass log file name prefix", "prefix" },
            { "deinterlace",  OPT_VIDEO | OPT_BOOL | OPT_EXPERT,                         { &do_deinterlace },
                "this option is deprecated, use the yadif filter instead" },
            { "psnr",         OPT_VIDEO | OPT_BOOL | OPT_EXPERT,                         { &do_psnr },
                "calculate PSNR of compressed frames" },
            { "vstats",       OPT_VIDEO | OPT_EXPERT ,                                   {.func_arg = std::bind(&FFmpeg::opt_vstats, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "dump video coding statistics to file" },
            { "vstats_file",  OPT_VIDEO | HAS_ARG | OPT_EXPERT ,                         {.func_arg = std::bind(&FFmpeg::opt_vstats_file, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "dump video coding statistics to file", "file" },
            { "vstats_version",  OPT_VIDEO | OPT_INT | HAS_ARG | OPT_EXPERT ,            { &vstats_version },
                "Version of the vstats format to use."},
            { "vf",           OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,           {.func_arg = std::bind(&FFmpeg::opt_video_filters, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set video filters", "filter_graph" },
            { "intra_matrix", OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_STRING | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(intra_matrices) },
                "specify intra matrix coeffs", "matrix" },
            { "inter_matrix", OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_STRING | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(inter_matrices) },
                "specify inter matrix coeffs", "matrix" },
            { "chroma_intra_matrix", OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_STRING | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(chroma_intra_matrices) },
                "specify intra matrix coeffs", "matrix" },
            { "top",          OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_INT | OPT_SPEC |
                              OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(top_field_first) },
                "top=1/bottom=0/auto=-1 field first", "" },
            { "vtag",         OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                              OPT_INPUT | OPT_OUTPUT,                                    {.func_arg = std::bind(&FFmpeg::opt_old2new, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force video tag/fourcc", "fourcc/tag" },
            { "qphist",       OPT_VIDEO | OPT_BOOL | OPT_EXPERT ,                        { &qp_hist },
                "show QP histogram" },
            { "force_fps",    OPT_VIDEO | OPT_BOOL | OPT_EXPERT | OPT_SPEC |
                              OPT_OUTPUT,                                                {.off = OFFSET(force_fps) },
                "force the selected framerate, disable the best supported framerate selection" },
            { "streamid",     OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                              OPT_OUTPUT,                                                {.func_arg = std::bind(&FFmpeg::opt_streamid, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the value of an outfile streamid", "streamIndex:value" },
            { "force_key_frames", OPT_VIDEO | OPT_STRING | HAS_ARG | OPT_EXPERT |
                                  OPT_SPEC | OPT_OUTPUT,                                 {.off = OFFSET(forced_key_frames) },
                "force key frames at specified timestamps", "timestamps" },
            { "ab",           OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,            {.func_arg = opt_bitrate },
                "audio bitrate (please use -b:a)", "bitrate" },
            { "b",            OPT_VIDEO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,            {.func_arg = opt_bitrate },
                "video bitrate (please use -b:v)", "bitrate" },
            { "hwaccel",          OPT_VIDEO | OPT_STRING | HAS_ARG | OPT_EXPERT |
                                  OPT_SPEC | OPT_INPUT,                                  {.off = OFFSET(hwaccels) },
                "use HW accelerated decoding", "hwaccel name" },
            { "hwaccel_device",   OPT_VIDEO | OPT_STRING | HAS_ARG | OPT_EXPERT |
                                  OPT_SPEC | OPT_INPUT,                                  {.off = OFFSET(hwaccel_devices) },
                "select a device for HW acceleration", "devicename" },
            { "hwaccel_output_format", OPT_VIDEO | OPT_STRING | HAS_ARG | OPT_EXPERT |
                                  OPT_SPEC | OPT_INPUT,                                  {.off = OFFSET(hwaccel_output_formats) },
                "select output format used with HW accelerated decoding", "format" },
        #if CONFIG_VIDEOTOOLBOX
            { "videotoolbox_pixfmt", HAS_ARG | OPT_STRING | OPT_EXPERT, { &videotoolbox_pixfmt}, "" },
        #endif
            { "hwaccels",         OPT_EXIT,                                              {.func_arg = std::bind(&FFmpeg::show_hwaccels, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "show available HW acceleration methods" },
            { "autorotate",       HAS_ARG | OPT_BOOL | OPT_SPEC |
                                  OPT_EXPERT | OPT_INPUT,                                {.off = OFFSET(autorotate) },
                "automatically insert correct rotate filters" },
            { "autoscale",        HAS_ARG | OPT_BOOL | OPT_SPEC |
                                  OPT_EXPERT | OPT_OUTPUT,                               {.off = OFFSET(autoscale) },
                "automatically insert a scale filter at the end of the filter graph" },

            /* audio options */
            { "aframes",        OPT_AUDIO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,           {.func_arg = std::bind(&FFmpeg::opt_audio_frames, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the number of audio frames to output", "number" },
            { "aq",             OPT_AUDIO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,           {.func_arg = std::bind(&FFmpeg::opt_audio_qscale, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set audio quality (codec-specific)", "quality", },
            { "ar",             OPT_AUDIO | HAS_ARG | OPT_INT | OPT_SPEC |
                                OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(audio_sample_rate) },
                "set audio sampling rate (in Hz)", "rate" },
            { "ac",             OPT_AUDIO | HAS_ARG | OPT_INT | OPT_SPEC |
                                OPT_INPUT | OPT_OUTPUT,                                    {.off = OFFSET(audio_channels) },
                "set number of audio channels", "channels" },
            { "an",             OPT_AUDIO | OPT_BOOL | OPT_OFFSET | OPT_INPUT | OPT_OUTPUT,{.off = OFFSET(audio_disable) },
                "disable audio" },
            { "acodec",         OPT_AUDIO | HAS_ARG | OPT_PERFILE |
                                OPT_INPUT | OPT_OUTPUT,                                    {.func_arg = std::bind(&FFmpeg::opt_audio_codec, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force audio codec ('copy' to copy stream)", "codec" },
            { "atag",           OPT_AUDIO | HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                                OPT_OUTPUT,                                                {.func_arg = std::bind(&FFmpeg::opt_old2new, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force audio tag/fourcc", "fourcc/tag" },
            { "vol",            OPT_AUDIO | HAS_ARG | OPT_INT,                            { &audio_volume },
                "change audio volume (256=normal)" , "volume" },
            { "sample_fmt",     OPT_AUDIO | HAS_ARG | OPT_EXPERT | OPT_SPEC |
                                OPT_STRING | OPT_INPUT | OPT_OUTPUT,                       {.off = OFFSET(sample_fmts) },
                "set sample format", "format" },
            { "channel_layout", OPT_AUDIO | HAS_ARG | OPT_EXPERT | OPT_PERFILE |
                                OPT_INPUT | OPT_OUTPUT,                                    {.func_arg = std::bind(&FFmpeg::opt_channel_layout, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set channel layout", "layout" },
            { "af",             OPT_AUDIO | HAS_ARG | OPT_PERFILE | OPT_OUTPUT,           {.func_arg = std::bind(&FFmpeg::opt_audio_filters, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set audio filters", "filter_graph" },
            { "guess_layout_max", OPT_AUDIO | HAS_ARG | OPT_INT | OPT_SPEC | OPT_EXPERT | OPT_INPUT, {.off = OFFSET(guess_layout_max) },
              "set the maximum number of channels to try to guess the channel layout" },

            /* subtitle options */
            { "sn",     OPT_SUBTITLE | OPT_BOOL | OPT_OFFSET | OPT_INPUT | OPT_OUTPUT, {.off = OFFSET(subtitle_disable) },
                "disable subtitle" },
            { "scodec", OPT_SUBTITLE | HAS_ARG | OPT_PERFILE | OPT_INPUT | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_subtitle_codec, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force subtitle codec ('copy' to copy stream)", "codec" },
            { "stag",   OPT_SUBTITLE | HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_old2new, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) }
                , "force subtitle tag/fourcc", "fourcc/tag" },
            { "fix_sub_duration", OPT_BOOL | OPT_EXPERT | OPT_SUBTITLE | OPT_SPEC | OPT_INPUT, {.off = OFFSET(fix_sub_duration) },
                "fix subtitles duration" },
            { "canvas_size", OPT_SUBTITLE | HAS_ARG | OPT_STRING | OPT_SPEC | OPT_INPUT, {.off = OFFSET(canvas_sizes) },
                "set canvas size (WxH or abbreviation)", "size" },

            /* grab options */
            { "vc", HAS_ARG | OPT_EXPERT | OPT_VIDEO, {.func_arg = std::bind(&FFmpeg::opt_video_channel, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "deprecated, use -channel", "channel" },
            { "tvstd", HAS_ARG | OPT_EXPERT | OPT_VIDEO, {.func_arg = std::bind(&FFmpeg::opt_video_standard, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "deprecated, use -standard", "standard" },
            { "isync", OPT_BOOL | OPT_EXPERT, { &input_sync }, "this option is deprecated and does nothing", "" },

            /* muxer options */
            { "muxdelay",   OPT_FLOAT | HAS_ARG | OPT_EXPERT | OPT_OFFSET | OPT_OUTPUT, {.off = OFFSET(mux_max_delay) },
                "set the maximum demux-decode delay", "seconds" },
            { "muxpreload", OPT_FLOAT | HAS_ARG | OPT_EXPERT | OPT_OFFSET | OPT_OUTPUT, {.off = OFFSET(mux_preload) },
                "set the initial demux-decode delay", "seconds" },
            { "sdp_file", HAS_ARG | OPT_EXPERT | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_sdp_file, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "specify a file in which to print sdp information", "file" },

            { "time_base", HAS_ARG | OPT_STRING | OPT_EXPERT | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(time_bases) },
                "set the desired time base hint for output stream (1:24, 1:48000 or 0.04166, 2.0833e-5)", "ratio" },
            { "enc_time_base", HAS_ARG | OPT_STRING | OPT_EXPERT | OPT_SPEC | OPT_OUTPUT, {.off = OFFSET(enc_time_bases) },
                "set the desired time base for the encoder (1:24, 1:48000 or 0.04166, 2.0833e-5). "
                "two special values are defined - "
                "0 = use frame rate (video) or sample rate (audio),"
                "-1 = match source time base", "ratio" },

            { "bsf", HAS_ARG | OPT_STRING | OPT_SPEC | OPT_EXPERT | OPT_OUTPUT, {.off = OFFSET(bitstream_filters) },
                "A comma-separated list of bitstream filters", "bitstream_filters" },
            { "absf", HAS_ARG | OPT_AUDIO | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_old2new, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "deprecated", "audio bitstream_filters" },
            { "vbsf", OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_old2new, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "deprecated", "video bitstream_filters" },

            { "apre", HAS_ARG | OPT_AUDIO | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT,    {.func_arg = std::bind(&FFmpeg::opt_preset, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the audio options to the indicated preset", "preset" },
            { "vpre", OPT_VIDEO | HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT,    {.func_arg = std::bind(&FFmpeg::opt_preset, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the video options to the indicated preset", "preset" },
            { "spre", HAS_ARG | OPT_SUBTITLE | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_preset, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set the subtitle options to the indicated preset", "preset" },
            { "fpre", HAS_ARG | OPT_EXPERT | OPT_PERFILE | OPT_OUTPUT,                {.func_arg = std::bind(&FFmpeg::opt_preset, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set options from indicated preset file", "filename" },

            { "max_muxing_queue_size", HAS_ARG | OPT_INT | OPT_SPEC | OPT_EXPERT | OPT_OUTPUT, {.off = OFFSET(max_muxing_queue_size) },
                "maximum number of packets that can be buffered while waiting for all streams to initialize", "packets" },
            { "muxing_queue_data_threshold", HAS_ARG | OPT_INT | OPT_SPEC | OPT_EXPERT | OPT_OUTPUT, {.off = OFFSET(muxing_queue_data_threshold) },
                "set the threshold after which max_muxing_queue_size is taken into account", "bytes" },

            /* data codec support */
            { "dcodec", HAS_ARG | OPT_DATA | OPT_PERFILE | OPT_EXPERT | OPT_INPUT | OPT_OUTPUT, {.func_arg = std::bind(&FFmpeg::opt_data_codec, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "force data codec ('copy' to copy stream)", "codec" },
            { "dn", OPT_BOOL | OPT_VIDEO | OPT_OFFSET | OPT_INPUT | OPT_OUTPUT, {.off = OFFSET(data_disable) },
                "disable data" },

        #if CONFIG_VAAPI
            { "vaapi_device", HAS_ARG | OPT_EXPERT, {.func_arg = opt_vaapi_device },
                "set VAAPI hardware device (DRM path or X11 display name)", "device" },
        #endif

        #if CONFIG_QSV
            { "qsv_device", HAS_ARG | OPT_STRING | OPT_EXPERT, { &qsv_device },
                "set QSV hardware device (DirectX adapter index, DRM path or X11 display name)", "device"},
        #endif

            { "init_hw_device", HAS_ARG | OPT_EXPERT, {.func_arg = std::bind(&FFmpeg::opt_init_hw_device, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "initialise hardware device", "args" },
            { "filter_hw_device", HAS_ARG | OPT_EXPERT, {.func_arg = std::bind(&FFmpeg::opt_filter_hw_device, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3) },
                "set hardware device used when filtering", "device" },

            { NULL, },
        };
    };
}
