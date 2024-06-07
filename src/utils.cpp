/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include "utils.hpp"

/********************
 * System utilities *
 ********************/

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

void _err_fatal_simple_core(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) 
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

int err_func_printf(const char *func, const char *format, ...)
{
	// 构造格式化时间字符串
	time_t rawtime;
    struct tm *info;
    char time_buf[80];
    time(&rawtime);
    info = localtime( &rawtime );
    strftime(time_buf,80,"%Y-%m-%d %X", info);

	// 打印函数名以及时间
    fprintf(stderr, "[%s - %s] ", func, time_buf);
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stderr, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stderr)", strerror(saveErrno));
	return done;
}

int err_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stderr, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stderr)", strerror(saveErrno));
	return done;
}

int stdout_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

void err_fgets(char *buff, size_t s, FILE *fp)
{
    if (fgets(buff, s, fp) == NULL) {
        err_fatal_simple("fgets error.\n");
    }
}

int err_puts(const char *s)
{
	int ret = puts(s);
	if (EOF == ret)
	{
		_err_fatal_simple("puts", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));
		
		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
	int ret = fclose(stream);
	if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}

/*********
 * alloc *
 *********/
void *err_malloc(const char *func, size_t s)
{
    void *ret = (void*)malloc(s);
    if (ret == NULL) err_fatal_core(func, "Malloc fail!\nSize: %lld\n", s);
    else return ret;
}

void *err_calloc(const char *func, size_t n, size_t s)
{
    void *ret = (void*)calloc(n, s);
    if (ret == NULL) err_fatal_core(func, "Calloc fail!\nN: %d\tSize: %lld\n", n, s);
    else return ret;
}

void *err_realloc(const char *func, void *p, size_t s)
{
    void *ret = (void*)realloc(p, s);
    if (ret == NULL) err_fatal_core(func, "Realloc fail!\nSize: %lld\n", s);
    else return ret;
}

/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void get_cur_time(const char *prefix)
{
    time_t now = time(0);
    struct tm ts; char buf[1024];
    ts = *localtime(&now);
    err_printf("[%s] ", prefix);
    strftime(buf, sizeof(buf), "%Y-%m-%d-%s", &ts);
}

void print_format_time(FILE *out)
{
    time_t rawtime;
    struct tm *info;
    char time_buf[80];

    time(&rawtime);
    info = localtime( &rawtime );
    strftime(time_buf,80,"%m-%d-%Y %X", info);
    fprintf(out, "=== %s === ", time_buf);
}

int err_func_format_printf(const char *func, const char *format, ...)
{
    print_format_time(stderr);
    fprintf(stderr, "[%s] ", func);
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stderr, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stderr)", strerror(saveErrno));
	return done;
}


/*********
 * String *
 *********/

#include <string>
using namespace std;
string uint64_to_str(uint64_t kmer, uint32_t k) {
	string a;
    char alphabet[4] = {'A','C','G','T'};
	for (uint32_t j = 0; j < k; ++j) {
		a += alphabet[(kmer >> ((k - 1 - j) * 2)) & 0x3];
	}
	return a;
}

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t str_to_uint64(string kmerStr, uint32_t k) {
	uint64_t a = 0;
	for (uint32_t i = 0; i < k; ++i) {
		a = (a << 2) | (uint32_t)(nst_nt4_table[(int)kmerStr[i]]);
	}
	return a;
}

std::string string_format(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    // 第一次调用 vsnprintf 获取格式化字符串的长度
    int len = vsnprintf(nullptr, 0, fmt, args);
    va_end(args);

    // 分配足够的内存来存储格式化的字符串
    char *buffer = new char[len + 1];

    va_start(args, fmt);
    // 第二次调用 vsnprintf 将格式化字符串写入缓冲区
    vsnprintf(buffer, len + 1, fmt, args);
    va_end(args);

    std::string result(buffer);
    delete[] buffer;
    return result;
}

#include <vector>
#include <sstream>
using namespace std;
vector<string> split(string s, char sep)
{
    istringstream iss(s);
    vector<string> res;
    string buffer;
    while(getline(iss, buffer, sep)){
        res.push_back(buffer);
    }
    return res;
}


//for kseq

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#define ks_err(ks) ((ks)->end == -1)
#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)


#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline kstream_t *ks_init(type_t f)						\
{																\
	kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
	ks->f = f;													\
	ks->buf = (unsigned char*)malloc(__bufsize);				\
	return ks;													\
}																\
static inline void ks_destroy(kstream_t *ks)					\
{																\
	if (ks) {													\
		free(ks->buf);											\
		free(ks);												\
	}															\
}

static inline int ks_getc(kstream_t *ks)				\
{														\
	if (ks_err(ks)) return -3;							\
	if (ks->is_eof && ks->begin >= ks->end) return -1;	\
	if (ks->begin >= ks->end) {							\
		ks->begin = 0;									\
		ks->end = __read(ks->f, ks->buf, __bufsize);	\
		if (ks->end == 0) { ks->is_eof = 1; return -1;}	\
		if (ks->end == -1) { ks->is_eof = 1; return -3;}\
	}													\
	return (int)ks->buf[ks->begin++];					\
}

static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append) \
{																	\
	int gotany = 0;													\
	if (dret) *dret = 0;											\
	str->l = append? str->l : 0;									\
	for (;;) {														\
		int i;														\
		if (ks_err(ks)) return -3;									\
		if (ks->begin >= ks->end) {									\
			if (!ks->is_eof) {										\
				ks->begin = 0;										\
				ks->end = __read(ks->f, ks->buf, __bufsize);		\
				if (ks->end == 0) { ks->is_eof = 1; break; }		\
				if (ks->end == -1) { ks->is_eof = 1; return -3; }	\
			} else break;											\
		}															\
		if (delimiter == KS_SEP_LINE) {								\
			unsigned char *sep = (unsigned char *)memchr(ks->buf + ks->begin, '\n', ks->end - ks->begin); \
			i = sep != NULL ? sep - ks->buf : ks->end;				\
		} else if (delimiter > KS_SEP_MAX) {						\
			unsigned char *sep = (unsigned char *)memchr(ks->buf + ks->begin, delimiter, ks->end - ks->begin); \
			i = sep != NULL ? sep - ks->buf : ks->end;				\
		} else if (delimiter == KS_SEP_SPACE) {						\
			for (i = ks->begin; i < ks->end; ++i)					\
				if (isspace(ks->buf[i])) break;						\
		} else if (delimiter == KS_SEP_TAB) {						\
			for (i = ks->begin; i < ks->end; ++i)					\
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
		} else i = 0; /* never come to here! */						\
		if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
			str->m = str->l + (i - ks->begin) + 1;					\
			kroundup32(str->m);										\
			str->s = (char*)realloc(str->s, str->m);				\
		}															\
		gotany = 1;													\
		memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
		str->l = str->l + (i - ks->begin);							\
		ks->begin = i + 1;											\
		if (i < ks->end) {											\
			if (dret) *dret = ks->buf[i];							\
			break;													\
		}															\
	}																\
	if (!gotany && ks_eof(ks)) return -1;							\
	if (str->s == 0) {												\
		str->m = 1;													\
		str->s = (char*)calloc(1, 1);								\
	} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
	str->s[str->l] = '\0';											\
	return str->l;													\
} \
static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

kseq_t *kseq_init(type_t fd)									\
{																	\
	kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
	s->f = ks_init(fd);												\
	return s;														\
}																	\
void kseq_destroy(kseq_t *ks)									\
{																	\
	if (!ks) return;												\
	free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
	ks_destroy(ks->f);												\
	free(ks);														\
}

int kseq_read(kseq_t *seq) \
{ \
	int c,r; \
	kstream_t *ks = seq->f; \
	if (seq->last_char == 0) { /* then jump to the next header line */ \
		while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '@'); \
		if (c < 0) return c; /* end of file or error*/ \
		seq->last_char = c; \
	} /* else: the first header char has been read in the previous call */ \
	seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
	if ((r=ks_getuntil(ks, 0, &seq->name, &c)) < 0) return r;  /* normal exit: EOF or error */ \
	if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
	if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
		seq->seq.m = 256; \
		seq->seq.s = (char*)malloc(seq->seq.m); \
	} \
	while ((c = ks_getc(ks)) >= 0 && c != '>' && c != '+' && c != '@') { \
		if (c == '\n') continue; /* skip empty lines */ \
		seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
		ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
	} \
	if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
	if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
		seq->seq.m = seq->seq.l + 2; \
		kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
		seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
	} \
	seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
	if (c != '+') return seq->seq.l; /* FASTA */ \
	if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
		seq->qual.m = seq->seq.m; \
		seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
	} \
	while ((c = ks_getc(ks)) >= 0 && c != '\n'); /* skip the rest of '+' line */ \
	if (c == -1) return -2; /* error: no quality string */ \
	while ((c = ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l)); \
	if (c == -3) return -3; /* stream error */ \
	seq->last_char = 0;	/* we have not come to the next header line */ \
	if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
	return seq->seq.l; \
}
