#include <stdio.h>
#include <stdint.h>
#include "umpriv.h"
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop

#include "khashl.h" // hash table
#define KC_BITS 12
#define KC_MAX ((1<<KC_BITS) - 1)
#define cnt_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 8 bits for counts; higher bits for k-mer
#define cnt_hash(a) ((a)>>KC_BITS)
KHASHL_SET_INIT(KH_LOCAL, cnthash_t, cnt_h, uint64_t, cnt_hash, cnt_eq)

#define XCALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define XMALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define XREALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

typedef struct {
	int p; // suffix length; at least 8
	uint64_t n_tot, n_dup;
	cnthash_t **h; // 1<<p hash tables
} cnthashx_t;

static cnthashx_t *chx_init(int p)
{
	int i;
	cnthashx_t *h;
	XCALLOC(h, 1);
	XMALLOC(h->h, 1<<p);
	h->p = p;
	for (i = 0; i < 1<<p; ++i)
		h->h[i] = cnt_h_init();
	return h;
}

typedef struct {
	int n, m;
	uint64_t *a;
} cntbuf_t;

static inline void chx_insert_buf(cntbuf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	cntbuf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		XREALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static void count_seq_buf(const mm_idx_t *mi, cntbuf_t *buf, int k, int p, int rid) // insert k-mers in $seq to linear buffer $buf
{
	int i, l, len = mi->seq[rid].len;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2, off = mi->seq[rid].offset;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = mm_seq4_get(mi->S, off + i);
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				chx_insert_buf(buf, p, um_hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	int k, n_thread;
	int rid0;
	uint64_t mini_batch_size;
	const mm_idx_t *mi;
	cnthashx_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n;
	uint64_t sum_len, nk;
	int *len;
	char **seq;
	cntbuf_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	cntbuf_t *b = &s->buf[i];
	cnthash_t *h = s->p->h->h[i];
	int j, p = s->p->h->p;
	for (j = 0; j < b->n; ++j) {
		khint_t k;
		int absent;
		k = cnt_h_put(h, b->a[j]>>p<<KC_BITS, &absent);
		if ((kh_key(h, k)&KC_MAX) < KC_MAX) ++kh_key(h, k);
	}
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int i;
		stepdat_t *s;
		XCALLOC(s, 1);
		s->p = p;
		for (i = p->rid0; i < (int)p->mi->n_seq; ++i) {
			s->sum_len += p->mi->seq[i].len;
			assert(p->mi->seq[i].len >= p->k);
			s->nk += p->mi->seq[i].len - p->k + 1;
			if (s->sum_len >= p->mini_batch_size)
				break;
		}
		s->n = i - p->rid0;
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->h->p, m;
		XCALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			XMALLOC(s->buf[i].a, m);
		}
		for (i = 0; i < s->n; ++i)
			count_seq_buf(p->mi, s->buf, p->k, p->h->p, p->rid0 + i);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->h->p;
		kt_for(p->n_thread, worker_for, s, n);
		p->rid0 += s->n;
		for (i = 0; i < n; ++i) free(s->buf[i].a);
		free(s->buf); free(s);
	}
	return 0;
}

// shrinking
typedef struct {
	int min, max;
	cnthashx_t *h;
} shrink_aux_t;

static void worker_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	cnthashx_t *h = a->h;
	cnthash_t *g = h->h[i], *f;
	khint_t k;
	int size = 0;
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int c = kh_key(g, k) & KC_MAX;
			if (c >= a->min && c <= a->max)
				++size;
		}
	}
	f = cnt_h_init();
	cnt_h_resize(f, size);
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & KC_MAX;
			if (c >= a->min && c <= a->max)
				cnt_h_put(f, kh_key(g, k), &absent);
		}
	}
	cnt_h_destroy(g);
	h->h[i] = f;
}

static void chx_shrink(cnthashx_t *h, int min, int max, int n_thread)
{
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_shrink, &a, 1<<h->p);
}

/*******************
 * Main interfaces *
 *******************/

void *um_didx_gen(const mm_idx_t *mi, int k, int pre, uint64_t mini_batch_size, int n_thread)
{
	int i;
	pldat_t pl;
	if (pre < KC_BITS) pre = KC_BITS;
	pl.mi = mi;
	pl.k = k;
	pl.n_thread = n_thread;
	pl.h = chx_init(pre);
	pl.mini_batch_size = mini_batch_size;
	kt_pipeline(2, worker_pipeline, &pl, 3);
	for (i = 0; i < 1<<pre; ++i)
		pl.h->n_tot += kh_size(pl.h->h[i]);
	chx_shrink(pl.h, 2, KC_MAX, n_thread);
	for (i = 0; i < 1<<pre; ++i)
		pl.h->n_dup += kh_size(pl.h->h[i]);
	return pl.h;
}

void um_didx_destroy(void *dh)
{
	int i;
	cnthashx_t *h = (cnthashx_t*)dh;
	for (i = 0; i < 1<<h->p; ++i) cnt_h_destroy(h->h[i]);
	free(h->h); free(h);
}
