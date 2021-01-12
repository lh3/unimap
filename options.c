#include <stdio.h>
#include <limits.h>
#include "umpriv.h"

void mm_idxopt_init(mm_idxopt_t *opt)
{
	memset(opt, 0, sizeof(mm_idxopt_t));
	opt->k = 21, opt->w = 11, opt->flag = 0;
	opt->bucket_bits = 14;
	opt->bf_bits = 35;
	opt->mini_batch_size = 1000000000;
	opt->adap_occ = 50;
}

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->seed = 11;
	opt->mid_occ = 50;
	opt->mid_occ_cap = 250;
	opt->mid_occ_frac = 2e-4f;

	opt->min_cnt = 3;
	opt->min_chain_score = 40;
	opt->bw = 100000;
	opt->max_gap = 5000;
	opt->max_gap_ref = -1;
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->adap_dist = -1;
	opt->chain_gap_scale = 1.0f;
	opt->rmq_inner_dist = 1000;
	opt->rmq_size_cap = 100000;

	opt->end_len_frac = 0.15f;
	opt->gap_flank_frac = 0.5f;

	opt->mask_level = 0.5f;
	opt->mask_len = INT_MAX;
	opt->pri_ratio = 0.8f;
	opt->best_n = 10;

	opt->alt_drop = 0.15f;

	opt->a = 1, opt->b = 3, opt->q = 5, opt->e = 2, opt->q2 = 25, opt->e2 = 1;
	opt->sc_ambi = 1;
	opt->zdrop = 800, opt->zdrop_inv = 200;
	opt->end_bonus = -1;
	opt->min_dp_max = 200;
	opt->min_ksw_len = 200;
	opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
	opt->max_clip_ratio = 1.0f;
	opt->mini_batch_size = 1000000000;
	opt->max_sw_mat = 100000000;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	int mid_occ;
	if ((opt->flag & MM_F_SPLICE_FOR) || (opt->flag & MM_F_SPLICE_REV))
		opt->flag |= MM_F_SPLICE;
	mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
	if (mid_occ > opt->mid_occ_cap) mid_occ = opt->mid_occ_cap;
	if (opt->mid_occ < mid_occ) opt->mid_occ = mid_occ;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), opt->mid_occ);
}

void mm_mapopt_max_intron_len(mm_mapopt_t *opt, int max_intron_len)
{
	if ((opt->flag & MM_F_SPLICE) && max_intron_len > 0)
		opt->max_gap_ref = opt->bw = max_intron_len;
}

int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	if (preset == 0) {
		mm_idxopt_init(io);
		mm_mapopt_init(mo);
	} else if (strcmp(preset, "clr") == 0 || strcmp(preset, "ont") == 0) {
		mo->a = 2, mo->b = 4, mo->q = 4, mo->e = 2, mo->q2 = 24, mo->e2 = 1, mo->zdrop = 400, mo->zdrop_inv = 200;
		mo->min_dp_max = 80;
		mo->bw = 10000;
		mo->best_n = 5;
		mo->flag |= MM_F_NO_RMQ;
		mo->mini_batch_size = 500000000;
	} else if (strcmp(preset, "hifi") == 0 || strcmp(preset, "ccs") == 0) {
		mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1, mo->zdrop = 800, mo->zdrop_inv = 200;
		mo->min_dp_max = 200;
		mo->bw = 10000;
		mo->best_n = 5;
		mo->flag |= MM_F_NO_RMQ;
		mo->mini_batch_size = 500000000;
	} else if (strcmp(preset, "asm5") == 0) {
		mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = 800, mo->zdrop_inv = 200;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "asm10") == 0) {
		mo->a = 1, mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = 800, mo->zdrop_inv = 200;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strcmp(preset, "asm20") == 0) {
		mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1, mo->zdrop = 800, mo->zdrop_inv = 200;
		mo->min_dp_max = 200;
		mo->best_n = 50;
	} else if (strncmp(preset, "splice", 6) == 0 || strcmp(preset, "cdna") == 0) {
		io->k = 15, io->w = 5;
		mo->flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK | MM_F_NO_RMQ;
		mo->max_gap = 2000, mo->max_gap_ref = mo->bw = 200000;
		mo->a = 1, mo->b = 2, mo->q = 2, mo->e = 1, mo->q2 = 32, mo->e2 = 0;
		mo->noncan = 9;
		mo->junc_bonus = 9;
		mo->zdrop = 200, mo->zdrop_inv = 100; // because mo->a is halved
		if (strcmp(preset, "splice:hq") == 0)
			mo->junc_bonus = 5, mo->b = 4, mo->q = 6, mo->q2 = 24;
	} else return -1;
	return 0;
}

int mm_check_opt(const mm_idxopt_t *io, const mm_mapopt_t *mo)
{
	if (io->k <= 0 || io->w <= 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -k and -w must be positive\033[0m\n");
		return -5;
	}
	if (mo->best_n < 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -N must be no less than 0\033[0m\n");
		return -4;
	}
	if (mo->best_n == 0 && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m '-N 0' reduces mapping accuracy. Please use '--secondary=no' instead.\033[0m\n");
	if (mo->pri_ratio < 0.0f || mo->pri_ratio > 1.0f) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -p must be within 0 and 1 (including 0 and 1)\033[0m\n");
		return -4;
	}
	if ((mo->flag & MM_F_FOR_ONLY) && (mo->flag & MM_F_REV_ONLY)) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m --for-only and --rev-only can't be applied at the same time\033[0m\n");
		return -3;
	}
	if (mo->e <= 0 || mo->q <= 0) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -O and -E must be positive\033[0m\n");
		return -1;
	}
	if ((mo->q != mo->q2 || mo->e != mo->e2) && !(mo->e > mo->e2 && mo->q + mo->e < mo->q2 + mo->e2)) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m dual gap penalties violating E1>E2 and O1+E1<O2+E2\033[0m\n");
		return -2;
	}
	if ((mo->q + mo->e) + (mo->q2 + mo->e2) > 127) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m scoring system violating ({-O}+{-E})+({-O2}+{-E2}) <= 127\033[0m\n");
		return -1;
	}
	if (mo->zdrop < mo->zdrop_inv) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m Z-drop should not be less than inversion-Z-drop\033[0m\n");
		return -5;
	}
	if ((mo->flag & MM_F_NO_PRINT_2ND) && (mo->flag & MM_F_ALL_CHAINS)) {
		if (mm_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m -X/-P and --secondary=no can't be applied at the same time\033[0m\n");
		return -5;
	}
	return 0;
}
