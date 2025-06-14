// Block Lanczos algorithm released "as it" into the public domain, without any warranty, express or implied.
// Results of operations are last AND non-const arguments.

// Kornél Lánczos (1893 - 1974) was a Hungarian-American and later Hungarian-Irish mathematician and physicist.

void lanczos_mul_MxN_Nx64(const QSSheet *qs, const uint64_t *X, const uint32 max_size, uint64_t *Y) {
	assert(Y != X);
	memset(Y, 0, max_size * sizeof(uint64_t));
	for (uint32 a = 0, b; a < qs->relations.length.now; ++a) {
		struct QSRelation *const rel = qs->relations.data[a];
		for (b = 0; b < rel->Y.length; ++b)
			Y[rel->Y.data[b]] ^= X[a];
	}
}

void lanczos_mul_trans_MxN_Nx64(const QSSheet *qs, const uint64_t *X, uint64_t *Y) {
	assert(Y != X);
	for (uint32 a = 0, b; a < qs->relations.length.now; ++a) {
		struct QSRelation *const rel = qs->relations.data[a];
		for (Y[a] = 0, b = 0; b < rel->Y.length; ++b)
			Y[a] ^= X[rel->Y.data[b]];
	}
}

void lanczos_mul_64xN_Nx64(const QSSheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	assert(X != Z && Y != Z);
	uint32 a, b, c, d;
	memset(Z, 0, 256 * 8 * sizeof(*Z));
	memset(T, 0, 64 * sizeof(*T));
	for (a = 0; a < qs->relations.length.now; ++a) {
		const uint64_t tmp = X[a]; // read while writing ?!
		for (b = 0, c = 0; c < 64; c += 8, b += 256)
			Z[b + (tmp >> c & 0xff)] ^= Y[a];
	}
	for (a = 0; a < 8; ++a, ++T) {
		uint64_t tmp[8] = {0};
		for (b = 0; b < 256; ++b)
			if (b >> a & 1)
				for (c = d = 0; c < 8; ++c, d += 256)
					tmp[c] ^= Z[b + d];
		for (b = 0, c = 0; b < 8; ++b, c += 8)
			T[c] = tmp[b];
	}
}

uint64_t lanczos_find_non_singular_sub(const uint64_t *t, const uint64_t *last_s, uint64_t *s, uint64_t last_dim, uint64_t *w) {
	uint64_t i, j, dim, cols[64];
	uint64_t M[64][2], mask, *row_i, *row_j, m_0, m_1;
	for (i = 0; i < 64; ++i)
		M[i][0] = t[i], M[i][1] = (uint64_t) 1 << i;
	mask = 0;
	for (i = 0; i < last_dim; ++i)
		mask |= (uint64_t) 1 << (cols[63 - i] = last_s[i]);
	for (i = j = 0; i < 64; ++i)
		if (!(mask & ((uint64_t) 1 << i)))
			cols[j++] = i;
	for (i = dim = 0; i < 64; ++i) {
		mask = (uint64_t) 1 << (cols[i]);
		row_i = M[cols[i]];
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break;
			}
		}
		if (j < 64) {
			for (j = 0; j < 64; ++j) {
				row_j = M[cols[j]];
				if (row_i != row_j && (row_j[0] & mask))
					row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
			}
			s[dim++] = cols[i];
			continue;
		}
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break;
			}
		}
		if (j == 64) return 0; // sub-matrix is not invertible

		for (j = 0; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_i != row_j && (row_j[1] & mask))
				row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
		}

		row_i[0] = row_i[1] = 0;
	}
	for (i = 0; i < 64; ++i)
		w[i] = M[i][1];
	mask = 0;
	for (i = 0; i < dim; ++i)
		mask |= (uint64_t) 1 << s[i];
	for (i = 0; i < last_dim; ++i)
		mask |= (uint64_t) 1 << last_s[i];
	dim *= mask == -(uint64_t) 1;
	return dim;
}

void lanczos_mul_Nx64_64x64_acc(QSSheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	uint32 a;
	uint64_t b, c, d, e;
	for (b = 0; b < 8; Y += 8, Z += 256, ++b)
		for (c = 0; c < 256; ++c)
			for (d = Z[c] = 0, e = c; e; e >>= 1, ++d)
				Z[c] ^= (e & 1) * Y[d];
	for (a = 0, Z -= 2048; a < qs->relations.length.now; ++a)
		for (c = d = 0; c < 64; c += 8, d += 256) {
			const uint64_t w = X[a];
			T[a] ^= Z[d + (w >> c & 0xff)];
		}
}

void lanczos_mul_64x64_64x64(const uint64_t *X, const uint64_t *Y, uint64_t *Z) {
	uint64_t a, b, c, d, tmp[64];
	for (a = 0; a < 64; tmp[a++] = c) {
		for (b = 0, c = 0, d = X[a]; d; d >>= 1, ++b)
			c ^= (d & 1) * Y[b];
	}
	memcpy(Z, tmp, sizeof(tmp));
}

void lanczos_transpose_vector(QSSheet *qs, const uint64_t *X, uint64_t **Y) {
	uint32 a; // Z will be zeroed automatically by the loop.
	uint64_t b, c, d, *Z;
	Z = memcpy(qs->mem.now, X, qs->relations.length.now * sizeof(*X));
	for (a = 0; a < qs->relations.length.now; ++a)
		for (b = 0, c = a >> 6, d = (uint64_t) 1 << (a % 64); Z[a]; Z[a] >>= 1, ++b)
			Y[b][c] |= (Z[a] & 1) * d;
}

void lanczos_combine_cols(QSSheet *qs, uint64_t *x, uint64_t *v, uint64_t *ax, uint64_t *av) {
	// Use Gaussian elimination on the columns of [ax | av] to find all linearly dependent columns.
	uint32 i, j, bit_pos, num_deps;
	uint64_t k, col, col_words;
	uint64_t mask, *mat_1[128], *mat_2[128], *tmp;
	num_deps = 64 << (v && av);
	col_words = (qs->relations.length.now + 63) / 64;
	for (i = 0; i < num_deps; ++i) {
		mat_1[i] = qs->mem.now;
		mat_2[i] = mat_1[i] + col_words;
		qs->mem.now = mat_2[i] + col_words;
	}
	lanczos_transpose_vector(qs, x, mat_1);
	lanczos_transpose_vector(qs, ax, mat_2);
	if (num_deps == 128) {
		lanczos_transpose_vector(qs, v, mat_1 + 64);
		lanczos_transpose_vector(qs, av, mat_2 + 64);
	}
	for (i = bit_pos = 0; i < num_deps && bit_pos < qs->relations.length.now; ++bit_pos) {
		mask = (uint64_t) 1 << (bit_pos % 64);
		col = bit_pos / 64;
		for (j = i; j < num_deps; ++j)
			if (mat_2[j][col] & mask) {
				tmp = mat_1[i];
				mat_1[i] = mat_1[j];
				mat_1[j] = tmp;
				tmp = mat_2[i];
				mat_2[i] = mat_2[j];
				mat_2[j] = tmp;
				break;
			}
		if (j == num_deps)
			continue;
		for (++j; j < num_deps; ++j)
			if (mat_2[j][col] & mask)
				for (k = 0; k < col_words; ++k) {
					mat_2[j][k] ^= mat_2[i][k];
					mat_1[j][k] ^= mat_1[i][k];
				}
		++i;
	}

	for (j = 0; j < qs->relations.length.now; ++j) {
		uint64_t word = 0;
		col = j / 64;
		mask = (uint64_t) 1 << (j % 64);
		for (k = i; k < 64; ++k)
			if (mat_1[k][col] & mask)
				word |= (uint64_t) 1 << k;
		x[j] = word;
	}
	char *open = (char *) mat_1[0], *close = qs->mem.now;
	qs->mem.now = memset(open, 0, close - open);
}

void lanczos_build_array(QSSheet *qs, uint64_t **target, const size_t quantity, const size_t size) {
	// ensure it remains memory for linear algebra
	const char *mem_needs = (char *) qs->mem.now + quantity * size * sizeof(uint64_t);
	const char *mem_ends = (char *) qs->mem.now + qs->mem.bytesAllocated;
	assert(mem_needs < mem_ends);
	for (size_t i = 0; i < quantity; ++i)
		target[i] = qs->mem.now, qs->mem.now = memAligned(target[i] + size);
}

uint64_t *lanczos_block_worker(QSSheet *qs) {
	const uint32 n_cols = qs->relations.length.now, v_size = n_cols > qs->base.length ? n_cols : qs->base.length;
	const uint64_t safe_size = qs->lanczos.safeLength;
	uint64_t *md[6], *xl[2], *sm[13], *tmp, *res, mask_0, mask_1;
	uint32 i, dim_0, dim_1;
	qs->mem.now = memAligned((uint64_t *) qs->mem.now + 1); // keep some padding.
	lanczos_build_array(qs, md, 6, safe_size);
	lanczos_build_array(qs, sm, 13, 64);
	lanczos_build_array(qs, xl, 2, safe_size < 2048 ? 2048 : safe_size);
	res = (*md) - 1; // simple "trick" to return mask + null_rows
	for (i = 0; i < 64; ++i)
		sm[12][i] = i;
	dim_0 = 0;
	dim_1 = 64;
	mask_1 = (uint64_t) -1;
	for (i = 0; i < qs->relations.length.now; ++i)
		md[1][i] = xorRandom(qs->seed);
	memcpy(md[0], md[1], v_size * sizeof(uint64_t));
	lanczos_mul_MxN_Nx64(qs, md[1], v_size, xl[1]);
	lanczos_mul_trans_MxN_Nx64(qs, xl[1], md[1]);
	memcpy(xl[0], md[1], v_size * sizeof(uint64_t));
	qs->lanczos.nIterations = 0 ;
	do {
		lanczos_mul_MxN_Nx64(qs, md[1], v_size, xl[1]);
		lanczos_mul_trans_MxN_Nx64(qs, xl[1], md[4]);
		lanczos_mul_64xN_Nx64(qs, md[1], md[4], xl[1], sm[3]);
		lanczos_mul_64xN_Nx64(qs, md[4], md[4], xl[1], sm[5]);
		for (i = 0; i < 64 && !(sm[3][i]); ++i);
		if (i != 64) {
			dim_0 = (uint32) lanczos_find_non_singular_sub(sm[3], sm[12], sm[11], dim_1, sm[0]);
			if (dim_0) {
				mask_0 = 0;
				for (i = 0; i < dim_0; ++i)
					mask_0 |= (uint64_t) 1 << sm[11][i];
				for (i = 0; i < 64; ++i)
					sm[7][i] = (sm[5][i] & mask_0) ^ sm[3][i];
				lanczos_mul_64x64_64x64(sm[0], sm[7], sm[7]);
				for (i = 0; i < 64; ++i)
					sm[7][i] ^= (uint64_t) 1 << i;
				lanczos_mul_64x64_64x64(sm[1], sm[3], sm[8]);
				for (i = 0; i < 64; ++i)
					sm[8][i] &= mask_0;
				lanczos_mul_64x64_64x64(sm[4], sm[1], sm[9]);
				for (i = 0; i < 64; ++i)
					sm[9][i] ^= (uint64_t) 1 << i;
				lanczos_mul_64x64_64x64(sm[2], sm[9], sm[9]);
				for (i = 0; i < 64; ++i)
					sm[10][i] = ((sm[6][i] & mask_1) ^ sm[4][i]) & mask_0;
				lanczos_mul_64x64_64x64(sm[9], sm[10], sm[9]);
				for (i = 0; i < qs->relations.length.now; ++i)
					md[4][i] &= mask_0;
				lanczos_mul_Nx64_64x64_acc(qs, md[1], sm[7], xl[1], md[4]);
				lanczos_mul_Nx64_64x64_acc(qs, md[3], sm[8], xl[1], md[4]);
				lanczos_mul_Nx64_64x64_acc(qs, md[2], sm[9], xl[1], md[4]);
				lanczos_mul_64xN_Nx64(qs, md[1], xl[0], xl[1], sm[7]);
				lanczos_mul_64x64_64x64(sm[0], sm[7], sm[7]);
				lanczos_mul_Nx64_64x64_acc(qs, md[1], sm[7], xl[1], md[0]);
				tmp = md[2], md[2] = md[3], md[3] = md[1], md[1] = md[4], md[4] = tmp;
				tmp = sm[2], sm[2] = sm[1], sm[1] = sm[0], sm[0] = tmp;
				tmp = sm[4], sm[4] = sm[3], sm[3] = tmp;
				tmp = sm[6], sm[6] = sm[5], sm[5] = tmp;
				memcpy(sm[12], sm[11], 64 * sizeof(uint64_t));
				mask_1 = mask_0;
				dim_1 = dim_0;
			}
		}
	} while (++qs->lanczos.nIterations < 2048 && dim_0 && i != 64);

	// it sometimes succeeds at 700+ iterations during the tests.

	// ===== answer finalization =====
	// result will be a simple array of the form [mask, null_rows...]
	// it's assumed that a null mask means "miss, no answer"

	*res = 0; // mask

	if (qs->lanczos.nIterations < 2048) {
		lanczos_mul_MxN_Nx64(qs, md[0], v_size, md[3]);
		lanczos_mul_MxN_Nx64(qs, md[1], v_size, md[2]);
		lanczos_combine_cols(qs, md[0], md[1], md[3], md[2]);
		lanczos_mul_MxN_Nx64(qs, md[0], v_size, md[1]);
		if (*md[1] == 0) // should hold (the buffer must contain only zero)
			if (memcmp(md[1], md[1] + 1, (v_size - 1) * sizeof(uint64_t)) == 0)
				for (i = 0; i < n_cols; *res |= (*md)[i++]);
	}

	// if no mask found : clears everything, otherwise leave [mask + null rows]
	char *open = (char *) md[*res != 0], *close = qs->mem.now;
	qs->mem.now = memset(open, 0, close - open);

	return res;
}

void lanczos_reduce_matrix(QSSheet *qs) {
	// a filtering is not always necessary to make "lanczos_block_worker" succeed :
	// - it writes to the relations [ Y lengths, relation counter ] will change
	uint32 a, b, c, row, col, reduced_rows = qs->base.length, *counts;
	counts = memset(qs->buffer[1], 0, qs->base.length * sizeof(*qs->buffer[1]));
	if (qs->sieveAgainPerms)
		for (a = 0; a < qs->relations.length.now; ++a) {
			// "snapshot" pointers, so they can be restored if "sieve again" is fired.
			qs->lanczos.snapshot[a].relation = qs->relations.data[a];
			qs->lanczos.snapshot[a].yLength = qs->relations.data[a]->Y.length;
		}
	for (a = 0; a < qs->relations.length.now; ++a)
		for (b = 0; b < qs->relations.data[a]->Y.length; ++b)
			++counts[qs->relations.data[a]->Y.data[b]];
	// the counter of relations will change, store the original counter for debug purpose.
	qs->relations.length.prev = qs->relations.length.now;
	//
	do {
		row = reduced_rows;
		do {
			col = qs->relations.length.now;
			for (a = b = 0; a < qs->relations.length.now; ++a) {
				struct QSRelation *restrict const rel = qs->relations.data[a];
				for (c = 0; c < rel->Y.length && counts[rel->Y.data[c]] > 1; ++c);
				if (c < rel->Y.length)
					for (; rel->Y.length;)
						--counts[rel->Y.data[--rel->Y.length]];
				else
					qs->relations.data[b++] = qs->relations.data[a]; // relation is thrown.
			}
		} while (qs->relations.length.now = b, b != col);
		for (a = reduced_rows = 0; a < qs->base.length; reduced_rows += counts[a++] != 0);
		if (b = reduced_rows + 64, qs->relations.length.now > b) { // 64 extra rows
			for (a = b; a < qs->relations.length.now; ++a)
				for (struct QSRelation *restrict const rel = qs->relations.data[a]; rel->Y.length;)
					--counts[rel->Y.data[--rel->Y.length]];
			qs->relations.length.now = b;
		}
	} while (row != reduced_rows);
	LOGGER(4, "[x] Maintenance of linear algebra reduces the relations from %u to %u.\n", qs->relations.length.prev, qs->relations.length.now);
}

uint64_t *blockLanczos(QSSheet *qs) {
	// the worker algorithm is probabilistic with high success rate.
	// it is interested in the Y field of the relations (to build its matrix).
	// it receives as input the relations (reduced or not, depending on the settings).

	if (qs->time.tick == qs->state->params.QStick_end)
		return qs->mem.now ; // used for testing during the development.

	if (qs->lanczos.safeLength < qs->relations.length.now)
		qs->lanczos.safeLength = qs->relations.length.now;
	if (qs->lanczos.safeLength < qs->base.length)
		qs->lanczos.safeLength = qs->base.length;
	qs->lanczos.safeLength += 64 - qs->lanczos.safeLength % 64;
	//
	const uint32 tries = 5, reduce_at = 3;
	uint64_t *res;
	//
	for(uint32 i = 0; i < tries; ++i) {
		if (i == reduce_at && !qs->lanczos.snapshot[0].relation)
			lanczos_reduce_matrix(qs);
		// This algorithm is likely to quickly report the success of the Quadratic Sieve,
		// just having to find in its smooth relations a subset of all the exponent vectors
		// such that the sum of their exponent vectors is the zero vector ...
		res = lanczos_block_worker(qs);
		if(*res) {
			const char *ord = i == 0 ? "st" : i == 1 ? "nd" : i == 2 ? "rd" : "th" ;
			LOGGER(4, "Linear algebra passed on %u%s attempt after %u iterations.\n", i + 1, ord, qs->lanczos.nIterations);
			break;
		}
	}
	// When N is not completely factored by this result, another call to this
	// function will provide another result which can sometimes finish the job.
	return res; // This can be due (rarely) to a less favorable PRNG sequence.
}
