//  Enjoy the Classic Self-Initializing Quadratic Sieve (SIQS) written in C,
//  released "as it" into the public domain, without any warranty, express or implied.

int quadraticSieve(state * state, int bits) {
	// The state contains a number N from the factorization manager and the Quadratic Sieve must :
	// - use resources present in the state (temporary variables + parameters)
	// - register prime (as possible) factors of N using the state
	// - maintain the value of N consistent with the removed factors
	// - return a non-zero value if the factorization has progressed

	if (bits < 65)
		return 0; // For every additional 10 bits, the factorization duration roughly doubles.

	QSSheet qs = {0};
	initializeState(&qs, state, bits);
	adjustInputSize(&qs);
	selectMultiplier(&qs);
	parametrize(&qs);
	allocateMemory(&qs);
	generateFactorBase(&qs);
	setupPolynomialParameters(&qs);
	do {
		do {
			// Keep randomly trying various polynomials.
			getStartedIteration(&qs);
			iterationPart1(&qs, &qs.poly.D, &qs.poly.A);
			iterationPart2(&qs, &qs.poly.A, &qs.poly.B);
			iterationPart3(&qs, &qs.poly.A, &qs.poly.B);
			for (uint32 i = 0, addi, *corr; i < qs.poly.grayMax && qs.nBits != 1; ++i, ++qs.poly.curves) {
				addi = iterationPart4(&qs, i, &corr, &qs.poly.B);
				iterationPart5(&qs, &qs.constants.kN, &qs.poly.B);
				iterationPart6(&qs, &qs.constants.kN, &qs.poly.A, &qs.poly.B, &qs.poly.C);
				iterationPart7(&qs, addi, corr);
				iterationPart8(&qs, addi, corr);
				registerRelations(&qs, &qs.poly.A, &qs.poly.B, &qs.poly.C);
			}
		} while (innerContinuationCondition(&qs));
		// Analyzes all observations made by the algorithm.
		factorizeUsingNullVectors(&qs, blockLanczos(&qs));
	} while (outerContinuationCondition(&qs));
	const int res = processRemainingFactors(&qs);
	free(qs.mem.base);
	return res;
}

// Quadratic sieve main condition 1 (often) :
// - Returns 1 : to produce more polynomials and relations.
// - Returns 0 : to start linear algebra with Block Lanczos.
int innerContinuationCondition(QSSheet *qs) {
	int answer = 1 ;
	if ((qs->time.end || 2 < qs->state->params.verbose) && qs->time.tick % 16 == 1)
		qs->time.now = getTime();
	answer &= qs->nBits != 1 ; // the bit count of N may have changed.
	answer &= (qs->relations.length.peak = qs->relations.length.now) < qs->relations.length.needs; // the condition.
	answer &= !qs->time.end || qs->time.tick % 16 == 1 || qs->time.now < qs->time.end ;
	answer &= qs->time.tick != qs->state->params.QStick_end ;
	const double pct = .01 * (10000 * qs->relations.length.now / qs->relations.length.needs) ;
	if (pct != qs->time.prev_pct) {
		if (qs->state->params.tty && qs->state->params.verbose)
			displayProgress("Quadratic Sieve", pct); // progress isn't linear.
		else if ((uint32)qs->time.prev_pct != (uint32)pct)
			LOGGER(4, "%s", answer && !(qs->time.prev_pct < 50. && 50. <= pct) ? "." : ".\n");
		qs->time.prev_pct = pct;
	}
	return answer;
}

// Quadratic sieve main condition 2. Block Lanczos linear algebra has been completed :
// - Returns 1 : to get a new attempt at linear algebra with more relations.
// - Returns 0 : when N = 1 (meaning it is fully factored) or to give up.
int outerContinuationCondition(QSSheet *qs) {
	int answer = 1 ;
	answer &= qs->nBits != 1 ; // the bit count of N isn't 1.
	answer &= qs->sieveAgainPerms-- != 0 ; // avoid infinite loop.
	answer &= qs->divisors.totalPrimes == 0 ; // search prime factors.
	answer &= !qs->time.end || getTime() < qs->time.end ;
	answer &= qs->time.tick != qs->state->params.QStick_end ;
	if (answer) {
		uint32 new_needs = qs->relations.length.needs;
		new_needs += new_needs >> (2 + qs->sieveAgainPerms);
		LOGGER(3, "[x] Maintenance re-evaluates the needs for %u additional relations.\n", new_needs - qs->relations.length.needs);
		qs->relations.length.needs = new_needs ;
	}
	return answer;
}

void parametrize(QSSheet *qs) {

	const uint32 bits = qs->knBits; // N adjusted has at least 115-bit.
	qs->knBits = (uint32) cint_count_bits(qs->state->session.tmp); // kN may be slight larger.

	LOGGER(4, "N is a %u-bit number, and kN is a %u-bit number using %u words.\n", (uint32) cint_count_bits(&qs->state->session.num), qs->knBits, (unsigned)(qs->state->session.tmp->end - qs->state->session.tmp->mem));

	uint64 tmp ;
	// params as { bits, value } take the extremal value if bits exceed.
	static const double param_base_size [][2]= { {135, 1300}, {165, 4200}, {200, 10000}, {260, 20000}, {330, 55000}, {0} };
	qs->base.length = (tmp = qs->state->params.QSbase_size) ? tmp : linearParamResolution(param_base_size, bits);

	static const double param_laziness [][2]= {{150, 95}, {220, 100}, {0} };
	// collecting more/fewer relations than recommended (in percentage).
	qs->relations.length.needs = qs->base.length * ((tmp = qs->state->params.QSlaziness) ? tmp : linearParamResolution(param_laziness, bits)) / 100;
	LOGGER(4, "The algorithm use the seed %" PRIu64 " and targets %u relations.\n", qs->state->params.rand.custom, qs->relations.length.needs);

	static const double param_m_value [][2]= { {120, 1}, {330, 6}, {0} };
	qs->m.length = (qs->m.lengthHalf = (qs->state->params.QSsieve ? qs->state->params.QSsieve : 31744) * linearParamResolution(param_m_value, bits)) << 1;

	qs->m.cacheSize = 95232 ; // algorithm reaches "M length" by steps of "cache size".

	static const double param_error [][2]= { {120, 15}, {330, 35}, {0} };
	// Logarithms of primes are rounded and errors accumulate; this specifies the magnitude of the error.
	qs->errorBits = (tmp = qs->state->params.QSerrorBits) ? tmp : linearParamResolution(param_error, bits);

	static const double param_threshold [][2]= { {120, 60}, {330, 110}, {0} };
	// The threshold that the sieve value must exceed to be considered smooth.
	qs->threshold.value = (tmp = qs->state->params.QSthreshold) ? tmp : linearParamResolution(param_threshold, bits);

	// A good multiplier reduces memory usage up to twice.
	static const double param_alloc [][2]= { {130, 992}, {140, 1280}, {150, 2176}, {160, 3584}, {170, 7168}, {180, 12288}, {190, 14336}, {190, 14336}, {200, 24576}, {210, 30720}, {220, 40960}, {230, 49152}, {240, 57344}, {250, 67584}, {260, 81920}, {270, 98304}, {280, 114688}, {290, 122880}, {300, 139264}, {310, 163840}, {320, 196608}, {330, 229376}, {0} };
	qs->mem.bytesAllocated = (tmp = qs->state->params.QSalloc_mb) ? tmp << 20 : linearParamResolution(param_alloc, qs->knBits) << 10;

	qs->sieveAgainPerms = 3; // Sieve again up to 3 times before giving up

	// Iterative list
	qs->iterativeList[0] = 1; // one
	static const double param_first_prime [][2] = { {170, 8}, {210, 12}, {320, 40}, {0} };
	qs->iterativeList[1] = linearParamResolution(param_first_prime, bits); // first
	tmp = qs->state->params.QSsieve_cutoff ? qs->state->params.QSsieve_cutoff : 5120 ;
	const uint32 large_base = tmp < qs->base.length ? tmp : qs->base.length;
	qs->iterativeList[2] = large_base >> 2; // medium
	qs->iterativeList[3] = large_base >> 1; // mid
	qs->iterativeList[4] = large_base; // sec

	LOGGER(4, "The iterative list contains %u, %u, %u and %u.\n", qs->iterativeList[1], qs->iterativeList[2], qs->iterativeList[3], qs->iterativeList[4]);

	const uint64 last_prime_in_base = (uint64) (qs->base.length * 2.5 * log(qs->base.length));
	qs->relations.tooLargePrime = (tmp = qs->state->params.QSlarge_prime) ? tmp : last_prime_in_base << 4;

	LOGGER(4, "The single large-prime variation is being processed under %" PRIu64 ".\n", qs->relations.tooLargePrime);

	qs->s.values.doubleValue = (qs->s.values.defined = (qs->s.values.subtractOne = bits / 28) + 1) << 1;
	qs->poly.grayMax = 1 << (qs->s.values.defined - 3); // computing the roots of f(x) once for all these polynomials.

	LOGGER(4, "Other params include sieve=%u, errorBits=%u, threshold=%u and s=%u.\n", qs->m.lengthHalf, qs->errorBits, qs->threshold.value, qs->s.values.defined);

	// The algorithm itself completes its configuration during the last preparation part.
	assert(qs->s.values.defined >= 3);

}

// Quadratic sieve utility function for parameter extrapolation.
uint32 linearParamResolution(const double v[][2], const uint32 point) {
	uint32 res, i, j ;
	if (v[1][0] == 0)
		res = (uint32) v[0][1];
	else {
		for (j = 1; v[j + 1][0] && point > v[j][0]; ++j);
		i = j - 1;
		if (v[i][0] > point) res = (uint32) v[i][1];
		else if (v[j][0] < point) res = (uint32) v[j][1];
		else {
			const double a = (v[j][1] - v[i][1]) / (v[j][0] - v[i][0]);
			const double b = v[i][1] - a * v[i][0];
			res = (uint32) (a * point + b);
		}
	}
	return res + (res > 512) * (16 - res % 16) ;
}

// Quadratic sieve source (algorithm)
void initializeState(QSSheet *qs, state *state, int bits) {
	// initializing (until kN is computed) with manager resources.
	qs->state = state;
	LOGGER(4, "\nQuadratic Sieve on %s.\n", cintString(state, &state->session.num));
	qs->sheet = state->session.sheet;
	qs->seed = state->session.seed;
	qs->nBits = qs->knBits = bits;
	if (2 < state->params.verbose || state->params.timeout) {
		qs->time.start = getTime() ;
		if (state->params.timeout)
			qs->time.end = qs->time.start + 1000 * qs->state->params.timeout;
	}
}

void adjustInputSize(QSSheet *qs) {
	// The algorithm is suitable for numbers larger than 115-bit,
	// and may adjust kN by a prime number to reach this size.
	cint * N = &qs->state->session.num, * kN = qs->state->session.tmp, *ADJUSTOR = kN + 1 ;
	static const unsigned char prime_generator[] = {
			9, 7, 5, 3, 17, 27, 3, 1, 29, 3, 21, 7, 17, 15,
			9, 43, 35, 15, 29, 3, 11, 3, 11, 15, 17, 25, 53,
			31, 9, 7, 23, 15, 27, 15, 29, 7, 59, 15, 5, 21,
			69, 55, 21, 21, 5, 159, 3, 81, 9, 69, 131, 33, 15 };
	const uint32 bits = (uint32) qs->nBits;
	if (bits < 115) {
		qs->adjustor = (1LLU << (124 - bits)) + prime_generator[115 - bits] ;
		intToCint(ADJUSTOR, qs->adjustor);
		cint_mul(N, ADJUSTOR, kN);
		qs->knBits = (uint32) cint_count_bits(kN);
		LOGGER(4, "Input (%u bits) is multiplied by %" PRIu64 " to reach %u bits.\n", bits, qs->adjustor, qs->knBits);
	} else
		qs->adjustor = 1, cint_dup(kN, N);
}

void selectMultiplier(QSSheet *qs) {
	// Frequently select a small multiplier (under 8-bit) that will save time and memory.
	// After it, the algorithm will factor kN instead of N, where k is a constant named "multiplier".
	const uint32 mul = (uint32) qs->state->params.QSmultiplier ;
	if (mul){
		LOGGER(4, "The multiplier is %u.\n", mul);
		qs->multiplier = mul ;
	} else {
		const size_t total_best = 7;
		uint32 best[total_best];
		for (int i = qs->state->params.verbose < 2; i < 2; ++i) {
			if (i)
				scoreDefaultMultipliers(qs, best, total_best);
			else
				scoreAlternativeMultipliers(qs, best, total_best);
			LOGGER(4, "%s", "Suggested multipliers are [");
			for (size_t j = 0; j < total_best - 1; ++j)
				LOGGER(4, "%u, ", best[j]);
			LOGGER(4, "%u]%s", best[total_best - 1], i ? "" : ".\n");
		}
		qs->multiplier = *best;
		LOGGER(4, ", so use %u.\n", *best);
	}
	if (1 < qs->multiplier) {
		cint *kN = qs->state->session.tmp, *MUL = kN + 1, *N = kN + 2 ;
		intToCint(MUL, qs->multiplier);
		cint_dup(N, kN);
		cint_mul(MUL, N, kN);
	}
}

void scoreDefaultMultipliers(QSSheet *qs, uint32 *caller_res, const size_t caller_res_len) {
	// Choose a multiplier that make the input more favorable for smoothness
	// over the future factor base, and lead to faster relation gathering.
	struct {
		uint32 mul;
		double score;
	} res[128] ; // 127 is the greatest multiplier.

	cint *N = qs->state->session.tmp, *PRIME = N + 1, *Q = N + 2, *R = N + 3;
	const double log_2 = 0.6931471805599453;
	const size_t len = sizeof(res) / sizeof(*res) - 1;
	for (uint32 i = 0; i < len ; ++i) {
		res[i].mul = i + 1;
		res[i].score = -0.5 * log(res[i].mul);
		switch (*N->mem * res[i].mul % 8) {
			// Special case against 2, the first prime number.
			case 3 : case 7 : res[i].score += 0.5 * log_2; break;
			case 5 : res[i].score += 1.0 * log_2; break;
			case 1 : res[i].score += 3.0 * log_2; break;
		}
	}

	const int limit = qs->nBits * qs->nBits >> 5 ;
	for (uint32 prime = 3; prime < limit; prime += 2)
		if (isTinyPrime(prime)) {
			// Normal case against the odd primes.
			intToCint(PRIME, prime);
			cint_div(qs->sheet, N, PRIME, Q, R);
			const uint32 n_mod_prime = (uint32) cintToInt(R);
			const double intake = 2.0 / (prime - 1) * log(prime);
			const int kronecker = kroneckerSymbol(n_mod_prime, prime);
			for (uint32 i = 0; i < len; ++i)
				if (kronecker * kroneckerSymbol(res[i].mul, prime) == 1)
					res[i].score += intake;
		}

	// Sort the results.
	for (int i = 0; i < len; ++i)
		for (int j = 1 + i; j < len; ++j)
			if (res[i].score < res[j].score)
				res[len] = res[i], res[i] = res[j], res[j] = res[len];

	for(int i = 0; i < caller_res_len; ++i)
		caller_res[i] = res[i].mul ;
}

void scoreAlternativeMultipliers(QSSheet *qs, uint32 *caller_res, const size_t caller_res_len) {
	// Choose a multiplier that make the input more favorable for smoothness
	// over the future factor base, and lead to faster relation gathering.
	struct {
		uint32 mul;
		double score;
	} res[] = {{1, 0}, {2, 0}, {3, 0}, {5, 0}, {6, 0}, {7, 0}, {10, 0}, {11, 0}, {13, 0}, {14, 0}, {15, 0}, {17, 0}, {19, 0}, {21, 0}, {22, 0}, {23, 0}, {26, 0}, {29, 0}, {30, 0}, {31, 0}, {33, 0}, {34, 0}, {35, 0}, {37, 0}, {38, 0}, {39, 0}, {41, 0}, {42, 0}, {43, 0}, {46, 0}, {47, 0}, {51, 0}, {53, 0}, {55, 0}, {57, 0}, {58, 0}, {59, 0}, {61, 0}, {62, 0}, {65, 0}, {66, 0}, {67, 0}, {69, 0}, {70, 0}, {71, 0}, {73, 0}, {79, 0}, {83, 0}, {89, 0}, {97, 0}, {101, 0}, {103, 0}, {107, 0}, {109, 0}, {0, 0}};

	cint *N = qs->state->session.tmp, *TMP = N + 1, *Q = N + 2, *R = N + 3;
	const double log_2 = 0.6931471805599453;
	const size_t len = sizeof(res) / sizeof(*res) - 1;
	for (uint32 i = 0; i < len; ++i) {
		res[i].score = -0.5 * log(res[i].mul);
		switch (*N->mem * res[i].mul % 8) {
			// Special case against 2, the first prime number.
			case 3 : case 7 : res[i].score += 0.5 * log_2; break;
			case 5 : res[i].score += 1.0 * log_2; break;
			case 1 : res[i].score += 3.0 * log_2; break;
		}
	}

	for (uint32 prime = 3; prime < 504; prime += 2)
		if (isTinyPrime(prime)) {
			// Normal case against the odd primes.
			intToCint(TMP, prime);
			cint_div(qs->sheet, N, TMP, Q, R);
			const uint32 n_mod_prime = (uint32) cintToInt(R);
			const double intake = log(prime) / (prime - 1);
			for (uint32 i = 0; i < len; ++i) {
				const uint32 kn_mod_prime = n_mod_prime * res[i].mul % prime;
				if (kn_mod_prime == 0)
					res[i].score += intake;
				else if (kroneckerSymbol(kn_mod_prime, prime) == 1)
					res[i].score += 2.0 * intake;
			}
		}

	// Sort the results.
	for (int i = 0; i < len; ++i)
		for (int j = 1 + i; j < len; ++j)
			if (res[i].score < res[j].score)
				res[len] = res[i], res[i] = res[j], res[j] = res[len];

	for(int i = 0; i < caller_res_len; ++i)
		caller_res[i] = res[i].mul ;
}

void allocateMemory(QSSheet *qs) {
	void *mem;
	mem = qs->mem.base = calloc(1, qs->mem.bytesAllocated);
	assert(mem);

	// kN has been prepared in the manager memory, now the QS has been parameterized and allocated.
	const size_t kn_size = qs->state->session.tmp[0].end - qs->state->session.tmp[0].mem + 1 ;
	// the Quadratic Sieve variables can store at most kN ^ 2 in terms of bits
	const size_t vars_size = kn_size << 1 ;

	LOGGER(4, "The underlying calculations use %u-bit variables.\n", (unsigned)(vars_size * cint_exponent));

	const size_t buffers_size = qs->base.length + (qs->iterativeList[1] << 1);
	// more relation pointers than "guessed" are available (sieve again feature).
	const size_t relations_size = (qs->base.length < qs->relations.length.needs ? qs->relations.length.needs : qs->base.length) * 7 / 4 ;

	{
		// list of the numbers used by the algorithm
		cint * const n[] = {
				&qs->vars.n,
				// polynomial
				&qs->poly.A,
				&qs->poly.B,
				&qs->poly.C,
				&qs->poly.D,
				// temporary
				&qs->vars.temp[0], &qs->vars.temp[1], &qs->vars.temp[2], &qs->vars.temp[3], &qs->vars.temp[4],
				// relations finder
				&qs->vars.x,
				&qs->vars.key,
				&qs->vars.value,
				&qs->vars.cycle,
				// a factor of N
				&qs->vars.factor,
				// constants
				&qs->constants.kN,
				&qs->constants.one,
				&qs->constants.mHalf,
				&qs->constants.tooLargePrime,
				&qs->constants.multiplier,
				0,
		};
		for (int i = 0; n[i]; ++i) {
			n[i]->mem = n[i]->end = memAligned(mem) ;
			mem = n[i]->mem + (n[i]->size = vars_size);
		}
	}

	cint_dup(&qs->vars.n, &qs->state->session.num);
	cint_dup(&qs->constants.kN, qs->state->session.tmp);

	intToCint(&qs->constants.one, 1);
	intToCint(&qs->constants.mHalf, qs->m.lengthHalf);
	intToCint(&qs->constants.tooLargePrime, qs->relations.tooLargePrime);
	intToCint(&qs->constants.multiplier, qs->multiplier);

	// Allocates "s" rows
	qs->s.data = memAligned(mem);
	mem = memAligned(qs->s.data + qs->s.values.defined);
	for (uint32 i = 0; i < qs->s.values.defined; ++i) {
		inlineCint(&qs->s.data[i].bTerms, kn_size, &mem); // also "s" more cint
		qs->s.data[i].aInvDoubleValueBTerms = mem;
		mem = memAligned(qs->s.data[i].aInvDoubleValueBTerms + qs->base.length);
	}
	qs->s.aIndexes = memAligned(mem); // the indexes of the prime numbers that compose A.

	// Allocates "base length" rows
	qs->base.data = memAligned(qs->s.aIndexes + qs->s.values.doubleValue);
	qs->m.positions[0] = memAligned(qs->base.data + qs->base.length);
	qs->m.positions[1] = memAligned(qs->m.positions[0] + qs->base.length);
	qs->m.sieve = memAligned(qs->m.positions[1] + qs->base.length);
	qs->m.sieve[qs->m.length] = 0xFF ; // the end of the sieve evaluates to "true" under any "truthy" mask.
	qs->m.flags = memAligned(qs->m.sieve + qs->m.length + sizeof(uint64_t));
	// Usage: buffer[0] is zeroed after use, buffer[1] isn't supposed zeroed after use.
	qs->buffer[0] = memAligned(qs->m.flags + qs->base.length);
	qs->buffer[1] = memAligned(qs->buffer[0] + buffers_size);

	// Other allocations
	qs->relations.length.capacity = (uint32) relations_size ;
	// Block Lanczos has a part of memory, it takes a "lite" snapshot before throwing relations.
	qs->lanczos.snapshot = memAligned(qs->buffer[1] + buffers_size) ;
	qs->relations.data = memAligned(qs->lanczos.snapshot + relations_size);
	qs->divisors.data = memAligned(qs->relations.data + relations_size);
	// A lot of divisors isn't needed because the algorithm calculate their GCD to reduce N.
	qs->mem.now = memAligned(qs->divisors.data + (qs->nBits * qs->nBits >> 8));

	const uint32 n_trees = (uint32) (sizeof(qs->uniqueness) / sizeof(struct avl_manager));
	for (uint32 i = 0; i < n_trees; ++i) {
		// the trees are used to identify duplicates (relations, partials, divisors of N)
		qs->uniqueness[i].inserter_argument = &qs->mem.now;
		qs->uniqueness[i].inserter = &avl_cint_inserter;
		qs->uniqueness[i].comparator = (int (*)(const void *, const void *)) &h_cint_compare;
		// they use default sign-less comparator.
	}
	avl_at(&qs->uniqueness[2], &qs->constants.one); // Ignore 1 as a divisor of N.
	LOGGER(4, "Allocated %u MB of memory with a %u KB structure.\n", qs->mem.bytesAllocated >> 20, (unsigned)((char*)qs->mem.now - (char*)qs->mem.base) >> 10);
}

void generateFactorBase(QSSheet *qs) {
	// Prepare the factor base (a set of small prime numbers used to find smooth numbers).
	static const double inv_ln_2 = 1.4426950408889634;
	cint *A = qs->vars.temp, *B = A + 1, *C = A + 2;
	uint32 i = 0, prime;

	// the factor base contain the multiplier if different from 2.
	if (qs->multiplier != 2)
		qs->base.data[i].size = (uint32) (.35 + inv_ln_2 * log(qs->base.data[i].num = qs->multiplier)), ++i;

	// then it contains the number 2.
	qs->base.data[i].num = 2, qs->base.data[i].size = 1;
	qs->base.data[i].sqrtKNModPrime = *qs->constants.kN.mem % 8 == 1 || *qs->constants.kN.mem % 8 == 7, ++i;

	// then the prime numbers for which kN is a quadratic residue modulo.
	for (prime = 3; i < qs->base.length; prime += 2)
		if (isTinyPrime(prime)) {
			intToCint(A, prime);
			cint_div(qs->sheet, &qs->constants.kN, A, B, C);
			const uint32 kn_mod_prime = (uint32) cintToInt(C);
			qs->base.data[i].sqrtKNModPrime = tonelliShanks(kn_mod_prime, prime);
			if (qs->base.data[i].sqrtKNModPrime) {
				qs->base.data[i].num = prime;
				qs->base.data[i].size = (uint32) (.35 + inv_ln_2 * log(prime)), ++i;
			}
		}
	// 2.5 * (base size) * ln(base size) is close to the largest prime number in factor base.
	qs->base.largest = qs->base.data[i - 1].num ;

	LOGGER(4, "The factor base of %u suitable primes ends with %u.\n", qs->base.length, qs->base.largest);
}

void setupPolynomialParameters(QSSheet *qs) {
	// completes the configuration by the algorithm itself.
	// computes D : a template (optimal value of hypercube) for the A polynomial coefficient.
	uint32 i, min, span;
	const uint32 s = qs->s.values.defined ;
	qs->poly.span.half = (span = qs->base.length / (s * (s + s))) >> 1;
	cint *kN = qs->vars.temp, *TMP = kN + 1, *R = kN + 2;
	cint_dup(kN, &qs->constants.kN);
	cint_left_shifti(kN, 1);
	cint_sqrt(qs->sheet, kN, TMP, R);
	cint_div(qs->sheet, TMP, &qs->constants.mHalf, &qs->poly.D, R);
	qs->poly.dBits = (uint32) cint_count_bits(&qs->poly.D);
	cint_nthRoot(qs->sheet, &qs->poly.D, s, R); // use the s-th root of D.
	const uint32 root = (uint32) cintToInt(R) ;
	for (i = 1; qs->base.data[i].num <= root; ++i, assert(i < qs->base.length));
	if (i < span) {
		LOGGER(3, "[x] Maintenance adjusts the span value from %u to %u.\n", span, i);
		span = i; // Avoids a rare case of failure when factoring small numbers (graceful degradation).
	}
	assert(i >= span);
	for (min = i - qs->poly.span.half, i *= i; i / min < span + min; --min);
	qs->poly.span.x1 = min ;
	qs->poly.span.x2 = min + qs->poly.span.half ;
	qs->poly.span.x3 = qs->poly.span.x2 * qs->poly.span.x2 ;
	assert(qs->poly.span.x2 < qs->base.length);
}

void getStartedIteration(QSSheet *qs) {
	if (qs->lanczos.snapshot[0].relation) {
		// the operation is fast, it shouldn't happen in average case.
		// it restores the relations reduced by the linear algebra step that failed.
		uint32 i ;
		for(i = 0; qs->lanczos.snapshot[i].relation; ++i) {
			qs->relations.data[i] = qs->lanczos.snapshot[i].relation;
			qs->relations.data[i]->Y.length = qs->lanczos.snapshot[i].yLength;
			qs->lanczos.snapshot[i].relation = 0 ;
		}
		assert(qs->relations.length.prev == i) ;
		qs->relations.length.now = i;
		LOGGER(4, "[x] Maintenance restores the relations to a size of %u.\n", i);
	}
	//  Increase the tick value in each iteration of the algorithm.
	if (++qs->time.tick % 32 == 0) {
		if (qs->relations.length.prev == qs->relations.length.now) {
			// The algorithm notices that no relations accumulates, and reacts to unblock the situation.
			LOGGER(3, "[x] Maintenance randomizes D because the relation counter remains at %u.\n", qs->relations.length.now);
			cint_random_bits(&qs->poly.D, qs->poly.dBits, qs->seed);
			*qs->poly.D.mem |= 1; // Shouldn't happen, D becomes a randomized odd number.
		}
		qs->relations.length.prev = qs->relations.length.now;
	}
}

void iterationPart1(QSSheet * qs, const cint * D, cint * A) {
	uint32 n_tries = 0 ; // several attempts may rarely be necessary.
	retry:;
	// A is a "random" product of "s" distinct prime numbers from the factor base.
	cint * X = qs->vars.temp, * Y = X + 1, *TMP, *_A = A ;
	uint32 a, b, i = 0, j;
	if (qs->s.values.defined & 1) TMP = A, A = X, X = TMP ;
	// swap pointers so the last multiplication completes inside the A variable.
	intToCint(A, 1);
	for (a = 0, b = qs->poly.span.x3; a < qs->s.values.subtractOne; ++a) {
				if (a & 1) i = b / (i + qs->poly.span.x1) - (uint32) xorRandint(qs->seed, 0, 9);		else i = qs->poly.span.x2 + (uint32) xorRandint(qs->seed, 0, qs->poly.span.half);		for (j = 0; j < a; j = i == qs->s.data[j].primeIndex ? ++i, 0 : j + 1);		qs->s.data[a].primeIndex = i; // the selected divisor of A wasn't already present in the product.
		intToCint(Y, qs->base.data[i].num);
		cint_mul(A, Y, X), TMP = A, A = X, X = TMP;
	}
	// a prime number from the factor base completes A, which must be close to D.
	cint_div(qs->sheet, D, A, X, Y);
	const uint32 d_over_a = (uint32) cintToInt(X);
	for (i = qs->base.data[0].num != 2 ; qs->base.data[i].num <= d_over_a; ++i);
	for (j = 0; j < qs->s.values.subtractOne; j = i == qs->s.data[j].primeIndex ? ++i, 0 : j + 1);
	if (qs->base.length <= i) {
		const char *ord = n_tries == 0 ? "st" : n_tries == 1 ? "nd" : n_tries == 2 ? "rd" : "th" ;
		LOGGER(3, "[x] Maintenance discards A=%s on %u%s attempt.\n", cintString(qs->state, A), n_tries + 1, ord);
		assert(++n_tries <= 16); // clearly shouldn't happen 16 times, review the algorithm parameters otherwise.
		A = _A ;
		goto retry;
	}
	qs->s.data[qs->s.values.subtractOne].primeIndex = i ;
	intToCint(Y, qs->base.data[i].num);
	cint_mul(A, Y, X); // generated A values should always be distinct, A no longer change.
	assert(X == &qs->poly.A);
}

void iterationPart2(QSSheet * qs, const cint * A, cint * B) {
	cint *X = qs->vars.temp, *PRIME = X + 1, *Y = X + 2, *R = X + 3;
	uint32 i, *pen = qs->s.aIndexes;
	cint_erase(B);
	for (i = 0; i < qs->s.values.defined; ++i) {
		const uint32 idx = qs->s.data[i].primeIndex, prime = qs->base.data[idx].num;
		if (idx >= qs->iterativeList[3])
			qs->iterativeList[3] = 8 + idx - idx % 8 ;
		// write [index of prime number, power] of the A factors into buffer.
		*pen++ = qs->s.data[i].primeIndex, *pen++ = 1;
		qs->s.data[i].primeSquared = (uint64)prime * (uint64)prime ;
		intToCint(PRIME, prime);
		cint_div(qs->sheet, A, PRIME, X, R), assert(R->mem == R->end); // div exact.
		cint_div(qs->sheet, X, PRIME, Y, R);
		qs->s.data[i].aOverPrimeModPrime = (uint32) cintToInt(R);
		uint64 x = modularInverse(qs->s.data[i].aOverPrimeModPrime, prime);
		x = x * qs->base.data[qs->s.data[i].primeIndex].sqrtKNModPrime % prime;
		intToCint(X, x > prime >> 1 ? prime - x : x);
		cint_mul(A, X, Y);
		cint_div(qs->sheet, Y, PRIME, &qs->s.data[i].bTerms, R), assert(R->mem == R->end); // div exact.
		cint_addi(B, &qs->s.data[i].bTerms);
	}
}

void iterationPart3(QSSheet * qs, const cint * A, const cint * B) {
	cint *Q = qs->vars.temp, *R = Q + 1, *PRIME = Q + 2;
	uint64 i, j, x, y;
	for (i = 0; i < qs->base.length; ++i) {
		// prepare the "roots" and "aInvDoubleValueBTerms". The algorithm will
		// fill 2 ** (s - 3) sieves by using these values and adding "prime sizes".
		const uint32 prime = qs->base.data[i].num;
		intToCint(PRIME, prime);
		cint_div(qs->sheet, A, PRIME, Q, R);
		const uint32 a_mod_prime = (uint32) cintToInt(R) ;
		cint_div(qs->sheet, B, PRIME, Q, R) ;
		const uint32 b_mod_prime = (uint32) cintToInt(R) ;
		const uint32 a_inv_doubleValue = modularInverse(a_mod_prime, prime) << 1 ;
		// Arithmetic shifts "<<" and ">>" performs multiplication or division by powers of two.
		x = y = prime;
		x += qs->base.data[i].sqrtKNModPrime;
		y -= qs->base.data[i].sqrtKNModPrime;
		x -= b_mod_prime;
		x *= a_inv_doubleValue >> 1;
		y *= a_inv_doubleValue ;
		x += qs->m.lengthHalf ;
		x %= prime ;
		y += x ;
		y %= prime ;
		qs->base.data[i].root[0] = (uint32) x ; // First root of the polynomial mod prime.
		qs->base.data[i].root[1] = (uint32) y ; // Second root of the polynomial mod prime.
		for (j = 0; j < qs->s.values.defined; ++j) {
			// compute the roots update value for all "s".
			cint_div(qs->sheet, &qs->s.data[j].bTerms, PRIME, Q, R);
			const uint64 b_term = cintToInt(R);
			qs->s.data[j].aInvDoubleValueBTerms[i] = (uint32)(a_inv_doubleValue * b_term % prime);
		}
	}
	// The next function operates over "bTerms" multiplied by 2.
	for (i = 0; i < qs->s.values.defined; cint_left_shifti(&qs->s.data[i++].bTerms, 1));
}

uint32 iterationPart4(const QSSheet * qs, const uint32 nth_curve, uint32 ** corr, cint *B) {
	uint32 i, gray_act; // the Gray code in "nth_curve" indicates which "bTerms" to consider.
	for (i = 0; nth_curve >> i & 1; ++i);
	if (gray_act = (nth_curve >> i & 2) != 0, gray_act)
		cint_addi(B, &qs->s.data[i].bTerms) ;
	else // and which action to perform.
		cint_subi(B, &qs->s.data[i].bTerms) ;
	*corr = qs->s.data[i].aInvDoubleValueBTerms;
	return gray_act; // B values generated here should always be distinct.
}

void iterationPart5(QSSheet *  qs, const cint * kN, const cint * B) {
	cint *P = qs->vars.temp, *Q = P + 1, *R_kN = P + 2, *R_B = P + 3, *TMP = P + 4;
	for (uint32 a = 0; a < qs->s.values.defined; ++a) {
		const uint32 i = qs->s.data[a].primeIndex;
		const int64 prime = qs->base.data[i].num ;
		intToCint(P, qs->s.data[a].primeSquared);
		cint_div(qs->sheet, B, P, Q, R_B);
		cint_div(qs->sheet, kN, P, Q, R_kN);
		if (B->nat < 0) cint_addi(R_B, P); // if B is negative.
		const int64 rem_b = (int64) cintToInt(R_B);
		const int64 rem_kn = (int64) cintToInt(R_kN);
		int64 s ; // both remainders are modulo the prime number squared.
		if (rem_b < 0xb504f334) {
			// the multiplication is straightforward.
			s = rem_b * rem_b - rem_kn;
			s /= prime ;
		} else {
			// the common multiplication would overflow.
			cint_mul(R_B, R_B, TMP);
			cint_subi(TMP, R_kN);
			intToCint(P, (uint64) prime);
			cint_div(qs->sheet, TMP, P, Q, R_B);
			s = (int64) cintToInt(Q);
			if (Q->nat < 0) s = -s ;
		}
		//
		int64 bezout = (rem_b % prime) * (int64) qs->s.data[a].aOverPrimeModPrime ;
		bezout = (int64) modularInverse((uint32) (bezout % prime), (uint32) prime);
		//
		s = (int64) qs->m.lengthHalf - s * bezout ;
		s %= prime ;
		s += (s < 0) * prime ;
		qs->base.data[i].root[0] = (uint32) s;
		qs->base.data[i].root[1] = (uint32) -1; // Zero out roots corresponding to the factors of A.
	}
}

void iterationPart6(QSSheet *qs, const cint *kN, const cint *A, const cint *B, cint *C) {
	cint *TMP = qs->vars.temp, *R = TMP + 1;
	cint_mul(B, B, TMP); // (B * B) % A = kN % A
	cint_subi(TMP, kN); // C = (B * B - kN) / A
	cint_div(qs->sheet, TMP, A, C, R), assert(R->mem == R->end); // div exact.
}

void iterationPart7(QSSheet * qs, const uint32 gray_addi, const uint32 * restrict corr) {
	// Sieve for larger prime numbers.
	memset(qs->m.sieve, 0, qs->m.length * sizeof(*qs->m.sieve));
	memset(qs->m.flags, 0, qs->base.length * sizeof(*qs->m.flags));
	uint8 * restrict end = qs->m.sieve + qs->m.length, *p_0, *p_1;
	for(uint32 i = qs->iterativeList[3], j = qs->iterativeList[4]; i < j; ++i) {
		const uint32 prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i];
		qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
		qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
		p_0 = qs->m.sieve + qs->base.data[i].root[0];
		p_1 = qs->m.sieve + qs->base.data[i].root[1];
		for (; end > p_0 && end > p_1;)
			*p_0 += size, p_0 += prime, *p_1 += size, p_1 += prime;
		*p_0 += (end > p_0) * size, *p_1 += (end > p_1) * size;
	}
	for(uint32 i = qs->iterativeList[4], j = qs->base.length; i < j; ++i){
		const uint32 prime = qs->base.data[i].num, size = qs->base.data[i].size, co = gray_addi ? prime - corr[i] : corr[i] ;
		qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
		qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
		for(p_0 = qs->m.sieve + qs->base.data[i].root[0]; end > p_0; qs->m.flags[i] |= 1 << ((p_0 - qs->m.sieve) & 7), *p_0 += size, p_0 += prime);
		for(p_1 = qs->m.sieve + qs->base.data[i].root[1]; end > p_1; qs->m.flags[i] |= 1 << ((p_1 - qs->m.sieve) & 7), *p_1 += size, p_1 += prime);
	}
}

void iterationPart8(QSSheet * qs, const uint32 gray_addi, const uint32 *  corr) {
	// Sieving means taking an interval [−M/2, +M/2] and determining for
	// which X in [−M/2, +M/2] a given prime number divides AX^2 + 2BX + C.
	uint8 * chunk_begin = qs->m.sieve, *chunk_end = chunk_begin;
	uint8 * sieve_end = chunk_begin + qs->m.length ;
	uint32 *buffer = qs->buffer[0], walk_idx, * walk = buffer;
	// Since the previous function, the check is performed for the prime numbers of the factor base.
	for(uint32 i = 0; i < qs->iterativeList[3]; ++i)
		if (qs->base.data[i].root[1] != (uint32) -1) {
			*walk++ = i ; // the current prime number isn't a factor of A.
			const uint32 prime = qs->base.data[i].num, co = gray_addi ? prime - corr[i] : corr[i] ;
			qs->base.data[i].root[0] += co; if (qs->base.data[i].root[0] >= prime) qs->base.data[i].root[0] -= prime;
			qs->base.data[i].root[1] += co; if (qs->base.data[i].root[1] >= prime) qs->base.data[i].root[1] -= prime;
			qs->m.positions[0][i] = chunk_begin + qs->base.data[i].root[0];
			qs->m.positions[1][i] = chunk_begin + qs->base.data[i].root[1];
		}
	for(walk_idx = 0; buffer[walk_idx] < qs->iterativeList[1]; ++walk_idx);
	do{ // iterates step by step until the entire sieve is filled.
		walk = buffer + walk_idx ;
		chunk_end = chunk_end + qs->m.cacheSize < sieve_end ? chunk_end + qs->m.cacheSize : sieve_end;
		do{
			const uint32 size = qs->base.data[*walk].size, prime = qs->base.data[*walk].num, times = 4 >> (*walk > qs->iterativeList[2]) ;
			uint8 ** const p_0 = qs->m.positions[0] + *walk, ** const p_1 = qs->m.positions[1] + *walk;
			const int64 diff = *p_1 - *p_0 ;
			for(const uint8 * const bound = chunk_end - prime * times; bound > *p_0;)
				for(uint32 i = 0; i < times; ++i)
					**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			for(; *p_0 < chunk_end && *p_0 + diff < chunk_end;)
				**p_0 += size, *(*p_0 + diff) += size, *p_0 += prime;
			*p_1 = *p_0 + diff ;
			if (*p_0 < chunk_end) **p_0 += size, *p_0 += prime;
			if (*p_1 < chunk_end) **p_1 += size, *p_1 += prime;
		} while(*++walk);
	} while(chunk_begin = chunk_end, chunk_begin < sieve_end);
	memset(qs->buffer[0], 0, (walk - qs->buffer[0]) * sizeof(*walk));
}

cint * divisorsUniquenessHelper(QSSheet * qs, const cint * num) {
	// Helper for uniqueness within the divisors of N.
	struct avl_node *node;
	node = avl_at(&qs->uniqueness[2], num) ;
	return qs->uniqueness[2].affected ? node->key : 0 ;
}

int registerDivisor(QSSheet *qs) {
	// Register a divisor of N, combine them with GCD and identify the perfect powers.
	// Returns 0 when the factorization is completed, 1 otherwise.
#define IN_RANGE(F) (h_cint_compare(&qs->constants.one, F) < 0 && h_cint_compare(F, &qs->vars.n) < 0)
	cint *F = &qs->vars.factor, *tmp;
	F->nat = 1 ; // Absolute value.
	if (!(IN_RANGE(F) && (tmp = divisorsUniquenessHelper(qs, F))))
		return 1; // Duplicates are ignored.
	struct task {
		cint * num ;
		char origin ;
	} tasks[63] ; // Implements a stack-based recursion.
	cint *curr, **divisors = qs->divisors.data, *Q = qs->vars.temp + 3, *R = Q + 1;
	int i = 0, pow;
	tasks[i++] = (struct task) {tmp, 0};
	LOGGER(4, "- New divisor %s shown.\n", cintString(qs->state, tmp));
	do {
		curr = tasks[--i].num; // Retrieve the top element.
		if (cint_is_prime(qs->sheet, curr, -1, qs->seed)) {
			pow = (int) cint_remove(qs->sheet, &qs->vars.n, curr);
			assert(pow); // Prime factors are removed from N.
			++qs->divisors.totalPrimes;
			qs->nBits = (uint32) cint_count_bits(&qs->vars.n);
			// Register this prime factor in the manager's routine.
			manager_add_factor(qs->state, curr, pow, 1);
			if (tasks[i].origin) {
				char * msg = 0 ; // Explain succinctly how this prime factor was found.
				switch(tasks[i].origin) {
					case 1 : msg = "And allows us for N"; break;
					case 2 : msg = "Prunes the divisors"; break;
					case 3 : msg = "Divides N"; break;
					case 4 : msg = "Notes a perfect power"; break;
					case 5 : msg = "Performs GCD within the divisors"; break;
				}
				LOGGER(4, "%*s- %s to get %s.\n", (i + 1) << 1, "", msg, cintString(qs->state, curr));
			}
			if (qs->nBits != 1) {
				LOGGER(4, "%*s- This prime factor reduces N to %d-bit.\n", (i + 1) << 1, "", qs->nBits);
				if ((tmp = divisorsUniquenessHelper(qs, &qs->vars.n)))
					tasks[i++] = (struct task){tmp, 1}; // 1. And allows us.
				for (uint32 j = 0; j < qs->divisors.length; ++j) {
					cint_dup(F, divisors[j]);
					pow = cint_remove(qs->sheet, F, curr);
					if (pow) {
						divisors[j--] = divisors[--qs->divisors.length];
						if ((tmp = divisorsUniquenessHelper(qs, F)) && IN_RANGE(tmp))
							tasks[i++] = (struct task){tmp, 2}; // 2. Prunes the divisors.
					}
				}
			} else
				LOGGER(4, "%*s- The factorization is complete since it's a prime.\n", (i + 1) << 1, "");
		} else {
			cint_div(qs->sheet, &qs->vars.n, curr, Q, R) ;
			if (R->mem == R->end && IN_RANGE(Q) && (tmp = divisorsUniquenessHelper(qs, Q)))
				tasks[i++] = (struct task){tmp, 3}; // 3. Divides N.
			pow = anyRootCheck(qs->state, curr, Q, R) ;
			if (pow && IN_RANGE(Q) && (tmp = divisorsUniquenessHelper(qs, Q)))
				tasks[i++] = (struct task){tmp, 4}; // 4. Notes a perfect power.
			for (uint32 j = 0; j < qs->divisors.length; ++j) {
				cint_gcd(qs->sheet, curr, divisors[j], Q);
				if (IN_RANGE(Q) && (tmp = divisorsUniquenessHelper(qs, Q)))
					tasks[i++] = (struct task){tmp, 5}; // 5. Performs GCD within the divisors.
			}
			divisors[qs->divisors.length++] = curr;
		}
	} while (i && qs->nBits != 1);
	return qs->nBits != 1;
#undef IN_RANGE
}

void registerRelations(QSSheet * qs, const cint * A, const cint * B, const cint * C) {
	cint *  TMP = qs->vars.temp, * K = &qs->vars.key, * V = &qs->vars.value ;
	uint32 m_idx, idx, mod;
	// iterates the values of X in [-M/2, +M/2].
	for (m_idx = 0; m_idx < qs->m.length; ++m_idx)
		if (qs->m.sieve[m_idx] >= qs->threshold.value) {
			// over the threshold, compute f(X) and check candidate for smoothness.
			intToCint(&qs->vars.x, m_idx);
			cint_subi(&qs->vars.x, &qs->constants.mHalf); // X = "current index" - M/2
			cint_mul(A, &qs->vars.x, TMP); // TMP = AX
			cint_addi(TMP, B); // TMP = AX + B
			cint_dup(K, TMP); // Key = copy of AX + B
			cint_addi(TMP, B); // TMP = AX + 2B
			cint_mul(TMP, &qs->vars.x, V); // V = AX^2 + 2BX
			cint_addi(V, C); // Value = f(X) = AX^2 + 2BX + C
			// We can inject X in the equation A * C + kN = B * B
			// So it should hold that (A * X + B)^2 - kN = A * f(X)
			V->nat = 1 ; // absolute value
			uint32 target_bits = (uint32) cint_count_bits(V) - qs->errorBits;
			uint32 removedBits = 0, * restrict pen = qs->buffer[1];
			//  writes the pairs [index of the prime number, power found in V].
			if (qs->base.data[0].num != 1) {
				intToCint(TMP, qs->base.data[0].num);
				*pen++ = 0; // remove powers of the multiplier.
				*pen = (uint32) cint_remove(qs->sheet, V, TMP);
				if (*pen) removedBits += *pen++ * qs->base.data[0].size; else --pen;
			}
			for (idx = 1; idx < qs->iterativeList[1]; ++idx)
				if (qs->base.data[idx].root[1] == (uint32) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
					intToCint(TMP, qs->base.data[idx].num);
					// for a given prime number of the factor base, "remove" returns
					// the numbers of powers that was present in V, and V is updated.
					*pen++ = idx;
					*pen = (uint32) cint_remove(qs->sheet, V, TMP);
					if (*pen) removedBits += *pen++ * qs->base.data[idx].size; else --pen;
				}
			if (removedBits + qs->m.sieve[m_idx] >= target_bits) {
				// there is a chance to register a new relation.
				for (removedBits = 0, target_bits = qs->m.sieve[m_idx]; idx < qs->iterativeList[4] && removedBits < target_bits; ++idx)
					if (qs->base.data[idx].root[1] == (uint32) -1 || (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1])) {
						intToCint(TMP, qs->base.data[idx].num);
						*pen++ = idx;
						*pen = (uint32) cint_remove(qs->sheet, V, TMP);
						if (*pen) removedBits += *pen++ * qs->base.data[idx].size; else --pen;
					}
				for (const uint8 mask = 1 << (m_idx & 7); idx < qs->base.length && removedBits < target_bits; ++idx)
					if (qs->m.flags[idx] & mask)
						if (mod = m_idx % qs->base.data[idx].num, mod == qs->base.data[idx].root[0] || mod == qs->base.data[idx].root[1]) {
							intToCint(TMP, qs->base.data[idx].num);
							*pen++ = idx;
							*pen = (uint32) cint_remove(qs->sheet, V, TMP);
							if (*pen) removedBits += *pen++ * qs->base.data[idx].size; else --pen;
						}
				const uint32 * restrict const primeIndexes_and_powers[4] = {
						qs->s.aIndexes, // really factoring A * f(X), commit outstanding A factors.
						qs->s.aIndexes + qs->s.values.doubleValue,
						qs->buffer[1],
						pen,
				};
				if (h_cint_compare(V, &qs->constants.one) == 0)
					registerRegularRelation(qs, K, primeIndexes_and_powers);
				else if (155 < qs->knBits && h_cint_compare(V, &qs->constants.tooLargePrime) < 0)
					//  Store it until another partial share the same variation (also called large prime, cofactor).
					registerPartialRelation(qs, K, V, primeIndexes_and_powers);
			}
		}
}

void registerRegularRelation(QSSheet * qs, const cint * KEY, const uint32 * const restrict args[4]) {
	struct avl_node *node = avl_at(&qs->uniqueness[0], KEY);
	if (node->value)
		return; // duplicates at this stage are ignored.
	struct QSRelation * rel = qs->mem.now;
	uint32 i, j ;
	qs->mem.now = rel + 1 ; // a relation must be swappable for Block Lanczos reducing.
	rel->X = node->key; // constant X is stored by the node key.
	rel->Y.data = qs->mem.now; // Y data length only decreases.
	const size_t yLength = (args[1] - args[0] + args[3] - args[2]) >> 1 ;
	rel->axis.Z.data = rel->Y.data + yLength; // writes Z ahead.
	for (i = 0; i < 4;) {
		// processes the given column arrays.
		const uint32 * restrict idx = args[i++], * restrict const end_index = args[i++];
		for (; idx < end_index; idx += 2) {
			const uint32 power = *(idx + 1) ;
			if (power & 1) {
				// remove from Y the indexes of the prime numbers that are already listed (factors of A).
				for (j = 0; j < rel->Y.length && rel->Y.data[j] != *idx; ++j);
				if (j == rel->Y.length) // add, the index wasn't present.
					rel->Y.data[rel->Y.length++] = *idx;
				else // or remove.
					memmove(rel->Y.data + j, rel->Y.data + j + 1, (--rel->Y.length - j) * sizeof(*rel->Y.data));
			}
			for (j = 0; j < power; ++j)
				rel->axis.Z.data[rel->axis.Z.length++] = *idx;
		}
	}
	qs->mem.now = rel->axis.Z.data + rel->axis.Z.length; // Z length is known.
	int verified = 0 ;
	if (rel->Y.length > qs->s.values.defined) {
		// it often passes.
		cint *A = qs->vars.temp, *B = A + 1;
		intToCint(A, 1);
		for (j = 0; j < rel->axis.Z.length; ++j) {
			intToCint(B, qs->base.data[rel->axis.Z.data[j]].num);
			cint_mulModi(qs->sheet, A, B, &qs->constants.kN);
		}
		cint_mulMod(qs->sheet, rel->X, rel->X, &qs->constants.kN, B);
		verified = !cint_compare(A, B) || (cint_addi(A, B), !cint_compare(A, &qs->constants.kN));
	}
	if (verified){
		node->value = qs->relations.data[qs->relations.length.now] = rel;
		qs->mem.now = rel->axis.Z.data + rel->axis.Z.length;
		rel->id = ++qs->relations.length.now; // Keep the relation.
	} else {
		LOGGER(3, "[x] Maintenance discards the relation at index %u.\n", qs->relations.length.now);
		char * open = (char*) rel, * close = qs->mem.now ;
		qs->mem.now = memset(open, 0, close - open); // Throw.
	}
}

void registerPartialRelation(QSSheet * qs, const cint * KEY, const cint * VALUE, const uint32 * const restrict args[4]) {
	// Process the single large-prime variation.
	// Searches 2 different KEY sharing the same VALUE.
	struct avl_node *node = avl_at(&qs->uniqueness[1], VALUE);
	struct QSRelation *old, *new;
	cint * BEZOUT = 0;
	old = node->value;
	if (old) {
		if (old->X == 0) return; // the value is already marked as "ignored".
		if (old->axis.next) return; // accepting all "next" without caring compromise the linear algebra.
		for (; old && h_cint_compare(KEY, old->X); old = old->axis.next);
		if (old) return; // same KEY already registered.
		old = node->value;
		if (old->axis.next == 0) {
			// there is an "old" using the same VALUE, and it has no "next" yet.
			cint *A = qs->vars.temp, *B = A + 1;
			if (qs->multiplier != 1)
				if (cint_gcd(qs->sheet, VALUE, &qs->constants.multiplier, A), cint_compare_char(A, 1)){
					old->X = 0; // VALUE shouldn't be related so close to the multiplier.
					return;
				}
			// so compute BEZOUT.
			cint_modularInverse(qs->sheet, VALUE, &qs->constants.kN, A);
			if (A->mem == A->end) {
				old->X = 0; // no solution to the linear congruence.
				cint_gcd(qs->sheet, VALUE, &qs->constants.kN, &qs->vars.factor);
				cint_div(qs->sheet, &qs->vars.n, &qs->vars.factor, A, B);
				if (B->mem == B->end) // found a small factor of N ?
					registerDivisor(qs);
				return; // nothing.
			} else
				BEZOUT = A;
		}
	}

	new = memAligned(qs->mem.now);
	qs->mem.now = new + 1;
	new->X = qs->mem.now;

	if (BEZOUT) {
		// BEZOUT is stored directly after the new X, like in an array.
		qs->mem.now = new->X + 2;
		dupCint(new->X, KEY, &qs->mem.now);
		dupCint(new->X + 1, BEZOUT, &qs->mem.now);
		// The 2nd newcomer become the root of the linked list.
		new->axis.next = old, node->value = new = old = new;
	} else {
		// All but the 2nd have no special treatment.
		qs->mem.now = new->X + 1; // they come at the end of the linked list.
		dupCint(new->X, KEY, &qs->mem.now);
		if (old) {
			for (; old->axis.next; old = old->axis.next);
			old->axis.next = new, old = node->value;
		} else node->value = new;
	}

	// data buffered isn't persistent, it may be needed, so it's copied.
	uint32 * data = new->Y.data = memAligned(qs->mem.now);
	new->Y.length = (uint32) (args[1] - args[0]);
	memcpy(data, args[0], new->Y.length * sizeof(*data));
	memcpy(data + new->Y.length, args[2], (args[3] - args[2]) * sizeof(*data));
	new->Y.length += (uint32) (args[3] - args[2]);
	qs->mem.now = new->Y.data + new->Y.length;

	if (old) {
		BEZOUT = old->X + 1 ; // the modular inverse was stored here.
		cint_mulMod(qs->sheet, new->X, BEZOUT, &qs->constants.kN, &qs->vars.cycle);
		do {
			if (old != new) {
				// combines, it registers a smooth relation using the 2 buffers.
				cint_mulMod(qs->sheet, &qs->vars.cycle, old->X, &qs->constants.kN, &qs->vars.key);
				uint32 * restrict begin = qs->buffer[0], * restrict pen = begin;
				data = memset(qs->buffer[1], 0, qs->base.length * sizeof(*data));
				for (uint32 i = 0; i < new->Y.length; i += 2)
					data[new->Y.data[i]] += new->Y.data[i + 1];
				for (uint32 i = 0; i < old->Y.length; i += 2)
					data[old->Y.data[i]] += old->Y.data[i + 1];
				for (uint32 i = 0; i < qs->base.length; ++i)
					if (data[i]) // writes [index of the prime number, power]
						*pen++ = i, *pen++ = data[i];
				args = (const uint32 * restrict const[4]){ begin, pen, 0, 0, };
				registerRegularRelation(qs, &qs->vars.key, args);
				++qs->relations.length.byPartial;
				memset(begin, 0, (char *) pen - (char *) begin); // zeroed.
			}
		} while ((old = old->axis.next));
		// the linked list can handle 2+ entries, but their more complex combinations isn't implemented.
	}
}

void factorizeUsingNullVectors(QSSheet * qs, const uint64_t * restrict const lanczos_answer) {
	// Block Lanczos linear algebra answer is simply "mask followed by null_rows", with read-only.
	if (*lanczos_answer == 0)
		return;
	const uint64_t mask = *lanczos_answer, * restrict null_rows = lanczos_answer + 1;
	cint *X = &qs->vars.x, *Y = qs->vars.temp, *TMP = Y + 1, *POW = Y + 2;
	uint32 * restrict power_of_primes;
	for(uint32 row = 0; row < 64 && qs->nBits != 1; ++row)
		if (mask >> row & 1){
			// use the Fermat's (1607 - 1665) method to compute a factorization of N.
			intToCint(X, 1), intToCint(TMP, 1), intToCint(Y, 1);
			power_of_primes = memset(qs->buffer[1], 0, qs->base.length * sizeof(*power_of_primes));
			for (uint32 i = 0; i < qs->relations.length.now; ++i)
				if (null_rows[i] >> row & 1) {
					const struct QSRelation * restrict const rel = qs->relations.data[i];
					// The algorithm must retrieve the X and Z relation fields
					// related to the Y field initially submitted to Block Lanczos.
					cint_mulModi(qs->sheet, X, rel->X, &qs->vars.n);
					for (uint32 j = 0; j < rel->axis.Z.length; ++j)
						++power_of_primes[rel->axis.Z.data[j]];
				}
			for (uint32 i = 0; i < qs->base.length; ++i)
				if (power_of_primes[i]){
					// powers are even ... square root ...
					intToCint(TMP, qs->base.data[i].num);
					intToCint(POW, power_of_primes[i] >> 1);
					cint_powModi(qs->sheet, TMP, POW, &qs->vars.n);
					cint_mulModi(qs->sheet, Y, TMP, &qs->vars.n);
				}
			h_cint_subi(Y, X);
			if (Y->mem != Y->end) {
				cint_gcd(qs->sheet, &qs->vars.n, Y, &qs->vars.factor);
				// 100 digits RSA number has been factored by the software in 2022.
				if (registerDivisor(qs) == 0)
					break;
			}
		}
}

int processRemainingFactors(QSSheet *qs) {
	if (qs->nBits != 1 && qs->divisors.length) {
		// In rare cases N must be partially factored.
		// Registers a divisor encountered by the algorithm in the manager routine.
		cint **divisors = qs->divisors.data, *tmp;
		for (uint32 i = 0; i < qs->divisors.length; ++i) {
			for (uint32 j = 1 + i; j < qs->divisors.length; ++j)
				if (h_cint_compare(divisors[j], divisors[i]) < 0)
					tmp = divisors[i], divisors[i] = divisors[j], divisors[j] = tmp; // Apply a selection sort.
			if (h_cint_compare(&qs->constants.one, divisors[i]) < 0 && h_cint_compare(divisors[i], &qs->vars.n) < 0) {
				const int power = (int) cint_remove(qs->sheet, &qs->vars.n, divisors[i]);
				if (power) {
					LOGGER(4, "Quadratic Sieve submits the composite divisor %s as a result.\n", cintString(qs->state, divisors[i]));
					manager_add_factor(qs->state, divisors[i], power, -1); // -1 marks this divisor as composite for recursive factorization.
					qs->nBits = (uint32) cint_count_bits(&qs->vars.n);
				}
			} else
				break; // No need to sort more.
		}
		// Recursively running the Quadratic Sieve will allow faster factorization than this oversized instance.
		LOGGER(3, "[x] Maintenance is to forward the remainder %s.\n", cintString(qs->state, &qs->vars.n));
	}

	if (qs->uniqueness[1].count)
		LOGGER(4, "The sieve reported %u partials which added %u smooth relations.\n", (unsigned) qs->uniqueness[1].count, qs->relations.length.byPartial);
	LOGGER(4, "The algorithm completed with %u polynomials and %u relations.\n", qs->poly.curves, qs->relations.length.peak);
	LOGGER(4, "Used %u MB of memory during %.02f second(s).\n", (unsigned) ((char *) qs->mem.now - (char *) qs->mem.base) >> 20, 0.001 * (getTime() - qs->time.start));

	// Tells the factorization manager whether its task has progressed using the Quadratic Sieve.
	const int has_fully_factored = qs->nBits == 1, has_partially_factored = qs->divisors.length != 0 ;
	const int res = has_fully_factored || has_partially_factored ;
	if (res) // Updates the input number accordingly, since the algorithm worked on a copy of it.
		cint_dup(&qs->state->session.num, &qs->vars.n);

	return res ;
}