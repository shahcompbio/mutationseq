#include <stdint.h>
#include "debug.h"
#include "base.h"

/* the nucleotide datatype */

char
base2_char(base2_t base) {
    char bases[4] = "ACGT";
    return bases[base];
}
base2_t
base4_base2(base4_t base) {
    switch(base) {
        case 0x1:
            return 0x0;
        case 0x2:
            return 0x1;
        case 0x4:
            return 0x2;
        case 0x8:
            return 0x3;
    }
    return 0x4;
}

base2_t
char_base2(char base) {
    switch(base) {
        case 'a':
        case 'A':
            return 0x0;
        case 'c':
        case 'C':
            return 0x1;
        case 'g':
        case 'G':
            return 0x2;
        case 't':
        case 'T':
            return 0x3;
    }
    return 0x4;
}

char
base4_char(base4_t b)
{
    switch(b) {
        case 0x01:
            return 'A';
        case 0x02:
            return 'C';
        case 0x04:
            return 'G';
        case 0x08:
            return 'T';
        case 0x0F:
            return 'N';
    }
    return '.';
}

void ngram16_push(ngram16_t *g, base2_t b) {
    *g = *g << 2;
    *g = *g | b;
}


void ngram32_push(ngram32_t *g, base2_t b) {
    *g = *g << 2;
    *g = *g | b;
}

void trace_ngram(ngram16_t *g) {
    int i;
    char n[8];
    for (i = 0; i < 8; i++) {
        n[i] = base2_char(g[i]);
    }
    trace("%s", n);
}
