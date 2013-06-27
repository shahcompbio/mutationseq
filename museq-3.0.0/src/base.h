#ifndef _base_h
#define _base_h

typedef uint8_t base4_t;
typedef uint8_t base2_t;

typedef uint16_t ngram16_t;
typedef uint32_t ngram32_t;

char base4_char(base4_t b);

base2_t char_base2(char base);
base2_t base4_base2(base4_t base);

void ngram32_push(ngram32_t *g, base2_t b);
void ngram16_push(ngram16_t *g, base2_t b);

ngram32_t str2ngram(char *s);
char* ngram16_str(ngram16_t *g);
char* ngram32_str(ngram32_t *g);

void trace_ngram(ngram16_t *g);
#endif
