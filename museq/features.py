from __future__ import division

name = "TCGA Benchmark 4 featureset with coverage info"
version = "5"

feature_set = (
("tumour_indels", lambda t, n, r: t[9] / t[5][0]),
("normal_indels", lambda t, n, r: n[9] / n[5][0]),

("tumour_ref_depth", lambda t, n, r: t[r[0] + 1][0] / t[5][0]),
("normal_ref_depth", lambda t, n, r: n[r[0] + 1][0] / n[5][0]),
("normal_mapq_ratio", lambda t, n, r: n[5][2] / n[5][0]),
("tumour_mapq_ratio", lambda t, n, r: t[5][2] / t[5][0]),
("normal_ref_quality", lambda t, n, r: n[r[0] + 1][1] / n[5][0]),
("tumour_ref_quality", lambda t, n, r: t[r[0] + 1][1] / t[5][0]),
("normal_quality_ratio", lambda t, n, r: n[5][1] / n[5][0]),
("tumour_quality_ratio", lambda t, n, r: t[5][1] / t[5][0]),

("normal_tumour_quality", lambda t, n, r: (t[5][1] / t[5][0]) / ((n[5][1] / n[5][0]) + 0.00001)),

("normal_tumour_mapq", lambda t, n, r: (t[5][2] / t[5][0]) / ((n[5][2] / n[5][0]) + 0.00001)),
("normal_tumour_ref_depth", lambda t, n, r: ((t[r[0] + 1][1] / t[5][0]) + 0.00001) / ((n[r[0] + 1][1] / n[5][0]) + 0.00001)),
("normal_tumour_direction", lambda t, n, r: (t[5][4] / t[5][0]) / ((n[5][4] / n[5][0]) + 0.00001)),
("normal_tumour_ref_direction", lambda t, n, r: (t[r[0] + 1][4] / (t[r[0] + 1][0] + 0.00001)) / ((n[r[0] + 1][4] / (n[r[0] + 1][0] + 0.00001)) + 0.00001)),
("normal_distance_ratio", lambda t, n, r: n[5][3] / n[5][0]),
("tumour_distance_ratio", lambda t, n, r: t[5][3] / t[5][0]),

("normal_direction_ratio", lambda t, n, r: n[5][4] / n[5][0]),
("tumour_direction_ratio", lambda t, n, r: t[5][4] / t[5][0]),

("normal_ref_direction_total", lambda t, n, r: n[r[0] + 1][4] / n[5][0]),
("tumour_ref_direction_total", lambda t, n, r: t[r[0] + 1][4] / t[5][0]),

("normal_ref_direction", lambda t, n, r: n[r[0] + 1][4] / (n[r[0] + 1][0] + 0.00001)),
("tumour_ref_direction", lambda t, n, r: t[r[0] + 1][4] / (t[r[0] + 1][0] + 0.00001)),


("region_entropy", lambda t, n, r: r[4]),
("region_gc_content", lambda t, n, r: r[3]),
("homopolymer_f", lambda t, n, r: r[1]),
("homopolymer_b", lambda t, n, r: r[2]),

("tumour_variant_quality_ratio", lambda t, n, r: ((t[5][1] - t[r[0] + 1][1]) / (t[5][0] - t[r[0] + 1][0] + 0.00001))),
("normal_variant_quality_ratio", lambda t, n, r: ((n[5][1] - n[r[0] + 1][1]) / (n[5][0] - n[r[0] + 1][0] + 0.00001))),
("tumour_variant_quality", lambda t, n, r: (t[5][1] - t[r[0] + 1][1])),
("normal_variant_quality", lambda t, n, r: (n[5][1] - n[r[0] + 1][1])),


#("normal_direction", lambda t, n, r: n[5][4]),
#("tumour_direction", lambda t, n, r: t[5][4]),

("normal_variant_direction_ratio", lambda t, n, r: (n[5][4] - n[r[0] + 1][4]) / (n[5][0] - n[r[0] + 1][0] + 0.00001)),
("tumour_variant_direction_ratio", lambda t, n, r: (t[5][4] - t[r[0] + 1][4]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),

#("normal_variant_direction", lambda t, n, r: n[5][4] - n[r[0] + 1][4]),
#("tumour_variant_direction", lambda t, n, r: t[5][4] - t[r[0] + 1][4]),

#("normal_tumour_variant_ratio", lambda t, n, r: (t[5][4] - t[r[0] + 1][4]) / (n[5][4] - n[r[0] + 1][4] + 0.00001)),

#("normal_tumour_variant_quality", lambda t, n, r: (t[5][1] - t[r[0] + 1][1]) / (n[5][1] - n[r[0] + 1][1] + 0.00001)),

#("normal_tumour_variant_distance", lambda t, n, r: (t[5][3] - t[r[0] + 1][3]) / (n[5][3] - n[r[0] + 1][3] + 0.00001)),

#("normal_quality", lambda t, n, r: n[5][1]),
#("tumour_quality", lambda t, n, r: t[5][1]),

("normal_tail_distance", lambda t, n, r: n[5][3] / n[5][0]),
("tumour_tail_distance", lambda t, n, r: t[5][3] / t[5][0]),


#("normal_tumour_variant_mapq", lambda t, n, r: (t[5][2] - t[r[0] + 1][2]) / (n[5][2] - n[r[0] + 1][2] + 0.00001)),
("normal_variant_distance", lambda t, n, r: (n[5][3] - n[r[0] + 1][3]) / (n[5][0] - n[r[0] + 1][0] + 0.00001)),
("tumour_variant_distance", lambda t, n, r: (t[5][3] - t[r[0] + 1][3]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),

#("normal_minor_allele", lambda t, n, r: (r[0] != t[7] and n[t[6] + 1][0] > 0) or (r[0] != t[6] and n[t[7] + 1][0] > 0)),

# Classic features

("normal_depth", lambda t, n, c: n[5][0]),
("tumour_depth", lambda t, n, c: t[5][0]),

#("normal_variant_depth", lambda t, n, r: (n[5][0] - n[r[0] + 1][0])),
#("tumour_variant_depth", lambda t, n, r: (t[5][0] - t[r[0] + 1][0])),
("normal_variant_depth_ratio", lambda t, n, r: ((n[5][0] - n[r[0] + 1][0]) / n[5][0])),
("tumour_variant_depth_ratio", lambda t, n, r: ((t[5][0] - t[r[0] + 1][0]) / t[5][0])),

("normal_tumour_depth", lambda t, n, r: (t[5][0] / n[5][0])),
#("normal_tumour_variant_depth", lambda t, n, r: (t[5][0] - t[r[0] + 1][0]) / (n[5][0] - n[r[0] + 1][0] + 0.00001)),
("normal_tumour_variant_depth_ratio", lambda t, n, r: ((t[5][0] - t[r[0] + 1][0]) / t[5][0]) / (((n[5][0] - n[r[0] + 1][0]) / n[5][0]) + 0.00001)),

("tumour_entropy", lambda t, n, r: t[10]),
("normal_entropy", lambda t, n, r: n[10]),
("normal_tumour_entropy", lambda t, n, r: n[10] / (t[10] + 0.00000001)),

("normal_variant_mapq_mean", lambda t, n, r: (n[5][2] - n[r[0] + 1][2]) / (n[5][0] - n[r[0] + 1][0] + 0.00001)),
("tumour_variant_mapq_mean", lambda t, n, r: (t[5][2] - t[r[0] + 1][2]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),

("normal_tumour_variant_direction_ratio", lambda t, n, r: ((t[5][4] - t[r[0] + 1][4]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)) / ((n[5][4] - n[r[0] + 1][4]) / (n[5][0] - n[r[0] + 1][0] + 0.00001) + 0.00001)),
("normal_tumour_variant_mapq_ratio", lambda t, n, r: ((t[5][2] - t[r[0] + 1][2]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)) / ((n[5][2] - n[r[0] + 1][2]) / (n[5][0] - n[r[0] + 1][0] + 0.00001) + 0.00001)),
("normal_tumour_variant_quality_ratio", lambda t, n, r: ((t[5][1] - t[r[0] + 1][1]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)) / ((n[5][1] - n[r[0] + 1][1]) / (n[5][0] - n[r[0] + 1][0] + 0.00001) + 0.00001)),
("normal_tumour_variant_distance_ratio", lambda t, n, r: ((t[5][3] - t[r[0] + 1][3]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)) / ((n[5][3] - n[r[0] + 1][3]) / (n[5][0] - n[r[0] + 1][0] + 0.00001) + 0.00001)),


("normal_tumour_direction_ratio", lambda t, n, r: (t[5][4] / t[5][0]) / ((n[5][4] / n[5][0]) + 0.00001)),
("normal_tumour_mapq_ratio", lambda t, n, r: (t[5][2] / t[5][0]) / ((n[5][2] / n[5][0]) + 0.00001)),
("normal_tumour_distance_ratio", lambda t, n, r: (t[5][3] / t[5][0]) / ((n[5][3] / n[5][0]) + 0.00001)),
("normal_tumour_quality_ratio", lambda t, n, r: (t[5][1] / t[5][0]) / ((n[5][1] / n[5][0]) + 0.00001)),


# Unscaled features
#("tumour_mapq", lambda t, n, r: t[5][2]),
#("normal_mapq", lambda t, n, r: n[5][2]),
#("normal_variant_mapq", lambda t, n, r: (n[5][2] - n[r[0] + 1][2])),
#("tumour_variant_mapq", lambda t, n, r: (t[5][2] - t[r[0] + 1][2])),
)

coverage_features = (
("normal_depth_coverage", lambda t, n, c: n[5][0] / c[0]),
("tumour_depth_coverage", lambda t, n, c: t[5][0] / c[1]),
("tumour_contamination", lambda t, n, c: c[2] / 100),
("whole_genome", lambda t, n, c: c[3])
)

