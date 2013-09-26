name = "TCGA Benchmark 4 featureset with coverage info"
version = "5_single"

feature_set = (
("tumour_indels", lambda t, r: t[9] / t[5][0]),
#("tumour_depth", lambda t, r: t[5][0]),
("tumour_entropy", lambda t, r: t[10]),


("tumour_mapq_ratio", lambda t, r: t[5][2] / t[5][0]),
("tumour_quality_ratio", lambda t, r: t[5][1] / t[5][0]),
("tumour_distance_ratio", lambda t, r: t[5][3] / t[5][0]),
("tumour_direction_ratio", lambda t, r: t[5][4] / t[5][0]),
#("tumour_tail_distance", lambda t, r: t[5][3] / t[5][0]),
        
                
("tumour_ref_depth", lambda t, r: t[r[0] + 1][0] / t[5][0]),
("tumour_ref_quality", lambda t, r: t[r[0] + 1][1] / (t[r[0]+1][0]+ 0.00001)),
("tumour_ref_mapq", lambda t, r: t[r[0] + 1][2] / (t[r[0]+1][0]+ 0.00001)), 
("tumour_ref_direction_total", lambda t, r: t[r[0] + 1][4] / t[5][0]),
("tumour_ref_direction", lambda t, r: t[r[0] + 1][4] / (t[r[0] + 1][0] + 0.00001)),


("region_entropy", lambda t, r: r[4]),
("region_gc_content", lambda t, r: r[3]),
("homopolymer_f", lambda t, r: r[1]),
("homopolymer_b", lambda t, r: r[2]),

("tumour_variant_quality_ratio", lambda t, r: ((t[5][1] - t[r[0] + 1][1]) / (t[5][0] - t[r[0] + 1][0] + 0.00001))),
#("tumour_variant_quality", lambda t, r: (t[5][1] - t[r[0] + 1][1])),
("tumour_variant_direction_ratio", lambda t,  r: (t[5][4] - t[r[0] + 1][4]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),
("tumour_variant_distance", lambda t, r: (t[5][3] - t[r[0] + 1][3]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),
("tumour_variant_depth_ratio", lambda t, r: ((t[5][0] - t[r[0] + 1][0]) / t[5][0])),
("tumour_variant_mapq_mean", lambda t, r: (t[5][2] - t[r[0] + 1][2]) / (t[5][0] - t[r[0] + 1][0] + 0.00001)),
)

coverage_features = (
#("normal_depth_coverage", lambda t, n, c: n[5][0] / c[0]),
("tumour_depth_coverage", lambda t, c: t[5][0] / c[1]),
("tumour_contamination", lambda t, c: c[2] / 100),
("whole_genome", lambda t, c: c[3])
)

