// Interfacing with R.

// TODO: write out plots?

pub const LOESS_R: &'static str = r#"
# CNVetti supporting R script
#
# LOESS-fit for coverage and mapability.
#
# Based on the code by Kuilman et al., (2015)
#
# - https://github.com/PeeperLab/CopywriteR/
#
# That code is GPLv3-licensed, so is this.

write("Starting LOESS", stderr());

# Load coverage data; header is:
#
# - use_loess
# - gc_content
# - mapability
# * foreach sample
#   - coverage

write("Loading data", stderr());

df = read.table(
    "{{input_file}}",
    header = TRUE,
    sep = "\t",
    na.strings = c("NA", "."),
    colClasses = c("logical", "numeric", "numeric", "numeric")
);

# First Step: GC-content normalization ----------------------------------------

write("Performing GC-based LOESS", stderr());

# Vector of logical values describing whether value is to be used in GC
# normalization.

use_gc = (
    df$use_loess &
    !is.na(df$mapability) &
    (df$mapability > 0.8) &
    !is.na(df$gc_content) & (df$gc_content > 0)
);

# Perform initial "rough" fit.

rough = loess(
    formula = count ~ gc_content,
    data = df,
    subset = use_gc,
    span = 0.03
);

# Refine rough fit into "final" fit.

i = seq(0, 1, by = 0.001);

final = loess(
    predict(rough, i) ~ i,
    span = 0.3
);

# Then, predict based on input GC content.

normv = predict(final, df$gc_content);

df$count_loess_gc = df$count / (normv / median(normv, na.rm = TRUE));

# DISABLED code for plotting
#
# plot(count ~ gc_content, data = df, subset = use_gc,
#      ylim = quantile(df$count[use_gc], c(1e-04, 0.999)), xlim = c(0, 1),
#      pch = ".")
# points(count ~ gc_content, data = df, subset = !use_gc, col = rgb(1, 0, 0, 0.3),
#        pch = ".")
# lines(i, predict(rough, i), col = "green")
# points(df$gc_content, normv, col = "red", pch = ".")

if ({{loess_mapability}}) {
    # Second Step: mapability-based normalization -----------------------------

    write("Performing mapability-based LOESS", stderr());

    # Logicals describing whether to use value in mapability normalization.

    use_map = !is.na(df$mapability);

    # Perform initial rough fit.

    rough = loess(
        count_loess_gc ~ mapability,
        data = df,
        subset = use_map,
        span = 0.03
    );

    # Refine rough fit into "final" fit.

    i = seq(0, 1, by = 0.001);

    final = loess(
        predict(rough, i) ~ i,
        span = 0.3
    );

    # Then, predict based on input GC content.

    normv = predict(final, df$mapability);

    df$count_loess_map = df$count_loess_gc / (
        normv / median(normv, na.rm = TRUE));

    # DISABLED code for plotting
    #
    # plot(countgcloess ~ mapability, data = df, subset = use_map,
    #      ylim = quantile(df$countgcloess[use_map], c(1e-04, 0.999),
    #                      na.rm = TRUE), xlim = c(0, 1), pch = ".")
    # points(countgcloess ~ mapability, data = df, subset = !use_map,
    #        col = rgb(1, 0, 0, 0.3), pch = ".")
    # lines(i, predict(rough, i), col = "green")
    # forplot = df;
    # forplot = forplot[!is.na(forplot$mapability) && !is.na(normv)];
    # points(forplot$mapability, normv[!is.na(forplot$mapability) && !is.na(normv)],
    #        col = "red", pch = ".")

    df$count_loess_map = df$count_loess_map / median(
        df$count_loess_map[use_gc], na.rm = TRUE);
    if ({{LOG2_TRANSFORM}}) {
        df$count_loess_map = log2(df$count_loess_map);
    }
}

df$count_loess_gc = df$count_loess_gc / median(
    df$count_loess_gc[use_gc], na.rm = TRUE);
if (!{{LOG2_TRANSFORM}}) {
    df$count_loess_gc = log2(df$count_loess_gc);
}

write(sprintf("MAD GC    : %.4f",
    mad(df$count_loess_gc[use_gc], na.rm = TRUE)), stderr());
if (!{{LOG2_TRANSFORM}}) {
    write(sprintf("MAD GC    : %.4f (log2)",
        mad(log2(df$count_loess_gc[use_gc]), na.rm = TRUE)), stderr());
}
if ({{loess_mapability}}) {
    write(
        sprintf("MAD GC+MAP: %.4f",
            mad(df$count_loess_map[use_gc], na.rm = TRUE)),
        stderr());
    if (!{{LOG2_TRANSFORM}}) {
        write(
            sprintf("MAD GC+MAP: %.4f (log2)",
                mad(log2(df$count_loess_map[use_gc]), na.rm = TRUE)),
            stderr());
    }
}

# Final Step: Write out results -----------------------------------------------

write.table(
    df,
    "{{output_file}}",
    sep = "\t",
    row.names = FALSE,
)

write("All done with loess.R", stderr());
"#;
