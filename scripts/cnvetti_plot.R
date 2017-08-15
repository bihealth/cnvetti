#!/usr/bin/env Rscript

library(argparser, quietly=TRUE);
library(futile.logger, quietly=TRUE);
library(tidyr, quietly=TRUE);
library(plotrix, quietly=TRUE);
library(maptools, quietly=TRUE);
gpclibPermit();
library(readr, quietly=TRUE);
library(plyr, quietly=TRUE);
suppressPackageStartupMessages(library(dplyr));
library(Cairo, quietly=TRUE);

# Default arguments for scripting.
argv = list(
    reference_fai = '/home/mholtgre/bih_cluster/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa.fai',
    samples = 'P33-F,P33-T',
    input_bed = '/vol/sshfs/mholtgre/bih_cluster/fast/groups/cubi/projects/2017-05-04_Schmuck/GRCh37/somatic_wgs_cnv_calling/covs.bed',
    out_chrom = '/tmp/cnv_chrom%s.png',
    out_genome = '/tmp/cnv_genome.png'
);

# Parse command line.

p = arg_parser('CNVetti plotting') %>%
    add_argument('--reference-fai', help='Path to reference FAI file.') %>%
    add_argument('--out-chrom', help='Path to chromsome output PNG file, "%s" will be replaced by chromosome name.') %>%
    add_argument('--out-genome', help='Path to genome output PNG file.') %>%
    add_argument('--input-bed', help='Path to input BED file.') %>%
    add_argument('--samples', help='Comma-separated list of samples to plot for.') %>%
    add_argument('--gene-beds', help='Comma-separated list of gene BED files.') %>%
    add_argument('--axis-y-min', help='Lower bound of y axis', default=-4) %>%
    add_argument('--axis-y-max', help='Lower bound of y axis', default=4);

argv <- parse_args(p);

flog.info('Arguments are as follows');
str(argv);

for (name in c('reference_fai', 'out_chrom', 'out_genome', 'input_bed', 'samples')) {
    if (is.na(argv[[name]])) {
        print(p);
        stop(sprintf('Argument %s is empty!', name));
    }
}

if (is.na(argv$samples)) {
    stop('Argument to --samples is empty');
} else {
    argv$samples = strsplit(argv$samples, ',')[[1]];
}
gene_files = c();
if (!is.na(argv$gene_beds)) {
    gene_files = strsplit(argv$gene_beds, ',')[[1]];
}

# Define chromosomes.
CHROMS = c(as.character(1:22), "X", "Y")

# Load chromosome lengths.
flog.info('Loading FAI file %s', argv$reference_fai);
chrom_lens = read_tsv(argv$reference_fai, col_names = c('chrom', 'length', 'offset', 'line_chars', 'line_bytes'), col_types = 'cdcii') %>%
    select(chrom, length);
chrom_offsets = data.frame(
    chrom=chrom_lens$chrom,
    length=chrom_lens$length,
    offset=c(0, head(cumsum(chrom_lens$length), -1)),
    end=c(tail(cumsum(chrom_lens$length), -1), 0)) %>%
    mutate(label_pos = offset + length / 2.0);
rownames(chrom_offsets) = chrom_offsets$chrom;


# Load gene BEDs.
genes = list();
for (gene_bed in gene_files) {
    genes[[gene_bed]] = read_tsv(gene_bed, col_names = c('chrom', 'begin', 'end', 'gene'), col_types = 'ciic');
}
genes = do.call(rbind, genes);
if (is.null(genes)) {
    genes = data.frame(x=c());
}
genes$offset = chrom_offsets[genes$chrom,]$offset;
genes$apos = (genes$offset + genes$begin + (genes$end - genes$begin) / 2) * 1e-6;
genes$pos = (genes$begin + (genes$end - genes$begin) / 2) * 1e-6;
genes$abegin = (genes$offset + genes$begin) * 1e-6;
genes$aend = (genes$offset + genes$end) * 1e-6;

# Prepare loading of mass data.
col_names = c('chrom', 'begin', 'end');
col_types = 'cdd';
for (sample in argv$samples) {
    col_names = c(col_names, paste(sample, 'cov', sep = '.'))
    col_names = c(col_names, paste(sample, 'seg', sep = '.'))
    col_types = paste0(col_types, 'dd');
}

flog.info('Loading TSV file %s', argv$input_bed);
raw_data = read_tsv(argv$input_bed, col_names = col_names, col_types = col_types, na = c('', 'NA', '.')) %>%
    mutate(pos = (begin + (end - begin) / 2));

log2_data = bind_cols(select(raw_data, chrom, begin, end, pos), log2(select(raw_data, -chrom, -begin, -end, -pos)));

data_slices = list();
for (sample in argv$samples) {
    slice = log2_data %>%
        select(
            chrom, pos,
            ncov = starts_with(paste(sample, 'cov', sep = '.')),
            seg = starts_with(paste(sample, 'seg', sep = '.'))) %>%
        mutate(
            sample = sample
        );
    slice$offset = chrom_offsets[slice$chrom,]$offset;
    slice = slice %>%
        mutate(apos = pos + offset) %>%  # absolute pos
        mutate(  # to Mbp
            pos = pos * 1e-6,
            apos = apos * 1e-6
        ) %>%
        select(-offset);
    data_slices[[sample]] = slice;
}
all_data = do.call(rbind, data_slices);
all_data$sample = factor(all_data$sample, levels = argv$samples);

nsamples = length(argv$samples);

cols = list(
    ncov = "lightgray",
    seg = "red",
    extremes = "magenta",
    gene_region = "#DFFF00AA",
    gene_name = "black",
    grid = "darkgray"
);
xlimits = list(
    min = min(chrom_offsets$offset) * 1e-6,
    max = max(chrom_offsets$offset + chrom_offsets$length) * 1e-6
);
ylimits = list(
    min = argv$axis_y_min,
    max = argv$axis_y_max
);

# genome-wide plotting

flog.info("Plotting genome to %s", argv$out_genome);

CairoPNG(
    sprintf(argv$out_genome),
    width = 1200 * 4,
    height = 300 * nsamples * 4,
    res = 64 * 4
);

par(omi = rep(1.0, 4), mar = c(0,0,0,0), mfrow = c(nsamples, 1));

for (s in argv$samples) {
    # Limit data to current sample for plotting.
    flog.info("  => %s", s);
    sub_df = all_data %>%
        filter(as.character(sample) == s & !is.na(ncov) & !is.na(seg)) %>%
        mutate(rounded_seg = as.integer(round(seg * 10))) %>%
        mutate(seg_diff = ifelse(rounded_seg - lag(rounded_seg) != 0, 1, 0)) %>%
        mutate(seg_diff = ifelse(is.na(seg_diff), 0, seg_diff)) %>%
        mutate(seg_no = cumsum(seg_diff)) %>%
        select(-seg_diff, -rounded_seg);
    # Compute segment medians and join with original value
    seg_medians = sub_df %>%
        group_by(seg_no) %>%
        summarise(seg_median = median(seg));
    sub_df = left_join(sub_df, seg_medians, by = "seg_no");

    # Draw coverage as point cloud
    flog.info("    => coverage");
    plot(
        ncov ~ apos,
        data = sub_df,
        col = cols$ncov,
        bg = cols$ncov,
        axes = FALSE,
        cex = 0.1,
        pch = 22,
        xlim = as.double(xlimits),
        ylim = as.double(ylimits),
        xaxs='i'
    );
    mtext(s, 2, 4);
    mtext(expression('log2(coverage ratio)'), 2, 2);
    box();
    axis(2);

    # Overlay coverage point cloud with segmentation
    flog.info("    => segmentation");
    points(
        seg_median ~ apos,
        data = sub_df,
        col = cols$seg,
        bg = cols$seg,
        cex = 0.5,
        pch = 22
    );

    # Draw marker for too small and too large values
    flog.info("    => markers");
    too_small = sub_df %>% filter(seg < ylimits$min);
    too_small$seg = rep(ylimits$min, length(too_small$seg));
    points(
        seg ~ apos,
        data = too_small,
        col = cols$extremes,
        bg = cols$extremes,
        cex = 1,
        pch = 25);
    too_small = sub_df %>% filter(seg > ylimits$max);
    too_small$seg = rep(ylimits$max, length(too_small$seg));
    points(
        seg ~ pos,
        data = too_small,
        col = cols$extremes,
        bg = cols$extremes,
        cex = 1,
        pch = 24);

    for (offset in tail(chrom_offsets$offset, -1)) {
        abline(v = offset * 1e-6);
    }

    if (nrow(genes) > 0) {
        flog.info('    => %d genes', nrow(genes));

        ys = head(rep(c(-2, -2.5, -3, -3.5, -4), length(genes$apos)), length(genes$apos));
        rect(
            xleft = genes$abegin,
            ybottom = ylimits$min,
            xright = genes$aend,
            ytop = ylimits$max,
            density = NA,
            col = cols$gene_region
        );
        pointLabel(
            genes$apos,
            ys,
            labels = genes$gene,
            col = cols$gene_name);
    } else {
        flog.info('    => no genes');
    }


    # title and bottom label
    if (s == argv$samples[[1]]) {
        mtext(sprintf("whole genome", c), 3, 2);
    }
    if (s == argv$samples[[nsamples]]) {
        axis(1, at=chrom_offsets$label_pos * 1e-6, labels=chrom_offsets$chrom, tick=FALSE, line=NA);
        mtext("chromosomes", 1, 2);
    }
}

dev.off();

# chromosome-wise plotting

for (c in CHROMS) {
    flog.info("Plotting chr%s to %s", c, sprintf(argv$out_chrom, c));

    CairoPNG(
        sprintf(argv$out_chrom, c),
        width = 1200 * 4,
        height = 300 * nsamples * 4,
        res = 64 * 4
    );

    par(omi = rep(1.0, 4), mar = c(0,0,0,0), mfrow = c(nsamples, 1));

    for (s in argv$samples) {
        # Limit data to current chromosome.
        flog.info("  => %s", s);
        sub_df = all_data %>%
            filter(chrom == c & as.character(sample) == s & !is.na(ncov) & !is.na(seg)) %>%
            mutate(rounded_seg = as.integer(round(seg * 10))) %>%
            mutate(seg_diff = ifelse(rounded_seg - lag(rounded_seg) != 0, 1, 0)) %>%
            mutate(seg_diff = ifelse(is.na(seg_diff), 0, seg_diff)) %>%
            mutate(seg_no = cumsum(seg_diff)) %>%
            select(-seg_diff, -rounded_seg);
        # Compute segment medians and join with original value
        seg_medians = sub_df %>%
            group_by(seg_no) %>%
            summarise(seg_median = median(seg));
        sub_df = left_join(sub_df, seg_medians, by = "seg_no");

        # Draw coverage as point cloud
        plot(
            ncov ~ pos,
            data = sub_df,
            col = cols$ncov,
            bg = cols$ncov,
            axes = FALSE,
            cex = 0.1,
            pch = 22,
            xlim = c(0, chrom_offsets[c,]$length * 1e-6),
            ylim = as.double(ylimits),
            xaxs='i'
        );
        mtext(s, 2, 4);
        mtext("log2(relative coverage)", 2, 2);
        box();
        axis(2);
        # Overlay coverage point cloud with segmentation
        points(
            seg ~ pos,
            data = sub_df,
            col = cols$seg,
            bg = cols$seg,
            cex = 0.2,
            pch = 22
        );

        # Draw marker for too small and too large values
        too_small = sub_df %>% filter(seg < ylimits$min);
        too_small$seg = rep(ylimits$min, length(too_small$seg));
        points(
            seg ~ pos,
            data = too_small,
            col = cols$extremes,
            bg = cols$extremes,
            cex = 1,
            pch = 25);
        too_small = sub_df %>% filter(seg > ylimits$max);
        too_small$seg = rep(ylimits$max, length(too_small$seg));
        points(
            seg ~ pos,
            data = too_small,
            col = cols$extremes,
            bg = cols$extremes,
            cex = 1,
            pch = 24);

        if (s == argv$samples[[1]]) {
            mtext(sprintf("chr%s", c), 3, 2);
        }
        if (s == argv$samples[[nsamples]]) {
            axis(1);
            mtext("position [Mbp]", 1, 2);
        }

        chrom_genes = genes %>% filter(chrom == c);
        if (nrow(chrom_genes) > 0) {
            flog.info('    => %d genes', nrow(chrom_genes));

            ys = head(rep(c(-2, -2.5, -3, -3.5, -4), length(chrom_genes$pos)), length(chrom_genes$pos));
            rect(
                xleft = chrom_genes$begin * 1e-6,
                ybottom = ylimits$min,
                xright = chrom_genes$end * 1e-6,
                ytop = ylimits$max,
                density = NA,
                col = cols$gene_region
            );
            pointLabel(
                chrom_genes$pos,
                ys,
                labels = chrom_genes$gene,
                col = cols$gene_name);
        } else {
            flog.info('    => no genes');
        }
    }

    dev.off();
}
