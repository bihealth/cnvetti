// Code for collecting summaries on the windows.

use std::collections::BTreeMap;
use std::iter::FromIterator;

use histogram::Histogram;


/// Summary of the coverage distribution of.
#[derive(Debug, Serialize, Deserialize)]
pub struct DistributionSummary {
    /// 5-number description of the distribution (q0, q1, q2, q3, q4).
    pub summary5: Vec<u64>,
    /// Coverage mean.
    pub mean: u64,
    /// Coverage standard deviation
    pub stddev: u64,
}

/// Metric summary (metric name and coverage summaries).
#[derive(Debug, Serialize, Deserialize)]
pub struct SummarisedMetric {
    /// Name of the metric.
    pub metric: String,
    /// Name of the sample.
    pub sample: String,
    /// Coverage summaries (for each bucket).
    pub summaries: BTreeMap<i32, DistributionSummary>,
}


// Coverage summary.
pub struct CoverageSummarizer {
    // The name of the metric.
    metric: String,
    // The bucket size in percent.
    bucket_size: i32,
    // The bucket size as float.
    bucket_size_float: f32,
    // Samples in original order.
    samples: Vec<String>,
    // Histograms to use for summarizing.
    histograms: BTreeMap<String, BTreeMap<i32, Histogram>>,
}


impl CoverageSummarizer {
    /// Construct with `Vec` of sample names.
    pub fn new(bucket_size: i32, samples: &Vec<String>, metric: &String) -> CoverageSummarizer {
        if 100 % bucket_size != 0 {
            panic!("100 must be divisible by bucket_size");
        }
        CoverageSummarizer {
            metric: metric.clone(),
            bucket_size: bucket_size,
            bucket_size_float: (bucket_size as f32 / 100.0),
            samples: samples.clone(),
            histograms: BTreeMap::from_iter(samples.iter().map(|sample| {
                (
                    sample.clone(),
                    BTreeMap::from_iter((0..101).step_by(bucket_size as usize).map(|bucket| {
                        (bucket, Histogram::new())
                    })),
                )
            })),
        }
    }

    /// Push count into histogram for the given GC content and sample.
    pub fn increment(&mut self, sample: &String, gc: f32, coverage: i32) {
        let bucket = (gc / self.bucket_size_float).floor() as i32 * self.bucket_size;
        let hists = &mut self.histograms;
        // println!(
        //     "increment(sample={}, gc={}, coverage={}); bucket={}",
        //     sample,
        //     gc,
        //     coverage,
        //     bucket,
        // );
        if coverage > 0 {
            // Limitation in `histogram.rs`
            hists
                .get_mut(sample)
                .unwrap()
                .get_mut(&bucket)
                .unwrap()
                .increment(coverage as u64)
                .expect("Could not increment counter in histogram");
        }
    }

    /// Create `Vec` of `CoverageSummary` objects summarising the distributions for the individual
    /// samples.
    pub fn summaries(&self) -> Vec<SummarisedMetric> {
        self.samples
            .iter()
            .map(|sample| {
                let hists = &self.histograms[sample];

                SummarisedMetric {
                    metric: self.metric.clone(),
                    sample: sample.clone(),
                    summaries: BTreeMap::from_iter(hists.iter().map(|(bucket, hist)| {
                        (
                            *bucket,
                            DistributionSummary {
                                summary5: vec![
                                    hist.minimum().unwrap_or(0),
                                    hist.percentile(25.0).unwrap_or(0),
                                    hist.percentile(50.0).unwrap_or(0),
                                    hist.percentile(75.0).unwrap_or(0),
                                    hist.maximum().unwrap_or(0),
                                ],
                                mean: hist.mean().unwrap_or(0),
                                stddev: hist.stddev().unwrap_or(0),
                            },
                        )
                    })),
                }
            })
            .collect()
    }
}
