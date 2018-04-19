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
    #[serde(with = "map_as_pairs")]
    pub summaries: BTreeMap<(i32, i32), DistributionSummary>,
}

// Coverage summary.
pub struct CoverageSummarizer {
    // The name of the metric.
    metric: String,
    // The GC bucket size in per mille.
    gc_bucket_size: i32,
    // The GC bucket size as float.
    gc_bucket_size_float: f32,
    // The mapability bucket size in per mille.
    map_bucket_size: i32,
    // The mapability bucket size as float.
    map_bucket_size_float: f32,
    // Samples in original order.
    samples: Vec<String>,
    // Histograms to use for summarizing, indexed by (gc_bucket, map_bucket), if no mapability
    // is given, mapability 1.0 is used.
    histograms: BTreeMap<String, BTreeMap<(i32, i32), Histogram>>,
}

impl CoverageSummarizer {
    /// Construct with `Vec` of sample names.
    pub fn new(
        gc_bucket_size: i32,
        map_bucket_size: i32,
        samples: &Vec<String>,
        metric: &String,
    ) -> CoverageSummarizer {
        if 1000 % gc_bucket_size != 0 {
            panic!("1000 must be divisible by gc_bucket_size");
        }
        if 1000 % map_bucket_size != 0 {
            panic!("1000 must be divisible by map_bucket_size");
        }
        CoverageSummarizer {
            metric: metric.clone(),
            gc_bucket_size: gc_bucket_size,
            gc_bucket_size_float: gc_bucket_size as f32 / 1000.0,
            map_bucket_size: map_bucket_size,
            map_bucket_size_float: map_bucket_size as f32 / 1000.0,
            samples: samples.clone(),
            histograms: BTreeMap::from_iter(samples.iter().map(|sample| {
                (
                    sample.clone(),
                    BTreeMap::from_iter(
                        iproduct!((0..1001).step_by(gc_bucket_size as usize),
                                  (0..1001).step_by(map_bucket_size as usize))
                            .map(|(gc, map)| ((gc, map), Histogram::new())),
                    ),
                )
            })),
        }
    }

    /// Push count into histogram for the given GC content and sample.
    pub fn increment(&mut self, sample: &String, gc: f32, map: f32, coverage: i32) {
        let gc_bucket = (gc / self.gc_bucket_size_float).floor() as i32 * self.gc_bucket_size;
        let map_bucket = (map / self.map_bucket_size_float).floor() as i32 * self.map_bucket_size;
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
                .get_mut(&(gc_bucket, map_bucket))
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


// From: https://github.com/serde-rs/json/issues/402#issuecomment-356039082
mod map_as_pairs {
    use std::fmt;
    use std::marker::PhantomData;
    use serde::ser::{Serialize, Serializer};
    use serde::de::{Deserialize, Deserializer, Visitor, SeqAccess};

    pub fn serialize<K, V, M, S>(map: M, serializer: S) -> Result<S::Ok, S::Error>
    where
        K: Serialize,
        V: Serialize,
        M: IntoIterator<Item = (K, V)>,
        S: Serializer,
    {
        serializer.collect_seq(map)
    }

    pub fn deserialize<'de, K, V, M, D>(deserializer: D) -> Result<M, D::Error>
    where
        K: Deserialize<'de>,
        V: Deserialize<'de>,
        M: Default + Extend<(K, V)>,
        D: Deserializer<'de>,
    {
        struct MapVisitor<K, V, M> {
            keys: PhantomData<K>,
            values: PhantomData<V>,
            map: PhantomData<M>,
        }

        impl<'de, K, V, M> Visitor<'de> for MapVisitor<K, V, M>
        where
            K: Deserialize<'de>,
            V: Deserialize<'de>,
            M: Default + Extend<(K, V)>,
        {
            type Value = M;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a sequence of key-value pairs")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut map = M::default();
                while let Some((k, v)) = seq.next_element()? {
                    map.extend(Some((k, v)));
                }
                Ok(map)
            }
        }

        deserializer.deserialize_seq(MapVisitor {
            keys: PhantomData,
            values: PhantomData,
            map: PhantomData,
        })
    }
}
